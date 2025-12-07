% run_thickness_cop_constrained_optimization.m
% Multi-objective: Minimize TEC thickness, maximize COP
% Constraint: Keep T_max below a practical limit (default 80 °C)
%
% ALL configuration (bounds, defaults, fixed parameters, BCs) is sourced from:
%   src/config/optimization_variables.m
% Edit that file to change bounds/defaults, or use the LOCAL_OVERRIDE block
% below for quick experiments (e.g., different heat flux or coolant temp).
%
% Uses: Global Optimization Toolbox (gamultiobj)

clear; clc;

%% Paths
script_dir = fileparts(mfilename('fullpath'));
project_root = fullfile(script_dir, '..', '..');
addpath(genpath(fullfile(project_root, 'src')));

fprintf('╔════════════════════════════════════════════════════════════╗\n');
fprintf('║  THICKNESS vs COP (TEMP-CONSTRAINED) MULTI-OBJECTIVE GA   ║\n');
fprintf('╚════════════════════════════════════════════════════════════╝\n\n');

%% Load centralized configuration
[var_names, lb, ub, x0, all_vars, CONFIG] = optimization_variables();
nvars = length(var_names);

%% Quick-change overrides (leave [] to use values from optimization_variables.m)
LOCAL_OVERRIDE = struct();
LOCAL_OVERRIDE.T_limit_C    = 80;     % Temperature ceiling (°C)
LOCAL_OVERRIDE.q_flux_W_m2  = [];     % e.g., 8e4 to test higher heat flux
LOCAL_OVERRIDE.h_conv_W_m2K = [];     % e.g., 5e5 for weaker convection
LOCAL_OVERRIDE.T_water_K    = [];     % e.g., 310 for warmer coolant

if ~isempty(LOCAL_OVERRIDE.q_flux_W_m2)
    CONFIG.q_flux_W_m2 = LOCAL_OVERRIDE.q_flux_W_m2;
end
if ~isempty(LOCAL_OVERRIDE.h_conv_W_m2K)
    CONFIG.h_conv_W_m2K = LOCAL_OVERRIDE.h_conv_W_m2K;
end
if ~isempty(LOCAL_OVERRIDE.T_water_K)
    CONFIG.T_water_K = LOCAL_OVERRIDE.T_water_K;
end
CONFIG.T_limit_C = LOCAL_OVERRIDE.T_limit_C;

%% Output directory
timestamp = datestr(now, 'yyyy-mm-dd_HH-MM-SS');
OUTPUT_DIR = fullfile(project_root, 'output', 'multiobjective_optimization', 'thickness_cop_constrained', timestamp);
if ~exist(OUTPUT_DIR, 'dir')
    mkdir(OUTPUT_DIR);
end

%% Display configuration
fprintf('Optimization Variables (%d enabled):\n', nvars);
fprintf('═══════════════════════════════════════════════════════════════\n');
fprintf('  %-35s  %10s  %10s  %10s\n', 'Variable', 'Lower', 'Upper', 'Initial');
fprintf('───────────────────────────────────────────────────────────────\n');
for i = 1:nvars
    fprintf('  %-35s  %10.4f  %10.4f  %10.4f\n', var_names{i}, lb(i), ub(i), x0(i));
end
fprintf('═══════════════════════════════════════════════════════════════\n');
fprintf('\nFixed/Integer Parameters (Thermal_Network_For_Radial_TEC.tex notation):\n');
fprintf('  N (stages):   %d\n', CONFIG.N);
fprintf('  M (wedges):   %d  →  θ = %.1f°\n', CONFIG.M, CONFIG.theta_deg);
fprintf('  T_limit:      %.0f °C (constraint)\n', CONFIG.T_limit_C);
fprintf('  q_flux:       %.2e W/m² (%.0f kW/m²)\n', CONFIG.q_flux_W_m2, CONFIG.q_flux_W_m2/1e3);
fprintf('  h_conv:       %.2e W/m²K\n', CONFIG.h_conv_W_m2K);
fprintf('  T_water:      %.1f K (%.1f °C)\n', CONFIG.T_water_K, CONFIG.T_water_K - 273.15);
fprintf('\nObjectives:\n');
fprintf('  1. Minimize TEC thickness (µm)\n');
fprintf('  2. Maximize COP = Q_in / P_elec (implemented as minimize -COP)\n');
fprintf('Constraint:\n');
fprintf('  • T_{max} ≤ %.0f °C at the chip nodes\n', CONFIG.T_limit_C);
fprintf('───────────────────────────────────────────────────────────────\n\n');

%% Base configuration
base_config = create_base_config(CONFIG);

%% Objective and constraint
objective_fn = @(x) objective_thickness_cop(x, base_config, CONFIG);
temp_constraint_fn = @(x) temperature_constraint(x, base_config, CONFIG);

%% GA options (live Pareto plot + live COP/Thickness tracker)
options = optimoptions('gamultiobj', ...
    'PopulationSize', 100, ...
    'MaxGenerations', 100, ...
    'ParetoFraction', 0.35, ...
    'Display', 'iter', ...
    'PlotFcn', {@gaplotpareto}, ...
    'OutputFcn', @live_plot_output, ...
    'UseParallel', false);

fprintf('=== RUNNING gamultiobj (Thickness ↓, COP ↑, T_max constrained) ===\n\n');

%% Run optimizer
opt_start = tic;
[x_pareto, fval_pareto, exitflag, output] = ...
    gamultiobj(objective_fn, nvars, [], [], [], [], lb, ub, temp_constraint_fn, options);
opt_time = toc(opt_start);

fprintf('\nOptimization completed in %.1f seconds\n', opt_time);
fprintf('Found %d Pareto-optimal solutions (feasible wrt T_max)\n\n', size(x_pareto, 1));

%% Post-process Pareto set
num_pts = size(x_pareto, 1);
metrics = repmat(struct('thickness_um', NaN, 'T_max_C', NaN, 'Power_W', NaN, ...
    'Q_cool_W', NaN, 'COP', NaN, 'T_profile', [], 'config', [], 'valid', false), num_pts, 1);
for i = 1:num_pts
    metrics(i) = evaluate_design(x_pareto(i, :), base_config, CONFIG, true); % skip cache reuse
end

thickness_um = fval_pareto(:, 1);
COP_vals = -fval_pareto(:, 2);
T_max_vals = [metrics.T_max_C]';
Power_vals = [metrics.Power_W]';
Q_cool_vals = [metrics.Q_cool_W]';

% Identify key solutions
[~, idx_min_thickness] = min(thickness_um);
[~, idx_max_COP] = max(COP_vals);

% Knee point using distance to utopia (min thickness, max COP)
utopia = [min(thickness_um), max(COP_vals)];
norm_thk = (thickness_um - utopia(1)) / max(1e-12, range(thickness_um));
norm_cop = (COP_vals - utopia(2)) / max(1e-12, range(COP_vals));
distances = sqrt(norm_thk.^2 + norm_cop.^2);
[~, idx_knee] = min(distances);

fprintf('Key Solutions on Pareto Front:\n');
fprintf('───────────────────────────────────────────────────────────────────────\n');
fprintf('%-12s | %10s | %10s | %10s | %10s\n', 'Solution', 't (µm)', 'COP', 'T_max (°C)', 'P_elec (W)');
fprintf('───────────────────────────────────────────────────────────────────────\n');
fprintf('%-12s | %10.1f | %10.3f | %10.2f | %10.4f\n', 'Min t', thickness_um(idx_min_thickness), COP_vals(idx_min_thickness), T_max_vals(idx_min_thickness), Power_vals(idx_min_thickness));
fprintf('%-12s | %10.1f | %10.3f | %10.2f | %10.4f\n', 'Max COP', thickness_um(idx_max_COP), COP_vals(idx_max_COP), T_max_vals(idx_max_COP), Power_vals(idx_max_COP));
fprintf('%-12s | %10.1f | %10.3f | %10.2f | %10.4f\n', 'Balanced', thickness_um(idx_knee), COP_vals(idx_knee), T_max_vals(idx_knee), Power_vals(idx_knee));
fprintf('───────────────────────────────────────────────────────────────────────\n\n');

% Plot temperature profile for balanced solution using ResultsManager
try
    best_metrics = metrics(idx_knee);
    rm = ResultsManager(best_metrics.config);
    rm.plot_temperature_profile(best_metrics.T_profile, best_metrics.config.geometry, ...
        'temp_profile_best.png', 'Balanced Pareto', '', best_metrics.config.boundary_conditions.T_water_K);
catch ME
    fprintf('⚠ Could not plot temperature profile: %s\n', ME.message);
end

%% Plots
figure('Position', [100, 100, 1200, 500], 'Name', 'Pareto: Thickness vs COP');
subplot(1,2,1);
scatter(thickness_um, COP_vals, 50, T_max_vals, 'filled', 'MarkerFaceAlpha', 0.75);
colorbar; ylabel(colorbar, 'T_{max} (°C)');
xlabel('TEC thickness (µm)');
ylabel('COP = Q_{in} / P_{elec}');
title('Pareto Front (feasible: T_{max} ≤ limit)');
grid on;
xline(thickness_um(idx_min_thickness), 'k--', 'Min thickness');

subplot(1,2,2);
scatter(thickness_um, T_max_vals, 50, COP_vals, 'filled', 'MarkerFaceAlpha', 0.75);
colorbar; ylabel(colorbar, 'COP');
xlabel('TEC thickness (µm)');
ylabel('T_{max} (°C)');
title(sprintf('Temperature vs Thickness (limit = %.0f °C)', CONFIG.T_limit_C));
yline(CONFIG.T_limit_C, 'r--', 'T_{limit}');
grid on;

saveas(gcf, fullfile(OUTPUT_DIR, 'pareto_thickness_cop.png'));
saveas(gcf, fullfile(OUTPUT_DIR, 'pareto_thickness_cop.fig'));

%% Save results
results = struct();
results.timestamp = timestamp;
results.config = CONFIG;
results.variables = var_names;
results.bounds.lb = lb;
results.bounds.ub = ub;
results.optimization_time = opt_time;
results.exitflag = exitflag;
results.output = output;

results.pareto.x = x_pareto;
results.pareto.fval = fval_pareto;
results.pareto.thickness_um = thickness_um;
results.pareto.COP = COP_vals;
results.pareto.T_max_C = T_max_vals;
results.pareto.Power_W = Power_vals;
results.pareto.Q_cool_W = Q_cool_vals;

results.key.min_thickness = struct('x', x_pareto(idx_min_thickness, :), 'metrics', metrics(idx_min_thickness));
results.key.max_COP = struct('x', x_pareto(idx_max_COP, :), 'metrics', metrics(idx_max_COP));
results.key.balanced = struct('x', x_pareto(idx_knee, :), 'metrics', metrics(idx_knee));

save(fullfile(OUTPUT_DIR, 'thickness_cop_results.mat'), 'results');

valid_var_names = matlab.lang.makeValidName(var_names);
pareto_table = array2table(x_pareto, 'VariableNames', valid_var_names(:)');
pareto_table.t_TEC_um = thickness_um;
pareto_table.COP = COP_vals;
pareto_table.T_max_C = T_max_vals;
pareto_table.Power_W = Power_vals;
pareto_table.Q_cool_W = Q_cool_vals;
writetable(pareto_table, fullfile(OUTPUT_DIR, 'pareto_front_thickness_cop.csv'));

fprintf('\nResults saved to: %s\n', OUTPUT_DIR);
fprintf('\n✓ Thickness vs COP optimization complete!\n');

%% ==================== HELPER FUNCTIONS ====================

function config = create_base_config(CONFIG)
config = struct();

% Geometry (paper notation)
config.geometry.N_stages = CONFIG.N;
config.geometry.M_wedges = CONFIG.M;
config.geometry.w_chip_um = CONFIG.W_chip_um;
config.geometry.t_chip_um = CONFIG.t_chip_um;
config.geometry.wedge_angle_deg = CONFIG.theta_deg;  % θ = 360°/M

% Defaults for optimizable parameters
config.geometry.R_cyl_um = 1000;
config.geometry.thickness_um = 200;
config.geometry.t_ins_um = 10;
config.geometry.radial_expansion_factor = 1.15;
config.geometry.azimuthal_gap_um = 20;
config.geometry.insulation_width_ratio = 0.05;
config.geometry.interconnect_ratio = 0.15;
config.geometry.interconnect_thickness_ratio = 1.0;
config.geometry.interconnect_angle_ratio = 0.16;
config.geometry.outerconnect_ratio = 0.15;
config.geometry.outerconnect_thickness_ratio = 1.0;
config.geometry.outerconnect_angle_ratio = 0.16;

% Boundary conditions
config.boundary_conditions.q_flux_W_m2 = CONFIG.q_flux_W_m2;
config.boundary_conditions.T_water_K = CONFIG.T_water_K;
config.boundary_conditions.h_conv_W_m2K = CONFIG.h_conv_W_m2K;

% Operating conditions
config.operating_conditions.I_current_A = 0.025;

% Materials
config.materials = CONFIG.materials;
end

function f = objective_thickness_cop(x, base_config, CONFIG)
metrics = evaluate_design(x, base_config, CONFIG, false);

if ~metrics.valid
    f = [1e6, 1e6];
    return;
end

f = [metrics.thickness_um, -metrics.COP];
end

function [c, ceq] = temperature_constraint(x, base_config, CONFIG)
metrics = evaluate_design(x, base_config, CONFIG, false);

if ~metrics.valid
    c = 1e3;  % Force infeasible
else
    c = metrics.T_max_C - CONFIG.T_limit_C;
end
ceq = [];
end

function metrics = evaluate_design(x, base_config, CONFIG, skip_cache)
persistent last_x last_metrics call_count
if isempty(call_count)
    call_count = 0;
end
if nargin < 4
    skip_cache = false;
end

% Simple cache to avoid duplicate objective/nonlcon evaluations
if ~skip_cache && ~isempty(last_x) && isequal(x, last_x)
    metrics = last_metrics;
    return;
end

call_count = call_count + 1;

try
    % Fresh config copy per call
    config = struct();
    config.geometry = base_config.geometry;
    config.boundary_conditions = base_config.boundary_conditions;
    config.operating_conditions = base_config.operating_conditions;
    config.materials = base_config.materials;

    all_vars = CONFIG.all_vars;
    enabled_mask = [all_vars{:, 5}];

    var_map = containers.Map();
    x_idx = 1;
    for i = 1:size(all_vars, 1)
        name = all_vars{i, 1};
        if enabled_mask(i)
            var_map(name) = x(x_idx);
            x_idx = x_idx + 1;
        else
            var_map(name) = all_vars{i, 4};
        end
    end

    get_val = @(name, default) get_var_or_default(var_map, name, default);

    % Apply variables
    config.operating_conditions.I_current_A = get_val('I', 0.025);
    config.geometry.thickness_um = get_val('t_TEC_um', 200);
    config.geometry.R_cyl_um = get_val('r_cyl_um', 1000);
    config.geometry.t_ins_um = get_val('t_ins_um', 10);
    config.geometry.radial_expansion_factor = get_val('f_L', 1.15);
    config.geometry.azimuthal_gap_um = get_val('W_az_um', 20);
    config.geometry.insulation_width_ratio = get_val('W_is_ratio', 0.05);
    config.geometry.interconnect_ratio = get_val('f_ic_W', 0.15);
    config.geometry.interconnect_thickness_ratio = get_val('f_ic_t', 1.0);
    config.geometry.interconnect_angle_ratio = get_val('f_ic_beta', 0.16);
    config.geometry.outerconnect_ratio = get_val('f_oc_W', 0.15);
    config.geometry.outerconnect_thickness_ratio = get_val('f_oc_t', 1.0);
    config.geometry.outerconnect_angle_ratio = get_val('f_oc_beta', 0.16);

    materials = MaterialProperties(config);
    geometry = TECGeometry(config);
    network = ThermalNetwork(geometry, materials, config);

    N = geometry.N_stages;
    T_water = config.boundary_conditions.T_water_K;
    T = ones(2*N + 1, 1) * (T_water + 50);  % Warm initial guess

    valid = true;
    Q_out = NaN; Q_in = NaN;
    for iter = 1:120
        T_old = T;
        try
            [T_new, Q_out, Q_in] = network.solve(T);
        catch
            valid = false;
            break;
        end

        if any(isnan(T_new)) || any(isinf(T_new)) || any(T_new < 0)
            valid = false;
            break;
        end

        T = 0.5 * T_new + 0.5 * T;  % Relaxation

        if max(abs(T - T_old)) < 1e-6
            break;
        end
    end

    if ~valid
        metrics = invalid_metrics();
    else
        T_max_C = max(T) - 273.15;
        Power_W = Q_out - Q_in;  % Electrical power input
        Q_cool = Q_in;           % Net cooling power from chip
        COP_val = Q_cool / Power_W;

        if Power_W <= 0 || COP_val <= 0 || isnan(COP_val) || isinf(COP_val)
            metrics = invalid_metrics();
        else
            metrics = struct();
            metrics.thickness_um = config.geometry.thickness_um;
            metrics.T_max_C = T_max_C;
            metrics.Power_W = Power_W;
            metrics.Q_cool_W = Q_cool;
            metrics.COP = COP_val;
            metrics.T_profile = T;
            metrics.config = config;
            metrics.valid = true;
        end
    end

    % Debug for early calls
    if call_count <= 5
        fprintf('  [Thickness-COP] Call %d: t=%.1f µm, T_max=%.1f°C, COP=%.3f, q_flux=%.0e\n', ...
            call_count, metrics.thickness_um, metrics.T_max_C, metrics.COP, config.boundary_conditions.q_flux_W_m2);
    end

catch ME
    if mod(call_count, 500) == 1
        fprintf('  [Thickness-COP] Exception at call %d: %s\n', call_count, ME.message);
    end
    metrics = invalid_metrics();
end

if ~skip_cache
    last_x = x;
    last_metrics = metrics;
end
end

function metrics = invalid_metrics()
metrics = struct('thickness_um', 1e6, 'T_max_C', 1e6, 'Power_W', 1e6, ...
    'Q_cool_W', 0, 'COP', -inf, 'T_profile', [], 'config', [], 'valid', false);
end

function val = get_var_or_default(var_map, name, default)
if var_map.isKey(name)
    val = var_map(name);
else
    val = default;
end
end

function [state, options, optchanged] = live_plot_output(options, state, flag)
persistent h_fig ax1 ax2 gen_hist best_cop_hist min_t_hist
optchanged = false;

switch flag
    case 'init'
        h_fig = figure('Name', 'Live: Thickness vs COP (objectives)', 'NumberTitle', 'off');
        ax1 = subplot(1,2,1, 'Parent', h_fig);
        ax2 = subplot(1,2,2, 'Parent', h_fig);
        gen_hist = [];
        best_cop_hist = [];
        min_t_hist = [];
    case 'iter'
        if isempty(state.Score)
            return;
        end
        thickness = state.Score(:, 1);
        COP_vals = -state.Score(:, 2);

        % Subplot 1: current population
        cla(ax1);
        scatter(ax1, thickness, COP_vals, 30, 'filled');
        xlabel(ax1, 'TEC thickness (µm)');
        ylabel(ax1, 'COP');
        title(ax1, sprintf('Generation %d', state.Generation));
        grid(ax1, 'on');

        % Subplot 2: history of best COP and min thickness
        gen_hist(end+1) = state.Generation;
        best_cop_hist(end+1) = max(COP_vals);
        min_t_hist(end+1) = min(thickness);

        cla(ax2);
        plot(ax2, gen_hist, best_cop_hist, '-ob', 'LineWidth', 1); hold(ax2, 'on');
        plot(ax2, gen_hist, min_t_hist, '-sr', 'LineWidth', 1);
        xlabel(ax2, 'Generation');
        ylabel(ax2, 'Value');
        legend(ax2, {'Best COP', 'Min thickness'}, 'Location', 'best');
        title(ax2, 'Progress');
        grid(ax2, 'on');
        hold(ax2, 'off');

        drawnow limitrate;
    otherwise
        % no-op
end
end
