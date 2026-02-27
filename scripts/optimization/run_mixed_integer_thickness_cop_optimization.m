% run_mixed_integer_thickness_cop_optimization.m
% Mixed-integer multi-objective optimization:
%   • Objective 1: Minimize TEC thickness (t_TEC)
%   • Objective 2: Maximize COP (implemented as minimize -COP)
%   • Constraint: T_max ≤ T_limit_C (default 80 °C)
%   • Integer design vars: N (stages), M (wedges → θ = 360°/M)
%
% Pulls all continuous variable bounds/defaults from:
%   src/config/optimization_variables.m
% Quick-change overrides are provided below for heat flux, convection, coolant temp, T_limit.
%
% Uses: Global Optimization Toolbox (gamultiobj)

clear; clc;

%% Paths
script_dir = fileparts(mfilename('fullpath'));
project_root = fullfile(script_dir, '..', '..');
addpath(genpath(fullfile(project_root, 'src')));

fprintf('╔════════════════════════════════════════════════════════════════════╗\n');
fprintf('║  MIXED-INTEGER: THICKNESS vs COP (T_max-CONSTRAINED) PARETO SEARCH ║\n');
fprintf('╚════════════════════════════════════════════════════════════════════╝\n\n');

%% Load centralized configuration
[var_names, lb, ub, x0, all_vars, CONFIG] = optimization_variables();

% Integer parameter ranges
N_range = [2, 12];                        % stages, matched to stage sweep
M_valid = CONFIG.M_valid;                 % allowed wedge counts (divisors of 360)
M_range = [min(M_valid), max(M_valid)];   % relaxed bounds (snapped to valid inside objective)

% Extend variable vector with integers (placed at the end)
var_names_ext = [var_names; {'N_stages'; 'M_wedges'}];
lb_ext = [lb; N_range(1); M_range(1)];
ub_ext = [ub; N_range(2); M_range(2)];
x0_ext = [x0; CONFIG.N; CONFIG.M];
nvars_ext = numel(var_names_ext);

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
OUTPUT_DIR = fullfile(project_root, 'output', 'multiobjective_optimization', 'mixed_integer_thickness_cop', timestamp);
if ~exist(OUTPUT_DIR, 'dir')
    mkdir(OUTPUT_DIR);
end

%% Display configuration
fprintf('Continuous Variables (%d enabled):\n', numel(var_names));
fprintf('═══════════════════════════════════════════════════════════════\n');
fprintf('  %-35s  %10s  %10s  %10s\n', 'Variable', 'Lower', 'Upper', 'Initial');
fprintf('───────────────────────────────────────────────────────────────\n');
for i = 1:numel(var_names)
    fprintf('  %-35s  %10.4f  %10.4f  %10.4f\n', var_names{i}, lb(i), ub(i), x0(i));
end
fprintf('═══════════════════════════════════════════════════════════════\n');
fprintf('\nInteger Variables:\n');
fprintf('  N_stages ∈ [%d, %d] (rounded)\n', N_range(1), N_range(2));
fprintf('  M_wedges ∈ {%s} (snapped to nearest valid)\n', num2str(M_valid));
fprintf('\nFixed/BC Parameters (Thermal_Network_For_Radial_TEC.tex notation):\n');
fprintf('  T_limit:      %.0f °C (constraint)\n', CONFIG.T_limit_C);
fprintf('  q_flux:       %.2e W/m² (%.0f kW/m²)\n', CONFIG.q_flux_W_m2, CONFIG.q_flux_W_m2/1e3);
fprintf('  h_conv:       %.2e W/m²K\n', CONFIG.h_conv_W_m2K);
fprintf('  T_water:      %.1f K (%.1f °C)\n', CONFIG.T_water_K, CONFIG.T_water_K - 273.15);
fprintf('Objectives:\n');
fprintf('  1) Minimize t_TEC (µm)\n');
fprintf('  2) Maximize COP (minimize -COP)\n');
fprintf('Constraint: T_{max} ≤ %.0f °C at chip nodes\n', CONFIG.T_limit_C);
fprintf('───────────────────────────────────────────────────────────────\n\n');

%% Base configuration
base_config = create_base_config(CONFIG);

%% Objective and constraint handles
objective_fn = @(x) objective_thickness_cop_mi(x, base_config, CONFIG, N_range, M_valid);
constraint_fn = @(x) temperature_constraint_mi(x, base_config, CONFIG, N_range, M_valid);

%% GA options (live plots)
options = optimoptions('gamultiobj', ...
    'PopulationSize', 120, ...
    'MaxGenerations', 120, ...
    'ParetoFraction', 0.35, ...
    'Display', 'iter', ...
    'PlotFcn', {@gaplotpareto}, ...
    'OutputFcn', @live_plot_output, ...
    'UseParallel', false);

fprintf('=== RUNNING gamultiobj (mixed-integer) ===\n\n');
opt_start = tic;
[x_pareto, fval_pareto, exitflag, output] = ...
    gamultiobj(objective_fn, nvars_ext, [], [], [], [], lb_ext, ub_ext, constraint_fn, options);
opt_time = toc(opt_start);

fprintf('\nOptimization completed in %.1f seconds\n', opt_time);
fprintf('Found %d Pareto-optimal solutions (feasible wrt T_max)\n\n', size(x_pareto, 1));

%% Post-process Pareto set
num_pts = size(x_pareto, 1);
metrics = repmat(struct('thickness_um', NaN, 'T_max_C', NaN, 'Power_W', NaN, ...
    'Q_cool_W', NaN, 'COP', NaN, 'N', NaN, 'M', NaN, 'T_profile', [], 'config', [], 'valid', false), num_pts, 1);
for i = 1:num_pts
    metrics(i) = evaluate_design_mi(x_pareto(i, :), base_config, CONFIG, N_range, M_valid, true);
end

thickness_um = fval_pareto(:, 1);
COP_vals = -fval_pareto(:, 2);
T_max_vals = [metrics.T_max_C]';
Power_vals = [metrics.Power_W]';
Q_cool_vals = [metrics.Q_cool_W]';
N_vals = [metrics.N]';
M_vals = [metrics.M]';

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
fprintf('───────────────────────────────────────────────────────────────────────────────────────────────\n');
fprintf('%-12s | %8s | %6s | %6s | %10s | %10s | %10s\n', 'Solution', 't (µm)', 'N', 'M', 'COP', 'T_max (°C)', 'P_elec (W)');
fprintf('───────────────────────────────────────────────────────────────────────────────────────────────\n');
print_key('Min t', idx_min_thickness, thickness_um, COP_vals, T_max_vals, Power_vals, N_vals, M_vals);
print_key('Max COP', idx_max_COP, thickness_um, COP_vals, T_max_vals, Power_vals, N_vals, M_vals);
print_key('Balanced', idx_knee, thickness_um, COP_vals, T_max_vals, Power_vals, N_vals, M_vals);
fprintf('───────────────────────────────────────────────────────────────────────────────────────────────\n\n');

%% Plots
figure('Position', [100, 100, 1200, 500], 'Name', 'Pareto: Thickness vs COP (Mixed-Integer)');
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

saveas(gcf, fullfile(OUTPUT_DIR, 'pareto_thickness_cop_mi.png'));
saveas(gcf, fullfile(OUTPUT_DIR, 'pareto_thickness_cop_mi.fig'));

%% Save results (MAT + CSV) for reporting
results = struct();
results.timestamp = timestamp;
results.config = CONFIG;
results.variables = var_names_ext;
results.bounds.lb = lb_ext;
results.bounds.ub = ub_ext;
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
results.pareto.N = N_vals;
results.pareto.M = M_vals;

results.key.min_thickness = struct('x', x_pareto(idx_min_thickness, :), 'metrics', metrics(idx_min_thickness));
results.key.max_COP = struct('x', x_pareto(idx_max_COP, :), 'metrics', metrics(idx_max_COP));
results.key.balanced = struct('x', x_pareto(idx_knee, :), 'metrics', metrics(idx_knee));

% Plot temperature profile for balanced solution
try
    best_metrics = metrics(idx_knee);
    rm = ResultsManager(best_metrics.config);
    rm.plot_temperature_profile(best_metrics.T_profile, best_metrics.config.geometry, ...
        'temp_profile_best.png', 'Balanced Pareto', '', best_metrics.config.boundary_conditions.T_water_K);
catch ME
    fprintf('⚠ Could not plot temperature profile: %s\n', ME.message);
end

save(fullfile(OUTPUT_DIR, 'mixed_integer_thickness_cop_results.mat'), 'results');

valid_var_names = matlab.lang.makeValidName(var_names_ext);
pareto_table = array2table(x_pareto, 'VariableNames', valid_var_names(:)');
pareto_table.t_TEC_um = thickness_um;
pareto_table.COP = COP_vals;
pareto_table.T_max_C = T_max_vals;
pareto_table.Power_W = Power_vals;
pareto_table.Q_cool_W = Q_cool_vals;
pareto_table.N = N_vals;
pareto_table.M = M_vals;
writetable(pareto_table, fullfile(OUTPUT_DIR, 'pareto_front_mixed_integer.csv'));

fprintf('\nResults saved to: %s\n', OUTPUT_DIR);
fprintf('\n✓ Mixed-integer thickness/COP optimization complete!\n');

%% ==================== HELPER FUNCTIONS ====================

function config = create_base_config(CONFIG)
config = struct();
% Geometry defaults
config.geometry.N_stages = CONFIG.N;
config.geometry.M_wedges = CONFIG.M;
config.geometry.w_chip_um = CONFIG.W_chip_um;
config.geometry.t_chip_um = CONFIG.t_chip_um;
config.geometry.wedge_angle_deg = CONFIG.theta_deg;
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

function f = objective_thickness_cop_mi(x, base_config, CONFIG, N_range, M_valid)
metrics = evaluate_design_mi(x, base_config, CONFIG, N_range, M_valid, false);
if ~metrics.valid
    f = [1e6, 1e6];
else
    f = [metrics.thickness_um, -metrics.COP];
end
end

function [c, ceq] = temperature_constraint_mi(x, base_config, CONFIG, N_range, M_valid)
metrics = evaluate_design_mi(x, base_config, CONFIG, N_range, M_valid, false);
if ~metrics.valid
    c = 1e3;  % force infeasible
else
    c = metrics.T_max_C - CONFIG.T_limit_C;
end
ceq = [];
end

function metrics = evaluate_design_mi(x, base_config, CONFIG, N_range, M_valid, skip_cache)
persistent last_x last_metrics call_count
if isempty(call_count)
    call_count = 0;
end
if nargin < 5
    skip_cache = false;
end

if ~skip_cache && ~isempty(last_x) && isequal(x, last_x)
    metrics = last_metrics;
    return;
end

call_count = call_count + 1;

try
    % Enforce integer parameters
    N_raw = x(end-1);
    M_raw = x(end);
    N_int = min(max(round(N_raw), N_range(1)), N_range(2));
    M_int = snap_to_valid(round(M_raw), M_valid);

    % Build config
    config = struct();
    config.geometry = base_config.geometry;
    config.boundary_conditions = base_config.boundary_conditions;
    config.operating_conditions = base_config.operating_conditions;
    config.materials = base_config.materials;

    % Apply integer geometry updates
    config.geometry.N_stages = N_int;
    config.geometry.M_wedges = M_int;
    config.geometry.wedge_angle_deg = 360 / M_int;

    % Map continuous variables
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
    T = ones(2*N + 1, 1) * (T_water + 50);  % warm initial guess

    valid = true;
    Q_out = NaN; Q_in = NaN;
    for iter = 1:150
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

        T = 0.5 * T_new + 0.5 * T;  % relaxation

        if max(abs(T - T_old)) < 1e-6
            break;
        end
    end

    if ~valid
        metrics = invalid_metrics();
    else
        T_max_C = max(T) - 273.15;
        Power_W = Q_out - Q_in;     % Electrical input power (balance)
        Q_cool = Q_in;              % Net cooling power drawn from chip (per TEC Physics note)
        COP_val = Q_cool / Power_W; % COP = Q_c / P

        if Power_W <= 0 || COP_val <= 0 || isnan(COP_val) || isinf(COP_val)
            metrics = invalid_metrics();
        else
            metrics = struct();
            metrics.thickness_um = config.geometry.thickness_um;
            metrics.T_max_C = T_max_C;
            metrics.Power_W = Power_W;
            metrics.Q_cool_W = Q_cool;
            metrics.COP = COP_val;
            metrics.N = N_int;
            metrics.M = M_int;
            metrics.T_profile = T;
            metrics.config = config;
            metrics.valid = true;
        end
    end

    if call_count <= 5
        fprintf('  [MI Pareto] Call %d: t=%.1f µm, COP=%.3f, T_max=%.1f°C, N=%d, M=%d, q_flux=%.0e\n', ...
            call_count, metrics.thickness_um, metrics.COP, metrics.T_max_C, metrics.N, metrics.M, config.boundary_conditions.q_flux_W_m2);
    end

catch ME
    if mod(call_count, 500) == 1
        fprintf('  [MI Pareto] Exception at call %d: %s\n', call_count, ME.message);
    end
    metrics = invalid_metrics();
end

if ~skip_cache
    last_x = x;
    last_metrics = metrics;
end
end

function M_int = snap_to_valid(M_raw, M_valid)
% Snap M to nearest allowed wedge count
[~, idx] = min(abs(M_valid - M_raw));
M_int = M_valid(idx);
end

function metrics = invalid_metrics()
metrics = struct('thickness_um', 1e6, 'T_max_C', 1e6, 'Power_W', 1e6, ...
    'Q_cool_W', 0, 'COP', -inf, 'N', NaN, 'M', NaN, 'T_profile', [], 'config', [], 'valid', false);
end

function val = get_var_or_default(var_map, name, default)
if var_map.isKey(name)
    val = var_map(name);
else
    val = default;
end
end

function print_key(label, idx, thickness_um, COP_vals, T_max_vals, Power_vals, N_vals, M_vals)
fprintf('%-12s | %8.1f | %6d | %6d | %10.3f | %10.2f | %10.4f\n', ...
    label, thickness_um(idx), N_vals(idx), M_vals(idx), COP_vals(idx), T_max_vals(idx), Power_vals(idx));
end

function [state, options, optchanged] = live_plot_output(options, state, flag)
persistent h_fig ax1 ax2 gen_hist best_cop_hist min_t_hist
optchanged = false;

switch flag
    case 'init'
        h_fig = figure('Name', 'Live: Thickness vs COP (mixed-integer)', 'NumberTitle', 'off');
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

        cla(ax1);
        scatter(ax1, thickness, COP_vals, 30, 'filled');
        xlabel(ax1, 'TEC thickness (µm)');
        ylabel(ax1, 'COP');
        title(ax1, sprintf('Generation %d', state.Generation));
        grid(ax1, 'on');

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
