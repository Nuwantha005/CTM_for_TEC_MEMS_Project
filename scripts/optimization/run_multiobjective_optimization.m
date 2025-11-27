% run_multiobjective_optimization.m
% Multi-objective optimization: Minimize T_max AND Power consumption
%
% This finds the Pareto-optimal front showing trade-offs between
% cooling performance and electrical power consumption.
%
% ALL CONFIGURATION is loaded from: src/config/optimization_variables.m
% Edit that file to change bounds, defaults, fixed parameters, or boundary conditions.
%
% Uses: Global Optimization Toolbox (gamultiobj)

clear; clc;

% Get the directory where this script is located (works both interactively and batch)
script_dir = fileparts(mfilename('fullpath'));
project_root = fullfile(script_dir, '..', '..');

addpath(genpath(fullfile(project_root, 'src')));

fprintf('╔════════════════════════════════════════════════════════════╗\n');
fprintf('║        MULTI-OBJECTIVE OPTIMIZATION FOR TEC DESIGN        ║\n');
fprintf('╚════════════════════════════════════════════════════════════╝\n\n');

%% Load ALL configuration from single control point
% ═══════════════════════════════════════════════════════════════════════════
% ALL parameters (variables, bounds, fixed params, BCs) come from:
%   src/config/optimization_variables.m
% ═══════════════════════════════════════════════════════════════════════════

[var_names, lb, ub, x0, all_vars, CONFIG] = optimization_variables();
nvars = length(var_names);

% Output directory (use absolute path based on script location)
timestamp = datestr(now, 'yyyy-mm-dd_HH-MM-SS');
OUTPUT_DIR = fullfile(project_root, 'output', 'multiobjective_optimization', timestamp);
if ~exist(OUTPUT_DIR, 'dir')
    mkdir(OUTPUT_DIR);
end

fprintf('Optimization Variables (%d enabled):\n', nvars);
fprintf('═══════════════════════════════════════════════════════════════\n');
fprintf('  %-35s  %10s  %10s  %10s\n', 'Variable', 'Lower', 'Upper', 'Initial');
fprintf('───────────────────────────────────────────────────────────────\n');
for i = 1:nvars
    fprintf('  %-35s  %10.4f  %10.4f  %10.4f\n', var_names{i}, lb(i), ub(i), x0(i));
end
fprintf('═══════════════════════════════════════════════════════════════\n');
fprintf('\nFixed Parameters (from optimization_variables.m):\n');
fprintf('  N_stages:     %d\n', CONFIG.N_stages);
fprintf('  T_target:     %.0f °C\n', CONFIG.T_target_C);
fprintf('  q_flux:       %.2e W/m² (%.0f kW/m²)\n', CONFIG.q_flux_W_m2, CONFIG.q_flux_W_m2/1e3);
fprintf('  h_conv:       %.2e W/m²K\n', CONFIG.h_conv_W_m2K);
fprintf('  T_water:      %.1f K (%.1f °C)\n', CONFIG.T_water_K, CONFIG.T_water_K - 273.15);
fprintf('\nObjectives:\n');
fprintf('  1. Minimize T_max (°C)\n');
fprintf('  2. Minimize Power consumption (W)\n');
fprintf('───────────────────────────────────────────────────────────────\n\n');

%% Create base configuration
base_config = create_base_config(CONFIG);

%% Define multi-objective function
% Returns [T_max, Power]

multiobjective = @(x) tec_multiobjective(x, base_config, CONFIG);

%% Run Multi-Objective Genetic Algorithm
fprintf('=== RUNNING gamultiobj ===\n\n');
fprintf('Finding Pareto-optimal front...\n');
fprintf('Population: 100, Generations: 100\n\n');

options = optimoptions('gamultiobj', ...
    'PopulationSize', 100, ...
    'MaxGenerations', 100, ...
    'ParetoFraction', 0.35, ...
    'Display', 'iter', ...
    'PlotFcn', {@gaplotpareto}, ...
    'UseParallel', false);

tic;
[x_pareto, fval_pareto, exitflag, output] = ...
    gamultiobj(multiobjective, nvars, [], [], [], [], lb, ub, options);
opt_time = toc;

fprintf('\nOptimization completed in %.1f seconds\n', opt_time);
fprintf('Found %d Pareto-optimal solutions\n\n', size(x_pareto, 1));

%% Analyze Pareto Front
fprintf('╔════════════════════════════════════════════════════════════╗\n');
fprintf('║                    PARETO FRONT ANALYSIS                   ║\n');
fprintf('╚════════════════════════════════════════════════════════════╝\n\n');

T_max_pareto = fval_pareto(:,1);
Power_pareto = fval_pareto(:,2);

% Find extreme points
[~, idx_min_T] = min(T_max_pareto);
[~, idx_min_P] = min(Power_pareto);

% Find knee point (balanced solution) using distance to utopia point
utopia = [min(T_max_pareto), min(Power_pareto)];
distances = sqrt(((T_max_pareto - utopia(1))/range(T_max_pareto)).^2 + ...
                 ((Power_pareto - utopia(2))/range(Power_pareto)).^2);
[~, idx_knee] = min(distances);

fprintf('Key Solutions on Pareto Front:\n');
fprintf('─────────────────────────────────────────────────────────────────────────────────\n');
fprintf('%-15s | %8s | %8s | %8s | %8s | %8s | %8s\n', ...
    'Solution', 'I (mA)', 't (µm)', 'k_r', 'ff', 'T_max(°C)', 'Power(W)');
fprintf('─────────────────────────────────────────────────────────────────────────────────\n');

% Min Temperature
fprintf('%-15s | %8.1f | %8.0f | %8.3f | %8.3f | %8.2f | %8.4f\n', ...
    'Min T_max', x_pareto(idx_min_T,1)*1000, x_pareto(idx_min_T,2), ...
    x_pareto(idx_min_T,3), x_pareto(idx_min_T,4), ...
    T_max_pareto(idx_min_T), Power_pareto(idx_min_T));

% Knee Point (Balanced)
fprintf('%-15s | %8.1f | %8.0f | %8.3f | %8.3f | %8.2f | %8.4f\n', ...
    'Balanced', x_pareto(idx_knee,1)*1000, x_pareto(idx_knee,2), ...
    x_pareto(idx_knee,3), x_pareto(idx_knee,4), ...
    T_max_pareto(idx_knee), Power_pareto(idx_knee));

% Min Power
fprintf('%-15s | %8.1f | %8.0f | %8.3f | %8.3f | %8.2f | %8.4f\n', ...
    'Min Power', x_pareto(idx_min_P,1)*1000, x_pareto(idx_min_P,2), ...
    x_pareto(idx_min_P,3), x_pareto(idx_min_P,4), ...
    T_max_pareto(idx_min_P), Power_pareto(idx_min_P));

fprintf('─────────────────────────────────────────────────────────────────────────────────\n\n');

%% Plot Pareto Front
figure('Position', [100, 100, 1200, 500], 'Name', 'Pareto Front Analysis');

% Pareto front with key points
subplot(1,2,1);
scatter(T_max_pareto, Power_pareto * 1000, 50, 'b', 'filled', 'MarkerFaceAlpha', 0.6);
hold on;
scatter(T_max_pareto(idx_min_T), Power_pareto(idx_min_T)*1000, 200, 'g', 'filled', 'MarkerEdgeColor', 'k', 'LineWidth', 2);
scatter(T_max_pareto(idx_knee), Power_pareto(idx_knee)*1000, 200, 'm', 'filled', 'MarkerEdgeColor', 'k', 'LineWidth', 2);
scatter(T_max_pareto(idx_min_P), Power_pareto(idx_min_P)*1000, 200, 'r', 'filled', 'MarkerEdgeColor', 'k', 'LineWidth', 2);

xlabel('Maximum Temperature (°C)', 'FontSize', 12);
ylabel('Power Consumption (mW)', 'FontSize', 12);
title('Pareto Front: T_{max} vs Power', 'FontSize', 14);
legend('Pareto solutions', 'Min T_{max}', 'Balanced (Knee)', 'Min Power', 'Location', 'northeast');
grid on;

% Temperature constraint line
xline(85, 'r--', 'T_{max} = 85°C', 'LineWidth', 1.5, 'LabelHorizontalAlignment', 'left');

% Parameter distributions
subplot(1,2,2);
param_labels = {'I (mA)', 't (µm)', 'k_r', 'fill'};
param_data = [x_pareto(:,1)*1000, x_pareto(:,2), x_pareto(:,3), x_pareto(:,4)];

% Normalize for visualization
param_norm = (param_data - min(param_data)) ./ (max(param_data) - min(param_data));

boxplot(param_norm, 'Labels', param_labels);
ylabel('Normalized Value', 'FontSize', 12);
title('Parameter Distributions on Pareto Front', 'FontSize', 14);
grid on;

saveas(gcf, fullfile(OUTPUT_DIR, 'pareto_front.png'));
saveas(gcf, fullfile(OUTPUT_DIR, 'pareto_front.fig'));

%% Additional Analysis: Trade-off Curves
figure('Position', [100, 100, 1200, 800], 'Name', 'Trade-off Analysis');

% T_max vs Current
subplot(2,2,1);
scatter(x_pareto(:,1)*1000, T_max_pareto, 50, Power_pareto*1000, 'filled');
xlabel('Current (mA)');
ylabel('T_{max} (°C)');
title('Temperature vs Current');
colorbar; ylabel(colorbar, 'Power (mW)');
grid on;

% T_max vs Thickness
subplot(2,2,2);
scatter(x_pareto(:,2), T_max_pareto, 50, Power_pareto*1000, 'filled');
xlabel('TEC Thickness (µm)');
ylabel('T_{max} (°C)');
title('Temperature vs Thickness');
colorbar; ylabel(colorbar, 'Power (mW)');
grid on;

% Power vs Current
subplot(2,2,3);
scatter(x_pareto(:,1)*1000, Power_pareto*1000, 50, T_max_pareto, 'filled');
xlabel('Current (mA)');
ylabel('Power (mW)');
title('Power vs Current');
colorbar; ylabel(colorbar, 'T_{max} (°C)');
grid on;

% 3D view
subplot(2,2,4);
scatter3(x_pareto(:,1)*1000, x_pareto(:,2), T_max_pareto, 50, Power_pareto*1000, 'filled');
xlabel('Current (mA)');
ylabel('Thickness (µm)');
zlabel('T_{max} (°C)');
title('3D Parameter Space');
colorbar; ylabel(colorbar, 'Power (mW)');
grid on;
view(45, 30);

saveas(gcf, fullfile(OUTPUT_DIR, 'tradeoff_analysis.png'));
saveas(gcf, fullfile(OUTPUT_DIR, 'tradeoff_analysis.fig'));

%% Save Results
results = struct();
results.timestamp = timestamp;
results.config = CONFIG;
results.variables = var_names;
results.bounds.lb = lb;
results.bounds.ub = ub;
results.optimization_time = opt_time;

results.pareto.x = x_pareto;
results.pareto.fval = fval_pareto;
results.pareto.T_max = T_max_pareto;
results.pareto.Power = Power_pareto;

results.key_solutions.min_temp.x = x_pareto(idx_min_T,:);
results.key_solutions.min_temp.T_max = T_max_pareto(idx_min_T);
results.key_solutions.min_temp.Power = Power_pareto(idx_min_T);

results.key_solutions.balanced.x = x_pareto(idx_knee,:);
results.key_solutions.balanced.T_max = T_max_pareto(idx_knee);
results.key_solutions.balanced.Power = Power_pareto(idx_knee);

results.key_solutions.min_power.x = x_pareto(idx_min_P,:);
results.key_solutions.min_power.T_max = T_max_pareto(idx_min_P);
results.key_solutions.min_power.Power = Power_pareto(idx_min_P);

save(fullfile(OUTPUT_DIR, 'multiobjective_results.mat'), 'results');

% Save Pareto front to CSV with dynamic variable names
valid_var_names = matlab.lang.makeValidName(var_names);
pareto_table = array2table(x_pareto, 'VariableNames', valid_var_names(:)');
pareto_table.T_max_C = T_max_pareto;
pareto_table.Power_mW = Power_pareto * 1000;
writetable(pareto_table, fullfile(OUTPUT_DIR, 'pareto_front.csv'));

fprintf('\nResults saved to: %s\n', OUTPUT_DIR);
fprintf('\n✓ Multi-objective optimization complete!\n');

%% ==================== HELPER FUNCTIONS ====================

function config = create_base_config(CONFIG)
    % Create base config from centralized CONFIG struct
    config = struct();
    
    % Geometry from CONFIG
    config.geometry.N_stages = CONFIG.N_stages;
    config.geometry.N_tsv_limit = CONFIG.N_tsv_limit;
    config.geometry.w_chip_um = CONFIG.w_chip_um;
    config.geometry.t_chip_um = CONFIG.t_chip_um;
    
    % TSV parameters from CONFIG
    config.geometry.tsv.R_TSV_um = CONFIG.tsv.R_TSV_um;
    config.geometry.tsv.P_TSV_um = CONFIG.tsv.P_TSV_um;
    config.geometry.tsv.g_rad_um = CONFIG.tsv.g_rad_um;
    config.geometry.tsv.t_SOI_um = 100;  % Default, can be optimized
    
    % Defaults for optimizable parameters (will be overwritten by optimizer)
    config.geometry.R_cyl_um = 1000;
    config.geometry.wedge_angle_deg = 30;
    config.geometry.thickness_um = 200;
    config.geometry.radial_expansion_factor = 1.15;
    config.geometry.fill_factor = 0.95;
    config.geometry.interconnect_ratio = 0.15;
    config.geometry.outerconnect_ratio = 0.15;
    config.geometry.insulation_width_ratio = 0.04;
    config.geometry.interconnect_angle_ratio = 0.16;
    config.geometry.outerconnect_angle_ratio = 0.16;
    config.geometry.interconnect_thickness_ratio = 1.0;
    config.geometry.outerconnect_thickness_ratio = 1.0;
    
    % Boundary conditions from CONFIG
    config.boundary_conditions.q_flux_W_m2 = CONFIG.q_flux_W_m2;
    config.boundary_conditions.T_water_K = CONFIG.T_water_K;
    config.boundary_conditions.h_conv_W_m2K = CONFIG.h_conv_W_m2K;
    
    % Operating conditions (default, will be optimized)
    config.operating_conditions.I_current_A = 0.025;
    
    % Materials from CONFIG
    config.materials = CONFIG.materials;
end

function f = tec_multiobjective(x, base_config, CONFIG)
    % Multi-objective function
    % Returns: [T_max (°C), Power (W)]
    %
    % Uses same variable names as optimization_variables.m
    % CRITICAL: Uses same solver approach as TECOptimizer:
    % 1. Warm initial guess (T_water + 50)
    % 2. Relaxation iteration (0.5 blend)
    % 3. Proper error handling
    
    try
        config = base_config;
        
        % Get variable names from CONFIG
        all_vars = CONFIG.all_vars;
        enabled_mask = [all_vars{:, 5}];
        
        % Create a map of variable values
        var_map = containers.Map();
        x_idx = 1;
        for i = 1:size(all_vars, 1)
            name = all_vars{i, 1};
            if enabled_mask(i)  % enabled
                var_map(name) = x(x_idx);
                x_idx = x_idx + 1;
            else  % disabled - use initial value
                var_map(name) = all_vars{i, 4};
            end
        end
        
        % Helper function
        get_val = @(name, default) get_var_or_default(var_map, name, default);
        
        % Apply variables to config (same names as optimization_variables.m)
        config.operating_conditions.I_current_A = get_val('current', 0.025);
        config.geometry.thickness_um = get_val('thickness_um', 200);
        config.geometry.wedge_angle_deg = get_val('wedge_angle_deg', 30);
        config.geometry.R_cyl_um = get_val('R_cyl_um', 1000);
        config.geometry.radial_expansion_factor = get_val('k_r', 1.15);
        config.geometry.fill_factor = get_val('fill_factor', 0.95);
        config.geometry.insulation_width_ratio = get_val('insulation_width_ratio', 0.04);
        config.geometry.interconnect_ratio = get_val('interconnect_ratio', 0.15);
        config.geometry.outerconnect_ratio = get_val('outerconnect_ratio', 0.15);
        config.geometry.interconnect_angle_ratio = get_val('interconnect_angle_ratio', 0.16);
        config.geometry.outerconnect_angle_ratio = get_val('outerconnect_angle_ratio', 0.16);
        config.geometry.interconnect_thickness_ratio = get_val('interconnect_thickness_ratio', 1.0);
        config.geometry.outerconnect_thickness_ratio = get_val('outerconnect_thickness_ratio', 1.0);
        if isfield(config.geometry, 'tsv')
            config.geometry.tsv.t_SOI_um = get_val('t_SOI_um', 100);
        end
        
        materials = MaterialProperties(config);
        geometry = TECGeometry(config);
        network = ThermalNetwork(geometry, materials, config);
        
        N = geometry.N_stages;
        T_water = config.boundary_conditions.T_water_K;
        
        % CRITICAL: Warm initial guess (like TECOptimizer)
        T = ones(2*N + 1, 1) * (T_water + 50);
        
        for iter = 1:100
            T_old = T;
            try
                [T_new, Q_out, Q_in] = network.solve(T);
            catch
                f = [1e6, 1e6];
                return;
            end
            
            % Check for invalid values
            if any(isnan(T_new)) || any(isinf(T_new)) || any(T_new < 0)
                f = [1e6, 1e6];
                return;
            end
            
            % CRITICAL: Relaxation (like TECOptimizer)
            T = 0.5 * T_new + 0.5 * T;
            
            if max(abs(T - T_old)) < 1e-6
                break;
            end
        end
        
        T_max_C = max(T) - 273.15;
        Power_W = Q_out - Q_in;  % Electrical power input
        
        % Penalize invalid solutions
        if T_max_C < 0 || T_max_C > 500 || Power_W < 0
            f = [1e6, 1e6];
        else
            f = [T_max_C, Power_W];
        end
        
    catch
        f = [1e6, 1e6];  % Penalty for failed simulations
    end
end

function val = get_var_or_default(var_map, name, default)
    if var_map.isKey(name)
        val = var_map(name);
    else
        val = default;
    end
end
