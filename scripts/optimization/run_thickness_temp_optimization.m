% run_thickness_temp_optimization.m
% Multi-objective optimization: Minimize T_max AND TEC thickness (t_TEC)
%
% Strategy:
%   - Outer loop: Grid search over integer parameters (N, M)
%   - Inner loop: gamultiobj over continuous variables
%
% Objectives:
%   1. Minimize T_max (chip temperature)
%   2. Minimize t_TEC (TEC layer thickness)
%
% No COP constraint - focus purely on thermal performance and compactness.
%
% Uses: Global Optimization Toolbox (gamultiobj)

clear; clc; close all;

%% Paths
script_dir = fileparts(mfilename('fullpath'));
project_root = fullfile(script_dir, '..', '..');
addpath(genpath(fullfile(project_root, 'src')));

fprintf('╔════════════════════════════════════════════════════════════════════╗\n');
fprintf('║   MULTI-OBJECTIVE: T_max vs THICKNESS (with N, M grid search)     ║\n');
fprintf('╚════════════════════════════════════════════════════════════════════╝\n\n');

%% Load centralized configuration
[var_names, lb, ub, x0, all_vars, CONFIG] = optimization_variables();
nvars = length(var_names);

%% Integer parameter grid
% Based on your sweep result showing N=3 is optimal, we'll focus around that
N_values = [2, 3, 4, 5, 6, 7, 8, 9, 10];           % Number of stages
M_values = [8, 10, 12, 15, 18];    % Number of wedges (θ = 360°/M)

fprintf('Integer Parameter Grid:\n');
fprintf('  N (stages): '); fprintf('%d ', N_values); fprintf('\n');
fprintf('  M (wedges): '); fprintf('%d ', M_values); fprintf('\n');
fprintf('  Total combinations: %d\n\n', length(N_values) * length(M_values));

%% Output directory
timestamp = datestr(now, 'yyyy-mm-dd_HH-MM-SS');
OUTPUT_DIR = fullfile(project_root, 'output', 'thickness_temp_optimization', timestamp);
if ~exist(OUTPUT_DIR, 'dir')
    mkdir(OUTPUT_DIR);
end

%% Display continuous variables
fprintf('Continuous Optimization Variables (%d):\n', nvars);
fprintf('═══════════════════════════════════════════════════════════════\n');
fprintf('  %-20s  %10s  %10s  %10s\n', 'Variable', 'Lower', 'Upper', 'Initial');
fprintf('───────────────────────────────────────────────────────────────\n');
for i = 1:nvars
    fprintf('  %-20s  %10.4f  %10.4f  %10.4f\n', var_names{i}, lb(i), ub(i), x0(i));
end
fprintf('═══════════════════════════════════════════════════════════════\n');
fprintf('\nBoundary Conditions:\n');
fprintf('  q_flux:   %.2e W/m² (%.0f kW/m²)\n', CONFIG.q_flux_W_m2, CONFIG.q_flux_W_m2/1e3);
fprintf('  T_water:  %.1f K (%.1f °C)\n', CONFIG.T_water_K, CONFIG.T_water_K - 273.15);
fprintf('  h_conv:   %.2e W/m²K\n', CONFIG.h_conv_W_m2K);
fprintf('\nObjectives:\n');
fprintf('  1. Minimize T_max (°C)\n');
fprintf('  2. Minimize t_TEC (µm)\n');
fprintf('───────────────────────────────────────────────────────────────\n\n');

%% GA options (reduced for speed since we have multiple (N,M) combinations)
ga_options = optimoptions('gamultiobj', ...
    'PopulationSize', 40, ...
    'MaxGenerations', 25, ...
    'ParetoFraction', 0.35, ...
    'Display', 'iter', ...      % Show progress
    'UseParallel', false);

% First, test the objective function works
fprintf('Testing objective function at initial point...\n');
try
    test_config = create_base_config(CONFIG, 3, 12);
    test_f = objective_temp_thickness(x0, test_config, CONFIG);
    fprintf('  Test result: T_max=%.1f °C, t_TEC=%.0f µm\n', test_f(1), test_f(2));
catch ME
    fprintf('  ERROR: %s\n', ME.message);
    return;
end
fprintf('\n');

%% Main grid search loop
all_results = struct();
all_results.N_values = N_values;
all_results.M_values = M_values;
all_results.grid = cell(length(N_values), length(M_values));

best_T_max = inf;
best_combo = [];

total_combos = length(N_values) * length(M_values);
combo_idx = 0;

for i_N = 1:length(N_values)
    for i_M = 1:length(M_values)
        combo_idx = combo_idx + 1;
        N = N_values(i_N);
        M = M_values(i_M);
        theta_deg = 360 / M;
        
        fprintf('═══════════════════════════════════════════════════════════════\n');
        fprintf('[%d/%d] N=%d, M=%d (θ=%.1f°)\n', combo_idx, total_combos, N, M, theta_deg);
        fprintf('═══════════════════════════════════════════════════════════════\n');
        
        % Create base config for this (N, M)
        base_config = create_base_config(CONFIG, N, M);
        
        % Define objective function
        objective_fn = @(x) objective_temp_thickness(x, base_config, CONFIG);
        
        % Run multi-objective GA
        try
            tic;
            [x_pareto, fval_pareto, exitflag, output] = ...
                gamultiobj(objective_fn, nvars, [], [], [], [], lb, ub, [], ga_options);
            opt_time = toc;
            
            % Extract results
            T_max_pareto = fval_pareto(:, 1);
            thickness_pareto = fval_pareto(:, 2);
            
            % Find key points
            [min_T, idx_min_T] = min(T_max_pareto);
            [min_t, idx_min_t] = min(thickness_pareto);
            
            % Knee point (balanced)
            if length(T_max_pareto) > 1
                utopia = [min(T_max_pareto), min(thickness_pareto)];
                norm_T = (T_max_pareto - utopia(1)) / max(1e-12, range(T_max_pareto));
                norm_t = (thickness_pareto - utopia(2)) / max(1e-12, range(thickness_pareto));
                distances = sqrt(norm_T.^2 + norm_t.^2);
                [~, idx_knee] = min(distances);
            else
                idx_knee = 1;
            end
            
            % Store results
            result = struct();
            result.N = N;
            result.M = M;
            result.theta_deg = theta_deg;
            result.x_pareto = x_pareto;
            result.fval_pareto = fval_pareto;
            result.T_max_pareto = T_max_pareto;
            result.thickness_pareto = thickness_pareto;
            result.min_T_max = min_T;
            result.min_thickness = min_t;
            result.x_min_T = x_pareto(idx_min_T, :);
            result.x_min_t = x_pareto(idx_min_t, :);
            result.x_knee = x_pareto(idx_knee, :);
            result.T_knee = T_max_pareto(idx_knee);
            result.t_knee = thickness_pareto(idx_knee);
            result.opt_time = opt_time;
            result.converged = true;
            
            all_results.grid{i_N, i_M} = result;
            
            fprintf('  Pareto points: %d\n', size(x_pareto, 1));
            fprintf('  Min T_max: %.1f °C (at t=%.0f µm)\n', min_T, thickness_pareto(idx_min_T));
            fprintf('  Min t_TEC: %.0f µm (at T=%.1f °C)\n', min_t, T_max_pareto(idx_min_t));
            fprintf('  Knee point: T=%.1f °C, t=%.0f µm\n', T_max_pareto(idx_knee), thickness_pareto(idx_knee));
            fprintf('  Time: %.1f s\n', opt_time);
            
            % Track global best
            if min_T < best_T_max
                best_T_max = min_T;
                best_combo = [N, M];
            end
            
        catch ME
            fprintf('  FAILED: %s\n', ME.message);
            result = struct();
            result.N = N;
            result.M = M;
            result.converged = false;
            result.error = ME.message;
            all_results.grid{i_N, i_M} = result;
        end
    end
end

%% Summary
fprintf('\n');
fprintf('╔════════════════════════════════════════════════════════════════════╗\n');
fprintf('║                         OPTIMIZATION SUMMARY                       ║\n');
fprintf('╚════════════════════════════════════════════════════════════════════╝\n\n');

fprintf('Best (N, M) combination for minimum T_max:\n');
fprintf('  N = %d, M = %d (θ = %.1f°)\n', best_combo(1), best_combo(2), 360/best_combo(2));
fprintf('  Minimum T_max achieved: %.1f °C\n\n', best_T_max);

% Create summary table
fprintf('Grid Summary (Min T_max for each N, M):\n');
fprintf('═══════════════════════════════════════════════════════════════\n');
fprintf('     M→  ');
for M = M_values
    fprintf('%8d ', M);
end
fprintf('\n');
fprintf('  N↓     ');
for M = M_values
    fprintf('%8.0f°', 360/M);
end
fprintf('\n');
fprintf('───────────────────────────────────────────────────────────────\n');

for i_N = 1:length(N_values)
    fprintf('  %d      ', N_values(i_N));
    for i_M = 1:length(M_values)
        result = all_results.grid{i_N, i_M};
        if result.converged
            fprintf('%8.1f ', result.min_T_max);
        else
            fprintf('%8s ', 'FAIL');
        end
    end
    fprintf('\n');
end
fprintf('═══════════════════════════════════════════════════════════════\n');

%% Plot results
% 1. Heatmap of minimum T_max across (N, M) grid
figure('Position', [100, 100, 1400, 500], 'Name', 'Thickness-Temperature Optimization');

subplot(1, 3, 1);
T_matrix = nan(length(N_values), length(M_values));
for i_N = 1:length(N_values)
    for i_M = 1:length(M_values)
        result = all_results.grid{i_N, i_M};
        if result.converged
            T_matrix(i_N, i_M) = result.min_T_max;
        end
    end
end
imagesc(M_values, N_values, T_matrix);
colorbar;
xlabel('M (Number of Wedges)');
ylabel('N (Number of Stages)');
title('Min T_{max} (°C) across (N, M) Grid');
set(gca, 'YDir', 'normal');
colormap(flipud(hot));

% Add text annotations
for i_N = 1:length(N_values)
    for i_M = 1:length(M_values)
        if ~isnan(T_matrix(i_N, i_M))
            text(M_values(i_M), N_values(i_N), sprintf('%.0f', T_matrix(i_N, i_M)), ...
                'HorizontalAlignment', 'center', 'Color', 'w', 'FontWeight', 'bold');
        end
    end
end

% 2. Pareto fronts for best (N, M)
subplot(1, 3, 2);
best_result = all_results.grid{N_values == best_combo(1), M_values == best_combo(2)};
if best_result.converged
    scatter(best_result.thickness_pareto, best_result.T_max_pareto, 50, 'filled');
    hold on;
    plot(best_result.t_knee, best_result.T_knee, 'rp', 'MarkerSize', 15, 'MarkerFaceColor', 'r');
    xlabel('TEC Thickness (µm)');
    ylabel('T_{max} (°C)');
    title(sprintf('Pareto Front (N=%d, M=%d)', best_combo(1), best_combo(2)));
    legend('Pareto Points', 'Knee Point', 'Location', 'best');
    grid on;
end

% 3. Compare Pareto fronts across different N values (for best M)
subplot(1, 3, 3);
best_M_idx = find(M_values == best_combo(2));
colors = lines(length(N_values));
hold on;
for i_N = 1:length(N_values)
    result = all_results.grid{i_N, best_M_idx};
    if result.converged
        scatter(result.thickness_pareto, result.T_max_pareto, 30, colors(i_N,:), 'filled', ...
            'DisplayName', sprintf('N=%d', N_values(i_N)));
    end
end
xlabel('TEC Thickness (µm)');
ylabel('T_{max} (°C)');
title(sprintf('Pareto Fronts Comparison (M=%d)', best_combo(2)));
legend('Location', 'best');
grid on;

saveas(gcf, fullfile(OUTPUT_DIR, 'optimization_summary.png'));
saveas(gcf, fullfile(OUTPUT_DIR, 'optimization_summary.fig'));

%% Extract and display best design parameters
fprintf('\n');
fprintf('╔════════════════════════════════════════════════════════════════════╗\n');
fprintf('║                    RECOMMENDED DESIGN PARAMETERS                   ║\n');
fprintf('╚════════════════════════════════════════════════════════════════════╝\n\n');

best_result = all_results.grid{N_values == best_combo(1), M_values == best_combo(2)};
if best_result.converged
    % Get knee point parameters
    x_best = best_result.x_knee;
    
    fprintf('Knee Point Solution (balanced T_max vs thickness):\n');
    fprintf('───────────────────────────────────────────────────────────────\n');
    fprintf('Integer Parameters:\n');
    fprintf('  N (stages):        %d\n', best_combo(1));
    fprintf('  M (wedges):        %d  →  θ = %.1f°\n', best_combo(2), 360/best_combo(2));
    fprintf('\nContinuous Parameters:\n');
    for i = 1:nvars
        fprintf('  %-18s = %.4f\n', var_names{i}, x_best(i));
    end
    fprintf('\nPerformance:\n');
    fprintf('  T_max:             %.1f °C\n', best_result.T_knee);
    fprintf('  t_TEC:             %.0f µm\n', best_result.t_knee);
    fprintf('───────────────────────────────────────────────────────────────\n');
    
    % Also show min-T solution
    fprintf('\nMinimum Temperature Solution:\n');
    fprintf('  T_max:             %.1f °C\n', best_result.min_T_max);
    fprintf('  t_TEC:             %.0f µm\n', best_result.thickness_pareto(best_result.T_max_pareto == best_result.min_T_max));
end

%% Save results
save(fullfile(OUTPUT_DIR, 'thickness_temp_results.mat'), 'all_results', 'best_combo', 'var_names');

% Export best design to CSV
if best_result.converged
    design_table = table();
    design_table.Parameter = ['N'; 'M'; 'theta_deg'; var_names; 'T_max_C'; 't_TEC_um'];
    design_table.Value = [best_combo(1); best_combo(2); 360/best_combo(2); ...
        x_best'; best_result.T_knee; best_result.t_knee];
    writetable(design_table, fullfile(OUTPUT_DIR, 'best_design_parameters.csv'));
end

fprintf('\nResults saved to: %s\n', OUTPUT_DIR);
fprintf('\n✓ Multi-objective optimization complete!\n');

%% ==================== HELPER FUNCTIONS ====================

function config = create_base_config(CONFIG, N, M)
    config = struct();
    
    % Integer parameters
    config.geometry.N_stages = N;
    config.geometry.M_wedges = M;
    config.geometry.wedge_angle_deg = 360 / M;
    
    % Fixed geometry
    config.geometry.w_chip_um = CONFIG.W_chip_um;
    config.geometry.t_chip_um = CONFIG.t_chip_um;
    
    % Defaults (will be overridden by optimization)
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

function f = objective_temp_thickness(x, base_config, CONFIG)
    % Returns [T_max (°C), t_TEC (µm)]
    
    persistent call_count
    if isempty(call_count)
        call_count = 0;
    end
    call_count = call_count + 1;
    
    try
        % Build config
        config = struct();
        config.geometry = base_config.geometry;
        config.boundary_conditions = base_config.boundary_conditions;
        config.operating_conditions = base_config.operating_conditions;
        config.materials = base_config.materials;
        
        % Map optimization variables
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
        
        % Apply continuous variables
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
        
        % Solve thermal network
        materials = MaterialProperties(config);
        geometry = TECGeometry(config);
        network = ThermalNetwork(geometry, materials, config);
        
        N = geometry.N_stages;
        T_water = config.boundary_conditions.T_water_K;
        T = ones(2*N + 1, 1) * (T_water + 50);
        
        for iter = 1:100
            T_old = T;
            T_new = network.solve(T);
            
            if any(isnan(T_new)) || any(isinf(T_new)) || any(T_new < 0)
                f = [1e6, 1e6];
                return;
            end
            
            T = 0.5 * T_new + 0.5 * T;  % Relaxation
            
            if max(abs(T - T_old)) < 1e-6
                break;
            end
        end
        
        T_max_C = max(T) - 273.15;
        t_TEC_um = config.geometry.thickness_um;
        
        f = [T_max_C, t_TEC_um];
        
    catch
        f = [1e6, 1e6];
    end
end

function val = get_var_or_default(var_map, name, default)
    if var_map.isKey(name)
        val = var_map(name);
    else
        val = default;
    end
end
