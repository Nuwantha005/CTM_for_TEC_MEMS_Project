% run_thickness_temp_optimization_parallel.m
% PARALLEL Multi-objective optimization: Minimize T_max AND TEC thickness (t_TEC)
%
% This version uses parfor for parallel grid search - MUCH faster!
%
% Strategy:
%   - Parallel loop over all (N, M) combinations
%   - Each worker runs gamultiobj independently
%
% Requirements:
%   - Parallel Computing Toolbox
%
% Uses: Global Optimization Toolbox (gamultiobj), Parallel Computing Toolbox

clear; clc; close all;

%% Paths
script_dir = fileparts(mfilename('fullpath'));
project_root = fullfile(script_dir, '..', '..');
addpath(genpath(fullfile(project_root, 'src')));

fprintf('╔════════════════════════════════════════════════════════════════════╗\n');
fprintf('║   PARALLEL MULTI-OBJECTIVE: T_max vs THICKNESS (N, M grid search) ║\n');
fprintf('╚════════════════════════════════════════════════════════════════════╝\n\n');

%% Load centralized configuration
[var_names, lb, ub, x0, all_vars, CONFIG] = optimization_variables();
nvars = length(var_names);

%% Integer parameter grid - FINER GRID
% Expanded ranges for more comprehensive search
N_values = 2:1:12;                          % Stages: 2 to 12
M_values = [4, 6, 8, 10, 12, 14, 16, 18, 20, 24];  % Wedges (more options)

fprintf('Integer Parameter Grid (FINE):\n');
fprintf('  N (stages): '); fprintf('%d ', N_values); fprintf('\n');
fprintf('  M (wedges): '); fprintf('%d ', M_values); fprintf('\n');
fprintf('  Total combinations: %d\n\n', length(N_values) * length(M_values));

%% Output directory
timestamp = datestr(now, 'yyyy-mm-dd_HH-MM-SS');
OUTPUT_DIR = fullfile(project_root, 'output', 'thickness_temp_optimization_parallel', timestamp);
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

%% Start parallel pool
fprintf('Setting up parallel pool...\n');
pool = gcp('nocreate');
if isempty(pool)
    pool = parpool('local');  % Start default pool
    fprintf('  Started new parallel pool with %d workers\n', pool.NumWorkers);
else
    fprintf('  Using existing pool with %d workers\n', pool.NumWorkers);
end
num_workers = pool.NumWorkers;

%% GA options - optimized for parallel execution
% Since we parallelize the outer loop, don't parallelize inside gamultiobj
ga_options = optimoptions('gamultiobj', ...
    'PopulationSize', 50, ...
    'MaxGenerations', 30, ...
    'ParetoFraction', 0.35, ...
    'Display', 'off', ...       % Quiet for parallel
    'UseParallel', false);      % Don't nest parallelism

%% Create all (N, M) combinations
combos = [];
for i_N = 1:length(N_values)
    for i_M = 1:length(M_values)
        combos = [combos; N_values(i_N), M_values(i_M), i_N, i_M];
    end
end
total_combos = size(combos, 1);
fprintf('\nTotal combinations to evaluate: %d\n', total_combos);
fprintf('Estimated time with %d workers: ~%.0f minutes\n', num_workers, ...
    total_combos * 45 / 60 / num_workers);  % ~45s per combo
fprintf('\n');

%% Pre-allocate results array (for parfor)
results_cell = cell(total_combos, 1);

%% Main PARALLEL grid search loop
fprintf('═══════════════════════════════════════════════════════════════\n');
fprintf('Starting parallel optimization...\n');
fprintf('═══════════════════════════════════════════════════════════════\n');

tic;
parfor combo_idx = 1:total_combos
    N = combos(combo_idx, 1);
    M = combos(combo_idx, 2);
    theta_deg = 360 / M;
    
    % Create base config for this (N, M)
    base_config = create_base_config_local(CONFIG, N, M);
    
    % Define objective function (must be local for parfor)
    objective_fn = @(x) objective_temp_thickness_local(x, base_config, CONFIG);
    
    % Run multi-objective GA
    try
        t_start = tic;
        [x_pareto, fval_pareto, ~, ~] = ...
            gamultiobj(objective_fn, nvars, [], [], [], [], lb, ub, [], ga_options);
        opt_time = toc(t_start);
        
        % Extract results
        T_max_pareto = fval_pareto(:, 1);
        thickness_pareto = fval_pareto(:, 2);
        
        % Find key points
        [min_T, idx_min_T] = min(T_max_pareto);
        [min_t, idx_min_t] = min(thickness_pareto);
        
        % Knee point (balanced)
        if length(T_max_pareto) > 1
            utopia = [min(T_max_pareto), min(thickness_pareto)];
            T_range = max(T_max_pareto) - min(T_max_pareto);
            t_range = max(thickness_pareto) - min(thickness_pareto);
            norm_T = (T_max_pareto - utopia(1)) / max(1e-12, T_range);
            norm_t = (thickness_pareto - utopia(2)) / max(1e-12, t_range);
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
        result.thickness_at_min_T = thickness_pareto(idx_min_T);
        result.T_at_min_thickness = T_max_pareto(idx_min_t);
        result.opt_time = opt_time;
        result.converged = true;
        
        results_cell{combo_idx} = result;
        
        fprintf('  [%3d/%d] N=%2d, M=%2d: T_min=%.1f°C (t=%.0fµm), Knee: T=%.1f°C, t=%.0fµm [%.1fs]\n', ...
            combo_idx, total_combos, N, M, min_T, thickness_pareto(idx_min_T), ...
            T_max_pareto(idx_knee), thickness_pareto(idx_knee), opt_time);
        
    catch ME
        result = struct();
        result.N = N;
        result.M = M;
        result.converged = false;
        result.error = ME.message;
        results_cell{combo_idx} = result;
        
        fprintf('  [%3d/%d] N=%2d, M=%2d: FAILED - %s\n', ...
            combo_idx, total_combos, N, M, ME.message);
    end
end
total_time = toc;

fprintf('\n═══════════════════════════════════════════════════════════════\n');
fprintf('Parallel optimization complete! Total time: %.1f minutes\n', total_time/60);
fprintf('═══════════════════════════════════════════════════════════════\n\n');

%% Reorganize results into grid format
all_results = struct();
all_results.N_values = N_values;
all_results.M_values = M_values;
all_results.grid = cell(length(N_values), length(M_values));

for combo_idx = 1:total_combos
    i_N = combos(combo_idx, 3);
    i_M = combos(combo_idx, 4);
    all_results.grid{i_N, i_M} = results_cell{combo_idx};
end

%% Find global best
best_T_max = inf;
best_combo = [];
best_knee_T = inf;
best_knee_combo = [];

for i_N = 1:length(N_values)
    for i_M = 1:length(M_values)
        result = all_results.grid{i_N, i_M};
        if ~isempty(result) && result.converged
            if result.min_T_max < best_T_max
                best_T_max = result.min_T_max;
                best_combo = [N_values(i_N), M_values(i_M)];
            end
            if result.T_knee < best_knee_T
                best_knee_T = result.T_knee;
                best_knee_combo = [N_values(i_N), M_values(i_M)];
            end
        end
    end
end

%% Summary
fprintf('╔════════════════════════════════════════════════════════════════════╗\n');
fprintf('║                         OPTIMIZATION SUMMARY                       ║\n');
fprintf('╚════════════════════════════════════════════════════════════════════╝\n\n');

fprintf('Best (N, M) for MINIMUM T_max:\n');
fprintf('  N = %d, M = %d (θ = %.1f°)\n', best_combo(1), best_combo(2), 360/best_combo(2));
fprintf('  Minimum T_max achieved: %.1f °C\n\n', best_T_max);

fprintf('Best (N, M) for KNEE POINT:\n');
fprintf('  N = %d, M = %d (θ = %.1f°)\n', best_knee_combo(1), best_knee_combo(2), 360/best_knee_combo(2));
fprintf('  Knee T_max: %.1f °C\n\n', best_knee_T);

% Create summary table
fprintf('Grid Summary (Min T_max for each N, M):\n');
fprintf('═══════════════════════════════════════════════════════════════════════════\n');
fprintf('     M→  ');
for M = M_values
    fprintf('%6d ', M);
end
fprintf('\n');
fprintf('  N↓     ');
for M = M_values
    fprintf('%6.0f°', 360/M);
end
fprintf('\n');
fprintf('───────────────────────────────────────────────────────────────────────────\n');

for i_N = 1:length(N_values)
    fprintf('  %2d     ', N_values(i_N));
    for i_M = 1:length(M_values)
        result = all_results.grid{i_N, i_M};
        if ~isempty(result) && result.converged
            fprintf('%6.0f ', result.min_T_max);
        else
            fprintf('%6s ', 'FAIL');
        end
    end
    fprintf('\n');
end
fprintf('═══════════════════════════════════════════════════════════════════════════\n');

%% Create second table for knee point temperatures
fprintf('\nGrid Summary (Knee Point T_max for each N, M):\n');
fprintf('═══════════════════════════════════════════════════════════════════════════\n');
fprintf('     M→  ');
for M = M_values
    fprintf('%6d ', M);
end
fprintf('\n');
fprintf('───────────────────────────────────────────────────────────────────────────\n');

for i_N = 1:length(N_values)
    fprintf('  %2d     ', N_values(i_N));
    for i_M = 1:length(M_values)
        result = all_results.grid{i_N, i_M};
        if ~isempty(result) && result.converged
            fprintf('%6.0f ', result.T_knee);
        else
            fprintf('%6s ', 'FAIL');
        end
    end
    fprintf('\n');
end
fprintf('═══════════════════════════════════════════════════════════════════════════\n');

%% Create third table for knee point thickness
fprintf('\nGrid Summary (Knee Point Thickness µm for each N, M):\n');
fprintf('═══════════════════════════════════════════════════════════════════════════\n');
fprintf('     M→  ');
for M = M_values
    fprintf('%6d ', M);
end
fprintf('\n');
fprintf('───────────────────────────────────────────────────────────────────────────\n');

for i_N = 1:length(N_values)
    fprintf('  %2d     ', N_values(i_N));
    for i_M = 1:length(M_values)
        result = all_results.grid{i_N, i_M};
        if ~isempty(result) && result.converged
            fprintf('%6.0f ', result.t_knee);
        else
            fprintf('%6s ', 'FAIL');
        end
    end
    fprintf('\n');
end
fprintf('═══════════════════════════════════════════════════════════════════════════\n');

%% Plot results
try
    % Use 'Visible', 'off' for batch mode compatibility
    fig = figure('Position', [50, 50, 1600, 900], 'Name', 'Parallel Thickness-Temperature Optimization', 'Visible', 'off');
    
    % 1. Heatmap of minimum T_max across (N, M) grid
    subplot(2, 3, 1);
    T_matrix = nan(length(N_values), length(M_values));
    for i_N = 1:length(N_values)
        for i_M = 1:length(M_values)
            result = all_results.grid{i_N, i_M};
            if ~isempty(result) && result.converged
                T_matrix(i_N, i_M) = result.min_T_max;
            end
        end
    end
    imagesc(1:length(M_values), 1:length(N_values), T_matrix);
    colorbar;
    xticks(1:length(M_values));
    xticklabels(arrayfun(@num2str, M_values, 'UniformOutput', false));
    yticks(1:length(N_values));
    yticklabels(arrayfun(@num2str, N_values, 'UniformOutput', false));
    xlabel('M (Number of Wedges)');
    ylabel('N (Number of Stages)');
    title('Min T_{max} (°C) - Lower is Better');
    set(gca, 'YDir', 'normal');
    colormap(gca, flipud(hot));

% 2. Heatmap of knee point T_max
subplot(2, 3, 2);
T_knee_matrix = nan(length(N_values), length(M_values));
for i_N = 1:length(N_values)
    for i_M = 1:length(M_values)
        result = all_results.grid{i_N, i_M};
        if ~isempty(result) && result.converged
            T_knee_matrix(i_N, i_M) = result.T_knee;
        end
    end
end
imagesc(1:length(M_values), 1:length(N_values), T_knee_matrix);
colorbar;
xticks(1:length(M_values));
xticklabels(arrayfun(@num2str, M_values, 'UniformOutput', false));
yticks(1:length(N_values));
yticklabels(arrayfun(@num2str, N_values, 'UniformOutput', false));
xlabel('M (Number of Wedges)');
ylabel('N (Number of Stages)');
title('Knee Point T_{max} (°C)');
set(gca, 'YDir', 'normal');
colormap(gca, flipud(hot));

% 3. Heatmap of knee point thickness
subplot(2, 3, 3);
t_knee_matrix = nan(length(N_values), length(M_values));
for i_N = 1:length(N_values)
    for i_M = 1:length(M_values)
        result = all_results.grid{i_N, i_M};
        if ~isempty(result) && result.converged
            t_knee_matrix(i_N, i_M) = result.t_knee;
        end
    end
end
imagesc(1:length(M_values), 1:length(N_values), t_knee_matrix);
colorbar;
xticks(1:length(M_values));
xticklabels(arrayfun(@num2str, M_values, 'UniformOutput', false));
yticks(1:length(N_values));
yticklabels(arrayfun(@num2str, N_values, 'UniformOutput', false));
xlabel('M (Number of Wedges)');
ylabel('N (Number of Stages)');
title('Knee Point Thickness (µm)');
set(gca, 'YDir', 'normal');
colormap(gca, flipud(parula));

% 4. Pareto fronts for best (N, M) - minimum T
subplot(2, 3, 4);
if ~isempty(best_combo)
    i_N_best = find(N_values == best_combo(1));
    i_M_best = find(M_values == best_combo(2));
    best_result = all_results.grid{i_N_best, i_M_best};
    if ~isempty(best_result) && best_result.converged
        scatter(best_result.thickness_pareto, best_result.T_max_pareto, 50, 'b', 'filled');
        hold on;
        plot(best_result.t_knee, best_result.T_knee, 'rp', 'MarkerSize', 15, 'MarkerFaceColor', 'r');
        plot(best_result.thickness_at_min_T, best_result.min_T_max, 'gd', 'MarkerSize', 12, 'MarkerFaceColor', 'g');
        xlabel('TEC Thickness (µm)');
        ylabel('T_{max} (°C)');
        title(sprintf('Pareto Front - Best Min T (N=%d, M=%d)', best_combo(1), best_combo(2)));
        legend('Pareto Points', 'Knee Point', 'Min T_{max}', 'Location', 'best');
        grid on;
    end
else
    i_N_best = 1;
    i_M_best = 1;
    best_result = struct('converged', false);
    title('No valid results');
end

% 5. Compare Pareto fronts across different N values (for best M)
subplot(2, 3, 5);
colors = lines(length(N_values));
hold on;
for i_N = 1:length(N_values)
    result = all_results.grid{i_N, i_M_best};
    if ~isempty(result) && result.converged
        scatter(result.thickness_pareto, result.T_max_pareto, 30, colors(i_N,:), 'filled', ...
            'DisplayName', sprintf('N=%d', N_values(i_N)));
    end
end
xlabel('TEC Thickness (µm)');
ylabel('T_{max} (°C)');
title(sprintf('Pareto Fronts Comparison (M=%d)', best_combo(2)));
legend('Location', 'best');
grid on;

% 6. Compare across different M values (for best N)
subplot(2, 3, 6);
colors = lines(length(M_values));
hold on;
for i_M = 1:length(M_values)
    result = all_results.grid{i_N_best, i_M};
    if ~isempty(result) && result.converged
        scatter(result.thickness_pareto, result.T_max_pareto, 30, colors(i_M,:), 'filled', ...
            'DisplayName', sprintf('M=%d', M_values(i_M)));
    end
end
xlabel('TEC Thickness (µm)');
ylabel('T_{max} (°C)');
title(sprintf('Pareto Fronts Comparison (N=%d)', best_combo(1)));
legend('Location', 'best');
grid on;

    saveas(fig, fullfile(OUTPUT_DIR, 'parallel_optimization_summary.png'));
    saveas(fig, fullfile(OUTPUT_DIR, 'parallel_optimization_summary.fig'));
    fprintf('  Figures saved successfully.\n');
    close(fig);
catch ME
    fprintf('  Warning: Could not create/save figures: %s\n', ME.message);
end

%% Extract and display best design parameters
fprintf('\n');
fprintf('╔════════════════════════════════════════════════════════════════════╗\n');
fprintf('║                    RECOMMENDED DESIGN PARAMETERS                   ║\n');
fprintf('╚════════════════════════════════════════════════════════════════════╝\n\n');

if ~isempty(best_combo) && ~isempty(best_result) && isfield(best_result, 'converged') && best_result.converged
    % Get knee point parameters
    x_best = best_result.x_knee;
    
    fprintf('KNEE POINT Solution (balanced T_max vs thickness):\n');
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
    x_min_T = best_result.x_min_T;
    fprintf('\nMINIMUM TEMPERATURE Solution:\n');
    fprintf('───────────────────────────────────────────────────────────────\n');
    fprintf('Integer Parameters:\n');
    fprintf('  N (stages):        %d\n', best_combo(1));
    fprintf('  M (wedges):        %d  →  θ = %.1f°\n', best_combo(2), 360/best_combo(2));
    fprintf('\nContinuous Parameters:\n');
    for i = 1:nvars
        fprintf('  %-18s = %.4f\n', var_names{i}, x_min_T(i));
    end
    fprintf('\nPerformance:\n');
    fprintf('  T_max:             %.1f °C\n', best_result.min_T_max);
    fprintf('  t_TEC:             %.0f µm\n', best_result.thickness_at_min_T);
    fprintf('───────────────────────────────────────────────────────────────\n');
else
    fprintf('  No valid optimization results found.\n');
end

%% Save results
save(fullfile(OUTPUT_DIR, 'parallel_thickness_temp_results.mat'), 'all_results', 'best_combo', 'best_knee_combo', 'var_names', 'combos', 'total_time');

% Export best design to CSV - Knee Point
if ~isempty(best_combo) && ~isempty(best_result) && isfield(best_result, 'converged') && best_result.converged
    design_table = table();
    design_table.Parameter = ['N'; 'M'; 'theta_deg'; var_names(:); 'T_max_C'; 't_TEC_um'];
    design_table.Knee_Value = [best_combo(1); best_combo(2); 360/best_combo(2); ...
        best_result.x_knee(:); best_result.T_knee; best_result.t_knee];
    design_table.MinT_Value = [best_combo(1); best_combo(2); 360/best_combo(2); ...
        best_result.x_min_T(:); best_result.min_T_max; best_result.thickness_at_min_T];
    writetable(design_table, fullfile(OUTPUT_DIR, 'best_design_parameters.csv'));
end

% Export full grid summary - pre-allocate arrays to avoid table warnings
gs_N = [];
gs_M = [];
gs_theta = [];
gs_min_T = [];
gs_min_t = [];
gs_T_knee = [];
gs_t_knee = [];
gs_time = [];

for i_N = 1:length(N_values)
    for i_M = 1:length(M_values)
        result = all_results.grid{i_N, i_M};
        if ~isempty(result) && result.converged
            gs_N(end+1) = result.N;
            gs_M(end+1) = result.M;
            gs_theta(end+1) = result.theta_deg;
            gs_min_T(end+1) = result.min_T_max;
            gs_min_t(end+1) = result.min_thickness;
            gs_T_knee(end+1) = result.T_knee;
            gs_t_knee(end+1) = result.t_knee;
            gs_time(end+1) = result.opt_time;
        end
    end
end

% Create table from arrays
grid_summary = table(gs_N(:), gs_M(:), gs_theta(:), gs_min_T(:), gs_min_t(:), ...
    gs_T_knee(:), gs_t_knee(:), gs_time(:), ...
    'VariableNames', {'N', 'M', 'theta_deg', 'min_T_max', 'min_thickness', ...
    'T_knee', 't_knee', 'opt_time_s'});
writetable(grid_summary, fullfile(OUTPUT_DIR, 'grid_summary.csv'));

fprintf('\nResults saved to: %s\n', OUTPUT_DIR);
fprintf('\n✓ Parallel multi-objective optimization complete!\n');
fprintf('  Total time: %.1f minutes (%.1f seconds)\n', total_time/60, total_time);
fprintf('  Speedup vs serial: ~%.1fx\n', total_combos / (total_time / 45));

%% ==================== LOCAL HELPER FUNCTIONS (for parfor) ====================
% Note: These must be defined as local functions, not nested, for parfor compatibility

function config = create_base_config_local(CONFIG, N, M)
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

function f = objective_temp_thickness_local(x, base_config, CONFIG)
    % Returns [T_max (°C), t_TEC (µm)]
    
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
        
        get_val = @(name, default) get_var_or_default_local(var_map, name, default);
        
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

function val = get_var_or_default_local(var_map, name, default)
    if var_map.isKey(name)
        val = var_map(name);
    else
        val = default;
    end
end
