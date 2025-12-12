% analyze_specific_nm_optimization.m
% Re-run or analyze multi-objective optimization for a SPECIFIC (N, M) pair
%
% This script:
%   1. Either re-runs gamultiobj OR uses existing results from grid search
%   2. Generates detailed Pareto front plots
%   3. Shows knee point, min-T, min-thickness solutions
%   4. Displays all dimensions for key solutions
%   5. Plots temperature distributions
%
% Use this when you want detailed analysis of a specific integer parameter pair
%
% ═══════════════════════════════════════════════════════════════════════════

clear; clc; close all;

%% ═══════════════════════════════════════════════════════════════════════════
%  USER CONFIGURATION - SET YOUR DESIRED (N, M) VALUES
%  ═══════════════════════════════════════════════════════════════════════════
SPECIFIC_N = 7;                              % Number of stages
SPECIFIC_M = 16;                             % Number of wedges

% MODE: Use existing grid search results or run new optimization
USE_EXISTING_RESULTS = false;                 % true = use existing Pareto data from grid search
                                             % false = run new gamultiobj optimization
EXISTING_RESULT_FOLDER = '2025-12-08_13-38-13';  % Grid search folder (only if USE_EXISTING_RESULTS=true)

% GA Options (only used if USE_EXISTING_RESULTS = false)
POPULATION_SIZE = 150;                       % Default: 150
MAX_GENERATIONS = 100;                       % Default: 100
PARETO_FRACTION = 0.35;                      % Fraction of population on Pareto front
% ═══════════════════════════════════════════════════════════════════════════

%% Setup paths
script_dir = fileparts(mfilename('fullpath'));
project_root = fullfile(script_dir, '..', '..');
addpath(genpath(fullfile(project_root, 'src')));

N = SPECIFIC_N;
M = SPECIFIC_M;
theta_deg = 360 / M;

fprintf('╔════════════════════════════════════════════════════════════════════╗\n');
fprintf('║    MULTI-OBJECTIVE OPTIMIZATION FOR SPECIFIC (N, M) PAIR          ║\n');
fprintf('╚════════════════════════════════════════════════════════════════════╝\n\n');

fprintf('Target Configuration:\n');
fprintf('  N (stages)      = %d\n', N);
fprintf('  M (wedges)      = %d\n', M);
fprintf('  θ (wedge angle) = %.2f°\n\n', theta_deg);

%% Load centralized configuration
[var_names, lb, ub, x0, all_vars, CONFIG] = optimization_variables();
nvars = length(var_names);

%% Output directory
timestamp = datestr(now, 'yyyy-mm-dd_HH-MM-SS');
OUTPUT_DIR = fullfile(project_root, 'output', 'specific_nm_optimization', ...
    sprintf('N%d_M%d_%s', N, M, timestamp));
if ~exist(OUTPUT_DIR, 'dir')
    mkdir(OUTPUT_DIR);
end

%% Display optimization setup
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

%% Create base configuration (same format as parallel optimization)
config = struct();

% Integer parameters
config.geometry.N_stages = N;
config.geometry.M_wedges = M;
config.geometry.wedge_angle_deg = theta_deg;

% Fixed geometry
config.geometry.w_chip_um = CONFIG.W_chip_um;
config.geometry.t_chip_um = CONFIG.t_chip_um;

% Defaults (will be overridden by optimization)
config.geometry.R_cyl_um = 1000;
config.geometry.thickness_um = 200;
config.geometry.t_ins_um = 10;
config.geometry.radial_expansion_factor = 1.15;
config.geometry.azimuthal_gap_um = 50;
config.geometry.insulation_width_um = 50;
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

%% Define objective function (for new optimization mode)
objective_fn = @(x) objective_temp_thickness_local(x, config, CONFIG);

%% Get Pareto Data (either from existing results or run new optimization)
if USE_EXISTING_RESULTS
    %% ═══════════════════════════════════════════════════════════════════════
    %  MODE: USE EXISTING GRID SEARCH RESULTS
    %  ═══════════════════════════════════════════════════════════════════════
    fprintf('Loading existing Pareto data from grid search...\n');
    result_file = fullfile(project_root, 'output', 'thickness_temp_optimization_parallel', ...
        EXISTING_RESULT_FOLDER, 'parallel_thickness_temp_results.mat');
    
    if ~exist(result_file, 'file')
        error('Result file not found: %s', result_file);
    end
    
    data = load(result_file);
    all_results = data.all_results;
    
    % Find N, M in results
    i_N = find(all_results.N_values == N);
    i_M = find(all_results.M_values == M);
    
    if isempty(i_N) || isempty(i_M)
        error('N=%d, M=%d not found in grid search results!', N, M);
    end
    
    existing_result = all_results.grid{i_N, i_M};
    if isempty(existing_result) || ~existing_result.converged
        error('No valid result for N=%d, M=%d in grid search!', N, M);
    end
    
    % Extract Pareto data
    x_pareto = existing_result.x_pareto;
    fval_pareto = existing_result.fval_pareto;
    T_max_pareto = existing_result.T_max_pareto;
    thickness_pareto = existing_result.thickness_pareto;
    
    % Key solutions from existing results
    x_min_T = existing_result.x_min_T;
    x_min_t = existing_result.x_min_t;
    x_knee = existing_result.x_knee;
    
    min_T = existing_result.min_T_max;
    t_at_min_T = existing_result.thickness_at_min_T;
    T_at_min_t = existing_result.T_at_min_thickness;
    min_t = existing_result.min_thickness;
    T_knee = existing_result.T_knee;
    t_knee = existing_result.t_knee;
    
    opt_time = existing_result.opt_time;
    
    fprintf('  Loaded %d Pareto-optimal solutions\n', size(x_pareto, 1));
    fprintf('  Original optimization time: %.1f seconds\n', opt_time);
    
    % Create output directory
    timestamp = datestr(now, 'yyyy-mm-dd_HH-MM-SS');
    OUTPUT_DIR = fullfile(project_root, 'output', 'specific_nm_optimization', ...
        sprintf('N%d_M%d_%s', N, M, timestamp));
    if ~exist(OUTPUT_DIR, 'dir')
        mkdir(OUTPUT_DIR);
    end
    
else
    %% ═══════════════════════════════════════════════════════════════════════
    %  MODE: RUN NEW OPTIMIZATION
    %  ═══════════════════════════════════════════════════════════════════════
    
    % Output directory
    timestamp = datestr(now, 'yyyy-mm-dd_HH-MM-SS');
    OUTPUT_DIR = fullfile(project_root, 'output', 'specific_nm_optimization', ...
        sprintf('N%d_M%d_%s', N, M, timestamp));
    if ~exist(OUTPUT_DIR, 'dir')
        mkdir(OUTPUT_DIR);
    end
    
    fprintf('GA Settings:\n');
    fprintf('  Population Size: %d\n', POPULATION_SIZE);
    fprintf('  Max Generations: %d\n', MAX_GENERATIONS);
    fprintf('  Pareto Fraction: %.2f\n\n', PARETO_FRACTION);
    
    ga_options = optimoptions('gamultiobj', ...
        'PopulationSize', POPULATION_SIZE, ...
        'MaxGenerations', MAX_GENERATIONS, ...
        'ParetoFraction', PARETO_FRACTION, ...
        'Display', 'iter', ...
        'PlotFcn', {@gaplotpareto}, ...
        'UseParallel', true);
    
    fprintf('═══════════════════════════════════════════════════════════════\n');
    fprintf('Starting multi-objective optimization...\n');
    fprintf('═══════════════════════════════════════════════════════════════\n\n');
    
    tic;
    [x_pareto, fval_pareto, ~, ~] = ...
        gamultiobj(objective_fn, nvars, [], [], [], [], lb, ub, [], ga_options);
    opt_time = toc;
    
    fprintf('\n═══════════════════════════════════════════════════════════════\n');
    fprintf('Optimization complete! Time: %.1f seconds (%.1f minutes)\n', opt_time, opt_time/60);
    fprintf('═══════════════════════════════════════════════════════════════\n\n');
    
    % Extract results
    T_max_pareto = fval_pareto(:, 1);
    thickness_pareto = fval_pareto(:, 2);
    
    % Find key points
    [min_T, idx_min_T] = min(T_max_pareto);
    x_min_T = x_pareto(idx_min_T, :);
    t_at_min_T = thickness_pareto(idx_min_T);
    
    [min_t, idx_min_t] = min(thickness_pareto);
    x_min_t = x_pareto(idx_min_t, :);
    T_at_min_t = T_max_pareto(idx_min_t);
    
    % Knee point
    utopia = [min(T_max_pareto), min(thickness_pareto)];
    T_range = max(T_max_pareto) - min(T_max_pareto);
    t_range = max(thickness_pareto) - min(thickness_pareto);
    if T_range > 0 && t_range > 0
        norm_T = (T_max_pareto - utopia(1)) / T_range;
        norm_t = (thickness_pareto - utopia(2)) / t_range;
        distances = sqrt(norm_T.^2 + norm_t.^2);
        [~, idx_knee] = min(distances);
    else
        idx_knee = 1;
        norm_T = zeros(size(T_max_pareto));
        norm_t = zeros(size(thickness_pareto));
        distances = zeros(size(T_max_pareto));
    end
    x_knee = x_pareto(idx_knee, :);
    T_knee = T_max_pareto(idx_knee);
    t_knee = thickness_pareto(idx_knee);
end

%% Statistics
num_pareto = size(x_pareto, 1);

% Compute normalized values and find indices for plotting
T_range = max(T_max_pareto) - min(T_max_pareto);
t_range = max(thickness_pareto) - min(thickness_pareto);
utopia = [min(T_max_pareto), min(thickness_pareto)];
if T_range > 0 && t_range > 0
    norm_T = (T_max_pareto - utopia(1)) / T_range;
    norm_t = (thickness_pareto - utopia(2)) / t_range;
    distances = sqrt(norm_T.^2 + norm_t.^2);
else
    norm_T = zeros(size(T_max_pareto));
    norm_t = zeros(size(thickness_pareto));
    distances = zeros(size(T_max_pareto));
end

% Find indices of key points in Pareto front
[~, idx_min_T] = min(T_max_pareto);
[~, idx_min_t] = min(thickness_pareto);
[~, idx_knee] = min(abs(T_max_pareto - T_knee) + abs(thickness_pareto - t_knee));

fprintf('\n═══════════════════════════════════════════════════════════════\n');
fprintf('Pareto Front Statistics:\n');
fprintf('  Number of Pareto-optimal solutions: %d\n', num_pareto);
fprintf('  T_max range:     [%.1f, %.1f] °C\n', min(T_max_pareto), max(T_max_pareto));
fprintf('  Thickness range: [%.0f, %.0f] µm\n\n', min(thickness_pareto), max(thickness_pareto));

%% Display key solutions
fprintf('╔════════════════════════════════════════════════════════════════════╗\n');
fprintf('║                         KEY PARETO SOLUTIONS                       ║\n');
fprintf('╚════════════════════════════════════════════════════════════════════╝\n\n');

fprintf('1. MINIMUM TEMPERATURE Solution:\n');
fprintf('───────────────────────────────────────────────────────────────\n');
fprintf('   T_max = %.2f °C, t_TEC = %.0f µm\n', min_T, t_at_min_T);
fprintf('   Parameters:\n');
for i = 1:nvars
    fprintf('     %-18s = %.4f\n', var_names{i}, x_min_T(i));
end

fprintf('\n2. MINIMUM THICKNESS Solution:\n');
fprintf('───────────────────────────────────────────────────────────────\n');
fprintf('   T_max = %.2f °C, t_TEC = %.0f µm\n', T_at_min_t, min_t);
fprintf('   Parameters:\n');
for i = 1:nvars
    fprintf('     %-18s = %.4f\n', var_names{i}, x_min_t(i));
end

fprintf('\n3. KNEE POINT Solution (Balanced):\n');
fprintf('───────────────────────────────────────────────────────────────\n');
fprintf('   T_max = %.2f °C, t_TEC = %.0f µm\n', T_knee, t_knee);
fprintf('   Parameters:\n');
for i = 1:nvars
    fprintf('     %-18s = %.4f\n', var_names{i}, x_knee(i));
end
fprintf('═══════════════════════════════════════════════════════════════\n\n');

%% Create comprehensive plots
fprintf('Generating plots...\n');

% Figure 1: Main Pareto Front Analysis
fig1 = figure('Position', [50, 50, 1400, 1000], 'Name', sprintf('Pareto Analysis N=%d M=%d', N, M));

% Subplot 1: Pareto Front with key points
subplot(2, 2, 1);
scatter(thickness_pareto, T_max_pareto, 60, 'b', 'filled', 'MarkerFaceAlpha', 0.6);
hold on;
plot(t_knee, T_knee, 'rp', 'MarkerSize', 20, 'MarkerFaceColor', 'r', 'LineWidth', 2);
plot(t_at_min_T, min_T, 'g^', 'MarkerSize', 15, 'MarkerFaceColor', 'g', 'LineWidth', 2);
plot(min_t, T_at_min_t, 'ms', 'MarkerSize', 15, 'MarkerFaceColor', 'm', 'LineWidth', 2);
xlabel('TEC Thickness (µm)', 'FontSize', 12);
ylabel('T_{max} (°C)', 'FontSize', 12);
title(sprintf('Pareto Front (N=%d, M=%d, θ=%.1f°)', N, M, theta_deg), 'FontSize', 14);
legend('Pareto Points', 'Knee Point', 'Min Temperature', 'Min Thickness', 'Location', 'best');
grid on;
set(gca, 'FontSize', 11);

% Subplot 2: Trade-off curve (sorted)
subplot(2, 2, 2);
[sorted_t, sort_idx] = sort(thickness_pareto);
sorted_T = T_max_pareto(sort_idx);
plot(sorted_t, sorted_T, 'b-', 'LineWidth', 2);
hold on;
scatter(sorted_t, sorted_T, 40, 'b', 'filled');
plot(t_knee, T_knee, 'rp', 'MarkerSize', 18, 'MarkerFaceColor', 'r');
xlabel('TEC Thickness (µm)', 'FontSize', 12);
ylabel('T_{max} (°C)', 'FontSize', 12);
title('Temperature vs Thickness Trade-off', 'FontSize', 14);
grid on;
set(gca, 'FontSize', 11);

% Add annotations for trade-off
if t_range > 0
    % Calculate trade-off slope at knee point
    delta_T = max(T_max_pareto) - min(T_max_pareto);
    delta_t = max(thickness_pareto) - min(thickness_pareto);
    text(t_knee + 10, T_knee + 2, sprintf('Knee: %.1f°C @ %.0fµm', T_knee, t_knee), ...
        'FontSize', 10, 'Color', 'r');
end

% Subplot 3: Normalized Pareto Front
subplot(2, 2, 3);
if T_range > 0 && t_range > 0
    scatter(norm_t, norm_T, 60, distances, 'filled');
    colorbar;
    hold on;
    plot(norm_t(idx_knee), norm_T(idx_knee), 'rp', 'MarkerSize', 18, 'MarkerFaceColor', 'r');
    plot([0 1], [0 1], 'k--', 'LineWidth', 1);  % Diagonal reference
    xlabel('Normalized Thickness', 'FontSize', 12);
    ylabel('Normalized T_{max}', 'FontSize', 12);
    title('Normalized Pareto Front (colored by distance to utopia)', 'FontSize', 14);
    colormap(gca, flipud(hot));
else
    text(0.5, 0.5, 'Insufficient range for normalization', 'HorizontalAlignment', 'center');
end
grid on;
set(gca, 'FontSize', 11);
axis equal;
xlim([0, 1.1]);
ylim([0, 1.1]);

% Subplot 4: Parameter sensitivity (key variables)
subplot(2, 2, 4);
% Plot how key parameters vary along Pareto front
key_vars = {'I', 't_TEC_um', 'r_cyl_um', 'fill_factor'};
key_indices = [];
for kv = key_vars
    idx = find(strcmp(var_names, kv{1}));
    if ~isempty(idx)
        key_indices(end+1) = idx;
    end
end

if ~isempty(key_indices) && T_range > 0
    colors = lines(length(key_indices));
    hold on;
    for ki = 1:length(key_indices)
        var_idx = key_indices(ki);
        % Normalize variable values
        var_vals = x_pareto(:, var_idx);
        var_norm = (var_vals - min(var_vals)) / max(1e-12, max(var_vals) - min(var_vals));
        scatter(thickness_pareto, var_norm, 30, colors(ki,:), 'filled', ...
            'DisplayName', var_names{var_idx}, 'MarkerFaceAlpha', 0.6);
    end
    xlabel('TEC Thickness (µm)', 'FontSize', 12);
    ylabel('Normalized Parameter Value', 'FontSize', 12);
    title('Parameter Variation Along Pareto Front', 'FontSize', 14);
    legend('Location', 'best');
    grid on;
end
set(gca, 'FontSize', 11);

sgtitle(sprintf('Multi-Objective Optimization Results: N=%d, M=%d', N, M), 'FontSize', 16, 'FontWeight', 'bold');

% Save figure 1
saveas(fig1, fullfile(OUTPUT_DIR, 'pareto_analysis.png'));
saveas(fig1, fullfile(OUTPUT_DIR, 'pareto_analysis.fig'));

%% Figure 2: Temperature Distribution for Key Solutions
fig2 = figure('Position', [100, 100, 1600, 600], 'Name', 'Temperature Distributions');

solutions = {x_min_T, x_knee, x_min_t};
solution_names = {'Min Temperature', 'Knee Point', 'Min Thickness'};
solution_colors = {'g', 'r', 'm'};
solution_T = [min_T, T_knee, T_at_min_t];
solution_t = [t_at_min_T, t_knee, min_t];

for sol_idx = 1:3
    subplot(1, 3, sol_idx);
    
    x = solutions{sol_idx};
    
    % Run solver and get temperature distribution
    [T_solution, converged] = run_solver_for_solution(x, config, CONFIG);
    
    if converged
        % Plot temperature
        stages = 0:N;  % 0 = chip, 1:N = stages
        T_plot = T_solution(1:N+1) - 273.15;  % Silicon layer temps in °C
        
        plot(stages, T_plot, '-o', 'Color', solution_colors{sol_idx}, ...
            'LineWidth', 2, 'MarkerSize', 10, 'MarkerFaceColor', solution_colors{sol_idx});
        hold on;
        
        % Add TEC layer temps
        T_TEC = T_solution(N+2:2*N+1) - 273.15;
        plot(1:N, T_TEC, '--s', 'Color', solution_colors{sol_idx}, ...
            'LineWidth', 1.5, 'MarkerSize', 8);
        
        xlabel('Stage', 'FontSize', 12);
        ylabel('Temperature (°C)', 'FontSize', 12);
        title(sprintf('%s\nT_{max}=%.1f°C, t=%.0fµm', solution_names{sol_idx}, ...
            solution_T(sol_idx), solution_t(sol_idx)), 'FontSize', 12);
        legend('Silicon Layer', 'TEC Layer', 'Location', 'best');
        grid on;
        xlim([-0.5, N+0.5]);
    else
        text(0.5, 0.5, 'Solver did not converge', 'HorizontalAlignment', 'center');
    end
    set(gca, 'FontSize', 11);
end

sgtitle(sprintf('Temperature Profiles for Key Solutions (N=%d, M=%d)', N, M), 'FontSize', 14, 'FontWeight', 'bold');

saveas(fig2, fullfile(OUTPUT_DIR, 'temperature_distributions.png'));
saveas(fig2, fullfile(OUTPUT_DIR, 'temperature_distributions.fig'));

%% Figure 3: Detailed Knee Point Analysis
fig3 = figure('Position', [150, 150, 1200, 800], 'Name', 'Knee Point Detailed Analysis');

% Run solver for knee point
[T_knee_sol, converged] = run_solver_for_solution(x_knee, config, CONFIG);

% Calculate radial positions
idx_t_TEC = find(strcmp(var_names, 't_TEC_um'));
idx_r_cyl = find(strcmp(var_names, 'r_cyl_um'));
idx_f_L = find(strcmp(var_names, 'f_L'));
t_TEC = x_knee(idx_t_TEC);
r_cyl = x_knee(idx_r_cyl);
f_L = x_knee(idx_f_L);

if converged
    % Subplot 1: Temperature vs radial distance
    subplot(2, 2, 1);
    r_positions = zeros(N+1, 1);
    r_positions(1) = 0;  % Chip center
    r_in = r_cyl;
    for stage = 1:N
        L_leg = t_TEC * f_L^(stage-1);
        r_positions(stage+1) = r_in;
        r_in = r_in + L_leg;
    end
    
    T_Si = T_knee_sol(1:N+1) - 273.15;
    plot(r_positions, T_Si, '-ob', 'LineWidth', 2, 'MarkerSize', 10, 'MarkerFaceColor', 'b');
    xlabel('Radial Position (µm)', 'FontSize', 12);
    ylabel('Temperature (°C)', 'FontSize', 12);
    title('Temperature vs Radial Distance', 'FontSize', 12);
    grid on;
    
    % Subplot 2: Stage-by-stage breakdown
    subplot(2, 2, 2);
    stages = 0:N;
    bar(stages, T_Si, 'FaceColor', [0.3, 0.6, 0.9]);
    xlabel('Stage (0 = Chip)', 'FontSize', 12);
    ylabel('Temperature (°C)', 'FontSize', 12);
    title('Temperature by Stage', 'FontSize', 12);
    grid on;
    
    % Subplot 3: Temperature drop per stage
    subplot(2, 2, 3);
    dT = -diff(T_Si);  % Temperature drop from center to periphery
    bar(1:N, dT, 'FaceColor', [0.9, 0.4, 0.3]);
    xlabel('Stage', 'FontSize', 12);
    ylabel('ΔT (°C)', 'FontSize', 12);
    title('Temperature Drop per Stage', 'FontSize', 12);
    grid on;
    
    % Subplot 4: Cumulative heat removal
    subplot(2, 2, 4);
    cumulative_dT = cumsum(dT);
    plot(1:N, cumulative_dT, '-sg', 'LineWidth', 2, 'MarkerSize', 10, 'MarkerFaceColor', 'g');
    xlabel('Stage', 'FontSize', 12);
    ylabel('Cumulative ΔT (°C)', 'FontSize', 12);
    title('Cumulative Temperature Reduction', 'FontSize', 12);
    grid on;
end

sgtitle(sprintf('Knee Point Detailed Analysis (T=%.1f°C, t=%.0fµm)', T_knee, t_knee), ...
    'FontSize', 14, 'FontWeight', 'bold');

saveas(fig3, fullfile(OUTPUT_DIR, 'knee_point_detailed.png'));
saveas(fig3, fullfile(OUTPUT_DIR, 'knee_point_detailed.fig'));

%% Print dimension table for knee point
fprintf('\n');
fprintf('╔════════════════════════════════════════════════════════════════════╗\n');
fprintf('║               KNEE POINT STAGE-BY-STAGE DIMENSIONS                 ║\n');
fprintf('╚════════════════════════════════════════════════════════════════════╝\n\n');

% Extract parameters
t_TEC = x_knee(find(strcmp(var_names, 't_TEC_um')));
r_cyl = x_knee(find(strcmp(var_names, 'r_cyl_um')));
f_L = x_knee(find(strcmp(var_names, 'f_L')));
fill_factor = x_knee(find(strcmp(var_names, 'fill_factor')));

fprintf('Main Geometry (all values in µm):\n');
fprintf('─────────────────────────────────────────────────────────────────────\n');
fprintf(' Stage   r_in      r_out     L_leg     W_az      W_leg\n');
fprintf('─────────────────────────────────────────────────────────────────────\n');

r_in = r_cyl;
for stage = 1:N
    L_leg = t_TEC * f_L^(stage-1);
    r_out = r_in + L_leg;
    r_mid = (r_in + r_out) / 2;
    W_leg = 2 * pi * r_mid / M;  % Full wedge width
    W_az = W_leg * fill_factor;  % Active material width
    
    fprintf('   %2d   %7.1f   %7.1f   %7.1f   %7.1f   %7.1f\n', ...
        stage, r_in, r_out, L_leg, W_az, W_leg);
    
    r_in = r_out;
end
fprintf('─────────────────────────────────────────────────────────────────────\n');
fprintf('\nTotal outer radius: %.1f µm (%.3f mm)\n', r_in, r_in/1000);

%% Save results
results = struct();
results.N = N;
results.M = M;
results.theta_deg = theta_deg;
results.x_pareto = x_pareto;
results.fval_pareto = fval_pareto;
results.T_max_pareto = T_max_pareto;
results.thickness_pareto = thickness_pareto;
results.num_pareto = num_pareto;

% Key solutions
results.knee.x = x_knee;
results.knee.T = T_knee;
results.knee.t = t_knee;
results.knee.idx = idx_knee;

results.min_T.x = x_min_T;
results.min_T.T = min_T;
results.min_T.t = t_at_min_T;
results.min_T.idx = idx_min_T;

results.min_t.x = x_min_t;
results.min_t.T = T_at_min_t;
results.min_t.t = min_t;
results.min_t.idx = idx_min_t;

results.opt_time = opt_time;
results.ga_options.PopulationSize = POPULATION_SIZE;
results.ga_options.MaxGenerations = MAX_GENERATIONS;
results.var_names = var_names;
results.CONFIG = CONFIG;

save(fullfile(OUTPUT_DIR, 'specific_nm_results.mat'), 'results');

% Export Pareto front to CSV (avoid duplicate column names)
pareto_var_names = [{'thickness_um', 'T_max_C'}, var_names(:)'];
pareto_table = array2table([thickness_pareto, T_max_pareto, x_pareto], ...
    'VariableNames', pareto_var_names);
writetable(pareto_table, fullfile(OUTPUT_DIR, 'pareto_front.csv'));

% Export key solutions
key_solutions = table();
key_solutions.Solution = {'MinTemperature'; 'KneePoint'; 'MinThickness'};
key_solutions.T_max_C = [min_T; T_knee; T_at_min_t];
key_solutions.t_TEC_um = [t_at_min_T; t_knee; min_t];
for i = 1:nvars
    key_solutions.(var_names{i}) = [x_min_T(i); x_knee(i); x_min_t(i)];
end
writetable(key_solutions, fullfile(OUTPUT_DIR, 'key_solutions.csv'));

fprintf('\n═══════════════════════════════════════════════════════════════\n');
fprintf('Results saved to: %s\n', OUTPUT_DIR);
fprintf('═══════════════════════════════════════════════════════════════\n');
fprintf('\nFiles generated:\n');
fprintf('  - pareto_analysis.png/fig     : Pareto front analysis\n');
fprintf('  - temperature_distributions   : Temperature profiles\n');
fprintf('  - knee_point_detailed        : Detailed knee analysis\n');
fprintf('  - specific_nm_results.mat    : All results\n');
fprintf('  - pareto_front.csv           : Pareto front data\n');
fprintf('  - key_solutions.csv          : Key solution parameters\n');
fprintf('\n✓ Analysis complete!\n');

%% ==================== HELPER FUNCTIONS ====================

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
        config.geometry.fill_factor = get_val('fill_factor', 0.7);
        config.geometry.insulation_width_um = get_val('W_is_um', 50);
        config.geometry.interconnect_ratio = get_val('f_ic_W', 0.15);
        config.geometry.interconnect_thickness_ratio = get_val('f_ic_t', 0.8);
        config.geometry.interconnect_angle_ratio = get_val('f_ic_beta', 0.16);
        config.geometry.outerconnect_ratio = get_val('f_oc_W', 0.15);
        config.geometry.outerconnect_thickness_ratio = get_val('f_oc_t', 0.8);
        config.geometry.outerconnect_angle_ratio = get_val('f_oc_beta', 0.16);
        
        % Manufacturability constraint check
        MIN_LEN = CONFIG.MIN_LENGTH_UM;
        N = config.geometry.N_stages;
        M = config.geometry.M_wedges;
        t_TEC = config.geometry.thickness_um;
        f_L = config.geometry.radial_expansion_factor;
        r_cyl = config.geometry.R_cyl_um;
        fill_factor = config.geometry.fill_factor;
        theta = 2 * pi / M;
        
        % Check azimuthal gap at innermost radius
        W_az_at_r_cyl = (1 - fill_factor) * r_cyl * theta;
        if W_az_at_r_cyl < MIN_LEN
            f = [1e6, 1e6];
            return;
        end
        
        % Check all stage leg lengths
        for i = 1:N
            L_i = t_TEC * f_L^(i-1);
            if L_i < MIN_LEN
                f = [1e6, 1e6];
                return;
            end
        end
        
        % Check interconnector widths
        f_ic_W = config.geometry.interconnect_ratio;
        for i = 1:N
            L_i = t_TEC * f_L^(i-1);
            W_ic = f_ic_W * L_i;
            if W_ic < MIN_LEN
                f = [1e6, 1e6];
                return;
            end
        end
        
        % Check outerconnector widths
        f_oc_W = config.geometry.outerconnect_ratio;
        for i = 1:N
            L_i = t_TEC * f_L^(i-1);
            W_oc = f_oc_W * L_i;
            if W_oc < MIN_LEN
                f = [1e6, 1e6];
                return;
            end
        end
        
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

function [T, converged] = run_solver_for_solution(x, base_config, CONFIG)
    % Run solver for a specific solution x and return temperature distribution
    
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
        config.geometry.fill_factor = get_val('fill_factor', 0.7);
        config.geometry.insulation_width_um = get_val('W_is_um', 50);
        config.geometry.interconnect_ratio = get_val('f_ic_W', 0.15);
        config.geometry.interconnect_thickness_ratio = get_val('f_ic_t', 0.8);
        config.geometry.interconnect_angle_ratio = get_val('f_ic_beta', 0.16);
        config.geometry.outerconnect_ratio = get_val('f_oc_W', 0.15);
        config.geometry.outerconnect_thickness_ratio = get_val('f_oc_t', 0.8);
        config.geometry.outerconnect_angle_ratio = get_val('f_oc_beta', 0.16);
        
        % Solve thermal network
        materials = MaterialProperties(config);
        geometry = TECGeometry(config);
        network = ThermalNetwork(geometry, materials, config);
        
        N = geometry.N_stages;
        T_water = config.boundary_conditions.T_water_K;
        T = ones(2*N + 1, 1) * (T_water + 50);
        
        converged = false;
        for iter = 1:100
            T_old = T;
            T_new = network.solve(T);
            
            if any(isnan(T_new)) || any(isinf(T_new)) || any(T_new < 0)
                break;
            end
            
            T = 0.5 * T_new + 0.5 * T;  % Relaxation
            
            if max(abs(T - T_old)) < 1e-6
                converged = true;
                break;
            end
        end
        
    catch
        T = [];
        converged = false;
    end
end
