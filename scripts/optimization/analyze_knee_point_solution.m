% analyze_knee_point_solution.m
% Post-optimization analysis script for knee point solution
%
% This script:
%   1. Loads optimization results from specified folder
%   2. Extracts either global best OR specific (N, M) knee point solution
%   3. Calculates and displays ALL stage dimensions in a detailed table
%   4. Runs the solver and shows temperature distribution
%
% Run AFTER: run_thickness_temp_optimization_parallel.m
%
% ═══════════════════════════════════════════════════════════════════════════
% HOW TO USE:
%   Option 1: Analyze LATEST results with GLOBAL BEST knee point
%       - Set USE_LATEST = true
%       - Set USE_GLOBAL_BEST = true
%
%   Option 2: Analyze SPECIFIC folder with GLOBAL BEST knee point
%       - Set USE_LATEST = false
%       - Set RESULT_FOLDER = 'yyyy-mm-dd_HH-MM-SS'
%       - Set USE_GLOBAL_BEST = true
%
%   Option 3: Analyze SPECIFIC (N, M) pair from a folder
%       - Set USE_LATEST = false (or true for latest folder)
%       - Set RESULT_FOLDER = 'yyyy-mm-dd_HH-MM-SS'
%       - Set USE_GLOBAL_BEST = false
%       - Set SPECIFIC_N = desired N (number of stages)
%       - Set SPECIFIC_M = desired M (number of wedges)
% ═══════════════════════════════════════════════════════════════════════════

clear; clc; close all;

%% ═══════════════════════════════════════════════════════════════════════════
%  USER CONFIGURATION - CHANGE THESE TO ANALYZE DIFFERENT RESULTS
%  ═══════════════════════════════════════════════════════════════════════════
USE_LATEST = true;                          % true = use most recent results
RESULT_FOLDER = '2025-12-08_11-41-15';      % Only used if USE_LATEST = false

% (N, M) Selection Options
USE_GLOBAL_BEST = false;                     % true = find best across all (N,M)
                                             % false = use specific (N,M) pair below
SPECIFIC_N = 7;                              % Number of stages (only if USE_GLOBAL_BEST = false)
SPECIFIC_M = 16;                             % Number of wedges (only if USE_GLOBAL_BEST = false)
% ═══════════════════════════════════════════════════════════════════════════

%% Setup paths
script_dir = fileparts(mfilename('fullpath'));
project_root = fullfile(script_dir, '..', '..');
addpath(genpath(fullfile(project_root, 'src')));

fprintf('╔════════════════════════════════════════════════════════════════════╗\n');
fprintf('║          KNEE POINT SOLUTION - DETAILED DIMENSION ANALYSIS         ║\n');
fprintf('╚════════════════════════════════════════════════════════════════════╝\n\n');

%% Find optimization results
results_base = fullfile(project_root, 'output', 'thickness_temp_optimization_parallel');
if ~exist(results_base, 'dir')
    error('No optimization results found! Run run_thickness_temp_optimization_parallel.m first.');
end

if USE_LATEST
    % Get latest timestamp folder
    folders = dir(results_base);
    folders = folders([folders.isdir] & ~ismember({folders.name}, {'.', '..'}));
    if isempty(folders)
        error('No optimization result folders found!');
    end
    
    % Sort by date
    [~, idx] = sort([folders.datenum], 'descend');
    selected_folder = folders(idx(1)).name;
    fprintf('Using LATEST results folder\n');
else
    % Use user-specified folder
    selected_folder = RESULT_FOLDER;
    if ~exist(fullfile(results_base, selected_folder), 'dir')
        % List available folders
        folders = dir(results_base);
        folders = folders([folders.isdir] & ~ismember({folders.name}, {'.', '..'}));
        fprintf('\nERROR: Folder "%s" not found!\n', selected_folder);
        fprintf('\nAvailable folders:\n');
        for i = 1:length(folders)
            fprintf('  - %s\n', folders(i).name);
        end
        error('Set RESULT_FOLDER to one of the folders listed above.');
    end
    fprintf('Using USER-SPECIFIED results folder\n');
end

results_path = fullfile(results_base, selected_folder);
fprintf('Loading results from: %s\n', selected_folder);

%% Load results
% Try different possible filenames
possible_files = {
    'optimization_results.mat',
    'parallel_thickness_temp_results.mat',
    'results.mat'
};

results_file = '';
for i = 1:length(possible_files)
    test_file = fullfile(results_path, possible_files{i});
    if exist(test_file, 'file')
        results_file = test_file;
        break;
    end
end

if isempty(results_file)
    % List available .mat files
    mat_files = dir(fullfile(results_path, '*.mat'));
    if ~isempty(mat_files)
        results_file = fullfile(results_path, mat_files(1).name);
    else
        error('No .mat results file found in: %s', results_path);
    end
end

fprintf('Loading: %s\n', results_file);
data = load(results_file);

% Handle different variable names
if isfield(data, 'all_results')
    all_results = data.all_results;
elseif isfield(data, 'results')
    all_results = data.results;
else
    error('Could not find results data in file');
end

if isfield(data, 'var_names')
    var_names = data.var_names;
else
    % Load from optimization_variables
    [var_names, ~, ~, ~, ~, ~] = optimization_variables();
end

if isfield(data, 'CONFIG')
    CONFIG = data.CONFIG;
else
    [~, ~, ~, ~, ~, CONFIG] = optimization_variables();
end

fprintf('Results loaded successfully!\n\n');

%% Find knee point solution (global best or specific N,M)
N_values = all_results.N_values;
M_values = all_results.M_values;

if USE_GLOBAL_BEST
    % Find global best knee point across all (N, M) combinations
    fprintf('Finding GLOBAL BEST knee point...\n');
    best_knee_T = inf;
    best_result = [];
    best_N = 0;
    best_M = 0;

    for i_N = 1:length(N_values)
        for i_M = 1:length(M_values)
            result = all_results.grid{i_N, i_M};
            if ~isempty(result) && isfield(result, 'converged') && result.converged
                if result.T_knee < best_knee_T
                    best_knee_T = result.T_knee;
                    best_result = result;
                    best_N = N_values(i_N);
                    best_M = M_values(i_M);
                end
            end
        end
    end
else
    % Use specific (N, M) pair
    fprintf('Looking for specific (N=%d, M=%d) pair...\n', SPECIFIC_N, SPECIFIC_M);
    
    % Find indices for specific N and M
    i_N = find(N_values == SPECIFIC_N);
    i_M = find(M_values == SPECIFIC_M);
    
    if isempty(i_N)
        fprintf('\nAvailable N values: %s\n', mat2str(N_values));
        error('N=%d not found in optimization results!', SPECIFIC_N);
    end
    if isempty(i_M)
        fprintf('\nAvailable M values: %s\n', mat2str(M_values));
        error('M=%d not found in optimization results!', SPECIFIC_M);
    end
    
    best_result = all_results.grid{i_N, i_M};
    if isempty(best_result) || ~isfield(best_result, 'converged') || ~best_result.converged
        error('No valid (converged) result for N=%d, M=%d!', SPECIFIC_N, SPECIFIC_M);
    end
    
    best_N = SPECIFIC_N;
    best_M = SPECIFIC_M;
    best_knee_T = best_result.T_knee;
end

if isempty(best_result)
    error('No valid results found!');
end

%% Extract knee point parameters
x_knee = best_result.x_knee;
N = best_N;
M = best_M;
theta_deg = 360 / M;
theta_rad = deg2rad(theta_deg);

fprintf('═══════════════════════════════════════════════════════════════════════\n');
fprintf('                    GLOBAL BEST KNEE POINT SOLUTION\n');
fprintf('═══════════════════════════════════════════════════════════════════════\n\n');
fprintf('Integer Parameters:\n');
fprintf('  N (stages)     = %d\n', N);
fprintf('  M (wedges)     = %d\n', M);
fprintf('  θ (wedge angle)= %.2f°\n', theta_deg);
fprintf('\nObjective Values:\n');
fprintf('  T_max          = %.2f °C\n', best_result.T_knee);
fprintf('  t_TEC          = %.1f µm\n', best_result.t_knee);
fprintf('\n');

%% Map optimization variables
var_map = containers.Map();
for i = 1:length(var_names)
    var_map(var_names{i}) = x_knee(i);
end

get_val = @(name, default) get_var_or_default(var_map, name, default);

% Extract key parameters
I = get_val('I', 0.025);
t_TEC = get_val('t_TEC_um', 200);
r_cyl = get_val('r_cyl_um', 1000);
t_ins = get_val('t_ins_um', 10);
f_L = get_val('f_L', 1.15);
fill_factor = get_val('fill_factor', 0.7);
W_is = get_val('W_is_um', 50);
f_ic_W = get_val('f_ic_W', 0.15);
f_ic_t = get_val('f_ic_t', 0.8);
f_ic_beta = get_val('f_ic_beta', 0.16);
f_oc_W = get_val('f_oc_W', 0.15);
f_oc_t = get_val('f_oc_t', 0.8);
f_oc_beta = get_val('f_oc_beta', 0.16);

%% Display Continuous Variables Table
fprintf('═══════════════════════════════════════════════════════════════════════\n');
fprintf('                       OPTIMIZED PARAMETERS\n');
fprintf('═══════════════════════════════════════════════════════════════════════\n');
fprintf('  %-20s  %12s  %s\n', 'Parameter', 'Value', 'Description');
fprintf('───────────────────────────────────────────────────────────────────────\n');
fprintf('  %-20s  %12.4f  %s\n', 'I [A]', I, 'Operating current');
fprintf('  %-20s  %12.1f  %s\n', 't_TEC [µm]', t_TEC, 'TEC layer thickness');
fprintf('  %-20s  %12.1f  %s\n', 'r_cyl [µm]', r_cyl, 'Central cylinder radius');
fprintf('  %-20s  %12.1f  %s\n', 't_ins [µm]', t_ins, 'Insulator thickness');
fprintf('  %-20s  %12.4f  %s\n', 'f_L', f_L, 'Radial expansion ratio');
fprintf('  %-20s  %12.4f  %s\n', 'fill_factor', fill_factor, 'Azimuthal fill factor');
fprintf('  %-20s  %12.1f  %s\n', 'W_is [µm]', W_is, 'Interstage insulator width');
fprintf('  %-20s  %12.4f  %s\n', 'f_ic_W', f_ic_W, 'IC width ratio');
fprintf('  %-20s  %12.4f  %s\n', 'f_ic_t', f_ic_t, 'IC thickness ratio');
fprintf('  %-20s  %12.4f  %s\n', 'f_ic_β', f_ic_beta, 'IC angle ratio');
fprintf('  %-20s  %12.4f  %s\n', 'f_oc_W', f_oc_W, 'OC width ratio');
fprintf('  %-20s  %12.4f  %s\n', 'f_oc_t', f_oc_t, 'OC thickness ratio');
fprintf('  %-20s  %12.4f  %s\n', 'f_oc_β', f_oc_beta, 'OC angle ratio');
fprintf('═══════════════════════════════════════════════════════════════════════\n\n');

%% Calculate Stage-by-Stage Dimensions
fprintf('═══════════════════════════════════════════════════════════════════════\n');
fprintf('                    STAGE-BY-STAGE DIMENSIONS\n');
fprintf('═══════════════════════════════════════════════════════════════════════\n\n');

% Pre-compute cumulative radii
r_in = zeros(N, 1);   % Inner radius of each stage
r_out = zeros(N, 1);  % Outer radius of each stage
L = zeros(N, 1);      % Leg length of each stage

r_in(1) = r_cyl;
for i = 1:N
    L(i) = t_TEC * f_L^(i-1);
    r_out(i) = r_in(i) + L(i);
    if i < N
        r_in(i+1) = r_out(i) + W_is;  % Add interstage insulator
    end
end

R_outer = r_out(N);  % Total outer radius

% Calculate derived dimensions for each stage
stage_data = struct();
for i = 1:N
    stage_data(i).stage = i;
    stage_data(i).r_in_um = r_in(i);
    stage_data(i).r_out_um = r_out(i);
    stage_data(i).L_um = L(i);
    
    % Azimuthal gap at mid-radius (derived from fill_factor)
    r_mid = (r_in(i) + r_out(i)) / 2;
    arc_length = r_mid * theta_rad;  % Arc length at mid-radius
    stage_data(i).W_az_um = (1 - fill_factor) * arc_length;
    stage_data(i).W_leg_um = fill_factor * arc_length;
    
    % Interconnector dimensions
    stage_data(i).W_ic_um = f_ic_W * L(i);
    stage_data(i).t_ic_um = f_ic_t * t_TEC;
    stage_data(i).beta_ic_deg = f_ic_beta * theta_deg;
    
    % Outerconnector dimensions
    stage_data(i).W_oc_um = f_oc_W * L(i);
    stage_data(i).t_oc_um = f_oc_t * t_TEC;
    stage_data(i).beta_oc_deg = f_oc_beta * theta_deg;
end

% Display main geometry table
fprintf('Main Geometry (all values in µm unless noted):\n');
fprintf('─────────────────────────────────────────────────────────────────────\n');
fprintf(' Stage   r_in      r_out     L_leg     W_az      W_leg     A_leg\n');
fprintf('─────────────────────────────────────────────────────────────────────\n');
for i = 1:N
    % Approximate leg area (trapezoidal)
    A_leg = 0.5 * (stage_data(i).r_in_um + stage_data(i).r_out_um) * theta_rad * ...
            fill_factor * t_TEC;  % in µm²
    fprintf('   %2d   %7.1f   %7.1f   %7.1f   %7.1f   %7.1f   %.2e\n', ...
        i, stage_data(i).r_in_um, stage_data(i).r_out_um, stage_data(i).L_um, ...
        stage_data(i).W_az_um, stage_data(i).W_leg_um, A_leg);
end
fprintf('─────────────────────────────────────────────────────────────────────\n\n');

% Display interconnector table
fprintf('Interconnector (IC) Dimensions - Cold Side Electrical Contact:\n');
fprintf('─────────────────────────────────────────────────────────────────────\n');
fprintf(' Stage   W_ic [µm]   t_ic [µm]   β_ic [°]   At radius [µm]\n');
fprintf('─────────────────────────────────────────────────────────────────────\n');
for i = 1:N
    fprintf('   %2d     %7.1f     %7.1f     %7.2f      %7.1f (inner)\n', ...
        i, stage_data(i).W_ic_um, stage_data(i).t_ic_um, ...
        stage_data(i).beta_ic_deg, stage_data(i).r_in_um);
end
fprintf('─────────────────────────────────────────────────────────────────────\n\n');

% Display outerconnector table
fprintf('Outerconnector (OC) Dimensions - Hot Side Electrical Contact:\n');
fprintf('─────────────────────────────────────────────────────────────────────\n');
fprintf(' Stage   W_oc [µm]   t_oc [µm]   β_oc [°]   At radius [µm]\n');
fprintf('─────────────────────────────────────────────────────────────────────\n');
for i = 1:N
    fprintf('   %2d     %7.1f     %7.1f     %7.2f      %7.1f (outer)\n', ...
        i, stage_data(i).W_oc_um, stage_data(i).t_oc_um, ...
        stage_data(i).beta_oc_deg, stage_data(i).r_out_um);
end
fprintf('─────────────────────────────────────────────────────────────────────\n\n');

%% Summary dimensions
fprintf('═══════════════════════════════════════════════════════════════════════\n');
fprintf('                       OVERALL DIMENSIONS\n');
fprintf('═══════════════════════════════════════════════════════════════════════\n');
fprintf('  Central cylinder radius (r_cyl)    : %8.1f µm = %.3f mm\n', r_cyl, r_cyl/1000);
fprintf('  Total outer radius (R_outer)       : %8.1f µm = %.3f mm\n', R_outer, R_outer/1000);
fprintf('  TEC annulus width (R_outer - r_cyl): %8.1f µm = %.3f mm\n', R_outer - r_cyl, (R_outer-r_cyl)/1000);
fprintf('  TEC layer thickness (t_TEC)        : %8.1f µm\n', t_TEC);
fprintf('  Interstage insulator width (W_is)  : %8.1f µm\n', W_is);
fprintf('  Total stages (N)                   : %8d\n', N);
fprintf('  Wedges per full circle (M)         : %8d\n', M);
fprintf('  Total TE legs (N × M)              : %8d\n', N * M);
fprintf('═══════════════════════════════════════════════════════════════════════\n\n');

%% Manufacturability Check
fprintf('═══════════════════════════════════════════════════════════════════════\n');
fprintf('                    MANUFACTURABILITY CHECK\n');
fprintf('═══════════════════════════════════════════════════════════════════════\n');
MIN_LEN = CONFIG.MIN_LENGTH_UM;
fprintf('  Minimum feature size constraint: %.0f µm\n\n', MIN_LEN);

issues = {};
for i = 1:N
    if stage_data(i).L_um < MIN_LEN
        issues{end+1} = sprintf('Stage %d: L_leg = %.1f µm < %.0f µm', i, stage_data(i).L_um, MIN_LEN);
    end
    if stage_data(i).W_az_um < MIN_LEN
        issues{end+1} = sprintf('Stage %d: W_az = %.1f µm < %.0f µm', i, stage_data(i).W_az_um, MIN_LEN);
    end
    if stage_data(i).W_ic_um < MIN_LEN
        issues{end+1} = sprintf('Stage %d: W_ic = %.1f µm < %.0f µm', i, stage_data(i).W_ic_um, MIN_LEN);
    end
    if stage_data(i).W_oc_um < MIN_LEN
        issues{end+1} = sprintf('Stage %d: W_oc = %.1f µm < %.0f µm', i, stage_data(i).W_oc_um, MIN_LEN);
    end
end

if isempty(issues)
    fprintf('  ✓ All dimensions meet minimum feature size constraint!\n');
else
    fprintf('  ✗ MANUFACTURABILITY ISSUES FOUND:\n');
    for i = 1:length(issues)
        fprintf('    - %s\n', issues{i});
    end
end
fprintf('═══════════════════════════════════════════════════════════════════════\n\n');

%% Run Solver for Temperature Distribution
fprintf('═══════════════════════════════════════════════════════════════════════\n');
fprintf('                    SOLVING THERMAL NETWORK\n');
fprintf('═══════════════════════════════════════════════════════════════════════\n\n');

% Build config for solver
config = struct();
config.geometry.N_stages = N;
config.geometry.M_wedges = M;
config.geometry.wedge_angle_deg = theta_deg;
config.geometry.thickness_um = t_TEC;
config.geometry.R_cyl_um = r_cyl;
config.geometry.t_ins_um = t_ins;
config.geometry.radial_expansion_factor = f_L;
config.geometry.fill_factor = fill_factor;
config.geometry.insulation_width_um = W_is;
config.geometry.interconnect_ratio = f_ic_W;
config.geometry.interconnect_thickness_ratio = f_ic_t;
config.geometry.interconnect_angle_ratio = f_ic_beta;
config.geometry.outerconnect_ratio = f_oc_W;
config.geometry.outerconnect_thickness_ratio = f_oc_t;
config.geometry.outerconnect_angle_ratio = f_oc_beta;
config.geometry.w_chip_um = CONFIG.W_chip_um;
config.geometry.t_chip_um = CONFIG.t_chip_um;

config.boundary_conditions.q_flux_W_m2 = CONFIG.q_flux_W_m2;
config.boundary_conditions.T_water_K = CONFIG.T_water_K;
config.boundary_conditions.h_conv_W_m2K = CONFIG.h_conv_W_m2K;

config.operating_conditions.I_current_A = I;
config.materials = CONFIG.materials;

% Initialize solver components
materials = MaterialProperties(config);
geometry = TECGeometry(config);
network = ThermalNetwork(geometry, materials, config);

% Initial guess
T_water = config.boundary_conditions.T_water_K;
T = ones(2*N + 1, 1) * (T_water + 50);

% Solve iteratively
fprintf('Running solver (relaxation method)...\n');
for iter = 1:200
    T_old = T;
    T_new = network.solve(T);
    T = 0.5 * T_new + 0.5 * T;  % Relaxation
    
    residual = max(abs(T - T_old));
    if residual < 1e-8
        fprintf('  Converged at iteration %d (residual = %.2e)\n', iter, residual);
        break;
    end
end

T_C = T - 273.15;  % Convert to Celsius
T_max_C = max(T_C);

fprintf('\nSolution Summary:\n');
fprintf('  T_max (chip center) = %.2f °C\n', T_C(1));
fprintf('  T_min (water side)  = %.2f °C\n', min(T_C));
fprintf('  ΔT across TEC       = %.2f °C\n', max(T_C) - min(T_C));

%% Display Temperature Distribution Table
fprintf('\n═══════════════════════════════════════════════════════════════════════\n');
fprintf('                    TEMPERATURE DISTRIBUTION\n');
fprintf('═══════════════════════════════════════════════════════════════════════\n');
fprintf('  Node    T [°C]\n');
fprintf('───────────────────────────────────────────────────────────────────────\n');
for i = 1:length(T_C)
    fprintf('  %3d     %7.2f\n', i, T_C(i));
end
fprintf('───────────────────────────────────────────────────────────────────────\n');

%% Plot Temperature Distribution
fig = figure('Position', [100, 100, 1400, 500], 'Name', 'Knee Point Solution Analysis');

% Get water temperature
T_water_C = config.boundary_conditions.T_water_K - 273.15;

% Extract in ResultsManager format
T_0 = T_C(1);                    % Chip center
T_Si = T_C(2:N+1);               % Silicon layer temps (nodes 2 to N+1)
T_c = T_C(N+2:2*N+1);            % TEC layer temps (nodes N+2 to 2N+1)

% ─────────────────────────────────────────────────────────────────────────────
% Plot 1: Radial Distance vs Temperature
% Silicon layer: r=0 (center), then r_in for each stage
% TEC layer: r_in for each stage, then r_out of last stage (water)
% ─────────────────────────────────────────────────────────────────────────────
subplot(1, 2, 1);

% Silicon/Chip layer: center (r=0) + r_in for each stage
chip_radii = [0];
chip_temps = [T_0];
for i = 1:N
    chip_radii = [chip_radii, stage_data(i).r_in_um];
    chip_temps = [chip_temps, T_Si(i)];
end

% TEC layer: r_in for each stage + r_out of last stage (water)
tec_radii = [];
tec_temps = [];
for i = 1:N
    tec_radii = [tec_radii, stage_data(i).r_in_um];
    tec_temps = [tec_temps, T_c(i)];
end
R_water = stage_data(N).r_out_um;
tec_radii = [tec_radii, R_water];
tec_temps = [tec_temps, T_water_C];

% Plot Silicon layer (blue circles)
plot(chip_radii/1000, chip_temps, '-ob', 'LineWidth', 1.8, 'MarkerFaceColor', 'b', 'DisplayName', 'Silicon Layer');
hold on;

% Plot TEC layer (red squares)
plot(tec_radii/1000, tec_temps, '-sr', 'LineWidth', 1.8, 'MarkerFaceColor', 'r', 'DisplayName', 'TEC Layer');

% Mark the water temperature point (green triangle)
plot(R_water/1000, T_water_C, 'g^', 'MarkerSize', 10, 'MarkerFaceColor', 'g', 'DisplayName', 'Coolant (T_{water})');

xlabel('Radial Distance (mm)');
ylabel('Temperature (°C)');
title('Temperature vs Distance');
legend('Location', 'best');
grid on;
hold off;

% ─────────────────────────────────────────────────────────────────────────────
% Plot 2: Stage Index vs Temperature (ResultsManager style)
% ─────────────────────────────────────────────────────────────────────────────
subplot(1, 2, 2);

% Silicon layer: stages 0 to N
r_chip = 0:N;
T_chip = [T_0; T_Si(:)];

% TEC layer: stages 1 to N+1 (includes water at N+1)
r_tec = 1:(N+1);
T_tec = [T_c(:); T_water_C];

% Plot Silicon layer (blue circles)
plot(r_chip, T_chip, '-ob', 'LineWidth', 1.8, 'MarkerFaceColor', 'b', 'DisplayName', 'Silicon Layer');
hold on;

% Plot TEC layer (red squares)
plot(r_tec, T_tec, '-sr', 'LineWidth', 1.8, 'MarkerFaceColor', 'r', 'DisplayName', 'TEC Layer');

% Mark the water temperature point distinctly (green triangle)
plot(N+1, T_water_C, 'g^', 'MarkerSize', 10, 'MarkerFaceColor', 'g', 'DisplayName', 'Coolant (T_{water})');

xlabel('Stage Index');
ylabel('Temperature (°C)');
title('Radial Temperature Profile');
legend('Location', 'best');
grid on;
hold off;

sgtitle(sprintf('Knee Point: T_{max} = %.1f°C, t_{TEC} = %.0fµm, I = %.4fA, N=%d, M=%d', ...
    T_max_C, t_TEC, I, N, M), 'FontWeight', 'bold');

%% Save figure
savefig(fig, fullfile(results_path, 'knee_point_analysis.fig'));
saveas(fig, fullfile(results_path, 'knee_point_analysis.png'));
fprintf('\nFigure saved to: %s\n', fullfile(results_path, 'knee_point_analysis.png'));

%% Export dimension table to CSV
dim_table = table();
dim_table.Stage = (1:N)';
dim_table.r_in_um = [stage_data.r_in_um]';
dim_table.r_out_um = [stage_data.r_out_um]';
dim_table.L_leg_um = [stage_data.L_um]';
dim_table.W_az_um = [stage_data.W_az_um]';
dim_table.W_leg_um = [stage_data.W_leg_um]';
dim_table.W_ic_um = [stage_data.W_ic_um]';
dim_table.t_ic_um = [stage_data.t_ic_um]';
dim_table.beta_ic_deg = [stage_data.beta_ic_deg]';
dim_table.W_oc_um = [stage_data.W_oc_um]';
dim_table.t_oc_um = [stage_data.t_oc_um]';
dim_table.beta_oc_deg = [stage_data.beta_oc_deg]';
dim_table.T_cold_C = T_C(2:2:2*N);
dim_table.T_hot_C = T_C(3:2:2*N+1);

writetable(dim_table, fullfile(results_path, 'knee_point_dimensions.csv'));
fprintf('Dimension table saved to: %s\n', fullfile(results_path, 'knee_point_dimensions.csv'));

fprintf('\n✓ Analysis complete!\n');

%% Helper function
function val = get_var_or_default(var_map, name, default)
    if var_map.isKey(name)
        val = var_map(name);
    else
        val = default;
    end
end
