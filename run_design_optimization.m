%% RADIAL TEC DESIGN OPTIMIZATION
% This script runs a comprehensive parameter sweep to find optimal TEC designs
% that can cool the chip below the target temperature (100°C by default).
%
% Outputs:
%   - Ranked list of feasible designs
%   - Publication-ready plots
%   - Best design configuration (JSON) for COMSOL import
%   - CSV of all feasible designs
%
% Author: Auto-generated
% Date: 2025-11-25

clear; clc; close all;

%% Add paths
addpath(genpath('src'));

%% Configuration
config_path = 'src/config/default_params.json';

% Modify config for your heat flux
config = jsondecode(fileread(config_path));
config.boundary_conditions.q_flux_W_m2 = 500;  % W/m² - adjust this!
fid = fopen(config_path, 'w');
fprintf(fid, '%s', jsonencode(config));
fclose(fid);

% Create optimizer
optimizer = DesignOptimizer(config_path);

%% Set optimization targets
optimizer.T_target_C = 100;  % Target: chip below 100°C

%% Define parameter ranges to sweep
% Adjust these based on your constraints

% Number of stages to test
N_stages_list = [2, 3, 4, 5, 6];

% Wedge angles (degrees) - affects number of wedges per chip
% 360/theta = number of wedges (e.g., 30° -> 12 wedges, 20° -> 18 wedges)
theta_list = [15, 20, 25, 30, 45];

% TEC thickness (μm) - constrained by 3D stacking requirements
t_TEC_list = [40, 50, 60, 80, 100];

% Radial expansion factor - how each stage length grows
k_r_list = [0.8, 0.9, 1.0, 1.1, 1.2];

%% Run optimization
fprintf('Starting optimization...\n');
fprintf('This may take several minutes depending on parameter ranges.\n\n');

results = optimizer.runFullOptimization(...
    'N_stages_list', N_stages_list, ...
    'Theta_list', theta_list, ...
    't_TEC_list', t_TEC_list, ...
    'k_r_list', k_r_list, ...
    'Verbose', true);

%% Display results location
fprintf('\n========================================\n');
fprintf('OPTIMIZATION COMPLETE\n');
fprintf('========================================\n');
fprintf('Results saved to: %s\n', optimizer.OutputDir);
fprintf('\nFiles generated:\n');
fprintf('  - optimization_results.mat (MATLAB data)\n');
fprintf('  - optimization_summary.txt (Human-readable summary)\n');
fprintf('  - best_design_config.json (For COMSOL import)\n');
fprintf('  - feasible_designs.csv (All feasible designs)\n');
fprintf('  - *.png, *.fig (Plots)\n');

%% Quick visualization of results
if isfield(results, 'best_feasible')
    fprintf('\n========================================\n');
    fprintf('BEST DESIGN FOR COMSOL VALIDATION\n');
    fprintf('========================================\n');
    best = results.best_feasible;
    fprintf('Copy these parameters to COMSOL:\n');
    fprintf('  N_stages = %d\n', best.N_stages);
    fprintf('  wedge_angle = %.1f degrees\n', best.theta_deg);
    fprintf('  t_TEC = %.0f μm\n', best.t_TEC_um);
    fprintf('  k_r = %.2f\n', best.k_r);
    fprintf('  I_current = %.4f A (%.2f mA)\n', best.I_opt, best.I_opt*1000);
    fprintf('  Expected T_max = %.1f°C\n', best.T_max_C);
    fprintf('  COP = %.3f\n', best.COP);
end

%% Open results folder
fprintf('\nOpening results folder...\n');
winopen(optimizer.OutputDir);
