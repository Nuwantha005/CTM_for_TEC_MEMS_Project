% Run Optimization Script
% Uses TECOptimizer with settings from src/config/optimization_variables.m

format longG;
clc;
clear all;

% Get the directory where this script is located (works both interactively and batch)
script_dir = fileparts(mfilename('fullpath'));
project_root = fullfile(script_dir, '..', '..');

addpath(genpath(fullfile(project_root, 'src')));

% Ensure output directory exists
output_dir = fullfile(project_root, 'output');
if ~exist(output_dir, 'dir')
    mkdir(output_dir);
end

% Load centralized configuration
[~, ~, ~, ~, ~, CONFIG] = optimization_variables();

% Create optimizer
config_path = fullfile(project_root, 'src', 'config', 'default_params.json');
optimizer = TECOptimizer(config_path);

% Apply configuration from centralized file
optimizer.Solver.Config.boundary_conditions.q_flux_W_m2 = CONFIG.q_flux_W_m2;
optimizer.Solver.Config.boundary_conditions.T_water_K = CONFIG.T_water_K;
optimizer.Solver.Config.boundary_conditions.h_conv_W_m2K = CONFIG.h_conv_W_m2K;

fprintf('=== Settings from optimization_variables.m ===\n');
fprintf('Target heat flux: %.0f W/m²\n', CONFIG.q_flux_W_m2);
fprintf('Coolant temperature: %.1f K\n', CONFIG.T_water_K);
fprintf('Convection coefficient: %.0e W/m²K\n', CONFIG.h_conv_W_m2K);
fprintf('Target max temperature: %d°C\n\n', CONFIG.T_target_C);

optimizer.run_optimization();
