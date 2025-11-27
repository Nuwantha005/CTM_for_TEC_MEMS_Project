% Run Optimization Script
% Uses TECOptimizer with settings from src/config/optimization_variables.m

format longG;
clc;
clear all;

addpath(genpath('src'));

% Ensure output directory exists
if ~exist('output', 'dir')
    mkdir('output');
end

% Load centralized configuration
[~, ~, ~, ~, ~, CONFIG] = optimization_variables();

% Create optimizer
optimizer = TECOptimizer('src/config/default_params.json');

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
