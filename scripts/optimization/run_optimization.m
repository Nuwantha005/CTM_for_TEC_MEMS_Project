% Run Optimization Script

format longG;
clc;
clear all;

addpath(genpath('src'));

% Ensure output directory exists
if ~exist('output', 'dir')
    mkdir('output');
end

% Create optimizer
optimizer = TECOptimizer('src/config/default_params.json');

% Set a feasible heat flux (from analysis: max ~2500 W/m² for T < 100°C with 200um TEC)
% Use 1000 W/m² for optimization with margin
target_heat_flux = 200000;  % W/m²
optimizer.Solver.Config.boundary_conditions.q_flux_W_m2 = target_heat_flux;
fprintf('Target heat flux set to: %.0f W/m² (%.2f W total chip)\n\n', ...
    target_heat_flux, target_heat_flux * 100e-6);

optimizer.run_optimization();
