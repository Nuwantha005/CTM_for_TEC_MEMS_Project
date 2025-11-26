%% PHYSICS-CONSTRAINED RADIAL TEC OPTIMIZATION
% This script uses a physics-constrained approach to find feasible TEC designs.
%
% Key improvements:
%   1. Enforces physical constraints (T >= T_water, realistic COP)
%   2. Validates solutions for physical consistency
%   3. Uses analytical estimates to guide search
%   4. Handles high heat flux scenarios better
%
% Author: Auto-generated
% Date: 2025-11-26

clear; clc; close all;

%% Add paths
addpath(genpath('src'));

%% Configuration
config_path = 'src/config/default_params.json';

% Load and modify config for your heat flux
config = jsondecode(fileread(config_path));

% =====================================================
% IMPORTANT: Set your heat flux here
% Typical values:
%   - Low power chip: 10,000 W/m² (10 kW/m²)
%   - Medium power: 50,000 W/m² (50 kW/m²)  
%   - High power: 100,000 W/m² (100 kW/m²)
%   - Extreme: 200,000+ W/m² (200+ kW/m²)
%
% Note: At very high heat flux, TEC may not provide enough cooling
% and you may need to reduce flux or accept higher temperatures.
% =====================================================
config.boundary_conditions.q_flux_W_m2 = 100000;  % 100 kW/m²

% Save modified config
fid = fopen(config_path, 'w');
fprintf(fid, '%s', jsonencode(config));
fclose(fid);

%% Create optimizer
optimizer = PhysicsConstrainedOptimizer(config_path);

%% Set target temperature (°C)
% Maximum allowable chip temperature
optimizer.T_max_K = 373.15;  % 100°C = 373.15 K

%% Run optimization
fprintf('Starting physics-constrained optimization...\n');
fprintf('This searches for designs that satisfy physical constraints.\n\n');

results = optimizer.runConstrainedOptimization(...
    'MaxDesigns', 1000, ...  % Max designs to evaluate
    'Verbose', true);

%% Display summary
if isfield(results, 'best_design')
    fprintf('\n\n========================================\n');
    fprintf('  READY FOR COMSOL SIMULATION\n');
    fprintf('========================================\n');
    fprintf('Best design configuration saved to:\n');
    fprintf('  %s/best_design_config.json\n', optimizer.OutputDir);
    fprintf('\nUse this configuration for high-fidelity COMSOL validation.\n');
else
    fprintf('\n\n========================================\n');
    fprintf('  NO FEASIBLE DESIGN FOUND\n');
    fprintf('========================================\n');
    fprintf('Suggestions:\n');
    fprintf('  1. Reduce heat flux (q_flux_W_m2)\n');
    fprintf('  2. Increase target temperature (T_max_K)\n');
    fprintf('  3. Use larger chip area (spreads heat over more wedges)\n');
    fprintf('  4. Consider alternative cooling approaches\n');
end

fprintf('\nOutput directory: %s\n', optimizer.OutputDir);
