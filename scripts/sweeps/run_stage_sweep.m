% Run Stage Optimization Sweep
% Uses TECOptimizer internally, which reads from optimization_variables.m

format longG;
clc;
clear all;

addpath(genpath('src'));

% Ensure output directory exists
if ~exist('output', 'dir')
    mkdir('output');
end

% Load centralized configuration to display settings
[~, ~, ~, ~, ~, CONFIG] = optimization_variables();
fprintf('=== Settings from optimization_variables.m ===\n');
fprintf('Heat flux: %.0e W/m²\n', CONFIG.q_flux_W_m2);
fprintf('Coolant temperature: %.1f K\n', CONFIG.T_water_K);
fprintf('Target max temperature: %d°C\n\n', CONFIG.T_target_C);

sweeper = StageSweeper('src/config/default_params.json');

% Sweep from 2 to 12 stages (or customize as needed)
sweeper.run_stage_sweep(2, 12);
