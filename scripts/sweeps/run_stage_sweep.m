% Run Stage Optimization Sweep
% Uses TECOptimizer internally, which reads from optimization_variables.m

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

% Load centralized configuration to display settings
[~, ~, ~, ~, ~, CONFIG] = optimization_variables();
fprintf('=== Settings from optimization_variables.m ===\n');
fprintf('Heat flux: %.0e W/m²\n', CONFIG.q_flux_W_m2);
fprintf('Coolant temperature: %.1f K\n', CONFIG.T_water_K);
fprintf('Target max temperature: %d°C\n\n', CONFIG.T_target_C);

config_path = fullfile(project_root, 'src', 'config', 'default_params.json');
sweeper = StageSweeper(config_path);

% Sweep from 2 to 12 stages (or customize as needed)
sweeper.run_stage_sweep(2, 12);
