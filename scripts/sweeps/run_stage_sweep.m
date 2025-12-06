% Run Stage Optimization Sweep
% Uses TECOptimizer internally, which reads from optimization_variables.m
% Notation follows: Thermal_Network_For_Radial_TEC.tex
%
% INTEGER PARAMETERS:
%   N: Number of stages (swept in this script)
%   M: Number of wedges (fixed, use run_wedge_sweep.m to sweep M)

format longG;
clc;
clear all;

% Get the directory where this script is located
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

fprintf('══════════════════════════════════════════════════════════════\n');
fprintf('         STAGE SWEEP (N) - Paper Notation Reference           \n');
fprintf('══════════════════════════════════════════════════════════════\n');
fprintf('=== Settings from optimization_variables.m ===\n');
fprintf('Heat flux (q_flux):     %.0e W/m²\n', CONFIG.q_flux_W_m2);
fprintf('Coolant temp (T_water): %.1f K\n', CONFIG.T_water_K);
fprintf('Target max temp:        %d°C\n', CONFIG.T_target_C);
fprintf('\nInteger Parameters:\n');
fprintf('  N (stages):   %d  (will sweep 2→12)\n', CONFIG.N);
fprintf('  M (wedges):   %d  →  θ = %.1f°  (fixed)\n', CONFIG.M, CONFIG.theta_deg);
fprintf('\n');

config_path = fullfile(project_root, 'src', 'config', 'default_params.json');
sweeper = StageSweeper(config_path);

% Sweep from 2 to 12 stages (N)
sweeper.run_stage_sweep(2, 12);
