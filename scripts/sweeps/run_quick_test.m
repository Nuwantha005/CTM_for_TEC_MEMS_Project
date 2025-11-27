%% QUICK TEST - Radial TEC Design Optimization
% A smaller test run to verify the optimizer works before full sweep
%
% This tests with fewer parameter combinations for quick feedback

clear; clc; close all;

%% Add paths
% Ensure src is on the MATLAB path and config is resolved from project root
config_path = get_config_path();
src_dir = fileparts(fileparts(config_path));
addpath(genpath(src_dir));

%% Configuration
% config_path already set from helper

% Create optimizer
optimizer = DesignOptimizer(config_path);

%% Set targets
optimizer.T_target_C = 100;  % Target: chip below 100°C

%% Quick test with minimal parameters
fprintf('Running quick test with minimal parameter set...\n\n');

results = optimizer.runFullOptimization(...
    'N_stages_list', [3, 4, 5], ...
    'Theta_list', [20, 30], ...
    't_TEC_list', [50, 80], ...
    'k_r_list', [1.0, 1.1], ...
    'Verbose', true);

%% Summary
fprintf('\n========================================\n');
fprintf('QUICK TEST COMPLETE\n');
fprintf('========================================\n');

if results.summary.feasible_count > 0
    fprintf('SUCCESS: Found %d feasible designs!\n', results.summary.feasible_count);
    fprintf('Results saved to: %s\n', optimizer.OutputDir);
    
    % Open output folder
    fprintf('\nOpening results folder...\n');
    winopen(optimizer.OutputDir);
else
    fprintf('WARNING: No feasible designs found in quick test.\n');
    fprintf('Try:\n');
    fprintf('  1. Reducing heat flux (q_flux_W_m2 in config)\n');
    fprintf('  2. Increasing TEC thickness range\n');
    fprintf('  3. Checking for solver errors\n');
end
