% main.m
% Main entry point for TEC optimization project
%
% Usage:
%   1. Run setup: run('main.m')
%   2. Then run any script from scripts/ folder

clear; clc;
fprintf('TEC MEMS Optimization Project\n');
fprintf('============================\n\n');

% Add paths
addpath(genpath('src'));
addpath(genpath('scripts'));

fprintf('Paths added. Available scripts:\n\n');

fprintf('OPTIMIZATION:\n');
fprintf('  run_global_optimization      - Global optimization (GA + PSO)\n');
fprintf('  run_multiobjective_optimization - Multi-objective (Pareto front)\n');
fprintf('  run_optimization             - Local optimization (fmincon)\n');
fprintf('\n');

fprintf('ANALYSIS:\n');
fprintf('  analyze_cooling_capacity     - What can be cooled\n');
fprintf('  diagnose_thermal_model       - Model diagnostics\n');
fprintf('  plot_comsol_comparison       - COMSOL vs MATLAB plots\n');
fprintf('\n');

fprintf('COMSOL:\n');
fprintf('  run_single_comsol            - Single COMSOL simulation\n');
fprintf('  run_overnight_comsol         - Batch COMSOL runs\n');
fprintf('  prepare_comsol_candidates    - Generate COMSOL candidates\n');
fprintf('\n');

fprintf('SWEEPS:\n');
fprintf('  run_solver                   - Run thermal solver\n');
fprintf('  run_sweep                    - Parameter sweep\n');
fprintf('\n');
