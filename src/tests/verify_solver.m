% verify_solver.m
% Script to verify the Radial TEC Solver implementation

clear all;
clc;

% 1. Setup Path
addpath(genpath('../core'));
addpath(genpath('../utils'));
addpath(genpath('../config'));

fprintf('Running Verification Script...\n');

% 2. Define Config Path
config_path = '../config/default_params.json';

% 3. Initialize Solver
try
    solver = RadialTECSolver(config_path);
    fprintf('Solver initialized successfully.\n');
catch e
    fprintf('Error initializing solver: %s\n', e.message);
    return;
end

% 4. Run Solver
try
    solver.run();
    fprintf('Solver run completed.\n');
catch e
    fprintf('Error running solver: %s\n', e.message);
    disp(e.stack(1));
    return;
end

% 5. Verify Results
% Find the latest output folder
results_mgr = solver.Results;
output_dir = results_mgr.OutputDir;
results_file = fullfile(output_dir, 'results.mat');

if exist(results_file, 'file')
    fprintf('Results file found: %s\n', results_file);
    load(results_file);
    
    % Check T dimensions
    N = config.geometry.N_stages;
    expected_dim = 2*N + 1;
    if length(T) == expected_dim
        fprintf('Temperature vector dimension correct: %d\n', length(T));
    else
        fprintf('ERROR: Temperature vector dimension mismatch. Expected %d, got %d\n', expected_dim, length(T));
    end
    
    % Check Energy Balance (Roughly)
    Q_gen = config.boundary_conditions.q_flux_W_m2;
    fprintf('Q_gen (Input): %.2f W\n', Q_gen);
    fprintf('Q_out (Output): %.2f W\n', Q_out);
    
    % Note: Q_out should be > Q_gen because of TEC power input
    if Q_out > Q_gen
        fprintf('Energy Balance Check: PASS (Q_out > Q_gen, implies active cooling/heating)\n');
    else
        fprintf('Energy Balance Check: WARNING (Q_out <= Q_gen). Check if I=0 or cooling is effective?\n');
    end
    
    % Check Node 0 Temp
    fprintf('Center Temperature (T0): %.2f K\n', T(1));
    
    % Check COMSOL Export
    comsol_file = fullfile(output_dir, 'comsol_params.txt');
    solver.Geometry.export_comsol_params(comsol_file);
    if exist(comsol_file, 'file')
        fprintf('COMSOL export successful: %s\n', comsol_file);
    else
        fprintf('ERROR: COMSOL export failed.\n');
    end
    
else
    fprintf('ERROR: Results file not found.\n');
end

fprintf('Verification Complete.\n');
