% reorganize_codebase.m
% Script to reorganize the codebase into a cleaner structure
%
% New structure:
%   scripts/
%       optimization/     - Optimization scripts (global, multi-objective, etc.)
%       analysis/         - Analysis and diagnostic scripts
%       comsol/          - COMSOL-related scripts
%   src/                 - Core library code (unchanged)
%   output/              - All output files
%   data/                - Input data files
%   notes/               - Documentation

fprintf('=== CODEBASE REORGANIZATION ===\n\n');

%% Define file mappings
% Format: {old_location, new_location, description}

moves = {
    % Optimization scripts -> scripts/optimization/
    'run_optimization.m', 'scripts/optimization/run_optimization.m', 'Main optimization';
    'run_design_optimization.m', 'scripts/optimization/run_design_optimization.m', 'Design optimization';
    'run_physics_constrained.m', 'scripts/optimization/run_physics_constrained.m', 'Physics constrained opt';
    'find_feasible_point.m', 'scripts/optimization/find_feasible_point.m', 'Feasibility finder';
    'check_feasibility.m', 'scripts/optimization/check_feasibility.m', 'Feasibility checker';
    
    % Analysis scripts -> scripts/analysis/
    'analyze_cooling_capacity.m', 'scripts/analysis/analyze_cooling_capacity.m', 'Cooling capacity analysis';
    'diagnose_thermal_model.m', 'scripts/analysis/diagnose_thermal_model.m', 'Thermal model diagnostics';
    'diagnose_network_debug.m', 'scripts/analysis/diagnose_network_debug.m', 'Network debug';
    'diagnose_optimal_current.m', 'scripts/analysis/diagnose_optimal_current.m', 'Optimal current analysis';
    'explore_design_improvements.m', 'scripts/analysis/explore_design_improvements.m', 'Design exploration';
    'test_various_conditions.m', 'scripts/analysis/test_various_conditions.m', 'Condition testing';
    'plot_comsol_comparison.m', 'scripts/analysis/plot_comsol_comparison.m', 'COMSOL comparison plots';
    
    % COMSOL scripts -> scripts/comsol/
    'run_single_comsol.m', 'scripts/comsol/run_single_comsol.m', 'Single COMSOL run';
    'run_overnight_comsol.m', 'scripts/comsol/run_overnight_comsol.m', 'Overnight COMSOL batch';
    'run_comsol_validation.m', 'scripts/comsol/run_comsol_validation.m', 'COMSOL validation';
    'run_comsol_test.m', 'scripts/comsol/run_comsol_test.m', 'COMSOL test';
    'test_comsol_connection.m', 'scripts/comsol/test_comsol_connection.m', 'COMSOL connection test';
    'extract_comsol_data.m', 'scripts/comsol/extract_comsol_data.m', 'COMSOL data extraction';
    'prepare_comsol_candidates.m', 'scripts/comsol/prepare_comsol_candidates.m', 'COMSOL candidate prep';
    'generate_comsol_design.m', 'scripts/comsol/generate_comsol_design.m', 'COMSOL design generator';
    
    % Sweep/solver scripts -> scripts/sweeps/
    'run_solver.m', 'scripts/sweeps/run_solver.m', 'Solver runner';
    'run_sweep.m', 'scripts/sweeps/run_sweep.m', 'Parameter sweep';
    'run_stage_sweep.m', 'scripts/sweeps/run_stage_sweep.m', 'Stage sweep';
    'run_quick_test.m', 'scripts/sweeps/run_quick_test.m', 'Quick test';
};

%% Create directories
dirs_to_create = {
    'scripts/optimization'
    'scripts/analysis'
    'scripts/comsol'
    'scripts/sweeps'
};

fprintf('Creating directories...\n');
for i = 1:length(dirs_to_create)
    if ~exist(dirs_to_create{i}, 'dir')
        mkdir(dirs_to_create{i});
        fprintf('  Created: %s\n', dirs_to_create{i});
    else
        fprintf('  Exists:  %s\n', dirs_to_create{i});
    end
end
fprintf('\n');

%% Move files
fprintf('Moving files...\n');
fprintf('─────────────────────────────────────────────────────────────\n');

moved = 0;
skipped = 0;

for i = 1:size(moves, 1)
    old_path = moves{i, 1};
    new_path = moves{i, 2};
    desc = moves{i, 3};
    
    if exist(old_path, 'file')
        % Check if already moved
        if exist(new_path, 'file')
            fprintf('  SKIP (exists): %s\n', old_path);
            skipped = skipped + 1;
        else
            % Move the file
            movefile(old_path, new_path);
            fprintf('  MOVED: %s -> %s\n', old_path, new_path);
            moved = moved + 1;
        end
    else
        % Check if already in new location
        if exist(new_path, 'file')
            fprintf('  OK (already moved): %s\n', new_path);
        else
            fprintf('  NOT FOUND: %s\n', old_path);
        end
    end
end

fprintf('─────────────────────────────────────────────────────────────\n');
fprintf('Moved: %d, Skipped: %d\n\n', moved, skipped);

%% Create main runner script
fprintf('Creating main runner script...\n');

main_script = sprintf([...
'%% main.m\n' ...
'%% Main entry point for TEC optimization project\n' ...
'%%\n' ...
'%% Usage:\n' ...
'%%   1. Run setup: run(''main.m'')\n' ...
'%%   2. Then run any script from scripts/ folder\n' ...
'\n' ...
'clear; clc;\n' ...
'fprintf(''TEC MEMS Optimization Project\\n'');\n' ...
'fprintf(''============================\\n\\n'');\n' ...
'\n' ...
'%% Add paths\n' ...
'addpath(genpath(''src''));\n' ...
'addpath(genpath(''scripts''));\n' ...
'\n' ...
'fprintf(''Paths added. Available scripts:\\n\\n'');\n' ...
'\n' ...
'fprintf(''OPTIMIZATION:\\n'');\n' ...
'fprintf(''  run_global_optimization      - Global optimization (GA + PSO)\\n'');\n' ...
'fprintf(''  run_multiobjective_optimization - Multi-objective (Pareto front)\\n'');\n' ...
'fprintf(''  run_optimization             - Local optimization (fmincon)\\n'');\n' ...
'fprintf(''\\n'');\n' ...
'\n' ...
'fprintf(''ANALYSIS:\\n'');\n' ...
'fprintf(''  analyze_cooling_capacity     - What can be cooled\\n'');\n' ...
'fprintf(''  diagnose_thermal_model       - Model diagnostics\\n'');\n' ...
'fprintf(''  plot_comsol_comparison       - COMSOL vs MATLAB plots\\n'');\n' ...
'fprintf(''\\n'');\n' ...
'\n' ...
'fprintf(''COMSOL:\\n'');\n' ...
'fprintf(''  run_single_comsol            - Single COMSOL simulation\\n'');\n' ...
'fprintf(''  run_overnight_comsol         - Batch COMSOL runs\\n'');\n' ...
'fprintf(''  prepare_comsol_candidates    - Generate COMSOL candidates\\n'');\n' ...
'fprintf(''\\n'');\n' ...
'\n' ...
'fprintf(''SWEEPS:\\n'');\n' ...
'fprintf(''  run_solver                   - Run thermal solver\\n'');\n' ...
'fprintf(''  run_sweep                    - Parameter sweep\\n'');\n' ...
'fprintf(''\\n'');\n' ...
]);

fid = fopen('main.m', 'w');
fprintf(fid, '%s', main_script);
fclose(fid);
fprintf('  Created: main.m\n\n');

%% Summary
fprintf('╔════════════════════════════════════════════════════════════╗\n');
fprintf('║              REORGANIZATION COMPLETE                       ║\n');
fprintf('╚════════════════════════════════════════════════════════════╝\n\n');

fprintf('New structure:\n');
fprintf('  scripts/\n');
fprintf('    optimization/   - Optimization scripts\n');
fprintf('    analysis/       - Analysis scripts\n');
fprintf('    comsol/         - COMSOL scripts\n');
fprintf('    sweeps/         - Sweep/solver scripts\n');
fprintf('  src/              - Core library\n');
fprintf('  output/           - Results\n');
fprintf('  data/             - Input data\n');
fprintf('  notes/            - Documentation\n');
fprintf('\n');
fprintf('To get started: run(''main.m'')\n');
