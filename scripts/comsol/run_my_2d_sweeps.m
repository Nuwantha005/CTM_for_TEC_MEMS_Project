%% CUSTOMIZABLE 2D SWEEP RUNNER (1D vs 3D Validation)
% This script sets up a 2D parameter sweep grid and calls the generalized function.

clear; clc; close all;

%% ============ PATH FIX & SETUP ============
root_dir = fullfile(fileparts(mfilename('fullpath')), '..', '..');
cd(root_dir);
addpath(genpath('src'));
addpath(genpath('scripts')); % Add scripts folder to path so helper functions are found

% Make sure COMSOL mli path is added
addpath('F:\EngineeringSoftware\COMSOL\COMSOL63\Multiphysics\mli');

COMSOL_PORT = 2036;
COMSOL_MODEL_PATH = 'E:\Semester 7\ME4311 - MicroNano Electro Mechanical Systems and Nanotechnology\Project\COMSOL\Test_Wedge\asym2.mph';      
COMSOL_SERVER_EXE = 'F:\EngineeringSoftware\COMSOL\COMSOL63\Multiphysics\bin\win64\comsolmphserver.exe';

timestamp = datestr(now, 'yyyy-mm-dd_HH-MM-SS');
OUTPUT_DIR = fullfile('output', 'validity_sweeps', timestamp);
if ~exist(OUTPUT_DIR, 'dir'), mkdir(OUTPUT_DIR); end

%% ============ BASELINE SETUP (Case 2 Anchor) ============
base.base_M = 12;
base.base_q = 500;
base.base_t_TEC = 150;

base.t_chip_um = 50; base.t_SOI_um = 20; base.R_cyl_um = 1000;
base.f_L = 1.15; base.fill_factor = 0.9; base.w_is_um = 50;
base.I_A = 0.10; base.T_water = 293.15;
base.ic_w_r = 0.1; base.ic_t_r = 0.6; base.ic_angle_r = 0.5;
base.oc_w_r = 0.1; base.oc_t_r = 0.6; base.oc_angle_r = 0.5;

%% ============ CONFIGURE YOUR 2D SWEEP PAIRS HERE ============
% Add as many pairs as you want; they will run sequentially in one session.
% Fields: name, var1, values1, var2, values2
vals_fL = linspace(0.8, 1.2, 10);
vals_fill = linspace(0.85, 0.99, 10);
vals_wis = linspace(30, 50, 10);
vals_q = linspace(200, 600, 10);
vals_ttec = linspace(100, 300, 10);

sweep_pairs = {
    %struct('name', 'LengthRatio_vs_FillFactor',        'var1', 'f_L',         'values1', vals_fL,   'var2', 'fill_factor', 'values2', vals_fill);
    %struct('name', 'LengthRatio_vs_RadialInsWidth',    'var1', 'f_L',         'values1', vals_fL,   'var2', 'w_is_um',     'values2', vals_wis);
    %struct('name', 'LengthRatio_vs_HeatFlux',          'var1', 'f_L',         'values1', vals_fL,   'var2', 'q_Wm2',       'values2', vals_q);
    struct('name', 'LengthRatio_vs_Thickness',         'var1', 'f_L',         'values1', vals_fL,   'var2', 't_TEC_um',    'values2', vals_ttec);
    struct('name', 'FillFactor_vs_RadialInsWidth',     'var1', 'fill_factor', 'values1', vals_fill, 'var2', 'w_is_um',     'values2', vals_wis);
    struct('name', 'FillFactor_vs_HeatFlux',           'var1', 'fill_factor', 'values1', vals_fill, 'var2', 'q_Wm2',       'values2', vals_q);
    struct('name', 'FillFactor_vs_Thickness',          'var1', 'fill_factor', 'values1', vals_fill, 'var2', 't_TEC_um',    'values2', vals_ttec);
    struct('name', 'RadialInsWidth_vs_HeatFlux',       'var1', 'w_is_um',     'values1', vals_wis,  'var2', 'q_Wm2',       'values2', vals_q);
    struct('name', 'RadialInsWidth_vs_Thickness',      'var1', 'w_is_um',     'values1', vals_wis,  'var2', 't_TEC_um',    'values2', vals_ttec);
    struct('name', 'HeatFlux_vs_Thickness',            'var1', 'q_Wm2',       'values1', vals_q,    'var2', 't_TEC_um',    'values2', vals_ttec);
};

%% Ensure COMSOL is ready
fprintf('Killing old COMSOL server...\n');
system('taskkill /F /IM comsolmphserver.exe 2>nul', '-echo');
pause(2);
fprintf('Starting COMSOL server...\n');
system(sprintf('start "" "%s" -port %d', COMSOL_SERVER_EXE, COMSOL_PORT));
pause(15);
try
    mphstart(COMSOL_PORT);
    import com.comsol.model.util.*
    fprintf('Connected to COMSOL!\n\n');
catch
    fprintf('Failed to connect to COMSOL on port %d.\n', COMSOL_PORT);
    return;
end

%% ============ EXECUTE ALL 2D SWEEP PAIRS ============
num_pairs = length(sweep_pairs);
total_runs_all = 0;
for k = 1:num_pairs
    total_runs_all = total_runs_all + length(sweep_pairs{k}.values1) * length(sweep_pairs{k}.values2);
end

fprintf('Starting batch 2D sweeps (%d pairs, %d total runs)\n', num_pairs, total_runs_all);

summary_rows = table();
for k = 1:num_pairs
    pair_cfg = sweep_pairs{k};
    pair_runs = length(pair_cfg.values1) * length(pair_cfg.values2);
    pair_output_dir = fullfile(OUTPUT_DIR, pair_cfg.name);
    if ~exist(pair_output_dir, 'dir'), mkdir(pair_output_dir); end

    fprintf('\n============================================================\n');
    fprintf('Pair %d/%d: %s\n', k, num_pairs, pair_cfg.name);
    fprintf('%s (%d pts)  vs  %s (%d pts)  =>  %d runs\n', ...
        pair_cfg.var1, length(pair_cfg.values1), pair_cfg.var2, length(pair_cfg.values2), pair_runs);
    fprintf('Output: %s\n', pair_output_dir);

    t_start = tic;
    results = run_2d_validity_sweep(pair_cfg.name, pair_cfg.var1, pair_cfg.values1, ...
                                    pair_cfg.var2, pair_cfg.values2, ...
                                    base, pair_output_dir, COMSOL_MODEL_PATH);
    elapsed_s = toc(t_start);

    summary_rows = [summary_rows; table(categorical({pair_cfg.name}), ...
                    categorical({pair_cfg.var1}), categorical({pair_cfg.var2}), ...
                    length(pair_cfg.values1), length(pair_cfg.values2), pair_runs, ...
                    height(results), elapsed_s, ...
                    'VariableNames', {'SweepName', 'Var1', 'Var2', 'N1', 'N2', 'PlannedRuns', 'CompletedRuns', 'Elapsed_s'})];

    writetable(summary_rows, fullfile(OUTPUT_DIR, 'batch_summary.csv'));
    fprintf('Completed pair %d/%d in %.1f min\n', k, num_pairs, elapsed_s/60);
end

fprintf('\n\n=== ALL 2D SWEEP PAIRS COMPLETE ===\n');
fprintf('Batch root output: %s\n', OUTPUT_DIR);
fprintf('Summary CSV: %s\n', fullfile(OUTPUT_DIR, 'batch_summary.csv'));