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

%% ============ CONFIGURE YOUR 2D SWEEP HERE ============
sweep_name = 'HeatFlux_vs_Thickness';
var_name1  = 'q_Wm2';
values1    = linspace(100, 1000, 10);

var_name2  = 't_TEC_um';
values2    = linspace(50, 500, 10);

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

%% ============ EXECUTE 2D SWEEP ============
fprintf('Starting 2D Sweep: %s\n', sweep_name);
fprintf('Grid size: %d x %d = %d runs total\n', length(values1), length(values2), length(values1)*length(values2));

results = run_2d_validity_sweep(sweep_name, var_name1, values1, var_name2, values2, ...
                                base, OUTPUT_DIR, COMSOL_MODEL_PATH);

fprintf('\n\n=== 2D SWEEP COMPLETE ===\n');
fprintf('Results and 2D plots saved to %s\n', OUTPUT_DIR);