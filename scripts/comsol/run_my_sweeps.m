%% CUSTOMIZABLE SWEEP RUNNER (1D vs 3D Validation)
% This script sets up flexible sweeps and calls the generalized function.
% You can easily modify ranges, add new sweep parameters, or change baselines.

clear; clc; close all;

%% ============ PATH FIX & SETUP ============
% Ensure we are in the workspace root so MATLAB finds all `src` classes
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
base.base_M = 8;
base.base_q = 500;
base.base_t_TEC = 150;

base.t_chip_um = 50; base.t_SOI_um = 20; base.R_cyl_um = 1000;
base.f_L = 1.15; base.fill_factor = 0.9; base.w_is_um = 50;
base.I_A = 0.10; base.T_water = 293.15;
base.ic_w_r = 0.1; base.ic_t_r = 0.6; base.ic_angle_r = 0.5;
base.oc_w_r = 0.1; base.oc_t_r = 0.6; base.oc_angle_r = 0.5;

%% ============ CONFIGURE YOUR SWEEPS HERE ============
% Format: {'Sweep Name Label', 'MATLAB_Variable_Name', [array_of_values]}
% Supported Variables: 'q_Wm2', 'M', 't_TEC_um', 'f_L', etc.
sweeps_to_run = {
    {'Current',  'I_A', linspace(0,1.0,11)},
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

%% ============ EXECUTE SWEEPS ============
all_results = table();
run_idx = 1;

for i = 1:length(sweeps_to_run)
    sweep_name = sweeps_to_run{i}{1};
    var_name = sweeps_to_run{i}{2};
    values = sweeps_to_run{i}{3};
    
    [res, next_idx] = run_single_validity_sweep(sweep_name, var_name, values, ...
                                                base, OUTPUT_DIR, run_idx, COMSOL_MODEL_PATH);
    
    all_results = [all_results; res];
    run_idx = next_idx;
end

% Save the master compilation
writetable(all_results, fullfile(OUTPUT_DIR, 'compiled_all_sweeps.csv'));
fprintf('\n\n=== ALL SWEEPS COMPLETE ===\n');
fprintf('Compiled results saved to %s\n', fullfile(OUTPUT_DIR, 'compiled_all_sweeps.csv'));