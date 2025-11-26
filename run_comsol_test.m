%% TEST COMSOL WITH DEFAULT PARAMETERS
% This script tests COMSOL with the default parameters from COMSOL_parameters.txt
% Use this to verify COMSOL LiveLink is working before running overnight simulations
%
% What this does:
%   1. Connects to COMSOL server
%   2. Loads your template model
%   3. Displays current parameters (should match COMSOL_parameters.txt)
%   4. Runs one simulation with defaults
%   5. Extracts and displays results
%
% Usage:
%   1. Start COMSOL server first:
%      - Option A: Run in command prompt: comsolmphserver -port 2036
%      - Option B: From COMSOL Desktop: File > Client/Server > Start Server
%      - Option C: Run MATLAB through COMSOL: comsol mphserver matlab
%   2. Run this script in MATLAB

clear; clc;
addpath(genpath('src'));

%% ==================== CONFIGURATION ====================

% Set your COMSOL model path here
MODEL_PATH = 'E:\Semester 7\ME4311 - MicroNano Electro Mechanical Systems and Nanotechnology\Project\COMSOL\MATLAB_API\template_3_stage.mph';

% Default parameters from COMSOL_parameters.txt (for reference)
DEFAULT_PARAMS = struct(...
    'I0',           0.15625, ...   % [A] - Current
    'q',            1000000, ...   % [W/m^2] - Heat flux
    'LL_beat_oc',   25, ...        % [deg] - Outer conductor beat angle
    'LL_beta_ic',   25, ...        % [deg] - Inner conductor beta angle
    'LL_k_r',       1.2, ...       % [-] - Radial length ratio
    'LL_L_1',       1, ...         % [mm] - Base length
    'LL_P_TSV',     30, ...        % [um] - TSV pitch
    'LL_r_chip',    10, ...        % [mm] - Chip radius
    'LL_R_cyl',     1, ...         % [mm] - Cylinder radius
    'LL_R_TSV',     10, ...        % [um] - TSV radius
    'LL_t_chip',    50, ...        % [um] - Chip thickness
    'LL_t_ic',      150, ...       % [um] - Inner conductor thickness
    'LL_t_oc',      150, ...       % [um] - Outer conductor thickness
    'LL_t_SOI',     100, ...       % [um] - SOI thickness
    'LL_t_TEC',     300, ...       % [um] - TEC thickness
    'LL_theta',     30, ...        % [deg] - Wedge angle
    'LL_w_az',      30, ...        % [um] - Azimuthal width
    'LL_w_ic',      50, ...        % [um] - Inner conductor width
    'LL_w_is',      40, ...        % [um] - Insulator width
    'LL_w_oc',      50 ...         % [um] - Outer conductor width
);

%% ==================== CHECK MODEL PATH ====================

fprintf('╔════════════════════════════════════════════════════════╗\n');
fprintf('║         COMSOL QUICK TEST - DEFAULT PARAMETERS        ║\n');
fprintf('╚════════════════════════════════════════════════════════╝\n\n');

if ~exist(MODEL_PATH, 'file')
    fprintf('ERROR: Model file not found:\n  %s\n\n', MODEL_PATH);
    fprintf('Please update MODEL_PATH in this script.\n');
    return;
end

fprintf('Model: %s\n\n', MODEL_PATH);

%% ==================== DISPLAY DEFAULT PARAMETERS ====================

fprintf('Default parameters (from COMSOL_parameters.txt):\n');
fprintf('┌────────────────┬────────────────┬────────────────────────┐\n');
fprintf('│ Parameter      │ Value          │ Description            │\n');
fprintf('├────────────────┼────────────────┼────────────────────────┤\n');
fprintf('│ I0             │ %.5f [A]     │ TEC Current            │\n', DEFAULT_PARAMS.I0);
fprintf('│ q              │ %.0e [W/m²] │ Heat Flux              │\n', DEFAULT_PARAMS.q);
fprintf('│ LL_theta       │ %d [deg]        │ Wedge Angle            │\n', DEFAULT_PARAMS.LL_theta);
fprintf('│ LL_t_TEC       │ %d [um]        │ TEC Thickness          │\n', DEFAULT_PARAMS.LL_t_TEC);
fprintf('│ LL_k_r         │ %.1f            │ Radial Length Ratio    │\n', DEFAULT_PARAMS.LL_k_r);
fprintf('│ LL_r_chip      │ %d [mm]        │ Chip Radius            │\n', DEFAULT_PARAMS.LL_r_chip);
fprintf('└────────────────┴────────────────┴────────────────────────┘\n\n');

%% ==================== CONNECT TO COMSOL ====================

fprintf('Connecting to COMSOL server...\n');

try
    % Try to connect
    if ~exist('mphstart', 'file')
        fprintf('\nERROR: COMSOL LiveLink not found in MATLAB path.\n');
        fprintf('Add LiveLink to path:\n');
        fprintf('  addpath(''C:\\Program Files\\COMSOL\\COMSOL63\\Multiphysics\\mli'')\n\n');
        fprintf('Or start MATLAB through COMSOL:\n');
        fprintf('  comsol mphserver matlab\n');
        return;
    end
    
    mphstart(2036);
    import com.comsol.model.*
    import com.comsol.model.util.*
    
    fprintf('✓ Connected to COMSOL server.\n\n');
    
catch ME
    fprintf('\nERROR: Could not connect to COMSOL server.\n');
    fprintf('  %s\n\n', ME.message);
    fprintf('Make sure COMSOL server is running:\n');
    fprintf('  comsolmphserver -port 2036\n');
    return;
end

%% ==================== LOAD MODEL ====================

fprintf('Loading model...\n');

try
    model = mphload(MODEL_PATH);
    fprintf('✓ Model loaded successfully.\n\n');
catch ME
    fprintf('\nERROR: Could not load model.\n');
    fprintf('  %s\n\n', ME.message);
    return;
end

%% ==================== DISPLAY CURRENT MODEL PARAMETERS ====================

fprintf('Current model parameters:\n');

try
    % Get parameter names and values
    paramNames = model.param.varnames();
    
    fprintf('┌────────────────┬────────────────────────┐\n');
    fprintf('│ Parameter      │ Value                  │\n');
    fprintf('├────────────────┼────────────────────────┤\n');
    
    for i = 1:length(paramNames)
        name = char(paramNames(i));
        value = char(model.param.get(name));
        fprintf('│ %-14s │ %-22s │\n', name, value);
    end
    fprintf('└────────────────┴────────────────────────┘\n\n');
    
catch ME
    fprintf('Warning: Could not read parameters: %s\n\n', ME.message);
end

%% ==================== RUN SIMULATION ====================

fprintf('Running simulation with default parameters...\n');
fprintf('(This may take a few minutes)\n\n');

simStart = tic;

try
    model.study('std1').run();
    simTime = toc(simStart);
    fprintf('✓ Simulation completed in %.1f seconds.\n\n', simTime);
catch ME
    fprintf('\nERROR: Simulation failed.\n');
    fprintf('  %s\n\n', ME.message);
    return;
end

%% ==================== EXTRACT RESULTS ====================

fprintf('Extracting results...\n\n');

try
    % Try to get maximum temperature
    % This depends on your model's derived values/datasets
    % Common approaches:
    
    % Method 1: Point evaluation at center
    T_center = mphinterp(model, 'T', 'coord', [0; 0; 0], 'unit', 'degC');
    fprintf('Temperature at center: %.2f °C\n', T_center);
    
    % Method 2: Global maximum (if available)
    try
        T_max = mphglobal(model, 'maxop1(T)', 'unit', 'degC');
        fprintf('Maximum temperature: %.2f °C\n', T_max);
    catch
        fprintf('(Could not extract T_max via maxop1 - check derived values)\n');
    end
    
    % Method 3: Dataset min/max
    try
        T_data = mpheval(model, 'T', 'dataset', 'dset1');
        T_min_K = min(T_data.d1);
        T_max_K = max(T_data.d1);
        fprintf('Temperature range: %.2f - %.2f K (%.2f - %.2f °C)\n', ...
            T_min_K, T_max_K, T_min_K-273.15, T_max_K-273.15);
    catch
        fprintf('(Could not extract temperature range from dset1)\n');
    end
    
catch ME
    fprintf('Warning: Result extraction had issues: %s\n', ME.message);
    fprintf('The simulation completed, but you may need to check results manually.\n');
end

%% ==================== SUMMARY ====================

fprintf('\n');
fprintf('════════════════════════════════════════════════════════\n');
fprintf('TEST COMPLETE\n');
fprintf('════════════════════════════════════════════════════════\n');
fprintf('Model loaded: ✓\n');
fprintf('Simulation ran: ✓\n');
fprintf('Time: %.1f seconds\n', simTime);
fprintf('\nCOMSOL LiveLink is working!\n');
fprintf('You can now run the overnight batch simulations.\n');
fprintf('════════════════════════════════════════════════════════\n');

%% ==================== CLEANUP ====================

% Optionally save results
% mphsave(model, 'output/comsol_results/test_result.mph');

% Clear model from memory (optional)
% ModelUtil.remove(model.tag);
