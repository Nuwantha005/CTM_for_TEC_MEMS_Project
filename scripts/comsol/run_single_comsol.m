% run_single_comsol.m
% Run a single COMSOL simulation and extract results
% Designed to work with server that stops after each run
%
% Usage: Set TEST_CASE_ID below and run. Restart COMSOL server between runs.

clear; clc;
addpath(genpath('src'));
addpath('F:\EngineeringSoftware\COMSOL\COMSOL63\Multiphysics\mli');

%% ============ CONFIGURATION ============
% Change this to run different test cases (1-5)
TEST_CASE_ID = 5;

COMSOL_PORT = 2036;
COMSOL_MODEL_PATH = 'E:\Semester 7\ME4311 - MicroNano Electro Mechanical Systems and Nanotechnology\Project\COMSOL\MATLAB_API\template_3_stage.mph';

% Output directory (shared across runs)
OUTPUT_DIR = fullfile('output', 'comsol_single_runs');
if ~exist(OUTPUT_DIR, 'dir')
    mkdir(OUTPUT_DIR);
end

%% ============ DEFINE TEST CASES ============
% Based on MATLAB optimization results for 3-stage TEC

test_cases = struct();

% Test Case 1: Baseline optimal (low flux, moderate current)
test_cases(1).name = 'Baseline_Optimal';
test_cases(1).I_A = 0.075;      % 75 mA
test_cases(1).t_TEC_um = 150;   % 150 um
test_cases(1).k_r = 1.15;
test_cases(1).q_Wm2 = 500;      % 500 W/m²
test_cases(1).theta_deg = 30;

% Test Case 2: Higher current
test_cases(2).name = 'High_Current';
test_cases(2).I_A = 0.100;      % 100 mA
test_cases(2).t_TEC_um = 150;
test_cases(2).k_r = 1.15;
test_cases(2).q_Wm2 = 500;
test_cases(2).theta_deg = 30;

% Test Case 3: Thicker TEC
test_cases(3).name = 'Thick_TEC';
test_cases(3).I_A = 0.100;
test_cases(3).t_TEC_um = 250;   % 250 um
test_cases(3).k_r = 1.0;
test_cases(3).q_Wm2 = 500;
test_cases(3).theta_deg = 30;

% Test Case 4: Higher heat flux (1000 W/m²)
test_cases(4).name = 'Medium_Flux';
test_cases(4).I_A = 0.100;
test_cases(4).t_TEC_um = 200;
test_cases(4).k_r = 1.15;
test_cases(4).q_Wm2 = 1000;     % 1000 W/m²
test_cases(4).theta_deg = 30;

% Test Case 5: Maximum feasible flux
test_cases(5).name = 'High_Flux';
test_cases(5).I_A = 0.100;
test_cases(5).t_TEC_um = 250;
test_cases(5).k_r = 1.0;
test_cases(5).q_Wm2 = 1500;     % 1500 W/m²
test_cases(5).theta_deg = 30;

%% ============ GET SELECTED TEST CASE ============
tc = test_cases(TEST_CASE_ID);

fprintf('╔══════════════════════════════════════════════════════════╗\n');
fprintf('║         SINGLE COMSOL SIMULATION RUN                     ║\n');
fprintf('╚══════════════════════════════════════════════════════════╝\n\n');

fprintf('Test Case %d: %s\n', TEST_CASE_ID, tc.name);
fprintf('─────────────────────────────────────────\n');
fprintf('  Current (I):     %.0f mA\n', tc.I_A * 1000);
fprintf('  TEC Thickness:   %.0f µm\n', tc.t_TEC_um);
fprintf('  Radial Factor:   %.2f\n', tc.k_r);
fprintf('  Heat Flux:       %.0f W/m²\n', tc.q_Wm2);
fprintf('  Wedge Angle:     %.0f°\n', tc.theta_deg);
fprintf('─────────────────────────────────────────\n\n');

%% ============ RUN MATLAB MODEL FOR COMPARISON ============
fprintf('Running MATLAB compact model...\n');

try
    % Create config for MATLAB model
    config = struct();
    config.geometry.N_stages = 3;
    config.geometry.wedge_angle_deg = tc.theta_deg;
    config.geometry.thickness_um = tc.t_TEC_um;
    config.geometry.radial_expansion_factor = tc.k_r;
    config.geometry.fill_factor = 0.90;
    config.geometry.w_chip_um = 10000;
    config.geometry.R_cyl_um = 1000;
    config.geometry.t_chip_um = 50;
    config.geometry.interconnect_ratio = 0.15;
    config.geometry.outerconnect_ratio = 0.15;
    config.geometry.insulation_width_ratio = 0.04;
    config.geometry.interconnect_angle_ratio = 0.16;
    config.geometry.outerconnect_angle_ratio = 0.16;
    config.geometry.interconnect_thickness_ratio = 1.0;
    config.geometry.outerconnect_thickness_ratio = 1.0;
    config.geometry.tsv.R_TSV_um = 10;
    config.geometry.tsv.P_TSV_um = 20;
    config.geometry.tsv.g_rad_um = 10;
    config.geometry.tsv.t_SOI_um = 100;
    
    config.operating_conditions.I_current_A = tc.I_A;
    config.boundary_conditions.q_flux_W_m2 = tc.q_Wm2;
    config.boundary_conditions.T_water_K = 300;
    config.boundary_conditions.h_conv_W_m2K = 1e6;
    
    config.materials.Bi2Te3 = struct('k', 1.2, 'rho', 1e-5, 'S', 0.0002);
    config.materials.Cu = struct('k', 400, 'rho', 1.7e-8);
    config.materials.Si = struct('k', 150, 'rho', 0.01);
    config.materials.AlN = struct('k', 170, 'rho', 1e10);
    config.materials.SiO2 = struct('k', 1.4, 'rho', 1e14);
    config.materials.Al2O3 = struct('k', 30, 'rho', 1e12);
    
    % Run MATLAB model
    materials = MaterialProperties(config);
    geometry = TECGeometry(config);
    network = ThermalNetwork(geometry, materials, config);
    
    N = geometry.N_stages;
    dim = 2*N + 1;
    T = ones(dim, 1) * 300;
    
    for iter = 1:100
        T_old = T;
        [T, Q_out, Q_in] = network.solve(T);
        if max(abs(T - T_old)) < 1e-6
            break;
        end
    end
    
    matlab_T_max = max(T) - 273.15;
    matlab_T_center = T(1) - 273.15;
    
    fprintf('  MATLAB T_max:    %.2f °C\n', matlab_T_max);
    fprintf('  MATLAB T_center: %.2f °C\n\n', matlab_T_center);
    
catch ME
    fprintf('  MATLAB model failed: %s\n\n', ME.message);
    matlab_T_max = NaN;
    matlab_T_center = NaN;
end

%% ============ CONNECT TO COMSOL ============
fprintf('Connecting to COMSOL...\n');

% Check if already connected (e.g., when launched from "COMSOL with MATLAB")
already_connected = false;
try
    % Try to access ModelUtil - if it works, we're already connected
    import com.comsol.model.*
    import com.comsol.model.util.*
    ModelUtil.tags();  % This will throw if not connected
    already_connected = true;
    fprintf('  Already connected to COMSOL (COMSOL with MATLAB mode)\n\n');
catch
    % Not connected yet, try to connect
    try
        fprintf('  Connecting to COMSOL server on port %d...\n', COMSOL_PORT);
        mphstart(COMSOL_PORT);
        import com.comsol.model.*
        import com.comsol.model.util.*
        fprintf('  Connected successfully!\n\n');
    catch ME
        fprintf('\nERROR: Could not connect to COMSOL server.\n');
        fprintf('  %s\n\n', ME.message);
        fprintf('Options:\n');
        fprintf('  1. Start COMSOL server: comsolmphserver -port 2036\n');
        fprintf('  2. Use "COMSOL with MATLAB" shortcut\n');
        return;
    end
end

%% ============ LOAD MODEL ============
fprintf('Loading COMSOL model...\n');
fprintf('  %s\n', COMSOL_MODEL_PATH);

try
    model = mphload(COMSOL_MODEL_PATH);
    fprintf('  Model loaded successfully!\n\n');
catch ME
    fprintf('\nERROR: Could not load model.\n');
    fprintf('  %s\n\n', ME.message);
    return;
end

%% ============ SET PARAMETERS ============
fprintf('Setting parameters in COMSOL...\n');

try
    % Set parameters - adjust names to match your model
    model.param.set('I0', sprintf('%g[A]', tc.I_A));
    model.param.set('LL_t_TEC', sprintf('%g[um]', tc.t_TEC_um));
    model.param.set('LL_k_r', sprintf('%g', tc.k_r));
    model.param.set('q', sprintf('%g[W/m^2]', tc.q_Wm2));
    model.param.set('LL_theta', sprintf('%g[deg]', tc.theta_deg));
    
    fprintf('  Parameters set successfully!\n\n');
catch ME
    fprintf('  Warning: Could not set some parameters: %s\n\n', ME.message);
end

%% ============ RUN SIMULATION ============
fprintf('Running COMSOL simulation...\n');
fprintf('  (This may take a few minutes...)\n');

simStart = tic;

try
    model.study('std1').run();
    simTime = toc(simStart);
    fprintf('  Simulation completed in %.1f seconds!\n\n', simTime);
    sim_success = true;
catch ME
    simTime = toc(simStart);
    fprintf('\nERROR: Simulation failed after %.1f seconds.\n', simTime);
    fprintf('  %s\n\n', ME.message);
    sim_success = false;
end

%% ============ EXTRACT RESULTS ============
comsol_T_center = NaN;
comsol_T_max = NaN;

if sim_success
    fprintf('Extracting results...\n');
    
    % Try to get temperature from "Node 0 Temp" probe
    try
        % Method 1: Try probe evaluation
        probeResult = mphglobal(model, 'comp1.ppb1', 'dataset', 'dset1');
        comsol_T_center = probeResult - 273.15;  % Convert to Celsius
        fprintf('  Node 0 Temp (probe): %.2f °C\n', comsol_T_center);
    catch
        try
            % Method 2: Evaluate at center point
            T_data = mphinterp(model, 'T', 'coord', [0; 0; 0], 'dataset', 'dset1');
            comsol_T_center = T_data - 273.15;
            fprintf('  Center temperature: %.2f °C\n', comsol_T_center);
        catch ME2
            fprintf('  Warning: Could not extract center temp: %s\n', ME2.message);
        end
    end
    
    % Get maximum temperature
    try
        T_all = mpheval(model, 'T', 'dataset', 'dset1');
        comsol_T_max = max(T_all.d1) - 273.15;
        fprintf('  Maximum temperature: %.2f °C\n', comsol_T_max);
    catch ME3
        try
            comsol_T_max = mphglobal(model, 'maxop1(T)', 'unit', 'degC');
            fprintf('  Maximum temperature: %.2f °C\n', comsol_T_max);
        catch
            fprintf('  Warning: Could not extract max temp: %s\n', ME3.message);
        end
    end
end

%% ============ SAVE RESULTS ============
fprintf('\nSaving results...\n');

result = struct();
result.test_case_id = TEST_CASE_ID;
result.test_case_name = tc.name;
result.params = tc;
result.matlab_T_max = matlab_T_max;
result.matlab_T_center = matlab_T_center;
result.comsol_T_max = comsol_T_max;
result.comsol_T_center = comsol_T_center;
result.error_T_max = comsol_T_max - matlab_T_max;
result.error_T_center = comsol_T_center - matlab_T_center;
result.sim_time = simTime;
result.sim_success = sim_success;
result.timestamp = datestr(now, 'yyyy-mm-dd HH:MM:SS');

% Save individual result
filename = sprintf('result_case%d_%s.mat', TEST_CASE_ID, tc.name);
save(fullfile(OUTPUT_DIR, filename), 'result');
fprintf('  Saved: %s\n', filename);

% Also append to summary CSV
csvFile = fullfile(OUTPUT_DIR, 'all_results.csv');
if exist(csvFile, 'file')
    % Append to existing
    fid = fopen(csvFile, 'a');
else
    % Create new with header
    fid = fopen(csvFile, 'w');
    fprintf(fid, 'ID,Name,I_mA,t_um,k_r,q_Wm2,MATLAB_Tmax,COMSOL_Tmax,Error,SimTime,Success,Timestamp\n');
end
fprintf(fid, '%d,%s,%.0f,%.0f,%.2f,%.0f,%.2f,%.2f,%.2f,%.1f,%d,%s\n', ...
    TEST_CASE_ID, tc.name, tc.I_A*1000, tc.t_TEC_um, tc.k_r, tc.q_Wm2, ...
    matlab_T_max, comsol_T_max, result.error_T_max, simTime, sim_success, result.timestamp);
fclose(fid);
fprintf('  Appended to: all_results.csv\n');

%% ============ SUMMARY ============
fprintf('\n');
fprintf('╔══════════════════════════════════════════════════════════╗\n');
fprintf('║                    RESULTS SUMMARY                       ║\n');
fprintf('╚══════════════════════════════════════════════════════════╝\n\n');

fprintf('Test Case %d: %s\n', TEST_CASE_ID, tc.name);
fprintf('─────────────────────────────────────────\n');
fprintf('  MATLAB  T_max:   %8.2f °C\n', matlab_T_max);
fprintf('  COMSOL  T_max:   %8.2f °C\n', comsol_T_max);
fprintf('  Error:           %8.2f °C (%.1f%%)\n', result.error_T_max, ...
    100*result.error_T_max/matlab_T_max);
fprintf('─────────────────────────────────────────\n');
fprintf('  MATLAB  T_center: %7.2f °C\n', matlab_T_center);
fprintf('  COMSOL  T_center: %7.2f °C\n', comsol_T_center);
fprintf('  Error:            %7.2f °C\n', result.error_T_center);
fprintf('─────────────────────────────────────────\n');
fprintf('  Simulation time: %.1f seconds\n', simTime);
fprintf('─────────────────────────────────────────\n\n');

if TEST_CASE_ID < 5
    fprintf('NEXT STEP:\n');
    fprintf('  1. Restart COMSOL server\n');
    fprintf('  2. Change TEST_CASE_ID to %d\n', TEST_CASE_ID + 1);
    fprintf('  3. Run this script again\n\n');
else
    fprintf('All test cases complete!\n');
    fprintf('Run: plot_comsol_comparison.m to generate plots\n\n');
end

fprintf('Output saved to: %s\n', OUTPUT_DIR);
