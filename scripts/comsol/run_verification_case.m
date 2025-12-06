%% RUN SINGLE VERIFICATION TEST CASE
% This script runs a single test case for model verification
% Designed for command-line execution with live logging
%
% Usage:
%   1. Start COMSOL server: comsolmphserver -port 2036
%   2. Set CASE_ID below (1-5)
%   3. Run this script
%   4. Restart server between test cases
%
% Reference: Thermal_Network_For_Radial_TEC.tex

clear; clc;

%% ============ SELECT TEST CASE ============
CASE_ID = 1;  % CHANGE THIS (1-5) FOR EACH RUN

%% ============ CONFIGURATION ============
% Get the script's directory and derive project root
script_path = mfilename('fullpath');
if isempty(script_path)
    % Running from command line with -batch, use pwd
    PROJECT_ROOT = pwd;
else
    % Running from MATLAB editor
    [script_dir, ~, ~] = fileparts(script_path);
    PROJECT_ROOT = fileparts(fileparts(script_dir));  % scripts/comsol -> scripts -> root
end

% Add paths using absolute paths
addpath(genpath(fullfile(PROJECT_ROOT, 'src')));
addpath('F:\EngineeringSoftware\COMSOL\COMSOL63\Multiphysics\mli');

COMSOL_PORT = 2036;
COMSOL_MODEL_PATH = 'E:\Semester 7\ME4311 - MicroNano Electro Mechanical Systems and Nanotechnology\Project\COMSOL\Test_Wedge\asym.mph';

% Save model after simulation for inspection
SAVE_MODEL = true;  % Set to true to save .mph file for review

OUTPUT_DIR = fullfile('output', 'model_verification');
if ~exist(OUTPUT_DIR, 'dir')
    mkdir(OUTPUT_DIR);
end

%% ============ DEFINE ALL TEST CASES ============
% Using realistic heat fluxes (~200 kW/m² typical for chip hotspots)

% Test Case 1: Baseline - Moderate heat flux, low current
tc(1).name = 'Baseline_Moderate_Flux';
tc(1).desc = 'Baseline validation with q=100 kW/m^2';
tc(1).M = 12; tc(1).t_TEC_um = 150; tc(1).t_chip_um = 50; tc(1).t_SOI_um = 20;
tc(1).R_cyl_um = 1000; tc(1).f_L = 1.15; tc(1).fill_factor = 0.9; tc(1).w_is_um = 50;
tc(1).I_A = 0.10; tc(1).q_Wm2 = 1e5;  % 100 kW/m²
tc(1).ic_w_r = 0.1; tc(1).ic_t_r = 0.6; tc(1).ic_angle_r = 0.5;
tc(1).oc_w_r = 0.1; tc(1).oc_t_r = 0.6; tc(1).oc_angle_r = 0.5;

% Test Case 2: High Heat Flux - Typical hotspot
tc(2).name = 'High_Heat_Flux';
tc(2).desc = 'High heat flux validation with q=200 kW/m^2';
tc(2).M = 12; tc(2).t_TEC_um = 150; tc(2).t_chip_um = 50; tc(2).t_SOI_um = 20;
tc(2).R_cyl_um = 1000; tc(2).f_L = 1.15; tc(2).fill_factor = 0.9; tc(2).w_is_um = 50;
tc(2).I_A = 0.15; tc(2).q_Wm2 = 2e5;  % 200 kW/m²
tc(2).ic_w_r = 0.1; tc(2).ic_t_r = 0.6; tc(2).ic_angle_r = 0.5;
tc(2).oc_w_r = 0.1; tc(2).oc_t_r = 0.6; tc(2).oc_angle_r = 0.5;

% Test Case 3: Very High Heat Flux - Stress test
tc(3).name = 'Very_High_Heat_Flux';
tc(3).desc = 'Very high heat flux validation with q=300 kW/m^2';
tc(3).M = 12; tc(3).t_TEC_um = 200; tc(3).t_chip_um = 50; tc(3).t_SOI_um = 20;
tc(3).R_cyl_um = 1000; tc(3).f_L = 1.15; tc(3).fill_factor = 0.9; tc(3).w_is_um = 50;
tc(3).I_A = 0.15; tc(3).q_Wm2 = 3e5;  % 300 kW/m²
tc(3).ic_w_r = 0.1; tc(3).ic_t_r = 0.6; tc(3).ic_angle_r = 0.5;
tc(3).oc_w_r = 0.1; tc(3).oc_t_r = 0.6; tc(3).oc_angle_r = 0.5;

% Test Case 4: Length Ratio Test - Different geometry
tc(4).name = 'Length_Ratio_Test';
tc(4).desc = 'Validates L_1 calculation with f_L=1.0 (equal lengths)';
tc(4).M = 12; tc(4).t_TEC_um = 150; tc(4).t_chip_um = 50; tc(4).t_SOI_um = 20;
tc(4).R_cyl_um = 1000; tc(4).f_L = 1.0; tc(4).fill_factor = 0.9; tc(4).w_is_um = 50;
tc(4).I_A = 0.12; tc(4).q_Wm2 = 2e5;  % 200 kW/m²
tc(4).ic_w_r = 0.1; tc(4).ic_t_r = 0.6; tc(4).ic_angle_r = 0.5;
tc(4).oc_w_r = 0.1; tc(4).oc_t_r = 0.6; tc(4).oc_angle_r = 0.5;

% Test Case 5: Wide Wedge Angle - Different number of wedges
tc(5).name = 'Wide_Wedge_Angle';
tc(5).desc = 'Validates angular scaling with theta=45 deg (M=8 wedges)';
tc(5).M = 8; tc(5).t_TEC_um = 150; tc(5).t_chip_um = 50; tc(5).t_SOI_um = 20;
tc(5).R_cyl_um = 1000; tc(5).f_L = 1.15; tc(5).fill_factor = 0.9; tc(5).w_is_um = 50;
tc(5).I_A = 0.12; tc(5).q_Wm2 = 2e5;  % 200 kW/m²
tc(5).ic_w_r = 0.1; tc(5).ic_t_r = 0.6; tc(5).ic_angle_r = 0.5;
tc(5).oc_w_r = 0.1; tc(5).oc_t_r = 0.6; tc(5).oc_angle_r = 0.5;

% Calculate derived parameters
for i = 1:5
    tc(i).theta_deg = 360 / tc(i).M;
    
    % Calculate L_1 using Eq. 238 from paper
    N_stages = 3;
    w_chip = 10000e-6;
    r_base = w_chip / sqrt(2);
    r_cyl = tc(i).R_cyl_um * 1e-6;
    w_is = tc(i).w_is_um * 1e-6;
    f_L = tc(i).f_L;
    
    L_total_active = r_base - r_cyl - (N_stages + 1) * w_is;
    if f_L == 1
        L_1 = L_total_active / N_stages;
    else
        L_1 = L_total_active * (1 - f_L) / (1 - f_L^N_stages);
    end
    tc(i).L_1_um = L_1 * 1e6;
end

%% ============ GET SELECTED TEST CASE ============
if CASE_ID < 1 || CASE_ID > 5
    error('CASE_ID must be between 1 and 5');
end
test = tc(CASE_ID);

%% ============ PRINT HEADER ============
fprintf('\n');
fprintf('========================================================================\n');
fprintf('     MODEL VERIFICATION - TEST CASE %d: %s\n', CASE_ID, test.name);
fprintf('========================================================================\n');
fprintf('     %s\n', test.desc);
fprintf('========================================================================\n\n');

fprintf('Parameters:\n');
fprintf('  theta     = %.1f deg (M = %d wedges)\n', test.theta_deg, test.M);
fprintf('  t_TEC     = %d um\n', test.t_TEC_um);
fprintf('  f_L       = %.2f (length ratio)\n', test.f_L);
fprintf('  L_1       = %.2f um (calculated)\n', test.L_1_um);
fprintf('  R_cyl     = %d um\n', test.R_cyl_um);
fprintf('  I         = %.0f mA\n', test.I_A * 1000);
fprintf('  q         = %.3g W/m^2 (%.0f kW/m^2)\n', test.q_Wm2, test.q_Wm2/1000);
fprintf('  fill_fact = %.2f\n', test.fill_factor);
fprintf('\n');

%% ============ RUN MATLAB MODEL ============
fprintf('Running MATLAB compact thermal model...\n');

try
    config = struct();
    config.geometry.N_stages = 3;
    config.geometry.M_wedges = test.M;
    config.geometry.wedge_angle_deg = test.theta_deg;
    config.geometry.thickness_um = test.t_TEC_um;
    config.geometry.radial_expansion_factor = test.f_L;
    config.geometry.w_chip_um = 10000;
    config.geometry.R_cyl_um = test.R_cyl_um;
    config.geometry.t_chip_um = test.t_chip_um;
    config.geometry.t_ins_um = test.t_SOI_um;
    config.geometry.interconnect_ratio = test.ic_w_r;
    config.geometry.outerconnect_ratio = test.oc_w_r;
    config.geometry.insulation_width_um = test.w_is_um;
    config.geometry.interconnect_angle_ratio = test.ic_angle_r;
    config.geometry.outerconnect_angle_ratio = test.oc_angle_r;
    config.geometry.interconnect_thickness_ratio = test.ic_t_r;
    config.geometry.outerconnect_thickness_ratio = test.oc_t_r;
    
    % Azimuthal gap from fill factor
    theta_rad = deg2rad(test.theta_deg);
    w_chip = 10000e-6;
    r_avg = w_chip / sqrt(2) / 2;
    arc_length = r_avg * theta_rad;
    w_az = (1 - test.fill_factor) * arc_length;
    config.geometry.azimuthal_gap_um = w_az * 1e6;
    
    config.operating_conditions.I_current_A = test.I_A;
    config.boundary_conditions.q_flux_W_m2 = test.q_Wm2;
    config.boundary_conditions.T_water_K = 300;
    config.boundary_conditions.h_conv_W_m2K = 1e6;
    
    config.materials.Bi2Te3 = struct('k', 1.2, 'rho', 1e-5, 'S', 0.0002);
    config.materials.Cu = struct('k', 400, 'rho', 1.7e-8);
    config.materials.Si = struct('k', 150, 'rho', 0.01);
    config.materials.AlN = struct('k', 170, 'rho', 1e10);
    config.materials.SiO2 = struct('k', 1.4, 'rho', 1e14);
    config.materials.Al2O3 = struct('k', 30, 'rho', 1e12);
    
    materials = MaterialProperties(config);
    geometry = TECGeometry(config);
    network = ThermalNetwork(geometry, materials, config);
    
    N_stages = 3;
    dim = 2*N_stages + 1;
    T = ones(dim, 1) * 300;
    
    for iter = 1:100
        T_old = T;
        [T, Q_out, Q_in] = network.solve(T);
        if max(abs(T - T_old)) < 1e-6
            break;
        end
    end
    
    matlab_T_center = T(1) - 273.15;
    matlab_T_max = max(T) - 273.15;
    
    fprintf('  MATLAB T_center = %.2f C\n', matlab_T_center);
    fprintf('  MATLAB T_max    = %.2f C\n', matlab_T_max);
    fprintf('  Q_out = %.4f W, Q_in = %.4f W\n', Q_out, Q_in);
    
catch ME
    fprintf('  ERROR: %s\n', ME.message);
    matlab_T_center = NaN;
    matlab_T_max = NaN;
end

fprintf('\n');

%% ============ CONNECT TO COMSOL ============
% Check if already connected (e.g., from previous run or test_load_model)
already_connected = false;
try
    import com.comsol.model.*
    import com.comsol.model.util.*
    ModelUtil.tags();  % This will throw if not connected
    already_connected = true;
    fprintf('Already connected to COMSOL server.\n\n');
catch
    % Not connected yet, need to connect
    fprintf('Connecting to COMSOL server on port %d...\n', COMSOL_PORT);
    try
        mphstart(COMSOL_PORT);
        import com.comsol.model.*
        import com.comsol.model.util.*
        fprintf('  Connected successfully!\n\n');
    catch ME
        fprintf('\n');
        fprintf('ERROR: Could not connect to COMSOL server.\n');
        fprintf('  %s\n\n', ME.message);
        fprintf('Please start the server with:\n');
        fprintf('  "F:\\EngineeringSoftware\\COMSOL\\COMSOL63\\Multiphysics\\bin\\win64\\comsolmphserver.exe" -port 2036\n\n');
        return;
    end
end

%% ============ LOAD MODEL ============
fprintf('Loading COMSOL model...\n');
fprintf('  %s\n', COMSOL_MODEL_PATH);

% Check if file exists first
if ~exist(COMSOL_MODEL_PATH, 'file')
    fprintf('ERROR: Model file does not exist!\n');
    fprintf('  Check the path: %s\n', COMSOL_MODEL_PATH);
    return;
end

fprintf('  Calling mphload (this may take a while with LiveLink)...\n');

% Flush output before potentially crashing call
drawnow;
pause(0.1);

model = [];
load_success = false;
load_error_msg = '';

try
    model = mphload(COMSOL_MODEL_PATH);
    load_success = true;
    fprintf('  Model loaded successfully!\n\n');
catch ME
    load_error_msg = ME.message;
    fprintf('\n');
    fprintf('========================================================================\n');
    fprintf('ERROR: Could not load COMSOL model\n');
    fprintf('========================================================================\n');
    fprintf('  Error ID: %s\n', ME.identifier);
    fprintf('  Message: %s\n', ME.message);
    if ~isempty(ME.cause)
        for ci = 1:length(ME.cause)
            fprintf('  Cause %d: %s\n', ci, ME.cause{ci}.message);
        end
    end
    fprintf('========================================================================\n');
    fprintf('\nPossible causes:\n');
    fprintf('  1. SolidWorks is not running or needs user interaction\n');
    fprintf('  2. The SolidWorks file linked to asym.mph is not open\n');
    fprintf('  3. COMSOL-SolidWorks connection timed out\n');
    fprintf('  4. License issues with LiveLink for SolidWorks\n\n');
end

if ~load_success
    fprintf('Model load failed. Exiting.\n');
    return;
end

%% ============ SET PARAMETERS ============
fprintf('Setting COMSOL parameters...\n');

try
    % Note: Parameters without units are in micrometers for COMSOL
    % The COMSOL model expects values in um, so we pass numerical values
    
    fprintf('  Setting parameters...\n');
    
    model.param.set('LL_k_r', sprintf('%g', test.f_L));
    fprintf('    LL_k_r      = %g\n', test.f_L);
    
    model.param.set('LL_L_1', sprintf('%g', test.L_1_um));
    fprintf('    LL_L_1      = %.2f (um)\n', test.L_1_um);
    
    model.param.set('LL_R_cyl', sprintf('%g', test.R_cyl_um));
    fprintf('    LL_R_cyl    = %g (um)\n', test.R_cyl_um);
    
    model.param.set('LL_t_chip', sprintf('%g', test.t_chip_um));
    fprintf('    LL_t_chip   = %g (um)\n', test.t_chip_um);
    
    model.param.set('LL_t_SOI', sprintf('%g', test.t_SOI_um));
    fprintf('    LL_t_SOI    = %g (um)\n', test.t_SOI_um);
    
    model.param.set('LL_t_TEC', sprintf('%g', test.t_TEC_um));
    fprintf('    LL_t_TEC    = %g (um)\n', test.t_TEC_um);
    
    model.param.set('LL_theta', sprintf('%g', test.theta_deg));
    fprintf('    LL_theta    = %g (deg)\n', test.theta_deg);
    
    model.param.set('LL_w_is', sprintf('%g', test.w_is_um));
    fprintf('    LL_w_is     = %g (um)\n', test.w_is_um);
    
    % Heat flux and current - these need units
    model.param.set('q_i', sprintf('%g[W/m^2]', test.q_Wm2));
    fprintf('    q_i         = %g W/m^2\n', test.q_Wm2);
    
    model.param.set('I_0', sprintf('%g[A]', test.I_A));
    fprintf('    I_0         = %g A\n', test.I_A);
    
    model.param.set('LL_fill_factor', sprintf('%g', test.fill_factor));
    fprintf('    fill_factor = %g\n', test.fill_factor);
    
    % Interconnect/outerconnect parameters (ratios - dimensionless)
    model.param.set('LL_ic_angle_r', sprintf('%g', test.ic_angle_r));
    model.param.set('LL_ic_t_r', sprintf('%g', test.ic_t_r));
    model.param.set('LL_ic_w_r', sprintf('%g', test.ic_w_r));
    model.param.set('LL_oc_angle_r', sprintf('%g', test.oc_angle_r));
    model.param.set('LL_oc_t_r', sprintf('%g', test.oc_t_r));
    model.param.set('LL_oc_w_r', sprintf('%g', test.oc_w_r));
    
    fprintf('  Parameters set.\n\n');
    
    % VERIFY: Read back critical parameters to confirm they were set
    fprintf('  Verifying parameters (readback)...\n');
    t_TEC_check = char(model.param.get('LL_t_TEC'));
    L_1_check = char(model.param.get('LL_L_1'));
    k_r_check = char(model.param.get('LL_k_r'));
    fprintf('    LL_t_TEC (readback) = %s\n', t_TEC_check);
    fprintf('    LL_L_1 (readback)   = %s\n', L_1_check);
    fprintf('    LL_k_r (readback)   = %s\n', k_r_check);
    fprintf('  Parameters verified.\n\n');
    
catch ME
    fprintf('  Warning: Could not set some parameters: %s\n\n', ME.message);
end

%% ============ REBUILD GEOMETRY & MESH ============
% This is CRITICAL for SolidWorks LiveLink to update the geometry
fprintf('Step 1: Rebuilding geometry from SolidWorks...\n');
fprintf('  (Watch SolidWorks - it should update to t_TEC = %d um)\n', test.t_TEC_um);

geom_start = tic;
try
    % Build geometry (this triggers SolidWorks update)
    model.component('comp1').geom('geom1').run();
    geom_time = toc(geom_start);
    fprintf('  Geometry rebuilt in %.1f seconds.\n', geom_time);
catch ME
    geom_time = toc(geom_start);
    fprintf('  Geometry rebuild FAILED after %.1f s: %s\n', geom_time, ME.message);
end

% VERIFY: Check parameters again after geometry rebuild
fprintf('\n  Verifying parameters AFTER geometry rebuild...\n');
try
    t_TEC_after = char(model.param.get('LL_t_TEC'));
    fprintf('    LL_t_TEC (after geom) = %s\n', t_TEC_after);
    if ~strcmp(t_TEC_after, sprintf('%g', test.t_TEC_um))
        fprintf('  WARNING: Parameter may have been reset!\n');
        fprintf('    Expected: %g, Got: %s\n', test.t_TEC_um, t_TEC_after);
    end
catch
    fprintf('    Could not verify parameters.\n');
end

fprintf('\nStep 2: Regenerating mesh...\n');
mesh_start = tic;
try
    model.component('comp1').mesh('mesh1').run();
    mesh_time = toc(mesh_start);
    fprintf('  Mesh regenerated in %.1f seconds.\n', mesh_time);
catch ME
    mesh_time = toc(mesh_start);
    fprintf('  Mesh regeneration FAILED after %.1f s: %s\n', mesh_time, ME.message);
end

% VERIFY: Check parameters again after meshing
fprintf('\n  Verifying parameters AFTER meshing...\n');
try
    t_TEC_after_mesh = char(model.param.get('LL_t_TEC'));
    fprintf('    LL_t_TEC (after mesh) = %s\n', t_TEC_after_mesh);
catch
    fprintf('    Could not verify parameters.\n');
end

%% ============ CHECK AND DISABLE PARAMETRIC SWEEPS ============
fprintf('\nStep 3: Checking for parametric sweeps...\n');

try
    % List all studies
    study_tags = model.study.tags();
    fprintf('  Found %d study/studies:\n', length(study_tags));
    for i = 1:length(study_tags)
        study_tag = char(study_tags(i));
        fprintf('    - %s\n', study_tag);
    end
    
    % Check study 'std1' for parametric sweep features
    std1_features = model.study('std1').feature.tags();
    fprintf('\n  Study std1 has %d feature(s):\n', length(std1_features));
    for i = 1:length(std1_features)
        feat_tag = char(std1_features(i));
        try
            feat_type = char(model.study('std1').feature(feat_tag).getType());
            is_active = model.study('std1').feature(feat_tag).isActive();
            fprintf('    - %s (type: %s, active: %d)\n', feat_tag, feat_type, is_active);
            
            % Disable any parametric sweep
            if contains(lower(feat_type), 'param') || contains(lower(feat_tag), 'param')
                fprintf('      >> DISABLING parametric sweep: %s\n', feat_tag);
                model.study('std1').feature(feat_tag).active(false);
            end
        catch
            fprintf('    - %s (could not get type)\n', feat_tag);
        end
    end
    
    % Also check for batch/parametric at study level
    try
        % Try to disable parametric sweep if it exists
        model.study('std1').feature('param').active(false);
        fprintf('  Disabled parametric sweep "param"\n');
    catch
        % No param feature or already disabled
    end
    
    try
        model.study('std1').feature('param1').active(false);
        fprintf('  Disabled parametric sweep "param1"\n');
    catch
    end
    
catch ME
    fprintf('  Warning checking studies: %s\n', ME.message);
end

%% ============ RUN SIMULATION ============
fprintf('\nStep 4: Running COMSOL simulation (stationary study only)...\n');
fprintf('  (This may take several minutes depending on mesh size)\n');

sim_start = tic;
sim_success = false;

try
    % Run only the stationary study, not any parametric sweeps
    model.study('std1').run();
    sim_time = toc(sim_start);
    fprintf('  Simulation completed in %.1f seconds.\n\n', sim_time);
    sim_success = true;
catch ME
    sim_time = toc(sim_start);
    fprintf('\n');
    fprintf('========================================================================\n');
    fprintf('ERROR: COMSOL Simulation Failed\n');
    fprintf('========================================================================\n');
    fprintf('  Time elapsed: %.1f seconds\n', sim_time);
    fprintf('  Error ID: %s\n', ME.identifier);
    fprintf('  Message: %s\n', ME.message);
    if ~isempty(ME.cause)
        for k = 1:length(ME.cause)
            fprintf('  Cause %d: %s\n', k, ME.cause{k}.message);
        end
    end
    fprintf('========================================================================\n\n');
    
    % Save model even on failure so user can debug
    fprintf('Saving model for debugging (simulation failed)...\n');
    debug_model_path = fullfile(OUTPUT_DIR, sprintf('DEBUG_case%d_%s.mph', CASE_ID, test.name));
    try
        if ~exist(OUTPUT_DIR, 'dir')
            mkdir(OUTPUT_DIR);
        end
        mphsave(model, debug_model_path);
        fprintf('  Debug model saved: %s\n', debug_model_path);
        fprintf('  Open this file in COMSOL GUI to investigate the error.\n\n');
    catch ME2
        fprintf('  Could not save debug model: %s\n\n', ME2.message);
    end
end

%% ============ EXTRACT RESULTS ============
comsol_T_center = NaN;
comsol_T_max = NaN;

if sim_success
    fprintf('Extracting results...\n');
    
    % Try probe first
    try
        T_probe = mphglobal(model, 'comp1.ppb1', 'dataset', 'dset1');
        comsol_T_center = T_probe - 273.15;
        fprintf('  T_center (probe): %.2f C\n', comsol_T_center);
    catch
        try
            T_data = mphinterp(model, 'T', 'coord', [0; 0; 0], 'dataset', 'dset1');
            comsol_T_center = T_data - 273.15;
            fprintf('  T_center (interp): %.2f C\n', comsol_T_center);
        catch ME2
            fprintf('  Warning: Could not get center temp: %s\n', ME2.message);
        end
    end
    
    % Get max temperature
    try
        T_all = mpheval(model, 'T', 'dataset', 'dset1');
        comsol_T_max = max(T_all.d1) - 273.15;
        fprintf('  T_max: %.2f C\n', comsol_T_max);
    catch ME3
        fprintf('  Warning: Could not get max temp: %s\n', ME3.message);
    end
end

%% ============ SAVE RESULTS ============
fprintf('\nSaving results...\n');

result = struct();
result.case_id = CASE_ID;
result.test = test;
result.matlab_T_center = matlab_T_center;
result.matlab_T_max = matlab_T_max;
result.comsol_T_center = comsol_T_center;
result.comsol_T_max = comsol_T_max;
result.sim_time = sim_time;
result.sim_success = sim_success;
result.timestamp = datestr(now);

filename = sprintf('verification_case%d_%s.mat', CASE_ID, test.name);
save(fullfile(OUTPUT_DIR, filename), 'result');
fprintf('  Saved: %s\n', fullfile(OUTPUT_DIR, filename));

% Append to CSV
csv_file = fullfile(OUTPUT_DIR, 'verification_results.csv');
if exist(csv_file, 'file')
    fid = fopen(csv_file, 'a');
else
    fid = fopen(csv_file, 'w');
    fprintf(fid, 'CaseID,Name,theta_deg,t_TEC_um,f_L,L_1_um,I_mA,q_Wm2,MATLAB_Tcenter,MATLAB_Tmax,COMSOL_Tcenter,COMSOL_Tmax,Error_Tmax,Error_pct,SimTime,Timestamp\n');
end

error_abs = comsol_T_max - matlab_T_max;
error_pct = 100 * error_abs / matlab_T_max;

fprintf(fid, '%d,%s,%.1f,%d,%.2f,%.2f,%.0f,%.0f,%.2f,%.2f,%.2f,%.2f,%.2f,%.1f,%.1f,%s\n', ...
    CASE_ID, test.name, test.theta_deg, test.t_TEC_um, test.f_L, test.L_1_um, ...
    test.I_A*1000, test.q_Wm2, matlab_T_center, matlab_T_max, ...
    comsol_T_center, comsol_T_max, error_abs, error_pct, sim_time, result.timestamp);
fclose(fid);
fprintf('  Appended to: %s\n', csv_file);

%% ============ SUMMARY ============
fprintf('\n');
fprintf('========================================================================\n');
fprintf('                         RESULTS SUMMARY\n');
fprintf('========================================================================\n');
fprintf('  Test Case %d: %s\n', CASE_ID, test.name);
fprintf('------------------------------------------------------------------------\n');
fprintf('                    MATLAB CTM        COMSOL FEM        Error\n');
fprintf('------------------------------------------------------------------------\n');
fprintf('  T_center:      %8.2f C        %8.2f C        %+6.2f C\n', ...
    matlab_T_center, comsol_T_center, comsol_T_center - matlab_T_center);
fprintf('  T_max:         %8.2f C        %8.2f C        %+6.2f C (%+.1f%%)\n', ...
    matlab_T_max, comsol_T_max, error_abs, error_pct);
fprintf('------------------------------------------------------------------------\n');
fprintf('  Simulation time: %.1f seconds\n', sim_time);
fprintf('========================================================================\n\n');

%% ============ NEXT STEPS ============
if CASE_ID < 5
    fprintf('NEXT STEPS:\n');
    fprintf('  1. The COMSOL server will disconnect.\n');
    fprintf('  2. Restart the server:\n');
    fprintf('     "F:\\EngineeringSoftware\\COMSOL\\COMSOL63\\Multiphysics\\bin\\win64\\comsolmphserver.exe" -port 2036\n');
    fprintf('  3. Change CASE_ID to %d in this script\n', CASE_ID + 1);
    fprintf('  4. Run this script again\n\n');
else
    fprintf('ALL TEST CASES COMPLETE!\n\n');
    fprintf('To generate comparison plots, run:\n');
    fprintf('  plot_verification_results\n\n');
end

fprintf('Results saved to: %s\n\n', OUTPUT_DIR);
