%% VERIFY MATHEMATICAL MODEL USING COMSOL
% This script verifies the compact thermal model (CTM) against COMSOL simulations
% Reference: Thermal_Network_For_Radial_TEC.tex
%
% Test Cases:
%   1. Low current, thin TEC - validates basic thermal resistance calculations
%   2. Moderate current - validates Joule heating effects
%   3. Different length ratio (f_L) - validates geometry calculations
%   4. Higher heat flux - validates heat generation terms
%   5. Different wedge angle - validates angular scaling
%
% Requirements:
%   - COMSOL 6.3 server running on port 2036
%   - Template model at specified path
%   - MATLAB compact thermal model classes in src/
%
% Usage:
%   1. Start COMSOL server: comsolmphserver -port 2036
%   2. Run: verify_model_comsol
%
% Author: Auto-generated verification script
% Date: 2025-12-06

clear; clc; close all;

%% ============ CONFIGURATION ============
% Add paths
addpath(genpath('src'));
addpath('F:\EngineeringSoftware\COMSOL\COMSOL63\Multiphysics\mli');

% COMSOL configuration
COMSOL_PORT = 2036;
COMSOL_MODEL_PATH = 'E:\Semester 7\ME4311 - MicroNano Electro Mechanical Systems and Nanotechnology\Project\COMSOL\MATLAB_API\template_3_stage.mph';
COMSOL_SERVER_PATH = 'F:\EngineeringSoftware\COMSOL\COMSOL63\Multiphysics\bin\win64\comsolmphserver.exe';

% Output directory
timestamp = datestr(now, 'yyyy-mm-dd_HH-MM-SS');
OUTPUT_DIR = fullfile('output', 'model_verification', timestamp);
if ~exist(OUTPUT_DIR, 'dir')
    mkdir(OUTPUT_DIR);
end

fprintf('╔══════════════════════════════════════════════════════════════════╗\n');
fprintf('║     MATHEMATICAL MODEL VERIFICATION USING COMSOL SIMULATIONS    ║\n');
fprintf('╠══════════════════════════════════════════════════════════════════╣\n');
fprintf('║  Reference: Thermal_Network_For_Radial_TEC.tex                  ║\n');
fprintf('║  Output: %s                               ║\n', OUTPUT_DIR);
fprintf('╚══════════════════════════════════════════════════════════════════╝\n\n');

%% ============ DEFINE 5 TEST CASES ============
% Using small thickness values to speed up meshing
% Parameters follow paper notation from COMSOL_parameters.txt

fprintf('Defining test cases...\n\n');

test_cases = struct();

% ----- Test Case 1: Baseline - Low current, thin TEC -----
% Purpose: Validate basic thermal resistance network
test_cases(1).name = 'Thin_TEC_Low_Current';
test_cases(1).description = 'Validates basic thermal resistance calculations';
test_cases(1).M = 12;                    % Number of wedges
test_cases(1).t_TEC_um = 100;            % t_TEC: 100 um (thin for fast meshing)
test_cases(1).t_chip_um = 50;            % t_chip
test_cases(1).t_SOI_um = 20;             % t_ins (SOI layer)
test_cases(1).R_cyl_um = 1000;           % r_cyl
test_cases(1).f_L = 1.0;                 % f_L: uniform stage lengths
test_cases(1).fill_factor = 0.9;         % Fill factor
test_cases(1).w_is_um = 50;              % W_is: radial insulator width
test_cases(1).I_A = 0.05;                % I: 50 mA (low current)
test_cases(1).q_Wm2 = 500;               % q: 500 W/m²
test_cases(1).ic_w_r = 0.1;              % W_ic ratio
test_cases(1).ic_t_r = 0.6;              % t_ic ratio  
test_cases(1).ic_angle_r = 0.5;          % β_ic ratio
test_cases(1).oc_w_r = 0.1;              % W_oc ratio
test_cases(1).oc_t_r = 0.6;              % t_oc ratio
test_cases(1).oc_angle_r = 0.5;          % β_oc ratio

% ----- Test Case 2: Moderate current -----
% Purpose: Validate Joule heating effects (I²R terms)
test_cases(2).name = 'Moderate_Current';
test_cases(2).description = 'Validates Joule heating effects (I²R terms)';
test_cases(2).M = 12;
test_cases(2).t_TEC_um = 150;
test_cases(2).t_chip_um = 50;
test_cases(2).t_SOI_um = 20;
test_cases(2).R_cyl_um = 1000;
test_cases(2).f_L = 1.15;
test_cases(2).fill_factor = 0.9;
test_cases(2).w_is_um = 50;
test_cases(2).I_A = 0.10;                % I: 100 mA (moderate)
test_cases(2).q_Wm2 = 500;
test_cases(2).ic_w_r = 0.1;
test_cases(2).ic_t_r = 0.6;
test_cases(2).ic_angle_r = 0.5;
test_cases(2).oc_w_r = 0.1;
test_cases(2).oc_t_r = 0.6;
test_cases(2).oc_angle_r = 0.5;

% ----- Test Case 3: Different length ratio -----
% Purpose: Validate geometry calculations (Eq. 238 for L_1)
test_cases(3).name = 'Length_Ratio_Test';
test_cases(3).description = 'Validates L_1 calculation with f_L=1.2 (Eq. 238)';
test_cases(3).M = 12;
test_cases(3).t_TEC_um = 150;
test_cases(3).t_chip_um = 50;
test_cases(3).t_SOI_um = 20;
test_cases(3).R_cyl_um = 1000;
test_cases(3).f_L = 1.2;                 % f_L: increasing stage lengths
test_cases(3).fill_factor = 0.9;
test_cases(3).w_is_um = 50;
test_cases(3).I_A = 0.08;
test_cases(3).q_Wm2 = 500;
test_cases(3).ic_w_r = 0.1;
test_cases(3).ic_t_r = 0.6;
test_cases(3).ic_angle_r = 0.5;
test_cases(3).oc_w_r = 0.1;
test_cases(3).oc_t_r = 0.6;
test_cases(3).oc_angle_r = 0.5;

% ----- Test Case 4: Higher heat flux -----
% Purpose: Validate heat generation terms (Q_gen)
test_cases(4).name = 'High_Heat_Flux';
test_cases(4).description = 'Validates heat generation terms with q=2000 W/m²';
test_cases(4).M = 12;
test_cases(4).t_TEC_um = 200;
test_cases(4).t_chip_um = 50;
test_cases(4).t_SOI_um = 20;
test_cases(4).R_cyl_um = 1000;
test_cases(4).f_L = 1.15;
test_cases(4).fill_factor = 0.9;
test_cases(4).w_is_um = 50;
test_cases(4).I_A = 0.12;                % Higher current for high flux
test_cases(4).q_Wm2 = 2000;              % q: 2000 W/m²
test_cases(4).ic_w_r = 0.1;
test_cases(4).ic_t_r = 0.6;
test_cases(4).ic_angle_r = 0.5;
test_cases(4).oc_w_r = 0.1;
test_cases(4).oc_t_r = 0.6;
test_cases(4).oc_angle_r = 0.5;

% ----- Test Case 5: Different wedge angle -----
% Purpose: Validate angular scaling (theta dependence)
test_cases(5).name = 'Wide_Wedge_Angle';
test_cases(5).description = 'Validates angular scaling with θ=45° (M=8 wedges)';
test_cases(5).M = 8;                     % 8 wedges -> θ = 45°
test_cases(5).t_TEC_um = 150;
test_cases(5).t_chip_um = 50;
test_cases(5).t_SOI_um = 20;
test_cases(5).R_cyl_um = 1000;
test_cases(5).f_L = 1.15;
test_cases(5).fill_factor = 0.9;
test_cases(5).w_is_um = 50;
test_cases(5).I_A = 0.08;
test_cases(5).q_Wm2 = 1000;
test_cases(5).ic_w_r = 0.1;
test_cases(5).ic_t_r = 0.6;
test_cases(5).ic_angle_r = 0.5;
test_cases(5).oc_w_r = 0.1;
test_cases(5).oc_t_r = 0.6;
test_cases(5).oc_angle_r = 0.5;

% Derive theta from M for all test cases
for i = 1:length(test_cases)
    test_cases(i).theta_deg = 360 / test_cases(i).M;
end

% Display test cases
fprintf('┌─────┬──────────────────────────┬───────────────┬───────────────┬────────────┬────────────────────────────────────┐\n');
fprintf('│  #  │ Name                     │ t_TEC [µm]    │ I [mA]        │ q [W/m²]   │ Description                        │\n');
fprintf('├─────┼──────────────────────────┼───────────────┼───────────────┼────────────┼────────────────────────────────────┤\n');
for i = 1:length(test_cases)
    fprintf('│  %d  │ %-24s │ %6.0f        │ %6.0f        │ %6.0f     │ %-34s │\n', ...
        i, test_cases(i).name, test_cases(i).t_TEC_um, test_cases(i).I_A*1000, ...
        test_cases(i).q_Wm2, test_cases(i).description(1:min(34,end)));
end
fprintf('└─────┴──────────────────────────┴───────────────┴───────────────┴────────────┴────────────────────────────────────┘\n\n');

%% ============ RUN MATLAB MODEL FOR ALL TEST CASES ============
fprintf('Running MATLAB compact thermal model for all test cases...\n\n');

matlab_results = struct();
for i = 1:length(test_cases)
    tc = test_cases(i);
    fprintf('  Test Case %d: %s... ', i, tc.name);
    
    try
        % Calculate L_1 using Eq. 238 from paper
        % L_1 = (r_base - r_cyl - (N+1)*W_is) * (1-f_L) / (1-f_L^N)
        N_stages = 3;  % Fixed for 3-stage model
        w_chip = 10000e-6;  % 10mm chip
        r_base = w_chip / sqrt(2);
        r_cyl = tc.R_cyl_um * 1e-6;
        w_is = tc.w_is_um * 1e-6;
        f_L = tc.f_L;
        
        L_total_active = r_base - r_cyl - (N_stages + 1) * w_is;
        if f_L == 1
            L_1 = L_total_active / N_stages;
        else
            L_1 = L_total_active * (1 - f_L) / (1 - f_L^N_stages);
        end
        
        % Build config for MATLAB model
        config = struct();
        config.geometry.N_stages = N_stages;
        config.geometry.M_wedges = tc.M;
        config.geometry.wedge_angle_deg = tc.theta_deg;
        config.geometry.thickness_um = tc.t_TEC_um;
        config.geometry.radial_expansion_factor = tc.f_L;
        config.geometry.w_chip_um = 10000;
        config.geometry.R_cyl_um = tc.R_cyl_um;
        config.geometry.t_chip_um = tc.t_chip_um;
        config.geometry.t_ins_um = tc.t_SOI_um;
        config.geometry.interconnect_ratio = tc.ic_w_r;
        config.geometry.outerconnect_ratio = tc.oc_w_r;
        config.geometry.insulation_width_um = tc.w_is_um;
        config.geometry.interconnect_angle_ratio = tc.ic_angle_r;
        config.geometry.outerconnect_angle_ratio = tc.oc_angle_r;
        config.geometry.interconnect_thickness_ratio = tc.ic_t_r;
        config.geometry.outerconnect_thickness_ratio = tc.oc_t_r;
        
        % Azimuthal gap from fill factor
        theta_rad = deg2rad(tc.theta_deg);
        r_avg = r_base / 2;  % Approximate average radius
        arc_length = r_avg * theta_rad;
        w_az = (1 - tc.fill_factor) * arc_length;
        config.geometry.azimuthal_gap_um = w_az * 1e6;
        
        config.operating_conditions.I_current_A = tc.I_A;
        config.boundary_conditions.q_flux_W_m2 = tc.q_Wm2;
        config.boundary_conditions.T_water_K = 300;
        config.boundary_conditions.h_conv_W_m2K = 1e6;
        
        % Material properties (Bi2Te3, Cu, Si, etc.)
        config.materials.Bi2Te3 = struct('k', 1.2, 'rho', 1e-5, 'S', 0.0002);
        config.materials.Cu = struct('k', 400, 'rho', 1.7e-8);
        config.materials.Si = struct('k', 150, 'rho', 0.01);
        config.materials.AlN = struct('k', 170, 'rho', 1e10);
        config.materials.SiO2 = struct('k', 1.4, 'rho', 1e14);
        config.materials.Al2O3 = struct('k', 30, 'rho', 1e12);
        
        config.simulation.max_iterations = 100;
        config.simulation.tolerance = 1e-6;
        
        % Run MATLAB model
        materials = MaterialProperties(config);
        geometry = TECGeometry(config);
        network = ThermalNetwork(geometry, materials, config);
        
        dim = 2*N_stages + 1;
        T = ones(dim, 1) * 300;
        
        for iter = 1:100
            T_old = T;
            [T, Q_out, Q_in] = network.solve(T);
            if max(abs(T - T_old)) < 1e-6
                break;
            end
        end
        
        % Store results
        matlab_results(i).T_center = T(1) - 273.15;  % Node 0 (center) in °C
        matlab_results(i).T_max = max(T) - 273.15;
        matlab_results(i).T_all = T - 273.15;
        matlab_results(i).Q_out = Q_out;
        matlab_results(i).Q_in = Q_in;
        matlab_results(i).L_1_um = L_1 * 1e6;
        matlab_results(i).success = true;
        
        fprintf('T_center = %.2f °C, T_max = %.2f °C, L_1 = %.1f µm\n', ...
            matlab_results(i).T_center, matlab_results(i).T_max, matlab_results(i).L_1_um);
        
    catch ME
        matlab_results(i).T_center = NaN;
        matlab_results(i).T_max = NaN;
        matlab_results(i).success = false;
        matlab_results(i).error = ME.message;
        fprintf('FAILED: %s\n', ME.message);
    end
end

%% ============ SAVE MATLAB RESULTS ============
save(fullfile(OUTPUT_DIR, 'matlab_results.mat'), 'matlab_results', 'test_cases');
fprintf('\nMATLAB results saved to: %s\n\n', fullfile(OUTPUT_DIR, 'matlab_results.mat'));

%% ============ GENERATE COMSOL PARAMETER FILES ============
fprintf('Generating COMSOL parameter files for each test case...\n\n');

for i = 1:length(test_cases)
    tc = test_cases(i);
    
    % Calculate L_1 (must match MATLAB calculation)
    N_stages = 3;
    w_chip = 10000e-6;
    r_base = w_chip / sqrt(2);
    r_cyl = tc.R_cyl_um * 1e-6;
    w_is = tc.w_is_um * 1e-6;
    f_L = tc.f_L;
    
    L_total_active = r_base - r_cyl - (N_stages + 1) * w_is;
    if f_L == 1
        L_1 = L_total_active / N_stages;
    else
        L_1 = L_total_active * (1 - f_L) / (1 - f_L^N_stages);
    end
    L_1_um = L_1 * 1e6;
    
    % Store L_1 in test_cases for later use
    test_cases(i).L_1_um = L_1_um;
    
    % Write parameter file
    param_file = fullfile(OUTPUT_DIR, sprintf('comsol_params_case%d.txt', i));
    fid = fopen(param_file, 'w');
    
    fprintf(fid, '%% COMSOL Parameters for Test Case %d: %s\n', i, tc.name);
    fprintf(fid, '%% %s\n', tc.description);
    fprintf(fid, '%% Generated: %s\n\n', datestr(now));
    
    % Write parameters matching COMSOL_parameters.txt format
    fprintf(fid, 'LL_k_r %.6f ""\n', tc.f_L);
    fprintf(fid, 'LL_L_1 %.6f ""\n', L_1_um);
    fprintf(fid, 'LL_R_cyl %.0f ""\n', tc.R_cyl_um);
    fprintf(fid, 'LL_t_chip %.0f ""\n', tc.t_chip_um);
    fprintf(fid, 'LL_t_SOI %.0f ""\n', tc.t_SOI_um);
    fprintf(fid, 'LL_t_TEC %.0f ""\n', tc.t_TEC_um);
    fprintf(fid, 'LL_theta %.0f ""\n', tc.theta_deg);
    fprintf(fid, 'LL_w_is %.0f ""\n', tc.w_is_um);
    fprintf(fid, 'q_i %.0e[W/m^2] ""\n', tc.q_Wm2);
    fprintf(fid, 'I_0 %.4f[A] ""\n', tc.I_A);
    fprintf(fid, 'LL_fill_factor %.2f ""\n', tc.fill_factor);
    fprintf(fid, 'LL_ic_angle_r %.2f ""\n', tc.ic_angle_r);
    fprintf(fid, 'LL_ic_t_r %.2f ""\n', tc.ic_t_r);
    fprintf(fid, 'LL_ic_w_r %.2f ""\n', tc.ic_w_r);
    fprintf(fid, 'LL_is_r_r 0.1 ""\n');  % Keep default
    fprintf(fid, 'LL_oc_angle_r %.2f ""\n', tc.oc_angle_r);
    fprintf(fid, 'LL_oc_t_r %.2f ""\n', tc.oc_t_r);
    fprintf(fid, 'LL_oc_w_r %.2f ""\n', tc.oc_w_r);
    
    fclose(fid);
    fprintf('  Case %d: %s\n', i, param_file);
end

fprintf('\n');

%% ============ INSTRUCTIONS FOR RUNNING COMSOL ============
fprintf('╔══════════════════════════════════════════════════════════════════╗\n');
fprintf('║                    COMSOL SIMULATION PHASE                       ║\n');
fprintf('╚══════════════════════════════════════════════════════════════════╝\n\n');

fprintf('The script will now run COMSOL simulations for each test case.\n');
fprintf('COMSOL server will be started/restarted between runs.\n\n');

fprintf('COMSOL Server Path: %s\n', COMSOL_SERVER_PATH);
fprintf('COMSOL Model Path:  %s\n', COMSOL_MODEL_PATH);
fprintf('Port: %d\n\n', COMSOL_PORT);

%% ============ INITIALIZE RESULTS STORAGE ============
comsol_results = struct();
for i = 1:length(test_cases)
    comsol_results(i).T_center = NaN;
    comsol_results(i).T_max = NaN;
    comsol_results(i).sim_time = NaN;
    comsol_results(i).success = false;
end

%% ============ RUN COMSOL SIMULATIONS ============
for tc_id = 1:length(test_cases)
    tc = test_cases(tc_id);
    
    fprintf('\n');
    fprintf('┌──────────────────────────────────────────────────────────────────┐\n');
    fprintf('│ TEST CASE %d: %-52s │\n', tc_id, tc.name);
    fprintf('│ %s │\n', sprintf('%-64s', tc.description));
    fprintf('└──────────────────────────────────────────────────────────────────┘\n');
    
    fprintf('\nParameters:\n');
    fprintf('  θ = %.1f° (M=%d), t_TEC = %d µm, f_L = %.2f\n', tc.theta_deg, tc.M, tc.t_TEC_um, tc.f_L);
    fprintf('  I = %.0f mA, q = %.0f W/m², L_1 = %.1f µm\n', tc.I_A*1000, tc.q_Wm2, tc.L_1_um);
    
    try
        %% Start COMSOL Server
        fprintf('\nStarting COMSOL server...\n');
        
        % Build the server start command
        server_cmd = sprintf('Start-Process -FilePath "%s" -ArgumentList "-port %d" -PassThru', ...
            COMSOL_SERVER_PATH, COMSOL_PORT);
        
        % Use system call to start server in background
        system(sprintf('powershell -Command "%s"', server_cmd));
        
        % Wait for server to initialize
        fprintf('  Waiting for server to initialize (15 seconds)...\n');
        pause(15);
        
        %% Connect to COMSOL
        fprintf('Connecting to COMSOL server...\n');
        mphstart(COMSOL_PORT);
        import com.comsol.model.*
        import com.comsol.model.util.*
        fprintf('  Connected successfully!\n');
        
        %% Load Model
        fprintf('Loading COMSOL model...\n');
        model = mphload(COMSOL_MODEL_PATH);
        fprintf('  Model loaded.\n');
        
        %% Set Parameters
        fprintf('Setting parameters...\n');
        model.param.set('LL_k_r', sprintf('%g', tc.f_L));
        model.param.set('LL_L_1', sprintf('%g', tc.L_1_um));
        model.param.set('LL_R_cyl', sprintf('%g', tc.R_cyl_um));
        model.param.set('LL_t_chip', sprintf('%g', tc.t_chip_um));
        model.param.set('LL_t_SOI', sprintf('%g', tc.t_SOI_um));
        model.param.set('LL_t_TEC', sprintf('%g', tc.t_TEC_um));
        model.param.set('LL_theta', sprintf('%g', tc.theta_deg));
        model.param.set('LL_w_is', sprintf('%g', tc.w_is_um));
        model.param.set('q_i', sprintf('%g[W/m^2]', tc.q_Wm2));
        model.param.set('I_0', sprintf('%g[A]', tc.I_A));
        model.param.set('LL_fill_factor', sprintf('%g', tc.fill_factor));
        model.param.set('LL_ic_angle_r', sprintf('%g', tc.ic_angle_r));
        model.param.set('LL_ic_t_r', sprintf('%g', tc.ic_t_r));
        model.param.set('LL_ic_w_r', sprintf('%g', tc.ic_w_r));
        model.param.set('LL_oc_angle_r', sprintf('%g', tc.oc_angle_r));
        model.param.set('LL_oc_t_r', sprintf('%g', tc.oc_t_r));
        model.param.set('LL_oc_w_r', sprintf('%g', tc.oc_w_r));
        fprintf('  Parameters set.\n');
        
        %% Run Simulation
        fprintf('Running simulation...\n');
        sim_start = tic;
        model.study('std1').run();
        sim_time = toc(sim_start);
        fprintf('  Simulation completed in %.1f seconds.\n', sim_time);
        
        %% Extract Results
        fprintf('Extracting results...\n');
        
        % Try to get center temperature from probe
        try
            T_center = mphglobal(model, 'comp1.ppb1', 'dataset', 'dset1');
            T_center_C = T_center - 273.15;
            fprintf('  T_center (probe): %.2f °C\n', T_center_C);
        catch
            try
                T_data = mphinterp(model, 'T', 'coord', [0; 0; 0], 'dataset', 'dset1');
                T_center_C = T_data - 273.15;
                fprintf('  T_center (interp): %.2f °C\n', T_center_C);
            catch
                T_center_C = NaN;
                fprintf('  Warning: Could not extract center temperature\n');
            end
        end
        
        % Get maximum temperature
        try
            T_all = mpheval(model, 'T', 'dataset', 'dset1');
            T_max_C = max(T_all.d1) - 273.15;
            fprintf('  T_max: %.2f °C\n', T_max_C);
        catch
            T_max_C = NaN;
            fprintf('  Warning: Could not extract max temperature\n');
        end
        
        % Store results
        comsol_results(tc_id).T_center = T_center_C;
        comsol_results(tc_id).T_max = T_max_C;
        comsol_results(tc_id).sim_time = sim_time;
        comsol_results(tc_id).success = true;
        
        %% Save intermediate results
        save(fullfile(OUTPUT_DIR, sprintf('comsol_result_case%d.mat', tc_id)), ...
            'tc', 'T_center_C', 'T_max_C', 'sim_time');
        
    catch ME
        fprintf('\nERROR in test case %d: %s\n', tc_id, ME.message);
        comsol_results(tc_id).error = ME.message;
        comsol_results(tc_id).success = false;
    end
    
    %% Clean up - Server disconnects automatically after each run
    fprintf('\nCOMSOL server will disconnect. Waiting before next run...\n');
    pause(5);
end

%% ============ SAVE ALL RESULTS ============
save(fullfile(OUTPUT_DIR, 'all_results.mat'), 'test_cases', 'matlab_results', 'comsol_results');

%% ============ GENERATE COMPARISON REPORT ============
fprintf('\n\n');
fprintf('╔══════════════════════════════════════════════════════════════════════════════════════╗\n');
fprintf('║                              VERIFICATION RESULTS SUMMARY                            ║\n');
fprintf('╚══════════════════════════════════════════════════════════════════════════════════════╝\n\n');

fprintf('┌─────┬──────────────────────────┬────────────────┬────────────────┬────────────────┬───────────┐\n');
fprintf('│  #  │ Test Case                │ MATLAB T_max   │ COMSOL T_max   │ Error          │ Status    │\n');
fprintf('├─────┼──────────────────────────┼────────────────┼────────────────┼────────────────┼───────────┤\n');

errors = [];
for i = 1:length(test_cases)
    matlab_T = matlab_results(i).T_max;
    comsol_T = comsol_results(i).T_max;
    
    if ~isnan(comsol_T) && ~isnan(matlab_T)
        error_abs = comsol_T - matlab_T;
        error_pct = 100 * error_abs / matlab_T;
        errors = [errors; abs(error_pct)];
        
        if abs(error_pct) < 5
            status = '  OK  ';
        elseif abs(error_pct) < 10
            status = ' WARN ';
        else
            status = ' FAIL ';
        end
        
        fprintf('│  %d  │ %-24s │ %10.2f °C  │ %10.2f °C  │ %+6.2f °C (%+5.1f%%) │   %s  │\n', ...
            i, test_cases(i).name, matlab_T, comsol_T, error_abs, error_pct, status);
    else
        if isnan(comsol_T)
            fprintf('│  %d  │ %-24s │ %10.2f °C  │      N/A       │       N/A      │ NO DATA │\n', ...
                i, test_cases(i).name, matlab_T);
        else
            fprintf('│  %d  │ %-24s │      N/A       │ %10.2f °C  │       N/A      │ NO DATA │\n', ...
                i, test_cases(i).name, comsol_T);
        end
    end
end
fprintf('└─────┴──────────────────────────┴────────────────┴────────────────┴────────────────┴───────────┘\n\n');

if ~isempty(errors)
    fprintf('Overall Statistics:\n');
    fprintf('  Mean Absolute Error: %.2f%%\n', mean(errors));
    fprintf('  Max Absolute Error:  %.2f%%\n', max(errors));
    fprintf('  Min Absolute Error:  %.2f%%\n', min(errors));
    fprintf('  Std Deviation:       %.2f%%\n', std(errors));
end

%% ============ CENTER TEMPERATURE COMPARISON ============
fprintf('\n\nCenter Temperature (Node 0) Comparison:\n');
fprintf('┌─────┬──────────────────────────┬────────────────┬────────────────┬────────────────┐\n');
fprintf('│  #  │ Test Case                │ MATLAB T_0     │ COMSOL T_0     │ Error          │\n');
fprintf('├─────┼──────────────────────────┼────────────────┼────────────────┼────────────────┤\n');

for i = 1:length(test_cases)
    matlab_T = matlab_results(i).T_center;
    comsol_T = comsol_results(i).T_center;
    
    if ~isnan(comsol_T) && ~isnan(matlab_T)
        error_abs = comsol_T - matlab_T;
        error_pct = 100 * error_abs / matlab_T;
        fprintf('│  %d  │ %-24s │ %10.2f °C  │ %10.2f °C  │ %+6.2f °C (%+5.1f%%) │\n', ...
            i, test_cases(i).name, matlab_T, comsol_T, error_abs, error_pct);
    else
        fprintf('│  %d  │ %-24s │ %10.2f °C  │      N/A       │       N/A      │\n', ...
            i, test_cases(i).name, matlab_T);
    end
end
fprintf('└─────┴──────────────────────────┴────────────────┴────────────────┴────────────────┘\n\n');

%% ============ GENERATE PLOTS ============
fprintf('Generating comparison plots...\n');

figure('Position', [100, 100, 1200, 500]);

% Plot 1: Max Temperature Comparison
subplot(1, 2, 1);
matlab_Tmax = [matlab_results.T_max];
comsol_Tmax = [comsol_results.T_max];
valid_idx = ~isnan(comsol_Tmax);

if any(valid_idx)
    bar_data = [matlab_Tmax(valid_idx); comsol_Tmax(valid_idx)]';
    bar(bar_data);
    set(gca, 'XTickLabel', {test_cases(valid_idx).name});
    xtickangle(45);
    legend('MATLAB CTM', 'COMSOL FEM', 'Location', 'best');
    ylabel('Maximum Temperature [°C]');
    title('Maximum Temperature Comparison');
    grid on;
end

% Plot 2: Error Analysis
subplot(1, 2, 2);
if any(valid_idx)
    errors_plot = comsol_Tmax(valid_idx) - matlab_Tmax(valid_idx);
    bar(errors_plot);
    set(gca, 'XTickLabel', {test_cases(valid_idx).name});
    xtickangle(45);
    ylabel('Error (COMSOL - MATLAB) [°C]');
    title('Temperature Prediction Error');
    grid on;
    
    % Add reference lines
    hold on;
    yline(0, 'k-', 'LineWidth', 1);
    yline(mean(errors_plot), 'r--', sprintf('Mean: %.2f°C', mean(errors_plot)));
    hold off;
end

% Save figure
savefig(fullfile(OUTPUT_DIR, 'verification_comparison.fig'));
saveas(gcf, fullfile(OUTPUT_DIR, 'verification_comparison.png'));
fprintf('  Plots saved.\n');

%% ============ WRITE SUMMARY REPORT ============
report_file = fullfile(OUTPUT_DIR, 'verification_report.txt');
fid = fopen(report_file, 'w');

fprintf(fid, '==============================================================================\n');
fprintf(fid, '         MATHEMATICAL MODEL VERIFICATION REPORT\n');
fprintf(fid, '         Reference: Thermal_Network_For_Radial_TEC.tex\n');
fprintf(fid, '==============================================================================\n\n');
fprintf(fid, 'Generated: %s\n', datestr(now));
fprintf(fid, 'COMSOL Model: %s\n\n', COMSOL_MODEL_PATH);

fprintf(fid, '------------------------------------------------------------------------------\n');
fprintf(fid, 'TEST CASE SUMMARY\n');
fprintf(fid, '------------------------------------------------------------------------------\n\n');

for i = 1:length(test_cases)
    tc = test_cases(i);
    fprintf(fid, 'Test Case %d: %s\n', i, tc.name);
    fprintf(fid, '  Description: %s\n', tc.description);
    fprintf(fid, '  Parameters:\n');
    fprintf(fid, '    theta = %.1f deg (M=%d wedges)\n', tc.theta_deg, tc.M);
    fprintf(fid, '    t_TEC = %d um\n', tc.t_TEC_um);
    fprintf(fid, '    f_L = %.2f\n', tc.f_L);
    fprintf(fid, '    L_1 = %.1f um (calculated)\n', tc.L_1_um);
    fprintf(fid, '    I = %.0f mA\n', tc.I_A * 1000);
    fprintf(fid, '    q = %.0f W/m^2\n', tc.q_Wm2);
    fprintf(fid, '  Results:\n');
    fprintf(fid, '    MATLAB T_max = %.2f C\n', matlab_results(i).T_max);
    fprintf(fid, '    COMSOL T_max = %.2f C\n', comsol_results(i).T_max);
    if ~isnan(comsol_results(i).T_max) && ~isnan(matlab_results(i).T_max)
        error_val = comsol_results(i).T_max - matlab_results(i).T_max;
        error_pct = 100 * error_val / matlab_results(i).T_max;
        fprintf(fid, '    Error = %.2f C (%.1f%%)\n', error_val, error_pct);
    end
    fprintf(fid, '\n');
end

fprintf(fid, '------------------------------------------------------------------------------\n');
fprintf(fid, 'CONCLUSIONS\n');
fprintf(fid, '------------------------------------------------------------------------------\n\n');

if ~isempty(errors)
    fprintf(fid, 'Mean Absolute Error: %.2f%%\n', mean(errors));
    fprintf(fid, 'Maximum Error: %.2f%%\n', max(errors));
    
    if mean(errors) < 5
        fprintf(fid, '\nThe compact thermal model shows GOOD agreement with COMSOL FEM simulations.\n');
        fprintf(fid, 'The model is suitable for preliminary design and optimization.\n');
    elseif mean(errors) < 10
        fprintf(fid, '\nThe compact thermal model shows ACCEPTABLE agreement with COMSOL FEM simulations.\n');
        fprintf(fid, 'Consider refinements for high-precision applications.\n');
    else
        fprintf(fid, '\nThe compact thermal model shows SIGNIFICANT deviation from COMSOL FEM simulations.\n');
        fprintf(fid, 'Model refinements are recommended before using for design.\n');
    end
end

fclose(fid);
fprintf('\nReport saved to: %s\n', report_file);

%% ============ FINAL MESSAGE ============
fprintf('\n');
fprintf('╔══════════════════════════════════════════════════════════════════╗\n');
fprintf('║                    VERIFICATION COMPLETE                         ║\n');
fprintf('╚══════════════════════════════════════════════════════════════════╝\n');
fprintf('\nAll results saved to: %s\n\n', OUTPUT_DIR);

% Display file list
fprintf('Output files:\n');
files = dir(OUTPUT_DIR);
for i = 1:length(files)
    if ~files(i).isdir
        fprintf('  - %s\n', files(i).name);
    end
end
fprintf('\n');
