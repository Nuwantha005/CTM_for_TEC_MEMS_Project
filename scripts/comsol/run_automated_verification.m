%% AUTOMATED MODEL VERIFICATION RUNNER
% This script automates the verification process:
%   1. Runs MATLAB model for all test cases
%   2. Starts COMSOL server
%   3. Runs each COMSOL simulation (with server restart between cases)
%   4. Compares results
%   5. Generates plots
%
% Usage:
%   Simply run this script. It will handle server restarts automatically.
%   Note: Uses command-line COMSOL execution for live logging.
%
% Reference: Thermal_Network_For_Radial_TEC.tex

clear; clc; close all;

%% ============ CONFIGURATION ============
addpath(genpath('src'));
addpath('F:\EngineeringSoftware\COMSOL\COMSOL63\Multiphysics\mli');

COMSOL_PORT = 2036;
COMSOL_MODEL_PATH = 'E:\Semester 7\ME4311 - MicroNano Electro Mechanical Systems and Nanotechnology\Project\COMSOL\Test_Wedge\asym2.mph';
COMSOL_SERVER_EXE = 'F:\EngineeringSoftware\COMSOL\COMSOL63\Multiphysics\bin\win64\comsolmphserver.exe';
SERVER_STARTUP_DELAY = 15;  % seconds to wait for server startup

% Output directory
timestamp = datestr(now, 'yyyy-mm-dd_HH-MM-SS');
OUTPUT_DIR = fullfile('output', 'model_verification', timestamp);
if ~exist(OUTPUT_DIR, 'dir')
    mkdir(OUTPUT_DIR);
end

%% ============ PRINT HEADER ============
fprintf('\n');
fprintf('================================================================================\n');
fprintf('           AUTOMATED MATHEMATICAL MODEL VERIFICATION\n');
fprintf('           Reference: Thermal_Network_For_Radial_TEC.tex\n');
fprintf('================================================================================\n');
fprintf('  COMSOL Server: %s\n', COMSOL_SERVER_EXE);
fprintf('  COMSOL Model:  %s\n', COMSOL_MODEL_PATH);
fprintf('  Output Dir:    %s\n', OUTPUT_DIR);
fprintf('================================================================================\n\n');

%% ============ DEFINE TEST CASES ============
fprintf('Defining 5 test cases...\n\n');

% Test cases with small thickness values (< 500 um) for fast meshing
tc = struct();

% Test Case 1: Thin TEC, Low Current
tc(1).name = 'Thin_TEC_Low_Current';
tc(1).desc = 'Validates basic thermal resistance calculations';
tc(1).M = 12; tc(1).t_TEC_um = 100; tc(1).t_chip_um = 50; tc(1).t_SOI_um = 20;
tc(1).R_cyl_um = 1000; tc(1).f_L = 1.0; tc(1).fill_factor = 0.9; tc(1).w_is_um = 50;
tc(1).I_A = 0.05; tc(1).q_Wm2 = 500;
tc(1).ic_w_r = 0.1; tc(1).ic_t_r = 0.6; tc(1).ic_angle_r = 0.5;
tc(1).oc_w_r = 0.1; tc(1).oc_t_r = 0.6; tc(1).oc_angle_r = 0.5;

% Test Case 2: Moderate Current
tc(2).name = 'Moderate_Current';
tc(2).desc = 'Validates Joule heating effects (I^2*R terms)';
tc(2).M = 12; tc(2).t_TEC_um = 150; tc(2).t_chip_um = 50; tc(2).t_SOI_um = 20;
tc(2).R_cyl_um = 1000; tc(2).f_L = 1.15; tc(2).fill_factor = 0.9; tc(2).w_is_um = 50;
tc(2).I_A = 0.10; tc(2).q_Wm2 = 500;
tc(2).ic_w_r = 0.1; tc(2).ic_t_r = 0.6; tc(2).ic_angle_r = 0.5;
tc(2).oc_w_r = 0.1; tc(2).oc_t_r = 0.6; tc(2).oc_angle_r = 0.5;

% Test Case 3: Length Ratio Test
tc(3).name = 'Length_Ratio_Test';
tc(3).desc = 'Validates L_1 calculation with f_L=1.2 (Eq. 238)';
tc(3).M = 12; tc(3).t_TEC_um = 150; tc(3).t_chip_um = 50; tc(3).t_SOI_um = 20;
tc(3).R_cyl_um = 1000; tc(3).f_L = 1.2; tc(3).fill_factor = 0.9; tc(3).w_is_um = 50;
tc(3).I_A = 0.08; tc(3).q_Wm2 = 500;
tc(3).ic_w_r = 0.1; tc(3).ic_t_r = 0.6; tc(3).ic_angle_r = 0.5;
tc(3).oc_w_r = 0.1; tc(3).oc_t_r = 0.6; tc(3).oc_angle_r = 0.5;

% Test Case 4: High Heat Flux
tc(4).name = 'High_Heat_Flux';
tc(4).desc = 'Validates heat generation terms with q=2000 W/m^2';
tc(4).M = 12; tc(4).t_TEC_um = 200; tc(4).t_chip_um = 50; tc(4).t_SOI_um = 20;
tc(4).R_cyl_um = 1000; tc(4).f_L = 1.15; tc(4).fill_factor = 0.9; tc(4).w_is_um = 50;
tc(4).I_A = 0.12; tc(4).q_Wm2 = 2000;
tc(4).ic_w_r = 0.1; tc(4).ic_t_r = 0.6; tc(4).ic_angle_r = 0.5;
tc(4).oc_w_r = 0.1; tc(4).oc_t_r = 0.6; tc(4).oc_angle_r = 0.5;

% Test Case 5: Wide Wedge Angle
tc(5).name = 'Wide_Wedge_Angle';
tc(5).desc = 'Validates angular scaling with theta=45 deg (M=8 wedges)';
tc(5).M = 8; tc(5).t_TEC_um = 150; tc(5).t_chip_um = 50; tc(5).t_SOI_um = 20;
tc(5).R_cyl_um = 1000; tc(5).f_L = 1.15; tc(5).fill_factor = 0.9; tc(5).w_is_um = 50;
tc(5).I_A = 0.08; tc(5).q_Wm2 = 1000;
tc(5).ic_w_r = 0.1; tc(5).ic_t_r = 0.6; tc(5).ic_angle_r = 0.5;
tc(5).oc_w_r = 0.1; tc(5).oc_t_r = 0.6; tc(5).oc_angle_r = 0.5;

% Calculate derived parameters
N_cases = length(tc);
for i = 1:N_cases
    tc(i).theta_deg = 360 / tc(i).M;
    
    % Calculate L_1 using Eq. 238
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

% Display test cases
fprintf('Test Cases:\n');
fprintf('  #   Name                        t_TEC   I       q        f_L     L_1\n');
fprintf('  --------------------------------------------------------------------------\n');
for i = 1:N_cases
    fprintf('  %d   %-25s %4d um  %3.0f mA  %4d W/m2  %.2f  %.1f um\n', ...
        i, tc(i).name, tc(i).t_TEC_um, tc(i).I_A*1000, tc(i).q_Wm2, tc(i).f_L, tc(i).L_1_um);
end
fprintf('\n');

%% ============ PHASE 1: MATLAB MODEL ============
fprintf('================================================================================\n');
fprintf('PHASE 1: Running MATLAB Compact Thermal Model\n');
fprintf('================================================================================\n\n');

matlab_results = struct();

for i = 1:N_cases
    test = tc(i);
    fprintf('Case %d: %s... ', i, test.name);
    
    try
        config = struct();
        config.geometry.N_stages = 3;
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
        config.geometry.fill_factor = test.fill_factor;  % For R_TE calculations
        
        config.operating_conditions.I_current_A = test.I_A;
        config.boundary_conditions.q_flux_W_m2 = test.q_Wm2;
        config.boundary_conditions.T_water_K = 300;
        config.boundary_conditions.h_conv_W_m2K = 1e6;
        
        % Material properties from COMSOL database (data/material_props/COMSOL.xml)
        % Bi2Te3: Temperature-dependent from data files, fallback values at ~300K
        config.materials.Bi2Te3 = struct('k', 1.6, 'rho', 1.15e-5, 'S', 210e-6);
        % Cu: k=400 W/mK, rho=1.667e-8 Ohm*m (from linearized resistivity at 293K)
        config.materials.Cu = struct('k', 400, 'rho', 1.667e-8);
        % Si: k=130 W/mK (single-crystal, isotropic)
        config.materials.Si = struct('k', 130, 'rho', 0.01);
        % AlN: Not in COMSOL.xml, keep existing value
        config.materials.AlN = struct('k', 170, 'rho', 1e10);
        % SiO2: k=1.4 W/mK
        config.materials.SiO2 = struct('k', 1.4, 'rho', 1e14);
        % Al2O3: k=35 W/mK
        config.materials.Al2O3 = struct('k', 35, 'rho', 1e12);
        
        materials = MaterialProperties(config);
        geometry = TECGeometry(config);
        network = ThermalNetwork(geometry, materials, config);
        
        dim = 7;  % 2*3 + 1
        T = ones(dim, 1) * 300;
        
        for iter = 1:100
            T_old = T;
            [T, Q_out, Q_in] = network.solve(T);
            if max(abs(T - T_old)) < 1e-6
                break;
            end
        end
        
        matlab_results(i).T_center = T(1) - 273.15;
        matlab_results(i).T_max = max(T) - 273.15;
        matlab_results(i).success = true;
        
        fprintf('T_max = %.2f C\n', matlab_results(i).T_max);
        
    catch ME
        matlab_results(i).T_center = NaN;
        matlab_results(i).T_max = NaN;
        matlab_results(i).success = false;
        fprintf('FAILED: %s\n', ME.message);
    end
end

fprintf('\nMATLAB phase complete.\n\n');

%% ============ PHASE 2: COMSOL SIMULATIONS ============
fprintf('================================================================================\n');
fprintf('PHASE 2: Running COMSOL FEM Simulations\n');
fprintf('================================================================================\n\n');

comsol_results = struct();

for case_id = 1:N_cases
    test = tc(case_id);
    
    fprintf('────────────────────────────────────────────────────────────────────────────────\n');
    fprintf('TEST CASE %d: %s\n', case_id, test.name);
    fprintf('  %s\n', test.desc);
    fprintf('────────────────────────────────────────────────────────────────────────────────\n');
    
    %% Start COMSOL Server
    fprintf('\nStarting COMSOL server...\n');
    
    % Kill any existing server first
    system('taskkill /F /IM comsolmphserver.exe 2>nul', '-echo');
    pause(2);
    
    % Start new server
    server_cmd = sprintf('start "" "%s" -port %d', COMSOL_SERVER_EXE, COMSOL_PORT);
    system(server_cmd);
    
    fprintf('  Waiting %d seconds for server initialization...\n', SERVER_STARTUP_DELAY);
    pause(SERVER_STARTUP_DELAY);
    
    %% Connect and Run
    try
        fprintf('Connecting to COMSOL...\n');
        mphstart(COMSOL_PORT);
        import com.comsol.model.*
        import com.comsol.model.util.*
        fprintf('  Connected!\n');
        
        fprintf('Loading model...\n');
        model = mphload(COMSOL_MODEL_PATH);
        fprintf('  Model loaded.\n');
        
        fprintf('Setting parameters...\n');
        model.param.set('LL_k_r', sprintf('%g', test.f_L));
        model.param.set('LL_L_1', sprintf('%g', test.L_1_um));
        model.param.set('LL_R_cyl', sprintf('%g', test.R_cyl_um));
        model.param.set('LL_t_chip', sprintf('%g', test.t_chip_um));
        model.param.set('LL_t_SOI', sprintf('%g', test.t_SOI_um));
        model.param.set('LL_t_TEC', sprintf('%g', test.t_TEC_um));
        model.param.set('LL_theta', sprintf('%g', test.theta_deg));
        model.param.set('LL_w_is', sprintf('%g', test.w_is_um));
        model.param.set('q_i', sprintf('%g[W/m^2]', test.q_Wm2));
        model.param.set('I_0', sprintf('%g[A]', test.I_A));
        model.param.set('LL_fill_factor', sprintf('%g', test.fill_factor));
        model.param.set('LL_ic_angle_r', sprintf('%g', test.ic_angle_r));
        model.param.set('LL_ic_t_r', sprintf('%g', test.ic_t_r));
        model.param.set('LL_ic_w_r', sprintf('%g', test.ic_w_r));
        model.param.set('LL_oc_angle_r', sprintf('%g', test.oc_angle_r));
        model.param.set('LL_oc_t_r', sprintf('%g', test.oc_t_r));
        model.param.set('LL_oc_w_r', sprintf('%g', test.oc_w_r));
        fprintf('  Parameters set.\n');
        
        fprintf('Running simulation...\n');
        sim_start = tic;
        model.study('std1').run();
        sim_time = toc(sim_start);
        fprintf('  Completed in %.1f seconds.\n', sim_time);
        
        fprintf('Extracting results...\n');
        
        % Get center temperature
        try
            T_center = mphglobal(model, 'comp1.ppb1', 'dataset', 'dset1');
            comsol_results(case_id).T_center = T_center - 273.15;
        catch
            try
                T_data = mphinterp(model, 'T', 'coord', [0; 0; 0], 'dataset', 'dset1');
                comsol_results(case_id).T_center = T_data - 273.15;
            catch
                comsol_results(case_id).T_center = NaN;
            end
        end
        
        % Get max temperature
        try
            T_all = mpheval(model, 'T', 'dataset', 'dset1');
            comsol_results(case_id).T_max = max(T_all.d1) - 273.15;
        catch
            comsol_results(case_id).T_max = NaN;
        end
        
        comsol_results(case_id).sim_time = sim_time;
        comsol_results(case_id).success = true;
        
        fprintf('  T_center = %.2f C, T_max = %.2f C\n', ...
            comsol_results(case_id).T_center, comsol_results(case_id).T_max);
        
    catch ME
        fprintf('\nERROR: %s\n', ME.message);
        comsol_results(case_id).T_center = NaN;
        comsol_results(case_id).T_max = NaN;
        comsol_results(case_id).sim_time = NaN;
        comsol_results(case_id).success = false;
    end
    
    % Clean up for next iteration
    fprintf('\nCleaning up...\n');
    clear model;
    pause(3);
end

%% ============ PHASE 3: RESULTS COMPARISON ============
fprintf('\n\n');
fprintf('================================================================================\n');
fprintf('PHASE 3: Results Comparison\n');
fprintf('================================================================================\n\n');

% Create results table
fprintf('VERIFICATION RESULTS:\n');
fprintf('┌─────┬──────────────────────────┬────────────┬────────────┬────────────┬──────────┐\n');
fprintf('│  #  │ Test Case                │ MATLAB [C] │ COMSOL [C] │ Error [C]  │ Error [%%] │\n');
fprintf('├─────┼──────────────────────────┼────────────┼────────────┼────────────┼──────────┤\n');

errors_pct = [];
for i = 1:N_cases
    m_T = matlab_results(i).T_max;
    c_T = comsol_results(i).T_max;
    
    if ~isnan(c_T) && ~isnan(m_T)
        err_abs = c_T - m_T;
        err_pct = 100 * err_abs / m_T;
        errors_pct = [errors_pct; abs(err_pct)];
        
        fprintf('│  %d  │ %-24s │ %10.2f │ %10.2f │ %+10.2f │ %+8.1f │\n', ...
            i, tc(i).name, m_T, c_T, err_abs, err_pct);
    else
        fprintf('│  %d  │ %-24s │ %10.2f │     N/A    │     N/A    │    N/A   │\n', ...
            i, tc(i).name, m_T);
    end
end
fprintf('└─────┴──────────────────────────┴────────────┴────────────┴────────────┴──────────┘\n\n');

% Statistics
if ~isempty(errors_pct)
    fprintf('Statistics:\n');
    fprintf('  Mean Absolute Error: %.2f%%\n', mean(errors_pct));
    fprintf('  Max Absolute Error:  %.2f%%\n', max(errors_pct));
    fprintf('  Min Absolute Error:  %.2f%%\n', min(errors_pct));
    fprintf('  Std Deviation:       %.2f%%\n\n', std(errors_pct));
    
    % Verdict
    if mean(errors_pct) < 5
        fprintf('VERDICT: EXCELLENT agreement - Model suitable for design optimization.\n');
    elseif mean(errors_pct) < 10
        fprintf('VERDICT: GOOD agreement - Model acceptable for preliminary design.\n');
    else
        fprintf('VERDICT: MODERATE agreement - Consider model refinements.\n');
    end
end

%% ============ SAVE ALL RESULTS ============
fprintf('\n\nSaving results...\n');

% Save to MAT file
save(fullfile(OUTPUT_DIR, 'verification_results.mat'), 'tc', 'matlab_results', 'comsol_results');
fprintf('  Saved: verification_results.mat\n');

% Save to CSV
csv_file = fullfile(OUTPUT_DIR, 'verification_results.csv');
fid = fopen(csv_file, 'w');
fprintf(fid, 'CaseID,Name,theta_deg,t_TEC_um,f_L,L_1_um,I_mA,q_Wm2,MATLAB_Tmax,COMSOL_Tmax,Error_C,Error_pct,SimTime\n');
for i = 1:N_cases
    m_T = matlab_results(i).T_max;
    c_T = comsol_results(i).T_max;
    err_abs = c_T - m_T;
    err_pct = 100 * err_abs / m_T;
    sim_t = comsol_results(i).sim_time;
    
    fprintf(fid, '%d,%s,%.1f,%d,%.2f,%.2f,%.0f,%.0f,%.2f,%.2f,%.2f,%.1f,%.1f\n', ...
        i, tc(i).name, tc(i).theta_deg, tc(i).t_TEC_um, tc(i).f_L, tc(i).L_1_um, ...
        tc(i).I_A*1000, tc(i).q_Wm2, m_T, c_T, err_abs, err_pct, sim_t);
end
fclose(fid);
fprintf('  Saved: verification_results.csv\n');

%% ============ GENERATE PLOTS ============
fprintf('\nGenerating comparison plots...\n');

fig = figure('Position', [100 100 1200 500], 'Color', 'w');

% Temperature comparison bar chart
subplot(1, 2, 1);
matlab_T = [matlab_results.T_max];
comsol_T = [comsol_results.T_max];
bar_data = [matlab_T; comsol_T]';
b = bar(bar_data);
b(1).FaceColor = [0.2 0.4 0.8];
b(2).FaceColor = [0.8 0.3 0.2];
set(gca, 'XTickLabel', strrep({tc.name}, '_', ' '));
xtickangle(30);
ylabel('Maximum Temperature [°C]');
title('Temperature Comparison');
legend('MATLAB CTM', 'COMSOL FEM', 'Location', 'northwest');
grid on;

% Error plot
subplot(1, 2, 2);
errors = comsol_T - matlab_T;
bar(errors, 'FaceColor', [0.4 0.7 0.4]);
set(gca, 'XTickLabel', strrep({tc.name}, '_', ' '));
xtickangle(30);
ylabel('Error (COMSOL - MATLAB) [°C]');
title('Prediction Error');
hold on;
yline(0, 'k-', 'LineWidth', 1.5);
yline(mean(errors(~isnan(errors))), 'r--', sprintf('Mean: %.2f°C', mean(errors(~isnan(errors)))));
hold off;
grid on;

savefig(fullfile(OUTPUT_DIR, 'verification_comparison.fig'));
saveas(fig, fullfile(OUTPUT_DIR, 'verification_comparison.png'));
fprintf('  Saved: verification_comparison.fig/png\n');

%% ============ FINAL SUMMARY ============
fprintf('\n');
fprintf('================================================================================\n');
fprintf('                         VERIFICATION COMPLETE\n');
fprintf('================================================================================\n');
fprintf('  Results saved to: %s\n', OUTPUT_DIR);
fprintf('  Files generated:\n');
fprintf('    - verification_results.mat\n');
fprintf('    - verification_results.csv\n');
fprintf('    - verification_comparison.fig\n');
fprintf('    - verification_comparison.png\n');
fprintf('================================================================================\n\n');
