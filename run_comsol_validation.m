%% COMSOL Validation of Candidate Designs
% This script validates the candidate TEC designs using COMSOL high-fidelity simulation
%
% Candidate Design Parameters (from optimization):
%   - TEC thickness: 6 mm
%   - Number of stages: 5
%   - Wedge angle: 30°
%   - Operating current: 5 A
%   - Heat flux: 20 W/cm² (200 kW/m²)
%   - Target: T_chip < 80°C

clear; clc; close all;

%% Setup paths
project_root = fileparts(mfilename('fullpath'));
addpath(genpath(fullfile(project_root, 'src')));

% COMSOL model path
comsol_model_path = 'E:\Semester 7\ME4311 - MicroNano Electro Mechanical Systems and Nanotechnology\Project\COMSOL\MATLAB_API\template_3_stage.mph';

% Output directory with timestamp
timestamp = datestr(now, 'yyyy-mm-dd_HH-MM-SS');
output_dir = fullfile(project_root, 'output', 'high_fidelity', timestamp);
if ~exist(output_dir, 'dir')
    mkdir(output_dir);
end

fprintf('==============================================\n');
fprintf('  COMSOL VALIDATION OF CANDIDATE DESIGNS\n');
fprintf('==============================================\n\n');
fprintf('Output directory: %s\n\n', output_dir);

%% Define Candidate Designs to Validate
% Based on our optimization results for 20 W/cm² cooling

candidates = struct([]);

% Candidate 1: Best design for 20 W/cm²
candidates(1).name = 'Design_20Wcm2_6mm_5stage';
candidates(1).description = 'Optimal design for 20 W/cm² heat flux';
candidates(1).N_stages = 5;
candidates(1).theta_deg = 30;
candidates(1).t_TEC_um = 6000;          % 6 mm
candidates(1).k_r = 1.2;
candidates(1).I_current_A = 5.0;
candidates(1).q_flux_W_m2 = 200000;     % 20 W/cm²
candidates(1).T_matlab_predicted_C = 40.3;  % From MATLAB model

% Candidate 2: Alternative design (5mm, 5 stages)
candidates(2).name = 'Design_15Wcm2_5mm_5stage';
candidates(2).description = 'Alternative for ~15 W/cm² heat flux';
candidates(2).N_stages = 5;
candidates(2).theta_deg = 30;
candidates(2).t_TEC_um = 5000;          % 5 mm
candidates(2).k_r = 1.2;
candidates(2).I_current_A = 5.0;
candidates(2).q_flux_W_m2 = 150000;     % 15 W/cm²
candidates(2).T_matlab_predicted_C = 16.9;

% Candidate 3: More moderate design (3mm, 4 stages, 10 W/cm²)
candidates(3).name = 'Design_10Wcm2_3mm_4stage';
candidates(3).description = 'Moderate design for 10 W/cm² heat flux';
candidates(3).N_stages = 4;
candidates(3).theta_deg = 30;
candidates(3).t_TEC_um = 3000;          % 3 mm
candidates(3).k_r = 1.2;
candidates(3).I_current_A = 3.0;
candidates(3).q_flux_W_m2 = 100000;     % 10 W/cm²
candidates(3).T_matlab_predicted_C = 60.0;

%% Geometry Parameters (from SW_equations.txt and config)
% These are common parameters for all candidates

common_params = struct();

% Chip geometry
common_params.t_chip_um = 50;           % Chip thickness
common_params.t_SOI_um = 100;           % SOI layer thickness
common_params.r_chip_mm = 10;           % Half chip width (10mm x 10mm chip)
common_params.R_cyl_um = 1000;          % Cylinder radius = 1 mm

% Interconnect parameters
common_params.w_ic_um = 50;             % Inner contact width (radial)
common_params.w_oc_um = 50;             % Outer contact width (radial)
common_params.t_ic_um = 20;             % Inner contact thickness (axial)
common_params.t_oc_um = 20;             % Outer contact thickness (axial)
common_params.beta_ic_deg = 5;          % Inner contact angle
common_params.beta_oc_deg = 10;         % Outer contact angle

% Insulation parameters
common_params.w_is_um = 40;             % Insulation width (radial)
common_params.w_az_um = 30;             % Azimuthal gap

% TSV parameters
common_params.R_TSV_um = 10;            % TSV radius
common_params.P_TSV_um = 30;            % TSV pitch

% Boundary conditions
common_params.T_water_K = 300;          % Water/heatsink temperature (27°C)

%% Display Candidates
fprintf('CANDIDATES TO VALIDATE:\n');
fprintf('%-5s | %-30s | %6s | %6s | %8s | %6s | %10s | %12s\n', ...
    'ID', 'Name', 'N', 'θ(°)', 't_TEC(μm)', 'I(A)', 'q(W/cm²)', 'T_matlab(°C)');
fprintf('%s\n', repmat('-', 1, 100));

for i = 1:length(candidates)
    c = candidates(i);
    fprintf('%-5d | %-30s | %6d | %6.0f | %8.0f | %6.1f | %10.0f | %12.1f\n', ...
        i, c.name, c.N_stages, c.theta_deg, c.t_TEC_um, c.I_current_A, ...
        c.q_flux_W_m2/10000, c.T_matlab_predicted_C);
end
fprintf('\n');

%% Try to connect to COMSOL server
fprintf('Attempting to connect to COMSOL server...\n');

try
    % Ensure mphstart is on the MATLAB path (LiveLink MLI)
    if ~exist('mphstart', 'file')
        comsol_mli_path = 'F:\EngineeringSoftware\COMSOL\COMSOL63\Multiphysics\mli';
        if exist(comsol_mli_path, 'dir')
            addpath(comsol_mli_path);
            fprintf('Added COMSOL LiveLink to path: %s\n', comsol_mli_path);
        end
    end

    % If mphstart still not available, warn now (later call will error and be caught)
    if ~exist('mphstart', 'file')
        warning('mphstart not found on MATLAB path. Ensure COMSOL LiveLink is installed and on the path.');
    end

    % Check whether a COMSOL server is already listening on localhost:2036
    server_port = 2036;
    server_host = '127.0.0.1';
    server_up = false;
    try
        addr = java.net.InetSocketAddress(server_host, server_port);
        sock = java.net.Socket();
        sock.connect(addr, 1000); % 1 second timeout
        sock.close();
        server_up = true;
    catch
        server_up = false;
    end

    if server_up
        fprintf('COMSOL server is reachable at %s:%d. Skipping new server start.\n', server_host, server_port);
        % Shadow mphstart with a no-op handle so the later unconditional call to mphstart(...) does nothing
        mphstart = @(varargin) fprintf('mphstart skipped: server already running at %s:%d\n', server_host, server_port);
    else
        fprintf('No COMSOL server detected at %s:%d. Will start server when mphstart(...) is called.\n', server_host, server_port);
        % Do not call mphstart here; the script later calls mphstart(2036) (kept as-is)
    end
    
    % Connect to COMSOL server
    mphstart(2036);
    fprintf('Successfully connected to COMSOL server on port 2036.\n');
    
    % Import COMSOL API
    import com.comsol.model.*
    import com.comsol.model.util.*
    
catch ME
    fprintf('\n');
    fprintf('!!! COMSOL CONNECTION FAILED !!!\n');
    fprintf('Error: %s\n\n', ME.message);
    fprintf('Please ensure COMSOL server is running:\n');
    fprintf('  1. Open Command Prompt as Administrator\n');
    fprintf('  2. Navigate to: F:\\EngineeringSoftware\\COMSOL\\COMSOL63\\Multiphysics\\bin\\win64\n');
    fprintf('  3. Run: comsolmphserver -port 2036\n');
    fprintf('  4. Wait for "COMSOL Multiphysics server started" message\n');
    fprintf('  5. Re-run this script\n\n');
    
    % Save parameters for later use
    save(fullfile(output_dir, 'candidates_to_validate.mat'), 'candidates', 'common_params');
    fprintf('Saved candidate parameters to: %s\n', fullfile(output_dir, 'candidates_to_validate.mat'));
    return;
end

%% Load COMSOL Model
fprintf('\nLoading COMSOL model: %s\n', comsol_model_path);

try
    model = mphload(comsol_model_path);
    fprintf('Model loaded successfully.\n');
    
    % Display current parameters
    fprintf('\nCurrent model parameters:\n');
    try
        params = mphgetexpressions(model.param);
        param_names = fieldnames(params);
        fprintf('%-20s | %s\n', 'Parameter', 'Value');
        fprintf('%s\n', repmat('-', 1, 50));
        for i = 1:length(param_names)
            fprintf('%-20s | %s\n', param_names{i}, params.(param_names{i}));
        end
    catch
        fprintf('Could not retrieve parameters.\n');
    end
    
catch ME
    fprintf('Failed to load model: %s\n', ME.message);
    return;
end

%% Run Simulations for Each Candidate
results = struct([]);

for idx = 1:length(candidates)
    c = candidates(idx);
    
    fprintf('\n==============================================\n');
    fprintf('  SIMULATING CANDIDATE %d: %s\n', idx, c.name);
    fprintf('==============================================\n');
    
    try
        %% Set Parameters
        fprintf('\nSetting parameters...\n');
        
        % Geometry parameters
        model.param.set('LL_theta', sprintf('%g[deg]', c.theta_deg));
        model.param.set('LL_t_TEC', sprintf('%g[um]', c.t_TEC_um));
        model.param.set('LL_k_r', sprintf('%g', c.k_r));
        
        % Common geometry
        model.param.set('LL_t_chip', sprintf('%g[um]', common_params.t_chip_um));
        model.param.set('LL_t_SOI', sprintf('%g[um]', common_params.t_SOI_um));
        model.param.set('LL_R_cyl', sprintf('%g[um]', common_params.R_cyl_um));
        model.param.set('LL_r_chip', sprintf('%g[mm]', common_params.r_chip_mm));
        
        % Interconnect parameters
        model.param.set('LL_w_ic', sprintf('%g[um]', common_params.w_ic_um));
        model.param.set('LL_w_oc', sprintf('%g[um]', common_params.w_oc_um));
        model.param.set('LL_t_ic', sprintf('%g[um]', common_params.t_ic_um));
        model.param.set('LL_t_oc', sprintf('%g[um]', common_params.t_oc_um));
        model.param.set('LL_beta_ic', sprintf('%g[deg]', common_params.beta_ic_deg));
        model.param.set('LL_beta_oc', sprintf('%g[deg]', common_params.beta_oc_deg));
        
        % Insulation parameters
        model.param.set('LL_w_is', sprintf('%g[um]', common_params.w_is_um));
        model.param.set('LL_w_az', sprintf('%g[um]', common_params.w_az_um));
        
        % TSV parameters
        model.param.set('LL_R_TSV', sprintf('%g[um]', common_params.R_TSV_um));
        model.param.set('LL_P_TSV', sprintf('%g[um]', common_params.P_TSV_um));
        
        % Operating conditions
        model.param.set('I0', sprintf('%g[A]', c.I_current_A));
        model.param.set('q', sprintf('%g[W/m^2]', c.q_flux_W_m2));
        model.param.set('T_water', sprintf('%g[K]', common_params.T_water_K));
        
        fprintf('Parameters set successfully.\n');
        
        %% Build Geometry (if needed)
        fprintf('\nBuilding geometry...\n');
        try
            model.component('comp1').geom('geom1').run();
            fprintf('Geometry built.\n');
        catch
            fprintf('Geometry already built or auto-build enabled.\n');
        end
        
        %% Run Mesh
        fprintf('\nGenerating mesh...\n');
        try
            model.component('comp1').mesh('mesh1').run();
            fprintf('Mesh generated.\n');
        catch ME_mesh
            fprintf('Mesh warning: %s\n', ME_mesh.message);
        end
        
        %% Run Study
        fprintf('\nRunning simulation...\n');
        tic;
        model.study('std1').run();
        sim_time = toc;
        fprintf('Simulation completed in %.1f seconds.\n', sim_time);
        
        %% Extract Results
        fprintf('\nExtracting results...\n');
        
        % Get temperature field statistics
        try
            % Try global operators first
            T_data = mphglobal(model, {'minop1(T)', 'maxop1(T)'}, 'dataset', 'dset1');
            T_min_K = T_data(1);
            T_max_K = T_data(2);
        catch
            % Alternative: evaluate on domain
            T_vals = mpheval(model, 'T', 'dataset', 'dset1');
            T_min_K = min(T_vals.d1(:));
            T_max_K = max(T_vals.d1(:));
        end
        
        T_min_C = T_min_K - 273.15;
        T_max_C = T_max_K - 273.15;
        
        fprintf('  T_min = %.2f K (%.2f °C)\n', T_min_K, T_min_C);
        fprintf('  T_max = %.2f K (%.2f °C)\n', T_max_K, T_max_C);
        
        % Compare with MATLAB prediction
        error_C = T_max_C - c.T_matlab_predicted_C;
        error_pct = 100 * error_C / c.T_matlab_predicted_C;
        
        fprintf('\n  MATLAB predicted: %.2f °C\n', c.T_matlab_predicted_C);
        fprintf('  COMSOL result:    %.2f °C\n', T_max_C);
        fprintf('  Error:            %.2f °C (%.1f%%)\n', error_C, error_pct);
        
        % Check if target is met
        if T_max_C <= 80
            fprintf('  ✓ Target T < 80°C: PASSED\n');
        else
            fprintf('  ✗ Target T < 80°C: FAILED (%.1f°C above target)\n', T_max_C - 80);
        end
        
        %% Store Results
        results(idx).name = c.name;
        results(idx).candidate = c;
        results(idx).T_min_K = T_min_K;
        results(idx).T_max_K = T_max_K;
        results(idx).T_min_C = T_min_C;
        results(idx).T_max_C = T_max_C;
        results(idx).T_matlab_C = c.T_matlab_predicted_C;
        results(idx).error_C = error_C;
        results(idx).error_pct = error_pct;
        results(idx).simulation_time = sim_time;
        results(idx).target_met = T_max_C <= 80;
        
        %% Export Temperature Plot
        fprintf('\nExporting temperature plot...\n');
        
        try
            % Plot using MATLAB figure
            figure('Position', [100, 100, 800, 600]);
            mphplot(model, 'pg1');
            title(sprintf('%s: T_{max} = %.1f°C', c.name, T_max_C));
            
            % Save figure
            plot_filename = sprintf('temperature_%s.png', c.name);
            saveas(gcf, fullfile(output_dir, plot_filename));
            fprintf('Plot saved: %s\n', plot_filename);
            
            results(idx).plot_file = plot_filename;
            
            close(gcf);
        catch ME_plot
            fprintf('Plot export failed: %s\n', ME_plot.message);
            results(idx).plot_file = '';
        end
        
        %% Save Model for This Candidate
        model_filename = sprintf('model_%s.mph', c.name);
        model_path = fullfile(output_dir, model_filename);
        fprintf('\nSaving model: %s\n', model_filename);
        mphsave(model, model_path);
        results(idx).model_file = model_filename;
        
    catch ME
        fprintf('\n!!! SIMULATION FAILED !!!\n');
        fprintf('Error: %s\n', ME.message);
        
        results(idx).name = c.name;
        results(idx).candidate = c;
        results(idx).error = ME.message;
        results(idx).T_max_C = NaN;
    end
end

%% Generate Summary Report
fprintf('\n==============================================\n');
fprintf('  VALIDATION SUMMARY\n');
fprintf('==============================================\n\n');

fprintf('%-30s | %10s | %10s | %10s | %10s | %8s\n', ...
    'Design', 'q(W/cm²)', 'T_matlab', 'T_COMSOL', 'Error', 'Target');
fprintf('%s\n', repmat('-', 1, 90));

for idx = 1:length(results)
    r = results(idx);
    if ~isnan(r.T_max_C)
        if r.target_met
            target_str = 'PASS';
        else
            target_str = 'FAIL';
        end
        fprintf('%-30s | %10.0f | %10.1f | %10.1f | %10.1f | %8s\n', ...
            r.name, r.candidate.q_flux_W_m2/10000, r.T_matlab_C, ...
            r.T_max_C, r.error_C, target_str);
    else
        fprintf('%-30s | %10.0f | %10.1f | %10s | %10s | %8s\n', ...
            r.name, r.candidate.q_flux_W_m2/10000, r.T_matlab_C, ...
            'FAILED', '-', '-');
    end
end

%% Save Results
results_file = fullfile(output_dir, 'validation_results.mat');
save(results_file, 'results', 'candidates', 'common_params');
fprintf('\nResults saved to: %s\n', results_file);

% Save summary as text
summary_file = fullfile(output_dir, 'validation_summary.txt');
fid = fopen(summary_file, 'w');
fprintf(fid, 'COMSOL VALIDATION SUMMARY\n');
fprintf(fid, '========================\n\n');
fprintf(fid, 'Date: %s\n', datestr(now));
fprintf(fid, 'Model: %s\n\n', comsol_model_path);

fprintf(fid, 'RESULTS:\n');
for idx = 1:length(results)
    r = results(idx);
    fprintf(fid, '\n%s:\n', r.name);
    fprintf(fid, '  Heat flux: %.0f W/cm2\n', r.candidate.q_flux_W_m2/10000);
    fprintf(fid, '  MATLAB prediction: %.1f C\n', r.T_matlab_C);
    if ~isnan(r.T_max_C)
        fprintf(fid, '  COMSOL result: %.1f C\n', r.T_max_C);
        fprintf(fid, '  Error: %.1f C (%.1f%%)\n', r.error_C, r.error_pct);
        if r.target_met
            fprintf(fid, '  Target (T<80C): PASSED\n');
        else
            fprintf(fid, '  Target (T<80C): FAILED\n');
        end
    else
        fprintf(fid, '  COMSOL: SIMULATION FAILED\n');
    end
end
fclose(fid);
fprintf('Summary saved to: %s\n', summary_file);

fprintf('\n==============================================\n');
fprintf('  VALIDATION COMPLETE\n');
fprintf('==============================================\n');
fprintf('Output directory: %s\n', output_dir);
