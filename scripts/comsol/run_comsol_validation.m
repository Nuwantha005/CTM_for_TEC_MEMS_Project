% run_comsol_validation.m
% Run COMSOL validation for candidate TEC designs
%
% Prerequisites:
%   1. Start COMSOL server: comsolmphserver -port 2036
%   2. Have the COMSOL template model ready
%   3. Model should have 3 stages with TSV in stage 1 only
%
% This script will:
%   1. Generate candidate designs compatible with 3-stage TSV-stage-1 template
%   2. Run each through COMSOL
%   3. Extract "Node 0 Temp" probe results
%   4. Run parametric sweeps
%   5. Export data and generate plots

clear; clc;
addpath(genpath('src'));

fprintf('=== COMSOL VALIDATION FOR TEC DESIGNS ===\n\n');

%% Configuration
COMSOL_PORT = 2036;
COMSOL_MODEL_PATH = '';  % Will prompt user if empty

% Output directory
timestamp = datestr(now, 'yyyy-mm-dd_HH-MM-SS');
OUTPUT_DIR = fullfile('output', 'comsol_validation', timestamp);
if ~exist(OUTPUT_DIR, 'dir')
    mkdir(OUTPUT_DIR);
end

fprintf('Output directory: %s\n\n', OUTPUT_DIR);

%% Step 1: Find or request COMSOL model path
if isempty(COMSOL_MODEL_PATH)
    fprintf('Please provide the path to your COMSOL template model (.mph file):\n');
    fprintf('Note: Model should have 3 stages with TSV in stage 1 only.\n\n');
    
    % Check common locations
    common_paths = {
        fullfile(pwd, '..', 'COMSOL', 'TEC_Template.mph'),
        fullfile(pwd, '..', 'COMSOL', 'radial_tec_3stage.mph'),
        'F:\Projects\TEC\comsol_model.mph'
    };
    
    for i = 1:length(common_paths)
        if exist(common_paths{i}, 'file')
            fprintf('Found model at: %s\n', common_paths{i});
            COMSOL_MODEL_PATH = common_paths{i};
            break;
        end
    end
    
    if isempty(COMSOL_MODEL_PATH)
        % Prompt user
        [file, path] = uigetfile('*.mph', 'Select COMSOL Template Model');
        if file == 0
            error('No model file selected. Exiting.');
        end
        COMSOL_MODEL_PATH = fullfile(path, file);
    end
end

fprintf('Using COMSOL model: %s\n\n', COMSOL_MODEL_PATH);

%% Step 2: Generate candidate designs for 3-stage, TSV-stage-1 configuration
fprintf('=== GENERATING CANDIDATE DESIGNS ===\n\n');

% Base configuration matching template constraints
base_config = struct();
base_config.geometry.N_stages = 3;  % FIXED - template has 3 stages
base_config.geometry.N_tsv_limit = 1;  % FIXED - TSV only in stage 1
base_config.geometry.w_chip_um = 10000;  % 10mm chip
base_config.geometry.R_cyl_um = 1000;
base_config.geometry.t_chip_um = 50;
base_config.boundary_conditions.T_water_K = 300;
base_config.boundary_conditions.h_conv_W_m2K = 1e6;

% Parameters to vary for candidate generation
% We'll create candidates by varying key parameters within feasible ranges
candidates = [];

% From our optimization results, these are the key parameters:
% Current: 0.02-0.15 A optimal range
% Thickness: 100-300 um optimal
% k_r: 0.8-1.2
% Wedge angle: 25-40 deg
% Fill factor: 0.9-0.99

currents = [0.02, 0.05, 0.1];
thicknesses = [100, 200, 300];  % um
k_r_values = [0.9, 1.0, 1.15];
wedge_angles = [30];  % Fixed to match template
fill_factors = [0.9, 0.95];

% Heat flux levels to test
q_flux_values = [500, 1000, 2000];  % W/m²

candidate_id = 0;
for I = currents
    for t = thicknesses
        for kr = k_r_values
            for ff = fill_factors
                for q = q_flux_values
                    candidate_id = candidate_id + 1;
                    
                    c = base_config;
                    c.operating_conditions.I_current_A = I;
                    c.geometry.thickness_um = t;
                    c.geometry.radial_expansion_factor = kr;
                    c.geometry.fill_factor = ff;
                    c.geometry.wedge_angle_deg = 30;  % Fixed
                    c.geometry.interconnect_ratio = 0.15;
                    c.geometry.outerconnect_ratio = 0.15;
                    c.geometry.insulation_width_ratio = 0.04;
                    c.geometry.interconnect_angle_ratio = 0.16;
                    c.geometry.outerconnect_angle_ratio = 0.16;
                    c.boundary_conditions.q_flux_W_m2 = q;
                    
                    c.id = candidate_id;
                    c.label = sprintf('I%.0fmA_t%dum_kr%.2f_ff%.2f_q%d', ...
                        I*1000, t, kr, ff, q);
                    
                    candidates = [candidates; c];
                end
            end
        end
    end
end

fprintf('Generated %d candidate designs\n\n', length(candidates));

% Run preliminary MATLAB simulation to filter obviously bad candidates
fprintf('Pre-filtering candidates with MATLAB model...\n');
materials = MaterialProperties(candidates(1));

feasible_candidates = [];
for i = 1:length(candidates)
    c = candidates(i);
    try
        config = struct();
        config.geometry = c.geometry;
        config.boundary_conditions = c.boundary_conditions;
        config.operating_conditions = c.operating_conditions;
        config.simulation = struct('max_iterations', 100, 'tolerance', 1e-6, ...
            'T_ambient', 300, 'T_initial_guess', 300);
        config.materials = struct('Bi2Te3', struct('k', 1.2, 'rho', 1e-5, 'S', 0.0002));
        
        geometry = TECGeometry(config);
        network = ThermalNetwork(geometry, materials, config);
        
        N = geometry.N_stages;
        dim = 2*N + 1;
        T = ones(dim, 1) * 300;
        
        for iter = 1:50
            T_old = T;
            [T, ~, ~] = network.solve(T);
            if max(abs(T - T_old)) < 1e-4
                break;
            end
        end
        
        T_max = max(T) - 273.15;  % Convert to Celsius
        
        if T_max > 0 && T_max < 200  % Reasonable range
            c.T_matlab = T_max;
            c.T_profile = T;
            feasible_candidates = [feasible_candidates; c];
            fprintf('  Candidate %d (%s): T_max = %.1f°C - FEASIBLE\n', ...
                c.id, c.label, T_max);
        else
            fprintf('  Candidate %d: T_max = %.1f°C - REJECTED\n', c.id, T_max);
        end
    catch ME
        fprintf('  Candidate %d: FAILED - %s\n', c.id, ME.message);
    end
end

fprintf('\n%d feasible candidates after MATLAB pre-filtering\n\n', length(feasible_candidates));

% Select top candidates for COMSOL (limit to save time)
MAX_COMSOL_RUNS = 20;
if length(feasible_candidates) > MAX_COMSOL_RUNS
    % Sort by temperature and select best
    [~, idx] = sort([feasible_candidates.T_matlab]);
    feasible_candidates = feasible_candidates(idx(1:MAX_COMSOL_RUNS));
    fprintf('Selected top %d candidates for COMSOL validation\n\n', MAX_COMSOL_RUNS);
end

% Save candidates
save(fullfile(OUTPUT_DIR, 'candidates.mat'), 'candidates', 'feasible_candidates');

%% Step 3: Connect to COMSOL
fprintf('=== CONNECTING TO COMSOL ===\n\n');

try
    % Add COMSOL LiveLink to path
    comsol_mli_path = 'F:\EngineeringSoftware\COMSOL\COMSOL63\Multiphysics\mli';
    if exist(comsol_mli_path, 'dir')
        addpath(comsol_mli_path);
        fprintf('Added COMSOL LiveLink to path: %s\n', comsol_mli_path);
    end
    
    % Create interface
    comsol = COMSOLInterface('Port', COMSOL_PORT, 'OutputDir', OUTPUT_DIR);
    
    % Connect
    if ~comsol.connect()
        fprintf('\n*** COMSOL Connection Failed ***\n');
        fprintf('Please ensure COMSOL server is running:\n');
        fprintf('  F:\\EngineeringSoftware\\COMSOL\\COMSOL63\\Multiphysics\\bin\\win64\\comsolmphserver -port %d\n\n', COMSOL_PORT);
        
        % Save what we have so far
        save(fullfile(OUTPUT_DIR, 'partial_results.mat'), 'feasible_candidates');
        fprintf('Saved MATLAB results to: %s\n', fullfile(OUTPUT_DIR, 'partial_results.mat'));
        
        % Generate MATLAB-only report
        generate_matlab_report(feasible_candidates, OUTPUT_DIR);
        return;
    end
    
    % Load model
    if ~comsol.loadModel(COMSOL_MODEL_PATH)
        error('Failed to load COMSOL model');
    end
    
catch ME
    fprintf('COMSOL setup error: %s\n', ME.message);
    fprintf('Continuing with MATLAB-only analysis...\n\n');
    
    generate_matlab_report(feasible_candidates, OUTPUT_DIR);
    return;
end

%% Step 4: Run COMSOL simulations for each candidate
fprintf('\n=== RUNNING COMSOL SIMULATIONS ===\n\n');

comsol_results = [];

for i = 1:length(feasible_candidates)
    c = feasible_candidates(i);
    fprintf('--- Candidate %d/%d: %s ---\n', i, length(feasible_candidates), c.label);
    
    try
        % Set parameters in COMSOL
        params = {
            'LL_theta', sprintf('%g[deg]', c.geometry.wedge_angle_deg), ...
            'LL_t_TEC', sprintf('%g[um]', c.geometry.thickness_um), ...
            'LL_k_r', sprintf('%g', c.geometry.radial_expansion_factor), ...
            'I0', sprintf('%g[A]', c.operating_conditions.I_current_A), ...
            'q', sprintf('%g[W/m^2]', c.boundary_conditions.q_flux_W_m2), ...
            'LL_R_cyl', sprintf('%g[um]', c.geometry.R_cyl_um), ...
            'LL_t_chip', sprintf('%g[um]', c.geometry.t_chip_um)
        };
        
        comsol.setParametersDirect(params);
        
        % Run study
        tic;
        if comsol.runStudy('std1')
            elapsed = toc;
            
            % Extract Node 0 temperature (central cylinder probe)
            try
                T_node0 = mphglobal(comsol.Model, 'comp1.ppb1', 'dataset', 'dset1');
                fprintf('  Node 0 Temp (COMSOL): %.2f K (%.2f°C)\n', T_node0, T_node0 - 273.15);
            catch
                % Try alternative extraction
                try
                    results = comsol.extractResults();
                    T_node0 = results.T_max_K;
                    fprintf('  T_max (COMSOL): %.2f K (%.2f°C)\n', T_node0, T_node0 - 273.15);
                catch
                    T_node0 = NaN;
                    fprintf('  Could not extract temperature\n');
                end
            end
            
            % Store result
            result = struct();
            result.id = c.id;
            result.label = c.label;
            result.params = c;
            result.T_comsol_K = T_node0;
            result.T_comsol_C = T_node0 - 273.15;
            result.T_matlab_C = c.T_matlab;
            result.error_percent = 100 * abs(result.T_comsol_C - c.T_matlab) / c.T_matlab;
            result.elapsed_time = elapsed;
            
            comsol_results = [comsol_results; result];
            
            fprintf('  MATLAB: %.2f°C, COMSOL: %.2f°C, Error: %.1f%%\n', ...
                c.T_matlab, result.T_comsol_C, result.error_percent);
            fprintf('  Simulation time: %.1f s\n\n', elapsed);
            
            % Save intermediate results
            save(fullfile(OUTPUT_DIR, 'comsol_results.mat'), 'comsol_results');
        else
            fprintf('  COMSOL simulation failed\n\n');
        end
        
    catch ME
        fprintf('  Error: %s\n\n', ME.message);
    end
end

%% Step 5: Run parametric sweeps for best design
fprintf('\n=== PARAMETRIC SWEEPS ===\n\n');

% Find best design
if ~isempty(comsol_results)
    [~, best_idx] = min([comsol_results.T_comsol_C]);
    best_design = comsol_results(best_idx);
    fprintf('Best design: %s (T = %.2f°C)\n\n', best_design.label, best_design.T_comsol_C);
    
    % Sweep 1: Current sweep
    fprintf('Running current sweep...\n');
    current_sweep = struct();
    current_sweep.I_values = linspace(0.01, 0.2, 10);
    current_sweep.T_values = zeros(size(current_sweep.I_values));
    
    for j = 1:length(current_sweep.I_values)
        try
            comsol.Model.param.set('I0', sprintf('%g[A]', current_sweep.I_values(j)));
            comsol.runStudy('std1');
            T = mphglobal(comsol.Model, 'comp1.ppb1', 'dataset', 'dset1');
            current_sweep.T_values(j) = T - 273.15;
            fprintf('  I = %.3f A: T = %.2f°C\n', current_sweep.I_values(j), current_sweep.T_values(j));
        catch
            current_sweep.T_values(j) = NaN;
        end
    end
    
    % Sweep 2: Heat flux sweep
    fprintf('\nRunning heat flux sweep...\n');
    % Reset to best current
    comsol.Model.param.set('I0', sprintf('%g[A]', best_design.params.operating_conditions.I_current_A));
    
    flux_sweep = struct();
    flux_sweep.q_values = linspace(100, 3000, 10);
    flux_sweep.T_values = zeros(size(flux_sweep.q_values));
    
    for j = 1:length(flux_sweep.q_values)
        try
            comsol.Model.param.set('q', sprintf('%g[W/m^2]', flux_sweep.q_values(j)));
            comsol.runStudy('std1');
            T = mphglobal(comsol.Model, 'comp1.ppb1', 'dataset', 'dset1');
            flux_sweep.T_values(j) = T - 273.15;
            fprintf('  q = %.0f W/m²: T = %.2f°C\n', flux_sweep.q_values(j), flux_sweep.T_values(j));
        catch
            flux_sweep.T_values(j) = NaN;
        end
    end
    
    % Sweep 3: TEC thickness sweep
    fprintf('\nRunning thickness sweep...\n');
    % Reset to best flux
    comsol.Model.param.set('q', sprintf('%g[W/m^2]', best_design.params.boundary_conditions.q_flux_W_m2));
    
    thickness_sweep = struct();
    thickness_sweep.t_values = [50, 100, 150, 200, 250, 300, 400, 500];
    thickness_sweep.T_values = zeros(size(thickness_sweep.t_values));
    
    for j = 1:length(thickness_sweep.t_values)
        try
            comsol.Model.param.set('LL_t_TEC', sprintf('%g[um]', thickness_sweep.t_values(j)));
            comsol.runStudy('std1');
            T = mphglobal(comsol.Model, 'comp1.ppb1', 'dataset', 'dset1');
            thickness_sweep.T_values(j) = T - 273.15;
            fprintf('  t = %d um: T = %.2f°C\n', thickness_sweep.t_values(j), thickness_sweep.T_values(j));
        catch
            thickness_sweep.T_values(j) = NaN;
        end
    end
    
    % Save sweep results
    sweeps = struct('current', current_sweep, 'flux', flux_sweep, 'thickness', thickness_sweep);
    save(fullfile(OUTPUT_DIR, 'parametric_sweeps.mat'), 'sweeps');
end

%% Step 6: Generate plots and report
fprintf('\n=== GENERATING PLOTS AND REPORT ===\n\n');

generate_validation_report(comsol_results, sweeps, OUTPUT_DIR);

% Disconnect from COMSOL
comsol.disconnect();

fprintf('\n=== VALIDATION COMPLETE ===\n');
fprintf('Results saved to: %s\n', OUTPUT_DIR);


%% Helper Functions

function generate_matlab_report(candidates, output_dir)
    % Generate report using only MATLAB results
    
    fprintf('Generating MATLAB-only report...\n');
    
    % Sort by temperature
    [temps, idx] = sort([candidates.T_matlab]);
    candidates = candidates(idx);
    
    % Create table
    fprintf('\nTop 10 Candidates (MATLAB prediction):\n');
    fprintf('%-5s | %-30s | %10s\n', 'Rank', 'Configuration', 'T_max (°C)');
    fprintf('%s\n', repmat('-', 1, 55));
    
    for i = 1:min(10, length(candidates))
        fprintf('%-5d | %-30s | %10.1f\n', i, candidates(i).label, candidates(i).T_matlab);
    end
    
    % Plot
    figure('Position', [100, 100, 800, 600]);
    
    subplot(2,2,1);
    bar(temps(1:min(20, length(temps))));
    xlabel('Candidate Rank');
    ylabel('T_{max} (°C)');
    title('MATLAB Temperature Predictions (Top 20)');
    yline(100, 'r--', 'Target: 100°C');
    
    subplot(2,2,2);
    scatter([candidates.T_matlab], [candidates(:).geometry.thickness_um], 50, 'filled');
    xlabel('T_{max} (°C)');
    ylabel('TEC Thickness (µm)');
    title('Temperature vs Thickness');
    
    subplot(2,2,3);
    scatter([candidates.T_matlab], [candidates(:).operating_conditions.I_current_A]*1000, 50, 'filled');
    xlabel('T_{max} (°C)');
    ylabel('Current (mA)');
    title('Temperature vs Current');
    
    subplot(2,2,4);
    scatter([candidates.T_matlab], [candidates(:).boundary_conditions.q_flux_W_m2], 50, 'filled');
    xlabel('T_{max} (°C)');
    ylabel('Heat Flux (W/m²)');
    title('Temperature vs Heat Flux');
    
    saveas(gcf, fullfile(output_dir, 'matlab_analysis.png'));
    saveas(gcf, fullfile(output_dir, 'matlab_analysis.fig'));
    
    % Save report
    fid = fopen(fullfile(output_dir, 'matlab_report.txt'), 'w');
    fprintf(fid, 'MATLAB Analysis Report\n');
    fprintf(fid, '======================\n\n');
    fprintf(fid, 'Date: %s\n\n', datestr(now));
    fprintf(fid, 'Total candidates evaluated: %d\n\n', length(candidates));
    fprintf(fid, 'Best design:\n');
    fprintf(fid, '  Configuration: %s\n', candidates(1).label);
    fprintf(fid, '  T_max: %.1f°C\n', candidates(1).T_matlab);
    fprintf(fid, '  Current: %.0f mA\n', candidates(1).operating_conditions.I_current_A * 1000);
    fprintf(fid, '  Thickness: %d µm\n', candidates(1).geometry.thickness_um);
    fprintf(fid, '  Heat flux: %d W/m²\n', candidates(1).boundary_conditions.q_flux_W_m2);
    fclose(fid);
    
    fprintf('Report saved to: %s\n', fullfile(output_dir, 'matlab_report.txt'));
end

function generate_validation_report(results, sweeps, output_dir)
    % Generate comprehensive validation report with COMSOL data
    
    fprintf('Generating validation report...\n');
    
    %% Figure 1: MATLAB vs COMSOL comparison
    figure('Position', [100, 100, 1200, 800], 'Name', 'COMSOL Validation Results');
    
    subplot(2,3,1);
    T_matlab = [results.T_matlab_C];
    T_comsol = [results.T_comsol_C];
    scatter(T_matlab, T_comsol, 80, 'filled');
    hold on;
    plot([min(T_matlab), max(T_matlab)], [min(T_matlab), max(T_matlab)], 'r--', 'LineWidth', 2);
    xlabel('MATLAB T_{max} (°C)');
    ylabel('COMSOL T_{max} (°C)');
    title('MATLAB vs COMSOL Correlation');
    legend('Data', 'Perfect Correlation', 'Location', 'best');
    grid on;
    
    subplot(2,3,2);
    errors = [results.error_percent];
    bar(errors);
    xlabel('Design Index');
    ylabel('Error (%)');
    title(sprintf('Prediction Error (Mean: %.1f%%)', mean(errors)));
    
    subplot(2,3,3);
    bar([T_matlab; T_comsol]');
    xlabel('Design Index');
    ylabel('Temperature (°C)');
    title('Temperature Comparison');
    legend('MATLAB', 'COMSOL');
    yline(100, 'r--', 'Target');
    
    %% Parametric sweep plots
    subplot(2,3,4);
    if isfield(sweeps, 'current') && ~isempty(sweeps.current.T_values)
        plot(sweeps.current.I_values * 1000, sweeps.current.T_values, '-o', 'LineWidth', 2);
        xlabel('Current (mA)');
        ylabel('T_{max} (°C)');
        title('Temperature vs Current (COMSOL)');
        grid on;
        
        % Find optimal current
        [T_min, idx] = min(sweeps.current.T_values);
        hold on;
        plot(sweeps.current.I_values(idx)*1000, T_min, 'r*', 'MarkerSize', 15);
        text(sweeps.current.I_values(idx)*1000 + 5, T_min, sprintf('Optimal: %.0f mA', sweeps.current.I_values(idx)*1000));
    end
    
    subplot(2,3,5);
    if isfield(sweeps, 'flux') && ~isempty(sweeps.flux.T_values)
        plot(sweeps.flux.q_values, sweeps.flux.T_values, '-o', 'LineWidth', 2);
        xlabel('Heat Flux (W/m²)');
        ylabel('T_{max} (°C)');
        title('Temperature vs Heat Flux (COMSOL)');
        grid on;
        yline(100, 'r--', 'Target: 100°C');
        
        % Find max feasible flux
        feasible_idx = find(sweeps.flux.T_values < 100, 1, 'last');
        if ~isempty(feasible_idx)
            hold on;
            plot(sweeps.flux.q_values(feasible_idx), sweeps.flux.T_values(feasible_idx), 'g*', 'MarkerSize', 15);
            text(sweeps.flux.q_values(feasible_idx) + 100, 95, sprintf('Max: %.0f W/m²', sweeps.flux.q_values(feasible_idx)));
        end
    end
    
    subplot(2,3,6);
    if isfield(sweeps, 'thickness') && ~isempty(sweeps.thickness.T_values)
        plot(sweeps.thickness.t_values, sweeps.thickness.T_values, '-o', 'LineWidth', 2);
        xlabel('TEC Thickness (µm)');
        ylabel('T_{max} (°C)');
        title('Temperature vs Thickness (COMSOL)');
        grid on;
        
        % Find optimal thickness
        [T_min, idx] = min(sweeps.thickness.T_values);
        hold on;
        plot(sweeps.thickness.t_values(idx), T_min, 'r*', 'MarkerSize', 15);
        text(sweeps.thickness.t_values(idx) + 20, T_min, sprintf('Optimal: %d µm', sweeps.thickness.t_values(idx)));
    end
    
    saveas(gcf, fullfile(output_dir, 'comsol_validation.png'));
    saveas(gcf, fullfile(output_dir, 'comsol_validation.fig'));
    
    %% Export data to CSV
    % Results table
    T = table([results.id]', {results.label}', [results.T_matlab_C]', [results.T_comsol_C]', ...
        [results.error_percent]', [results.elapsed_time]', ...
        'VariableNames', {'ID', 'Label', 'T_MATLAB_C', 'T_COMSOL_C', 'Error_Percent', 'Time_s'});
    writetable(T, fullfile(output_dir, 'validation_results.csv'));
    
    % Sweep data
    if isfield(sweeps, 'current')
        T_current = table(sweeps.current.I_values', sweeps.current.T_values', ...
            'VariableNames', {'Current_A', 'Temperature_C'});
        writetable(T_current, fullfile(output_dir, 'sweep_current.csv'));
    end
    
    if isfield(sweeps, 'flux')
        T_flux = table(sweeps.flux.q_values', sweeps.flux.T_values', ...
            'VariableNames', {'HeatFlux_Wm2', 'Temperature_C'});
        writetable(T_flux, fullfile(output_dir, 'sweep_flux.csv'));
    end
    
    if isfield(sweeps, 'thickness')
        T_thick = table(sweeps.thickness.t_values', sweeps.thickness.T_values', ...
            'VariableNames', {'Thickness_um', 'Temperature_C'});
        writetable(T_thick, fullfile(output_dir, 'sweep_thickness.csv'));
    end
    
    %% Text report
    fid = fopen(fullfile(output_dir, 'validation_report.txt'), 'w');
    fprintf(fid, 'COMSOL Validation Report\n');
    fprintf(fid, '========================\n\n');
    fprintf(fid, 'Date: %s\n\n', datestr(now));
    
    fprintf(fid, '--- SUMMARY ---\n');
    fprintf(fid, 'Candidates validated: %d\n', length(results));
    fprintf(fid, 'Mean prediction error: %.1f%%\n', mean([results.error_percent]));
    fprintf(fid, 'Max prediction error: %.1f%%\n', max([results.error_percent]));
    fprintf(fid, '\n');
    
    [~, best_idx] = min([results.T_comsol_C]);
    fprintf(fid, '--- BEST DESIGN (COMSOL) ---\n');
    fprintf(fid, 'Configuration: %s\n', results(best_idx).label);
    fprintf(fid, 'T_max (COMSOL): %.1f°C\n', results(best_idx).T_comsol_C);
    fprintf(fid, 'T_max (MATLAB): %.1f°C\n', results(best_idx).T_matlab_C);
    fprintf(fid, '\n');
    
    fprintf(fid, '--- PARAMETRIC SWEEP RESULTS ---\n');
    if isfield(sweeps, 'current')
        [T_min, idx] = min(sweeps.current.T_values);
        fprintf(fid, 'Optimal current: %.0f mA (T = %.1f°C)\n', sweeps.current.I_values(idx)*1000, T_min);
    end
    if isfield(sweeps, 'flux')
        feasible_idx = find(sweeps.flux.T_values < 100, 1, 'last');
        if ~isempty(feasible_idx)
            fprintf(fid, 'Max feasible flux (T<100°C): %.0f W/m²\n', sweeps.flux.q_values(feasible_idx));
        end
    end
    if isfield(sweeps, 'thickness')
        [T_min, idx] = min(sweeps.thickness.T_values);
        fprintf(fid, 'Optimal thickness: %d µm (T = %.1f°C)\n', sweeps.thickness.t_values(idx), T_min);
    end
    
    fclose(fid);
    
    fprintf('Report saved to: %s\n', fullfile(output_dir, 'validation_report.txt'));
    fprintf('Data exported to CSV files\n');
end
