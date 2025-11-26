%% RUN OVERNIGHT COMSOL HIGH-FIDELITY SIMULATIONS
% This script runs COMSOL simulations on candidate designs
% from the preliminary optimization. Run this overnight.
%
% Requirements:
%   1. COMSOL 6.3 with LiveLink for MATLAB
%   2. Completed preliminary optimization (feasible_designs.csv)
%   3. COMSOL template model (linked to SolidWorks)
%
% Usage:
%   1. Start COMSOL server: comsolmphserver -port 2036
%      OR start MATLAB through COMSOL: comsol mphserver matlab
%   2. Modify CONFIG section below with your paths
%   3. Run this script and leave computer overnight
%   4. Results saved to output/high_fidelity/

clear; clc;
addpath(genpath('src'));
addpath('F:\EngineeringSoftware\COMSOL\COMSOL63\Multiphysics\mli');

%% ==================== CONFIG (MODIFY THIS) ====================

% Path to preliminary optimization results directory
OPTIMIZATION_DIR = 'output/optimizations/2025-11-25_22-56-47';  % Change to your latest

% Path to COMSOL template model (linked to SolidWorks)
COMSOL_MODEL = 'E:\Semester 7\ME4311 - MicroNano Electro Mechanical Systems and Nanotechnology\Project\COMSOL\MATLAB_API\template_3_stage.mph';

% Maximum number of candidates to simulate (to limit overnight time)
MAX_CANDIDATES = 20;  % ~5-15 min each, 20 candidates = ~2-5 hours

% Heat flux for simulations (W/m^2) - should match your optimization
Q_FLUX = 500;  % Check your config file

% Output directory
OUTPUT_DIR = fullfile('output', 'high_fidelity', datestr(now, 'yyyy-mm-dd_HH-MM-SS'));

%% ==================== SETUP ====================

fprintf('╔════════════════════════════════════════════════════════╗\n');
fprintf('║     OVERNIGHT COMSOL HIGH-FIDELITY SIMULATIONS        ║\n');
fprintf('╚════════════════════════════════════════════════════════╝\n\n');

% Check for COMSOL model
if ~exist(COMSOL_MODEL, 'file')
    fprintf('ERROR: COMSOL model not found:\n  %s\n\n', COMSOL_MODEL);
    fprintf('Please set COMSOL_MODEL path in this script.\n');
    return;
end

fprintf('COMSOL Model: %s\n\n', COMSOL_MODEL);

% Create output directory
if ~exist(OUTPUT_DIR, 'dir')
    mkdir(OUTPUT_DIR);
end

% Create subfolder for saved COMSOL models
MODELS_DIR = fullfile(OUTPUT_DIR, 'comsol_models');
if ~exist(MODELS_DIR, 'dir')
    mkdir(MODELS_DIR);
end
fprintf('COMSOL models will be saved to: %s\n\n', MODELS_DIR);

%% ==================== LOAD CANDIDATE DESIGNS FROM CSV ====================

fprintf('Loading candidate designs from preliminary optimization...\n');

% Load from CSV file (feasible_designs.csv)
csvFile = fullfile(OPTIMIZATION_DIR, 'feasible_designs.csv');
if ~exist(csvFile, 'file')
    error('Results file not found: %s\nMake sure OPTIMIZATION_DIR points to a valid optimization output folder.', csvFile);
end

% Read CSV - columns: Rank,N_stages,theta_deg,t_TEC_um,k_r,I_opt_mA,T_max_C,COP,Q_in_W,Q_out_W
candidatesTable = readtable(csvFile);
fprintf('  Loaded %d feasible designs from CSV\n', height(candidatesTable));

% Convert table to struct array for easier handling
candidates = struct();
for i = 1:height(candidatesTable)
    candidates(i).rank = candidatesTable.Rank(i);
    candidates(i).N_stages = candidatesTable.N_stages(i);
    candidates(i).theta_deg = candidatesTable.theta_deg(i);
    candidates(i).t_TEC_um = candidatesTable.t_TEC_um(i);
    candidates(i).k_r = candidatesTable.k_r(i);
    candidates(i).I_opt_mA = candidatesTable.I_opt_mA(i);
    candidates(i).T_max_C = candidatesTable.T_max_C(i);
    candidates(i).COP = candidatesTable.COP(i);
end

% Limit number of candidates (already sorted by T_max in CSV)
nCandidates = min(length(candidates), MAX_CANDIDATES);
candidates = candidates(1:nCandidates);

fprintf('  Will simulate top %d designs\n\n', nCandidates);

% Display candidate summary
fprintf('Candidate designs to simulate:\n');
fprintf('┌──────┬───────┬───────────┬──────────┬──────────┬──────────┬───────────────┐\n');
fprintf('│  #   │   N   │ θ (deg)   │ t_TEC(µm)│  k_r     │ I (mA)   │ T_max (°C)    │\n');
fprintf('├──────┼───────┼───────────┼──────────┼──────────┼──────────┼───────────────┤\n');
for i = 1:min(10, nCandidates)
    c = candidates(i);
    fprintf('│ %4d │ %5d │ %9.1f │ %8.0f │ %8.2f │ %8.2f │ %13.2f │\n', ...
        i, c.N_stages, c.theta_deg, c.t_TEC_um, c.k_r, c.I_opt_mA, c.T_max_C);
end
if nCandidates > 10
    fprintf('│  ... │   ... │       ... │      ... │      ... │      ... │           ... │\n');
    c = candidates(end);
    fprintf('│ %4d │ %5d │ %9.1f │ %8.0f │ %8.2f │ %8.2f │ %13.2f │\n', ...
        nCandidates, c.N_stages, c.theta_deg, c.t_TEC_um, c.k_r, c.I_opt_mA, c.T_max_C);
end
fprintf('└──────┴───────┴───────────┴──────────┴──────────┴──────────┴───────────────┘\n\n');

%% ==================== ESTIMATE TIME ====================

% Estimate based on typical COMSOL thermal simulation times
EST_TIME_PER_SIM = 10;  % minutes (adjust based on your model complexity)

totalMinutes = nCandidates * EST_TIME_PER_SIM;
hours = floor(totalMinutes / 60);
mins = mod(totalMinutes, 60);

fprintf('Estimated runtime: %d hours %d minutes\n', hours, mins);
fprintf('Expected completion: %s\n\n', datestr(now + totalMinutes/1440, 'HH:MM on mmm dd'));

%% ==================== CHECK COMSOL CONNECTION ====================

fprintf('Checking COMSOL LiveLink...\n');

if ~exist('mphstart', 'file')
    fprintf('\nERROR: COMSOL LiveLink not found in MATLAB path.\n');
    fprintf('Either:\n');
    fprintf('  1. Add LiveLink to path:\n');
    fprintf('     addpath(''C:\\Program Files\\COMSOL\\COMSOL63\\Multiphysics\\mli'')\n\n');
    fprintf('  2. Or start MATLAB through COMSOL:\n');
    fprintf('     comsol mphserver matlab\n');
    return;
end

fprintf('✓ LiveLink found.\n\n');

%% ==================== CONNECT TO COMSOL ====================

fprintf('Connecting to COMSOL server...\n');

try
    mphstart(2036);
    import com.comsol.model.*
    import com.comsol.model.util.*
    fprintf('✓ Connected to COMSOL server.\n\n');
catch ME
    fprintf('\nERROR: Could not connect to COMSOL server.\n');
    fprintf('  %s\n\n', ME.message);
    fprintf('Start COMSOL server first:\n');
    fprintf('  comsolmphserver -port 2036\n');
    return;
end

%% ==================== LOAD MODEL ====================

fprintf('Loading COMSOL model...\n');

try
    model = mphload(COMSOL_MODEL);
    fprintf('✓ Model loaded successfully.\n\n');
catch ME
    fprintf('\nERROR: Could not load model.\n');
    fprintf('  %s\n\n', ME.message);
    return;
end

%% ==================== RUN COMSOL SIMULATIONS ====================

fprintf('Starting COMSOL simulations...\n');
fprintf('Press Ctrl+C to abort.\n\n');

% Initialize results storage
results = struct();
results.candidates = candidates;
results.comsol_T_max = zeros(nCandidates, 1);
results.matlab_T_max = zeros(nCandidates, 1);
results.error_C = zeros(nCandidates, 1);
results.sim_time = zeros(nCandidates, 1);
results.success = false(nCandidates, 1);

startTime = tic;

for i = 1:nCandidates
    c = candidates(i);
    
    fprintf('\n═══════════════════════════════════════════════════════\n');
    fprintf('Candidate %d/%d\n', i, nCandidates);
    fprintf('═══════════════════════════════════════════════════════\n');
    fprintf('  N_stages = %d\n', c.N_stages);
    fprintf('  theta    = %.1f deg\n', c.theta_deg);
    fprintf('  t_TEC    = %.0f um\n', c.t_TEC_um);
    fprintf('  k_r      = %.2f\n', c.k_r);
    fprintf('  I_opt    = %.2f mA\n', c.I_opt_mA);
    fprintf('  MATLAB prediction: T_max = %.2f °C\n\n', c.T_max_C);
    
    try
        % Set parameters in COMSOL
        fprintf('Setting COMSOL parameters...\n');
        model.param.set('LL_theta', sprintf('%g[deg]', c.theta_deg));
        model.param.set('LL_t_TEC', sprintf('%g[um]', c.t_TEC_um));
        model.param.set('LL_k_r', sprintf('%g', c.k_r));
        model.param.set('I0', sprintf('%g[A]', c.I_opt_mA / 1000));  % mA to A
        model.param.set('q', sprintf('%g[W/m^2]', Q_FLUX));
        
        % Run simulation
        fprintf('Running simulation...\n');
        simStart = tic;
        model.study('std1').run();
        simTime = toc(simStart);
        
        fprintf('✓ Simulation completed in %.1f seconds.\n', simTime);
        
        % Extract temperature
        try
            T_data = mpheval(model, 'T', 'dataset', 'dset1');
            T_max_K = max(T_data.d1);
            T_max_C = T_max_K - 273.15;
            fprintf('  COMSOL result: T_max = %.2f °C\n', T_max_C);
        catch
            % Try alternative method
            try
                T_max_C = mphglobal(model, 'maxop1(T)', 'unit', 'degC');
            catch
                T_max_C = NaN;
                fprintf('  Warning: Could not extract T_max\n');
            end
        end
        
        % Store results
        results.comsol_T_max(i) = T_max_C;
        results.matlab_T_max(i) = c.T_max_C;
        results.error_C(i) = T_max_C - c.T_max_C;
        results.sim_time(i) = simTime;
        results.success(i) = true;
        
        fprintf('  Error vs MATLAB: %.2f °C (%.1f%%)\n', ...
            results.error_C(i), 100*results.error_C(i)/c.T_max_C);
        
        % Save intermediate result (MATLAB data)
        resultFile = fullfile(OUTPUT_DIR, sprintf('result_%02d.mat', i));
        candidate = c;
        comsol_T_max = T_max_C;
        save(resultFile, 'candidate', 'comsol_T_max', 'simTime');
        
        % Save COMSOL model with solution for later analysis
        % Filename format: candidate_XX_N#_theta##_tTEC###.mph
        modelFileName = sprintf('candidate_%02d_N%d_theta%.0f_tTEC%.0f.mph', ...
            i, c.N_stages, c.theta_deg, c.t_TEC_um);
        modelFilePath = fullfile(MODELS_DIR, modelFileName);
        fprintf('  Saving COMSOL model: %s\n', modelFileName);
        try
            mphsave(model, modelFilePath);
            fprintf('  ✓ Model saved successfully.\n');
        catch ME_save
            fprintf('  Warning: Could not save model: %s\n', ME_save.message);
        end
        
    catch ME
        fprintf('✗ Simulation FAILED: %s\n', ME.message);
        results.success(i) = false;
    end
    
    % Progress update
    elapsed = toc(startTime);
    avgTime = elapsed / i;
    eta = avgTime * (nCandidates - i);
    fprintf('\nProgress: %d/%d complete - Elapsed: %.0f s - ETA: %.0f s\n', ...
        i, nCandidates, elapsed, eta);
end

%% ==================== SAVE FINAL RESULTS ====================

% Save all results
save(fullfile(OUTPUT_DIR, 'all_results.mat'), 'results');

% Save comparison table
compTable = table(...
    (1:nCandidates)', ...
    [candidates.N_stages]', ...
    [candidates.theta_deg]', ...
    [candidates.t_TEC_um]', ...
    [candidates.k_r]', ...
    results.matlab_T_max, ...
    results.comsol_T_max, ...
    results.error_C, ...
    results.success, ...
    'VariableNames', {'Rank', 'N_stages', 'theta_deg', 't_TEC_um', 'k_r', ...
                      'MATLAB_T_max_C', 'COMSOL_T_max_C', 'Error_C', 'Success'});
writetable(compTable, fullfile(OUTPUT_DIR, 'comparison_results.csv'));

%% ==================== SUMMARY ====================

elapsedTime = toc(startTime);
nSuccess = sum(results.success);

fprintf('\n');
fprintf('════════════════════════════════════════════════════════\n');
fprintf('SIMULATION COMPLETE\n');
fprintf('════════════════════════════════════════════════════════\n');
fprintf('Total runtime: %.1f minutes (%.2f hours)\n', elapsedTime/60, elapsedTime/3600);
fprintf('Simulations completed: %d / %d\n', nSuccess, nCandidates);
fprintf('Results saved to: %s\n', OUTPUT_DIR);
fprintf('COMSOL models saved to: %s\n\n', MODELS_DIR);

if nSuccess > 0
    successIdx = results.success;
    fprintf('Error Statistics (COMSOL - MATLAB):\n');
    fprintf('  Mean error: %.2f °C\n', mean(results.error_C(successIdx)));
    fprintf('  Std error:  %.2f °C\n', std(results.error_C(successIdx)));
    fprintf('  Max error:  %.2f °C\n', max(abs(results.error_C(successIdx))));
    
    % Find best design according to COMSOL
    [~, bestIdx] = min(results.comsol_T_max(successIdx));
    validIdx = find(successIdx);
    bestIdx = validIdx(bestIdx);
    
    fprintf('\nBest design (by COMSOL):\n');
    fprintf('  Candidate #%d: N=%d, θ=%.0f°, t_TEC=%.0fμm\n', ...
        bestIdx, candidates(bestIdx).N_stages, candidates(bestIdx).theta_deg, ...
        candidates(bestIdx).t_TEC_um);
    fprintf('  COMSOL T_max = %.2f °C\n', results.comsol_T_max(bestIdx));
    fprintf('  MATLAB T_max = %.2f °C\n', results.matlab_T_max(bestIdx));
end

fprintf('════════════════════════════════════════════════════════\n');
fprintf('Done!\n');
