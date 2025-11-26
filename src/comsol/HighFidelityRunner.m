classdef HighFidelityRunner < handle
    % HIGHFIDELITYRUNNER - Automated high-fidelity COMSOL simulations
    %
    % Takes candidates from preliminary optimization and runs them through
    % COMSOL for validation. Designed for overnight batch runs.
    %
    % Workflow:
    %   1. Load optimization results
    %   2. Select top N candidates
    %   3. For each candidate:
    %      a. Update COMSOL parameters
    %      b. Run simulation
    %      c. Extract and save results
    %   4. Compare MATLAB vs COMSOL results
    %   5. Generate validation report
    
    properties
        COMSOLInterface     % COMSOL connection
        Candidates          % Candidate designs to simulate
        COMSOLResults       % Results from COMSOL
        ComparisonReport    % Comparison with MATLAB predictions
        OutputDir           % Output directory
        ModelPath           % Path to COMSOL model
        
        % Default geometry parameters (from SW_equations.txt)
        DefaultParams
    end
    
    methods
        function obj = HighFidelityRunner(modelPath, varargin)
            % Constructor
            %
            % modelPath: path to COMSOL .mph file
            
            p = inputParser;
            addParameter(p, 'OutputDir', '');
            parse(p, varargin{:});
            
            obj.ModelPath = modelPath;
            
            % Create output directory
            if isempty(p.Results.OutputDir)
                timestamp = datestr(now, 'yyyy-mm-dd_HH-MM-SS');
                obj.OutputDir = fullfile('output', 'comsol_validation', timestamp);
            else
                obj.OutputDir = p.Results.OutputDir;
            end
            
            if ~exist(obj.OutputDir, 'dir')
                mkdir(obj.OutputDir);
            end
            
            % Set default parameters from SW_equations.txt
            obj.DefaultParams = struct(...
                't_chip_um', 50, ...
                't_SOI_um', 100, ...
                't_TEC_um', 50, ...
                'theta_deg', 30, ...
                'r_chip_mm', 10, ...
                'R_cyl_um', 1000, ...      % 1mm = 1000um
                'L_1_um', 1000, ...        % 1mm = 1000um
                'beta_oc_deg', 10, ...
                'beta_ic_deg', 5, ...
                'w_ic_um', 50, ...
                'w_oc_um', 50, ...
                't_oc_um', 20, ...
                't_ic_um', 20, ...
                'w_is_um', 40, ...
                'w_az_um', 30, ...
                'R_TSV_um', 10, ...
                'P_TSV_um', 30, ...
                'k_r', 1.2 ...
            );
            
            obj.Candidates = [];
            obj.COMSOLResults = [];
        end
        
        function loadCandidates(obj, source, varargin)
            % Load candidate designs from optimization results
            %
            % source can be:
            %   - Path to optimization_results.mat
            %   - Struct array of candidates
            %   - 'latest' to find most recent optimization
            
            p = inputParser;
            addParameter(p, 'TopN', 10);  % Number of top candidates
            addParameter(p, 'FilterFeasible', true);
            parse(p, varargin{:});
            
            if ischar(source)
                if strcmp(source, 'latest')
                    % Find most recent optimization
                    source = obj.findLatestOptimization();
                end
                
                fprintf('Loading candidates from: %s\n', source);
                data = load(source);
                
                if isfield(data, 'results')
                    results = data.results;
                else
                    results = data;
                end
                
                if isfield(results, 'all_results')
                    candidates = results.all_results;
                else
                    error('Could not find results in file.');
                end
            else
                candidates = source;
            end
            
            % Filter feasible designs
            if p.Results.FilterFeasible && ~isempty(candidates)
                meets_target = [candidates.meets_target];
                candidates = candidates(meets_target);
            end
            
            % Sort by T_max and take top N
            if ~isempty(candidates)
                T_max = [candidates.T_max_C];
                [~, sort_idx] = sort(T_max);
                n = min(p.Results.TopN, length(candidates));
                candidates = candidates(sort_idx(1:n));
            end
            
            obj.Candidates = candidates;
            fprintf('Loaded %d candidate designs.\n', length(candidates));
            
            % Display candidates
            obj.displayCandidates();
        end
        
        function path = findLatestOptimization(~)
            % Find the most recent optimization results
            
            optDir = 'output/optimizations';
            if ~exist(optDir, 'dir')
                error('No optimization results found.');
            end
            
            dirs = dir(optDir);
            dirs = dirs([dirs.isdir]);
            dirs = dirs(~ismember({dirs.name}, {'.', '..'}));
            
            if isempty(dirs)
                error('No optimization results found.');
            end
            
            % Sort by date (name is timestamp)
            [~, idx] = sort({dirs.name}, 'descend');
            latestDir = dirs(idx(1)).name;
            
            path = fullfile(optDir, latestDir, 'optimization_results.mat');
            
            if ~exist(path, 'file')
                error('Results file not found: %s', path);
            end
        end
        
        function displayCandidates(obj)
            % Display loaded candidates
            
            if isempty(obj.Candidates)
                fprintf('No candidates loaded.\n');
                return;
            end
            
            fprintf('\n=== CANDIDATE DESIGNS ===\n');
            fprintf('Rank | N  | θ(°) | t_TEC(μm) | k_r  | I(mA) | T_max(°C) | COP\n');
            fprintf('-----|----|----- |-----------|------|-------|-----------|------\n');
            
            for i = 1:length(obj.Candidates)
                c = obj.Candidates(i);
                fprintf('%4d | %2d | %4.0f | %9.0f | %4.2f | %5.1f | %9.1f | %.3f\n', ...
                    i, c.N_stages, c.theta_deg, c.t_TEC_um, c.k_r, ...
                    c.I_opt*1000, c.T_max_C, c.COP);
            end
            fprintf('\n');
        end
        
        function success = initialize(obj)
            % Initialize COMSOL connection and load model
            
            success = false;
            
            % Create COMSOL interface
            obj.COMSOLInterface = COMSOLInterface('OutputDir', obj.OutputDir);
            
            % Connect to COMSOL
            if ~obj.COMSOLInterface.connect()
                fprintf('\nCOMSOL connection failed.\n');
                fprintf('Please ensure COMSOL server is running.\n');
                return;
            end
            
            % Load model
            if ~obj.COMSOLInterface.loadModel(obj.ModelPath)
                fprintf('\nFailed to load COMSOL model.\n');
                return;
            end
            
            success = true;
        end
        
        function results = runAllCandidates(obj, varargin)
            % Run COMSOL simulation for all candidates
            
            p = inputParser;
            addParameter(p, 'SaveIntermediates', true);
            addParameter(p, 'ExportPlots', true);
            parse(p, varargin{:});
            
            if isempty(obj.Candidates)
                error('No candidates loaded. Call loadCandidates() first.');
            end
            
            if ~obj.COMSOLInterface.IsConnected
                error('COMSOL not connected. Call initialize() first.');
            end
            
            nCandidates = length(obj.Candidates);
            results = cell(nCandidates, 1);
            
            fprintf('\n========================================\n');
            fprintf('  HIGH-FIDELITY VALIDATION RUN\n');
            fprintf('========================================\n');
            fprintf('Candidates to simulate: %d\n', nCandidates);
            fprintf('Output directory: %s\n', obj.OutputDir);
            fprintf('----------------------------------------\n\n');
            
            startTime = tic;
            
            for i = 1:nCandidates
                candidate = obj.Candidates(i);
                
                fprintf('\n--- Candidate %d/%d ---\n', i, nCandidates);
                fprintf('N=%d, θ=%.0f°, t_TEC=%.0fμm, k_r=%.2f, I=%.1fmA\n', ...
                    candidate.N_stages, candidate.theta_deg, candidate.t_TEC_um, ...
                    candidate.k_r, candidate.I_opt*1000);
                fprintf('MATLAB prediction: T_max = %.1f°C\n', candidate.T_max_C);
                
                try
                    % Convert candidate to COMSOL parameters
                    params = obj.candidateToParams(candidate);
                    
                    % Set parameters in COMSOL
                    obj.COMSOLInterface.setParameters(params);
                    
                    % Run simulation
                    simStart = tic;
                    if obj.COMSOLInterface.runStudy()
                        simTime = toc(simStart);
                        
                        % Extract results
                        result = obj.COMSOLInterface.extractResults();
                        result.candidate_idx = i;
                        result.simulation_time = simTime;
                        result.matlab_T_max_C = candidate.T_max_C;
                        result.params = params;
                        
                        % Calculate error
                        if isfield(result, 'T_max_C')
                            result.error_C = result.T_max_C - candidate.T_max_C;
                            result.error_pct = 100 * result.error_C / candidate.T_max_C;
                            fprintf('COMSOL result: T_max = %.1f°C (error: %.1f°C, %.1f%%)\n', ...
                                result.T_max_C, result.error_C, result.error_pct);
                        end
                        
                        results{i} = result;
                        
                        % Export plot if requested
                        if p.Results.ExportPlots
                            plotFile = sprintf('candidate_%02d_temperature.png', i);
                            obj.COMSOLInterface.exportPlot('pg1', plotFile);
                        end
                        
                        % Save intermediate results
                        if p.Results.SaveIntermediates
                            obj.saveIntermediateResult(i, result, candidate);
                        end
                        
                    else
                        fprintf('Simulation FAILED for candidate %d\n', i);
                        results{i} = struct('error', 'Simulation failed', 'candidate_idx', i);
                    end
                    
                catch ME
                    fprintf('ERROR for candidate %d: %s\n', i, ME.message);
                    results{i} = struct('error', ME.message, 'candidate_idx', i);
                end
                
                % Progress update
                elapsed = toc(startTime);
                avgTime = elapsed / i;
                eta = avgTime * (nCandidates - i);
                fprintf('Progress: %d/%d - Elapsed: %.0fs - ETA: %.0fs\n', ...
                    i, nCandidates, elapsed, eta);
            end
            
            obj.COMSOLResults = results;
            
            % Generate comparison report
            obj.generateReport();
            
            totalTime = toc(startTime);
            fprintf('\n========================================\n');
            fprintf('VALIDATION COMPLETE\n');
            fprintf('Total time: %.1f minutes\n', totalTime/60);
            fprintf('Results saved to: %s\n', obj.OutputDir);
            fprintf('========================================\n');
        end
        
        function params = candidateToParams(obj, candidate)
            % Convert candidate design to COMSOL parameters
            
            % Start with defaults
            params = obj.DefaultParams;
            
            % Override with candidate values
            params.theta_deg = candidate.theta_deg;
            params.t_TEC_um = candidate.t_TEC_um;
            params.k_r = candidate.k_r;
            params.I_current_A = candidate.I_opt;
            
            % Get q_flux from candidate config if available
            if isfield(candidate, 'config') && isfield(candidate.config, 'boundary_conditions')
                params.q_flux_W_m2 = candidate.config.boundary_conditions.q_flux_W_m2;
            else
                params.q_flux_W_m2 = 500;  % Default
            end
            
            % Note: N_stages affects geometry but may need manual adjustment
            % in SolidWorks/COMSOL if the model is parametric
        end
        
        function saveIntermediateResult(obj, idx, result, candidate)
            % Save intermediate result to file
            
            filename = sprintf('result_candidate_%02d.mat', idx);
            filepath = fullfile(obj.OutputDir, filename);
            save(filepath, 'result', 'candidate');
        end
        
        function generateReport(obj)
            % Generate validation report
            
            if isempty(obj.COMSOLResults)
                fprintf('No results to report.\n');
                return;
            end
            
            % Create summary table
            fprintf('\n========================================\n');
            fprintf('  VALIDATION SUMMARY\n');
            fprintf('========================================\n');
            fprintf('Candidate | MATLAB T(°C) | COMSOL T(°C) | Error (°C) | Error (%%)\n');
            fprintf('----------|--------------|--------------|------------|----------\n');
            
            validResults = [];
            
            for i = 1:length(obj.COMSOLResults)
                r = obj.COMSOLResults{i};
                
                if isfield(r, 'T_max_C') && ~isfield(r, 'error')
                    fprintf('%9d | %12.1f | %12.1f | %10.1f | %9.1f\n', ...
                        i, r.matlab_T_max_C, r.T_max_C, r.error_C, r.error_pct);
                    validResults = [validResults; r];
                else
                    errMsg = 'Unknown';
                    if isfield(r, 'error')
                        errMsg = r.error;
                    end
                    fprintf('%9d | %12.1f | %12s | %10s | %9s\n', ...
                        i, obj.Candidates(i).T_max_C, 'FAILED', errMsg, '-');
                end
            end
            
            % Statistics
            if ~isempty(validResults)
                errors = [validResults.error_C];
                fprintf('\nStatistics:\n');
                fprintf('  Mean error: %.2f °C\n', mean(errors));
                fprintf('  Std error:  %.2f °C\n', std(errors));
                fprintf('  Max error:  %.2f °C\n', max(abs(errors)));
                fprintf('  Correlation: %.4f\n', corr([validResults.matlab_T_max_C]', [validResults.T_max_C]'));
            end
            
            % Save report
            obj.saveReport();
        end
        
        function saveReport(obj)
            % Save report to file
            
            filepath = fullfile(obj.OutputDir, 'validation_report.txt');
            fid = fopen(filepath, 'w');
            
            fprintf(fid, 'COMSOL VALIDATION REPORT\n');
            fprintf(fid, '========================\n\n');
            fprintf(fid, 'Date: %s\n', datestr(now));
            fprintf(fid, 'Model: %s\n\n', obj.ModelPath);
            
            fprintf(fid, 'CANDIDATES VALIDATED\n');
            fprintf(fid, '--------------------\n');
            
            for i = 1:length(obj.Candidates)
                c = obj.Candidates(i);
                fprintf(fid, 'Candidate %d:\n', i);
                fprintf(fid, '  N_stages: %d\n', c.N_stages);
                fprintf(fid, '  theta: %.1f deg\n', c.theta_deg);
                fprintf(fid, '  t_TEC: %.0f um\n', c.t_TEC_um);
                fprintf(fid, '  k_r: %.2f\n', c.k_r);
                fprintf(fid, '  I_opt: %.4f A\n', c.I_opt);
                fprintf(fid, '  MATLAB T_max: %.2f C\n', c.T_max_C);
                
                if i <= length(obj.COMSOLResults)
                    r = obj.COMSOLResults{i};
                    if isfield(r, 'T_max_C')
                        fprintf(fid, '  COMSOL T_max: %.2f C\n', r.T_max_C);
                        fprintf(fid, '  Error: %.2f C (%.1f%%)\n', r.error_C, r.error_pct);
                    else
                        fprintf(fid, '  COMSOL: FAILED\n');
                    end
                end
                fprintf(fid, '\n');
            end
            
            fclose(fid);
            fprintf('Report saved to: %s\n', filepath);
            
            % Also save as .mat
            results = obj.COMSOLResults;
            candidates = obj.Candidates;
            save(fullfile(obj.OutputDir, 'validation_results.mat'), 'results', 'candidates');
        end
        
        function cleanup(obj)
            % Cleanup COMSOL connection
            
            if ~isempty(obj.COMSOLInterface)
                obj.COMSOLInterface.disconnect();
            end
        end
    end
end
