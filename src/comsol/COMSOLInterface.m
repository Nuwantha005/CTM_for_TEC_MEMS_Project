classdef COMSOLInterface < handle
    % COMSOLINTERFACE - Interface for COMSOL LiveLink with MATLAB
    %
    % This class provides methods to:
    %   1. Connect to COMSOL server
    %   2. Load/modify models
    %   3. Update parameters from optimization results
    %   4. Run simulations
    %   5. Extract results
    %
    % Requirements:
    %   - COMSOL Multiphysics 6.3 with LiveLink for MATLAB
    %   - COMSOL server running (or will start automatically)
    %
    % Usage:
    %   comsol = COMSOLInterface();
    %   comsol.connect();
    %   comsol.loadModel('path/to/model.mph');
    %   comsol.setParameters(params);
    %   comsol.runStudy();
    %   results = comsol.extractResults();
    
    properties
        Model           % COMSOL model object
        ModelPath       % Path to .mph file
        IsConnected     % Connection status
        ServerPort      % COMSOL server port
        Results         % Simulation results
        OutputDir       % Output directory for results
    end
    
    properties (Constant)
        % Parameter name mapping: MATLAB -> COMSOL
        ParamMap = struct(...
            'theta_deg', 'LL_theta', ...
            't_TEC_um', 'LL_t_TEC', ...
            't_chip_um', 'LL_t_chip', ...
            't_SOI_um', 'LL_t_SOI', ...
            'R_cyl_um', 'LL_R_cyl', ...
            'L_1_um', 'LL_L_1', ...
            'k_r', 'LL_k_r', ...
            'w_ic_um', 'LL_w_ic', ...
            'w_oc_um', 'LL_w_oc', ...
            't_ic_um', 'LL_t_ic', ...
            't_oc_um', 'LL_t_oc', ...
            'w_is_um', 'LL_w_is', ...
            'w_az_um', 'LL_w_az', ...
            'R_TSV_um', 'LL_R_TSV', ...
            'P_TSV_um', 'LL_P_TSV', ...
            'beta_ic_deg', 'LL_beta_ic', ...
            'beta_oc_deg', 'LL_beat_oc', ...
            'r_chip_mm', 'LL_r_chip', ...
            'I_current_A', 'I0', ...
            'q_flux_W_m2', 'q' ...
        );
    end
    
    methods
        function obj = COMSOLInterface(varargin)
            % Constructor
            % Optional: COMSOLInterface('Port', 2036)
            
            p = inputParser;
            addParameter(p, 'Port', 2036);
            addParameter(p, 'OutputDir', 'output/comsol_results');
            parse(p, varargin{:});
            
            obj.ServerPort = p.Results.Port;
            obj.IsConnected = false;
            obj.OutputDir = p.Results.OutputDir;
            
            if ~exist(obj.OutputDir, 'dir')
                mkdir(obj.OutputDir);
            end
        end
        
        function success = connect(obj)
            % Connect to COMSOL server
            % Returns true if successful
            
            success = false;
            
            try
                fprintf('Connecting to COMSOL server on port %d...\n', obj.ServerPort);
                
                % Check if mphstart exists (LiveLink installed)
                if ~exist('mphstart', 'file')
                    error('COMSOL LiveLink for MATLAB not found. Please ensure it is installed and on the MATLAB path.');
                end
                
                % Try to connect
                mphstart(obj.ServerPort);
                
                % Import COMSOL API
                import com.comsol.model.*
                import com.comsol.model.util.*
                
                obj.IsConnected = true;
                fprintf('Successfully connected to COMSOL server.\n');
                success = true;
                
            catch ME
                fprintf('Failed to connect to COMSOL: %s\n', ME.message);
                fprintf('\nTroubleshooting:\n');
                fprintf('1. Start COMSOL server: comsolmphserver -port %d\n', obj.ServerPort);
                fprintf('2. Or start COMSOL with LiveLink: comsol mphserver matlab\n');
                fprintf('3. Ensure LiveLink for MATLAB is licensed\n');
            end
        end
        
        function success = loadModel(obj, modelPath)
            % Load a COMSOL model file
            
            success = false;
            
            if ~obj.IsConnected
                error('Not connected to COMSOL. Call connect() first.');
            end
            
            if ~exist(modelPath, 'file')
                error('Model file not found: %s', modelPath);
            end
            
            try
                fprintf('Loading model: %s\n', modelPath);
                obj.Model = mphload(modelPath);
                obj.ModelPath = modelPath;
                fprintf('Model loaded successfully.\n');
                
                % Display current parameters
                obj.displayParameters();
                
                success = true;
            catch ME
                fprintf('Failed to load model: %s\n', ME.message);
            end
        end
        
        function displayParameters(obj)
            % Display current model parameters
            
            if isempty(obj.Model)
                error('No model loaded.');
            end
            
            try
                params = mphgetexpressions(obj.Model.param);
                fprintf('\nCurrent Model Parameters:\n');
                fprintf('%-20s %-20s %s\n', 'Name', 'Value', 'Description');
                fprintf('%s\n', repmat('-', 1, 60));
                
                names = fieldnames(params);
                for i = 1:length(names)
                    name = names{i};
                    value = params.(name);
                    fprintf('%-20s %-20s\n', name, value);
                end
            catch ME
                fprintf('Could not retrieve parameters: %s\n', ME.message);
            end
        end
        
        function setParameters(obj, params)
            % Set model parameters from a struct
            %
            % params can be:
            %   - A struct with field names matching COMSOL parameters
            %   - Or field names matching MATLAB convention (will be mapped)
            
            if isempty(obj.Model)
                error('No model loaded.');
            end
            
            fprintf('\nUpdating parameters...\n');
            
            paramNames = fieldnames(params);
            for i = 1:length(paramNames)
                matlabName = paramNames{i};
                value = params.(matlabName);
                
                % Check if mapping exists
                if isfield(obj.ParamMap, matlabName)
                    comsolName = obj.ParamMap.(matlabName);
                else
                    comsolName = matlabName;  % Use as-is
                end
                
                % Format value with units
                valueStr = obj.formatValue(matlabName, value);
                
                try
                    obj.Model.param.set(comsolName, valueStr);
                    fprintf('  %s = %s\n', comsolName, valueStr);
                catch ME
                    fprintf('  Warning: Could not set %s: %s\n', comsolName, ME.message);
                end
            end
        end
        
        function setParametersDirect(obj, paramCell)
            % Set parameters directly using name-value cell array
            % paramCell = {'LL_theta', '30[deg]', 'LL_t_TEC', '50[um]', ...}
            
            if isempty(obj.Model)
                error('No model loaded.');
            end
            
            fprintf('\nUpdating parameters directly...\n');
            
            for i = 1:2:length(paramCell)
                name = paramCell{i};
                value = paramCell{i+1};
                
                try
                    obj.Model.param.set(name, value);
                    fprintf('  %s = %s\n', name, value);
                catch ME
                    fprintf('  Warning: Could not set %s: %s\n', name, ME.message);
                end
            end
        end
        
        function valueStr = formatValue(~, paramName, value)
            % Format value with appropriate units for COMSOL
            
            % Determine units based on parameter name
            if contains(paramName, '_deg') || contains(paramName, 'theta') || contains(paramName, 'beta')
                valueStr = sprintf('%g[deg]', value);
            elseif contains(paramName, '_um')
                valueStr = sprintf('%g[um]', value);
            elseif contains(paramName, '_mm')
                valueStr = sprintf('%g[mm]', value);
            elseif contains(paramName, '_A') || contains(paramName, 'current')
                valueStr = sprintf('%g[A]', value);
            elseif contains(paramName, '_W_m2') || contains(paramName, 'flux')
                valueStr = sprintf('%g[W/m^2]', value);
            elseif contains(paramName, 'k_r')
                valueStr = sprintf('%g', value);  % Dimensionless
            else
                valueStr = sprintf('%g', value);
            end
        end
        
        function success = runStudy(obj, studyTag)
            % Run a COMSOL study
            % studyTag: optional, defaults to 'std1'
            
            if nargin < 2
                studyTag = 'std1';
            end
            
            success = false;
            
            if isempty(obj.Model)
                error('No model loaded.');
            end
            
            try
                fprintf('\nRunning study "%s"...\n', studyTag);
                tic;
                
                obj.Model.study(studyTag).run();
                
                elapsed = toc;
                fprintf('Study completed in %.1f seconds.\n', elapsed);
                success = true;
                
            catch ME
                fprintf('Study failed: %s\n', ME.message);
            end
        end
        
        function results = extractResults(obj, varargin)
            % Extract results from the solved model
            %
            % Optional parameters:
            %   'Dataset': dataset tag (default: 'dset1')
            %   'Expressions': cell array of expressions to evaluate
            
            p = inputParser;
            addParameter(p, 'Dataset', 'dset1');
            addParameter(p, 'Expressions', {'T', 'minop1(T)', 'maxop1(T)'});
            parse(p, varargin{:});
            
            if isempty(obj.Model)
                error('No model loaded.');
            end
            
            results = struct();
            
            try
                fprintf('\nExtracting results...\n');
                
                % Get temperature field statistics
                % Note: These expressions depend on your COMSOL model setup
                
                % Try to get max/min temperatures
                try
                    % Global evaluation for min/max
                    T_data = mphglobal(obj.Model, {'minop1(T)', 'maxop1(T)'}, 'dataset', p.Results.Dataset);
                    results.T_min_K = T_data(1);
                    results.T_max_K = T_data(2);
                    results.T_min_C = T_data(1) - 273.15;
                    results.T_max_C = T_data(2) - 273.15;
                    
                    fprintf('  T_min = %.2f K (%.2f °C)\n', results.T_min_K, results.T_min_C);
                    fprintf('  T_max = %.2f K (%.2f °C)\n', results.T_max_K, results.T_max_C);
                catch
                    fprintf('  Could not extract T_min/T_max via global operators.\n');
                    fprintf('  Trying alternative method...\n');
                    
                    % Alternative: evaluate on domain
                    try
                        T_vals = mpheval(obj.Model, 'T', 'dataset', p.Results.Dataset);
                        results.T_min_K = min(T_vals.d1);
                        results.T_max_K = max(T_vals.d1);
                        results.T_min_C = results.T_min_K - 273.15;
                        results.T_max_C = results.T_max_K - 273.15;
                        
                        fprintf('  T_min = %.2f K (%.2f °C)\n', results.T_min_K, results.T_min_C);
                        fprintf('  T_max = %.2f K (%.2f °C)\n', results.T_max_K, results.T_max_C);
                    catch ME2
                        fprintf('  Alternative method also failed: %s\n', ME2.message);
                    end
                end
                
                % Try to get heat flux integrals
                try
                    % These depend on your model's boundary integration operators
                    Q_data = mphglobal(obj.Model, {'intop1(ht.ntflux)'}, 'dataset', p.Results.Dataset);
                    results.Q_total = Q_data(1);
                    fprintf('  Q_total = %.4f W\n', results.Q_total);
                catch
                    % Heat flux extraction failed, skip
                end
                
                obj.Results = results;
                
            catch ME
                fprintf('Result extraction failed: %s\n', ME.message);
                results.error = ME.message;
            end
        end
        
        function exportPlot(obj, plotTag, filename, varargin)
            % Export a plot to file
            %
            % plotTag: COMSOL plot tag (e.g., 'pg1')
            % filename: output filename (will be placed in OutputDir)
            
            p = inputParser;
            addParameter(p, 'Format', 'png');
            addParameter(p, 'Resolution', 300);
            parse(p, varargin{:});
            
            if isempty(obj.Model)
                error('No model loaded.');
            end
            
            filepath = fullfile(obj.OutputDir, filename);
            
            try
                fprintf('Exporting plot to: %s\n', filepath);
                mphplot(obj.Model, plotTag, 'rangenum', 1);
                
                % Save current figure
                saveas(gcf, filepath);
                fprintf('Plot saved.\n');
                
            catch ME
                fprintf('Plot export failed: %s\n', ME.message);
            end
        end
        
        function saveModel(obj, filename)
            % Save the current model to file
            
            if isempty(obj.Model)
                error('No model loaded.');
            end
            
            if nargin < 2
                [~, name, ~] = fileparts(obj.ModelPath);
                filename = fullfile(obj.OutputDir, [name '_modified.mph']);
            end
            
            try
                fprintf('Saving model to: %s\n', filename);
                mphsave(obj.Model, filename);
                fprintf('Model saved.\n');
            catch ME
                fprintf('Save failed: %s\n', ME.message);
            end
        end
        
        function disconnect(obj)
            % Disconnect from COMSOL server
            
            if obj.IsConnected
                try
                    % Clear model
                    if ~isempty(obj.Model)
                        ModelUtil.remove(obj.Model.tag);
                    end
                    obj.Model = [];
                    obj.IsConnected = false;
                    fprintf('Disconnected from COMSOL.\n');
                catch
                    % Ignore errors during disconnect
                end
            end
        end
    end
    
    methods (Static)
        function checkInstallation()
            % Check if COMSOL LiveLink is properly installed
            
            fprintf('Checking COMSOL LiveLink installation...\n\n');
            
            % Check for mphstart
            if exist('mphstart', 'file')
                fprintf('✓ mphstart found\n');
            else
                fprintf('✗ mphstart NOT found\n');
                fprintf('  Add COMSOL LiveLink to MATLAB path:\n');
                fprintf('  addpath(''C:\\Program Files\\COMSOL\\COMSOL63\\Multiphysics\\mli'')\n');
            end
            
            % Check for mphload
            if exist('mphload', 'file')
                fprintf('✓ mphload found\n');
            else
                fprintf('✗ mphload NOT found\n');
            end
            
            % Check COMSOL installation directory
            comsolPaths = {
                'C:\Program Files\COMSOL\COMSOL63\Multiphysics',
                'C:\Program Files\COMSOL\COMSOL62\Multiphysics',
                'C:\Program Files\COMSOL\COMSOL61\Multiphysics'
            };
            
            found = false;
            for i = 1:length(comsolPaths)
                if exist(comsolPaths{i}, 'dir')
                    fprintf('✓ COMSOL found at: %s\n', comsolPaths{i});
                    found = true;
                    break;
                end
            end
            
            if ~found
                fprintf('✗ COMSOL installation not found in standard locations\n');
            end
            
            fprintf('\nTo start COMSOL server:\n');
            fprintf('  1. Open Command Prompt as Administrator\n');
            fprintf('  2. Run: comsolmphserver -port 2036\n');
            fprintf('  OR\n');
            fprintf('  Start COMSOL Desktop with: File > Client/Server > Start Server\n');
        end
    end
end
