classdef ParametricSweeper < handle
    properties
        Solver
        BaseConfig
        OutputDir
        IndividualSubfolder
    end
    
    methods
        function obj = ParametricSweeper(config_path)
            % Initialize solver and config
            obj.Solver = RadialTECSolver(config_path);
            obj.BaseConfig = obj.Solver.Config;
            
            % Setup output directory
            timestamp = datestr(now, 'yyyy-mm-dd_HH-MM-SS');
            obj.OutputDir = fullfile('output', 'sweeps', timestamp);
            if ~exist(obj.OutputDir, 'dir')
                mkdir(obj.OutputDir);
            end
            
            obj.IndividualSubfolder = 'individual results';
            individual_path = fullfile(obj.OutputDir, obj.IndividualSubfolder);
            if ~exist(individual_path, 'dir')
                mkdir(individual_path);
            end

            % Redirect solver results to this folder temporarily or just use it for our plots
            % The solver's internal ResultsManager has its own folder. 
            % We can update it to point to our sweep folder so all logs go there.
            obj.Solver.Results.OutputDir = obj.OutputDir;
        end
        
        function run_sweep(obj, param_name, start_val, end_val, num_steps)
            fprintf('Starting Parametric Sweep for %s...\n', param_name);
            
            values = linspace(start_val, end_val, num_steps);
            results_T_cold = zeros(size(values));
            results_Q_in = zeros(size(values));
            
            N_stages = obj.Solver.Geometry.N_stages;
            chip_profiles = NaN(length(values), N_stages + 1);
            tec_profiles = NaN(length(values), N_stages);
            
            % Prepare CSV log
            log_file = fullfile(obj.OutputDir, 'sweep_data.csv');
            fid = fopen(log_file, 'w');
            fprintf(fid, '%s,T_cold_Si,Q_in\n', param_name);
            fclose(fid);
            
            for i = 1:length(values)
                val = values(i);
                fprintf('Step %d/%d: %s = %f\n', i, num_steps, param_name, val);
                
                % Update Config
                obj.update_config(param_name, val);
                
                try
                    max_iter = obj.Solver.Config.simulation.max_iterations;
                    tol = obj.Solver.Config.simulation.tolerance;
                    T_current = ones(2*obj.Solver.Geometry.N_stages+1, 1) * 300;
                    
                    for iter = 1:max_iter
                        [T_new, ~, Q_in] = obj.Solver.Network.solve(T_current);
                        if norm(T_new - T_current, inf) < tol
                            T_current = T_new;
                            break;
                        end
                        T_current = 0.5 * T_new + 0.5 * T_current;
                    end
                    
                    T_cold_Si = T_current(1);
                    results_T_cold(i) = T_cold_Si;
                    results_Q_in(i) = Q_in;
                    
                    % Log to CSV
                    fid = fopen(log_file, 'a');
                    fprintf(fid, '%f,%f,%f\n', val, T_cold_Si, Q_in);
                    fclose(fid);
                    
                    % Save Temperature Profile Plot
                    val_tag = regexprep(sprintf('%.4g', val), '\\.', 'p');
                    val_tag = regexprep(val_tag, '-', 'm');
                    plot_name = sprintf('%s_%s.png', param_name, val_tag);
                    title_suffix = sprintf('%s = %.4g, T_0 = %.2f K', param_name, val, T_cold_Si);
                    obj.Solver.Results.plot_temperature_profile(T_current, obj.Solver.Geometry, plot_name, title_suffix, obj.IndividualSubfolder);
                    
                    chip_profiles(i, :) = [T_current(1); T_current(2:N_stages+1)]';
                    tec_profiles(i, :) = T_current(N_stages+2:end)';
                catch ME
                    fprintf('Error at %s = %f: %s\n', param_name, val, ME.message);
                    results_T_cold(i) = NaN;
                end
            end
            
            % Plot Summary and comparisons
            obj.plot_summary(values, results_T_cold, param_name);
            obj.plot_layer_comparison(values, chip_profiles, param_name, 'Chip Layer', 0:N_stages, 'chip_layer_comparison.png');
            obj.plot_layer_comparison(values, tec_profiles, param_name, 'TEC Layer', 1:N_stages, 'tec_layer_comparison.png');
            
            fprintf('Sweep Complete. Results saved to %s\n', obj.OutputDir);
        end
        
        function update_config(obj, param_name, val)
            % Handle nested fields or specific known parameters
            if strcmp(param_name, 'current')
                obj.Solver.Config.operating_conditions.I_current_A = val;
            elseif strcmp(param_name, 'k_r')
                obj.Solver.Config.geometry.radial_expansion_factor = val;
            elseif isfield(obj.Solver.Config.geometry, param_name)
                obj.Solver.Config.geometry.(param_name) = val;
            elseif isfield(obj.Solver.Config.operating_conditions, param_name)
                obj.Solver.Config.operating_conditions.(param_name) = val;
            else
                % Try to find it recursively or error
                error('Parameter %s not found in config structure.', param_name);
            end
            
            % Re-init geometry and network
            obj.Solver.Geometry = TECGeometry(obj.Solver.Config);
            obj.Solver.Network = ThermalNetwork(obj.Solver.Geometry, obj.Solver.Materials, obj.Solver.Config);
        end
        
        function plot_summary(obj, x, y, param_name)
            h = figure('Visible', 'off');
            plot(x, y, '-o', 'LineWidth', 2);
            grid on;
            xlabel(param_name);
            ylabel('T_{cold} (Si Node 0) [K]');
            title(['Parametric Sweep: ' param_name]);
            
            saveas(h, fullfile(obj.OutputDir, 'sweep_summary.png'));
            close(h);
        end

        function plot_layer_comparison(obj, values, profiles, param_name, layer_label, x_values, filename)
            h = figure('Visible', 'off');
            hold on;
            colors = lines(size(profiles, 1));
            legend_entries = {};

            for idx = 1:size(profiles, 1)
                row = profiles(idx, :);
                if all(isnan(row))
                    continue;
                end
                cidx = mod(idx-1, size(colors, 1)) + 1;
                plot(x_values, row, '-o', 'LineWidth', 1.5, 'Color', colors(cidx, :));
                legend_entries{end+1} = sprintf('%s = %.3g', param_name, values(idx));
            end

            if ~isempty(legend_entries)
                legend(legend_entries, 'Interpreter', 'none', 'Location', 'best');
            end
            xlabel('Stage Index');
            ylabel('Temperature (K)');
            title(sprintf('%s Comparison (%s)', layer_label, param_name));
            grid on;

            saveas(h, fullfile(obj.OutputDir, filename));
            close(h);
        end
    end
end
