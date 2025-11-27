classdef TECOptimizer < handle
    properties
        Solver
        BaseConfig
        OutputDir
        % Live plot handles
        TempProfileFig
        TempProfileAxes
        LastT  % Store last temperature distribution
    end
    
    methods
        function obj = TECOptimizer(config_path, output_dir_override)
            obj.Solver = RadialTECSolver(config_path);
            obj.BaseConfig = obj.Solver.Config;
            obj.LastT = [];
            if nargin > 1 && ~isempty(output_dir_override)
                obj.OutputDir = output_dir_override;
            else
                timestamp = datestr(now, 'yyyy-mm-dd_HH-MM-SS');
                obj.OutputDir = fullfile('output', 'optimizations', timestamp);
            end
            if ~exist(obj.OutputDir, 'dir')
                mkdir(obj.OutputDir);
            end
        end
        
        function [x_opt, fval, T_final] = run_optimization(obj)
            fprintf('Starting Optimization...\n');
            
            % Show current heat flux setting
            q_flux = obj.Solver.Config.boundary_conditions.q_flux_W_m2;
            fprintf('Heat flux: %.0f W/m² (%.2f W total chip)\n', q_flux, q_flux * 100e-6);
            
            vars = {
            % name, lower_bound, upper_bound, initial_value
            'current', 0.0005, 1.0, 0.025;  % Lower current often better
            'k_r', 0.2, 3.0, 1.15;
            'interconnect_ratio', 0.1, 0.35, 0.15;
            'outerconnect_ratio', 0.1, 0.35, 0.15;
            'interconnect_angle_ratio', 0.1, 0.4, 0.16;
            'outerconnect_angle_ratio', 0.1, 0.4, 0.16;
            'fill_factor', 0.8, 0.99, 0.95;
            'thickness_um', 50, 1000, 200;  % Thicker TEC helps
            'wedge_angle_deg', 10, 90, 30;
            'insulation_width_ratio', 0.02, 0.1, 0.04;
            'interconnect_thickness_ratio', 0.5, 2.0, 1.0;
            'outerconnect_thickness_ratio', 0.5, 2.0, 1.0
            };
            var_names = vars(:, 1);
            lb = [vars{:, 2}];
            ub = [vars{:, 3}];
            x0 = [vars{:, 4}];
            
            % Initialize live plots
            obj.init_live_plots();
            
            % Use optimoptions with live plot callback
            options = optimoptions('fmincon', 'Display', 'iter', ...
            'MaxFunctionEvaluations', 500, ...
            'OptimalityTolerance', 1e-4, ...
            'StepTolerance', 1e-6, ...
            'PlotFcn', {@optimplotfval, @optimplotx}, ...
            'OutputFcn', @(x,optimValues,state) obj.update_live_plots(x, optimValues, state, var_names));
            
            if exist('fmincon', 'file')
            fprintf('Using fmincon...\n');
            [x_opt, fval] = fmincon(@(x) obj.evaluate_cost(x, var_names), x0, [], [], [], [], lb, ub, [], options);
            else
            fprintf('Using fminsearch (bounds handled by penalty)...\n');
            opts = optimset('Display', 'iter', 'MaxFunEvals', 500, 'PlotFcns', @optimplotfval);
            [x_opt, fval] = fminsearch(@(x) obj.evaluate_cost_penalty(x, var_names, lb, ub), x0, opts);
            end
            fprintf('Optimization Complete.\n');
            fprintf('Optimal T_max: %.4f K (%.1f °C)\n', fval, fval - 273.15);
            
            % Print optimal parameters
            fprintf('\nOptimal Parameters:\n==============================\n');
            for i = 1:length(var_names)
            fprintf('%s: %.4f\n', var_names{i}, x_opt(i));
            end
            fprintf('==============================\n');
            
            obj.update_config(x_opt, var_names);
            [T_final, ~, ~] = obj.Solver.run();
            obj.log_optimization_results(x_opt, var_names, fval, T_final);
            
            % Close live plot figure
            if ~isempty(obj.TempProfileFig) && isvalid(obj.TempProfileFig)
            close(obj.TempProfileFig);
            end
        end
        
        function cost = evaluate_cost(obj, x, var_names)
            obj.update_config(x, var_names);
            T_water = obj.Solver.Config.boundary_conditions.T_water_K;
            T_target = 373.15;  % 100°C target
            T_current = ones(2*obj.Solver.Geometry.N_stages+1, 1) * (T_water + 50);
            
            for i = 1:50
                try
                    [T_new, ~, ~] = obj.Solver.Network.solve(T_current);
                catch
                    cost = 1e8;  % Solver failed
                    obj.log_step(x, var_names, cost);
                    return;
                end
                
                % Check for invalid values
                if any(isnan(T_new)) || any(isinf(T_new)) || any(T_new < 0)
                    cost = 1e8;
                    obj.log_step(x, var_names, cost);
                    return;
                end
                
                T_current = 0.5 * T_new + 0.5 * T_current;
                
                % Early termination if converged
                if max(abs(T_new - T_current)) < 1e-4
                    break;
                end
            end
            
            % Store for live plotting
            obj.LastT = T_current;
            
            % Get temperatures
            T_max = max(T_current);
            T_min = min(T_current);
            
            % Primary objective: minimize max temperature
            cost = T_max;
            
            % Soft penalty for temperatures below water (shouldn't happen but don't reject)
            if T_min < T_water
                cost = cost + 100 * abs(T_min - T_water);
            end
            
            % Add penalty for exceeding target temperature (soft constraint)
            if T_max > T_target
                cost = cost + 10 * (T_max - T_target);
            end
            
            obj.log_step(x, var_names, cost);
        end
        
        function cost = evaluate_cost_penalty(obj, x, var_names, lb, ub)
            if any(x < lb) || any(x > ub)
                cost = 1e6 + sum(abs(x - lb) .* (x < lb)) + sum(abs(x - ub) .* (x > ub));
                return;
            end
            cost = obj.evaluate_cost(x, var_names);
        end
        
        function update_config(obj, x, var_names)
            for i = 1:length(var_names)
                name = var_names{i};
                val = x(i);
                if strcmp(name, 'current')
                    obj.Solver.Config.operating_conditions.I_current_A = val;
                elseif strcmp(name, 'k_r')
                    obj.Solver.Config.geometry.radial_expansion_factor = val;
                elseif strcmp(name, 'thickness_um')
                    obj.Solver.Config.geometry.thickness_um = val;
                elseif strcmp(name, 't_TEC')
                    obj.Solver.Config.geometry.thickness_um = val;
                elseif strcmp(name, 'wedge_angle_deg')
                    obj.Solver.Config.geometry.wedge_angle_deg = val;
                elseif strcmp(name, 'insulation_width_ratio')
                    obj.Solver.Config.geometry.insulation_width_ratio = val;
                elseif strcmp(name, 'interconnect_angle_ratio')
                    obj.Solver.Config.geometry.interconnect_angle_ratio = val;
                elseif strcmp(name, 'outerconnect_angle_ratio')
                    obj.Solver.Config.geometry.outerconnect_angle_ratio = val;
                elseif strcmp(name, 'interconnect_thickness_ratio')
                    obj.Solver.Config.geometry.interconnect_thickness_ratio = val;
                elseif strcmp(name, 'outerconnect_thickness_ratio')
                    obj.Solver.Config.geometry.outerconnect_thickness_ratio = val;
                else
                    obj.Solver.Config.geometry.(name) = val;
                end
            end
            obj.Solver.Geometry = TECGeometry(obj.Solver.Config);
            obj.Solver.Network = ThermalNetwork(obj.Solver.Geometry, obj.Solver.Materials, obj.Solver.Config);
        end
        
        function log_step(obj, x, var_names, cost)
            fname = fullfile(obj.OutputDir, 'optimization_log.csv');
            if ~exist(fname, 'file')
                fid = fopen(fname, 'w');
                fprintf(fid, 'Iter,Cost');
                for i = 1:length(var_names)
                    fprintf(fid, ',%s', var_names{i});
                end
                fprintf(fid, '\n');
                fclose(fid);
            end
            fid = fopen(fname, 'a');
            fprintf(fid, '0,%.4f', cost);
            for i = 1:length(x)
                fprintf(fid, ',%.4f', x(i));
            end
            fprintf(fid, '\n');
            fclose(fid);
        end
        
        function log_optimization_results(obj, x, var_names, cost, T_final)
            fid = fopen(fullfile(obj.OutputDir, 'optimal_params.txt'), 'w');
            fprintf(fid, 'Optimal Cost (T_max): %.4f K (%.1f °C)\n', cost, cost - 273.15);
            fprintf(fid, 'Parameters:\n');
            for i = 1:length(x)
                fprintf(fid, '%s: %.4f\n', var_names{i}, x(i));
            end
            fclose(fid);
            
            % Check if log file exists before trying to read it
            log_file = fullfile(obj.OutputDir, 'optimization_log.csv');
            if ~exist(log_file, 'file')
                fprintf('Warning: No optimization log file found. Skipping convergence plot.\n');
                return;
            end
            data = readtable(log_file);
            h = figure('Visible', 'off');
            plot(data.Cost, '-o');
            xlabel('Evaluation Index');
            ylabel('Objective (T_{cold})');
            title('Optimization Convergence');
            saveas(h, fullfile(obj.OutputDir, 'convergence.png'));
            close(h);
            for i = 1:length(var_names)
                name = var_names{i};
                h = figure('Visible', 'off');
                scatter(data.(name), data.Cost);
                xlabel(name);
                ylabel('Cost');
                title(['Impact of ' name]);
                saveas(h, fullfile(obj.OutputDir, ['scatter_' name '.png']));
                close(h);
            end
            
            % Save temperature profile plot to optimization output folder
            if nargin >= 5 && ~isempty(T_final)
                T_water = obj.Solver.Config.boundary_conditions.T_water_K;
                obj.plot_temperature_profile_to_folder(T_final, T_water);
            end
        end
        
        function plot_temperature_profile_to_folder(obj, T, T_water)
            % Plot and save temperature profile to the optimization output folder
            N = obj.Solver.Geometry.N_stages;
            T_0 = T(1);
            T_Si = T(2:N+1);
            T_c = T(N+2:end);

            r_chip = 0:N;
            T_chip = [T_0; T_Si];
            
            % TEC layer: stages 1 to N, plus water at "stage N+1"
            r_tec = 1:(N+1);
            T_tec = [T_c; T_water];

            h = figure('Visible', 'off');
            plot(r_chip, T_chip, '-ob', 'LineWidth', 1.8, 'MarkerFaceColor', 'b', 'DisplayName', 'Silicon Layer'); hold on;
            plot(r_tec, T_tec, '-sr', 'LineWidth', 1.8, 'MarkerFaceColor', 'r', 'DisplayName', 'TEC Layer');
            
            % Mark the water temperature point distinctly
            plot(N+1, T_water, 'g^', 'MarkerSize', 10, 'MarkerFaceColor', 'g', 'DisplayName', 'Coolant (T_{water})');
            
            xlabel('Stage Index');
            ylabel('Temperature (K)');
            title('Optimal Temperature Profile');
            legend('Location', 'best');
            grid on;

            saveas(h, fullfile(obj.OutputDir, 'optimal_temperature_profile.png'));
            close(h);
            fprintf('Temperature profile saved to %s\n', fullfile(obj.OutputDir, 'optimal_temperature_profile.png'));
        end
        
        function init_live_plots(obj)
            % Initialize live temperature profile plot
            obj.TempProfileFig = figure('Name', 'Live Temperature Profile', ...
                'Position', [100, 100, 600, 400]);
            obj.TempProfileAxes = axes(obj.TempProfileFig);
            title(obj.TempProfileAxes, 'Current Best Temperature Distribution');
            xlabel(obj.TempProfileAxes, 'Radial Position (Stage Index)');
            ylabel(obj.TempProfileAxes, 'Temperature (K)');
            grid(obj.TempProfileAxes, 'on');
            hold(obj.TempProfileAxes, 'on');
            
            % Add target temperature line
            yline(obj.TempProfileAxes, 373.15, 'r--', 'T_{target} = 100°C', 'LineWidth', 1.5);
            yline(obj.TempProfileAxes, 300, 'b--', 'T_{water} = 27°C', 'LineWidth', 1.5);
            
            drawnow;
        end
        
        function stop = update_live_plots(obj, x, optimValues, state, var_names)
            stop = false;
            
            if strcmp(state, 'done')
                return;
            end
            
            % Only update every few iterations to avoid slowdown
            if mod(optimValues.iteration, 2) ~= 0 && ~strcmp(state, 'init')
                return;
            end
            
            % Get current temperature distribution
            if ~isempty(obj.LastT)
                T = obj.LastT;
                N = obj.Solver.Geometry.N_stages;
                T_water = obj.Solver.Config.boundary_conditions.T_water_K;
                
                % Extract temperatures
                T_0 = T(1);
                T_Si = T(2:N+1);
                T_c = T(N+2:end);
                
                % Radial positions
                r_chip = 0:N;
                T_chip = [T_0; T_Si];
                r_tec = 1:(N+1);
                T_tec = [T_c; T_water];
                
                % Clear and replot
                if isvalid(obj.TempProfileAxes)
                    cla(obj.TempProfileAxes);
                    hold(obj.TempProfileAxes, 'on');
                    
                    % Plot temperature profiles
                    plot(obj.TempProfileAxes, r_chip, T_chip - 273.15, '-ob', 'LineWidth', 2, ...
                        'MarkerFaceColor', 'b', 'DisplayName', 'Silicon Layer');
                    plot(obj.TempProfileAxes, r_tec, T_tec - 273.15, '-sr', 'LineWidth', 2, ...
                        'MarkerFaceColor', 'r', 'DisplayName', 'TEC Cold Side');
                    plot(obj.TempProfileAxes, N+1, T_water - 273.15, 'g^', 'MarkerSize', 12, ...
                        'MarkerFaceColor', 'g', 'DisplayName', 'Coolant');
                    
                    % Reference lines
                    yline(obj.TempProfileAxes, 100, 'r--', 'LineWidth', 1.5);
                    yline(obj.TempProfileAxes, 27, 'b--', 'LineWidth', 1.5);
                    
                    % Labels
                    title(obj.TempProfileAxes, sprintf('Iter %d: T_{max} = %.1f°C (Target: 100°C)', ...
                        optimValues.iteration, max(T) - 273.15));
                    xlabel(obj.TempProfileAxes, 'Radial Position (Stage Index)');
                    ylabel(obj.TempProfileAxes, 'Temperature (°C)');
                    legend(obj.TempProfileAxes, 'Location', 'best');
                    grid(obj.TempProfileAxes, 'on');
                    
                    % Set axis limits
                    ylim(obj.TempProfileAxes, [20, max(150, max(T) - 273.15 + 10)]);
                    
                    drawnow;
                end
            end
        end
    end
end
