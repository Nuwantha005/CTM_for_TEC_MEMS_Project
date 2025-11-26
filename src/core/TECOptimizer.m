classdef TECOptimizer < handle
    properties
        Solver
        BaseConfig
        OutputDir
    end
    
    methods
        function obj = TECOptimizer(config_path, output_dir_override)
            obj.Solver = RadialTECSolver(config_path);
            obj.BaseConfig = obj.Solver.Config;
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
            vars = {
                'current', 0.01, 1.0, 0.15;
                'k_r', 0.5, 2.0, 1.2;
                'interconnect_ratio', 0.05, 0.45, 0.2;
                'outerconnect_ratio', 0.05, 0.45, 0.2;
                'interconnect_angle_ratio', 0.05, 0.5, 0.16;
                'outerconnect_angle_ratio', 0.05, 0.5, 0.16;
                'fill_factor', 0.5, 0.99, 0.9
            };
            var_names = vars(:, 1);
            lb = [vars{:, 2}];
            ub = [vars{:, 3}];
            x0 = [vars{:, 4}];
            options = optimset('Display', 'iter', 'PlotFcns', @optimplotfval);
            if exist('fmincon', 'file')
                fprintf('Using fmincon...\n');
                [x_opt, fval] = fmincon(@(x) obj.evaluate_cost(x, var_names), x0, [], [], [], [], lb, ub, [], options);
            else
                fprintf('Using fminsearch (bounds handled by penalty)...\n');
                [x_opt, fval] = fminsearch(@(x) obj.evaluate_cost_penalty(x, var_names, lb, ub), x0, options);
            end
            fprintf('Optimization Complete.\n');
            fprintf('Optimal Cost (T_cold): %.4f K\n', fval);
            obj.update_config(x_opt, var_names);
            [T_final, ~, ~] = obj.Solver.run();
            obj.log_optimization_results(x_opt, var_names, fval);
        end
        
        function cost = evaluate_cost(obj, x, var_names)
            obj.update_config(x, var_names);
            T_current = ones(2*obj.Solver.Geometry.N_stages+1, 1) * 300;
            for i = 1:30
                [T_new, ~, ~] = obj.Solver.Network.solve(T_current);
                T_current = 0.5 * T_new + 0.5 * T_current;
            end
            cost = T_current(1);
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
                elseif strcmp(name, 't_TEC')
                    obj.Solver.Config.geometry.thickness_um = val;
                elseif strcmp(name, 'interconnect_angle_ratio')
                    obj.Solver.Config.geometry.interconnect_angle_ratio = val;
                elseif strcmp(name, 'outerconnect_angle_ratio')
                    obj.Solver.Config.geometry.outerconnect_angle_ratio = val;
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
        
        function log_optimization_results(obj, x, var_names, cost)
            fid = fopen(fullfile(obj.OutputDir, 'optimal_params.txt'), 'w');
            fprintf(fid, 'Optimal Cost (T_cold): %.4f K\n', cost);
            fprintf(fid, 'Parameters:\n');
            for i = 1:length(x)
                fprintf(fid, '%s: %.4f\n', var_names{i}, x(i));
            end
            fclose(fid);
            data = readtable(fullfile(obj.OutputDir, 'optimization_log.csv'));
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
        end
    end
end
