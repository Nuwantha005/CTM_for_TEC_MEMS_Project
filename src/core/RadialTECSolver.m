classdef RadialTECSolver < handle
    properties
        Config
        Geometry
        Materials
        Network
        Results
    end

    methods
        function obj = RadialTECSolver(config_path)
            text = fileread(config_path);
            obj.Config = jsondecode(text);

            obj.Materials = MaterialProperties(obj.Config);
            obj.Geometry = TECGeometry(obj.Config);
            obj.Network = ThermalNetwork(obj.Geometry, obj.Materials, obj.Config);
            obj.Results = ResultsManager(obj.Config);
        end

        function [T_current, Q_out, Q_in] = run(obj)
            fprintf('Starting Radial TEC Simulation...\n');
            
            % Iteration parameters
            max_iter = obj.Config.simulation.max_iterations;
            tol = obj.Config.simulation.tolerance;
            
            N = obj.Geometry.N_stages;
            dim = 2*N + 1;
            
            % Initial guess
            T_current = ones(dim, 1) * obj.Config.simulation.T_initial_guess;
            
            convergence_history = [];
            
            fprintf('Iterative Solver Started (Max Iter: %d, Tol: %e)\n', max_iter, tol);
            
            for iter = 1:max_iter
                % Solve step
                [T_new, Q_out, Q_in] = obj.Network.solve(T_current);
                
                % Check convergence
                diff = norm(T_new - T_current, inf);
                convergence_history(end+1) = diff;
                
                fprintf('Iter %d: Max Diff = %e\n', iter, diff);
                
                if diff < tol
                    fprintf('Converged in %d iterations.\n', iter);
                    T_current = T_new; % Update final T
                    break;
                end
                
                % Relaxation (optional, but good for stability)
                alpha = 0.5; 
                T_current = alpha * T_new + (1 - alpha) * T_current;
            end
            
            if iter == max_iter
                fprintf('Warning: Maximum iterations reached without full convergence.\n');
            end

            obj.Results.log_results(T_current, Q_out, Q_in, obj.Config);
            obj.Results.save_results_text(T_current, Q_out, Q_in, obj.Config);
            obj.Results.save_results(T_current, Q_out, obj.Config);
            obj.Results.plot_temperature_profile(T_current, obj.Geometry, 'best_temperature_profile.png', 'Final Solution');
            
            % Save convergence plot
            obj.plot_convergence(convergence_history);

            fprintf('Simulation Complete.\n');
        end
        
        function plot_convergence(obj, history)
            h = figure('Visible', 'off');
            semilogy(1:length(history), history, '-o', 'LineWidth', 2);
            grid on;
            xlabel('Iteration');
            ylabel('Max Temperature Change (K)');
            title('Solver Convergence');
            
            filename = fullfile(obj.Results.OutputDir, 'convergence_plot.png');
            saveas(h, filename);
            close(h);
            fprintf('Convergence plot saved to %s\n', filename);
        end
    end
end
