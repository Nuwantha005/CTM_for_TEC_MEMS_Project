classdef StageSweeper < handle
    properties
        ConfigPath
        OutputDir
    end
    
    methods
        function obj = StageSweeper(config_path)
            obj.ConfigPath = config_path;
            
            % Setup output directory
            timestamp = datestr(now, 'yyyy-mm-dd_HH-MM-SS');
            obj.OutputDir = fullfile('output', 'stage_sweeps', timestamp);
            if ~exist(obj.OutputDir, 'dir')
                mkdir(obj.OutputDir);
            end
        end
        
        function run_stage_sweep(obj, min_stages, max_stages)
            fprintf('Starting Stage Optimization Sweep (%d to %d stages)...\n', min_stages, max_stages);
            
            stages = min_stages:max_stages;
            opt_temps_K = zeros(size(stages));
            opt_temps_C = zeros(size(stages));
            all_x_opt = cell(size(stages));
            actual_var_names = {};  % Will be populated from first successful optimization
            
            % Summary file path (header will be written after first successful run)
            summary_file = fullfile(obj.OutputDir, 'stage_sweep_summary.csv');
            header_written = false;
            
            for i = 1:length(stages)
                N = stages(i);
                fprintf('\n');
                fprintf('╔════════════════════════════════════════════════════════════╗\n');
                fprintf('║         OPTIMIZING FOR N = %d STAGES                        ║\n', N);
                fprintf('╚════════════════════════════════════════════════════════════╝\n');
                
                % 1. Read base config
                base_config = jsondecode(fileread(obj.ConfigPath));
                
                % 2. Modify N_stages
                base_config.geometry.N_stages = N;
                
                % 3. Create stage-specific directory and config file
                stage_output_dir = fullfile(obj.OutputDir, sprintf('stage_%d', N));
                if ~exist(stage_output_dir, 'dir')
                    mkdir(stage_output_dir);
                end
                temp_config_path = fullfile(stage_output_dir, sprintf('config_stage_%d.json', N));
                fid = fopen(temp_config_path, 'w');
                fprintf(fid, '%s', jsonencode(base_config));
                fclose(fid);
                
                % 4. Run Optimizer
                optimizer = TECOptimizer(temp_config_path, stage_output_dir);
                
                try
                    [x_opt, fval, T_final] = optimizer.run_optimization();
                    
                    % fval is T_max in Kelvin from the optimizer
                    % But we should verify using T_final which is the actual result
                    T_max_from_fval = fval;  % This is what optimizer returned
                    T_max_from_T_final = max(T_final);  % This is from final solver run
                    
                    % Use the T_final value as it's the actual simulated result
                    opt_temps_K(i) = T_max_from_T_final;
                    opt_temps_C(i) = T_max_from_T_final - 273.15;
                    all_x_opt{i} = x_opt;
                    
                    % Write CSV header on first successful run (uses actual var count)
                    if ~header_written
                        % Read var names from optimizer's log file
                        opt_log = fullfile(stage_output_dir, 'optimization_log.csv');
                        if exist(opt_log, 'file')
                            log_data = readtable(opt_log);
                            actual_var_names = log_data.Properties.VariableNames(3:end);  % Skip Iter, Cost
                        else
                            actual_var_names = arrayfun(@(k) sprintf('var%d', k), 1:length(x_opt), 'UniformOutput', false);
                        end
                        fid = fopen(summary_file, 'w');
                        fprintf(fid, 'Stages,T_max_K,T_max_C');
                        for v = 1:length(actual_var_names)
                            fprintf(fid, ',%s', actual_var_names{v});
                        end
                        fprintf(fid, '\n');
                        fclose(fid);
                        header_written = true;
                    end
                    
                    fprintf('\n--- Stage %d Results ---\n', N);
                    fprintf('  Optimizer fval:  %.2f K (%.1f °C)\n', T_max_from_fval, T_max_from_fval - 273.15);
                    fprintf('  T_final max:     %.2f K (%.1f °C)\n', T_max_from_T_final, T_max_from_T_final - 273.15);
                    if abs(T_max_from_fval - T_max_from_T_final) > 1
                        fprintf('  WARNING: Mismatch between optimizer fval and T_final!\n');
                    end
                    
                    % Log to summary CSV
                    fid = fopen(summary_file, 'a');
                    fprintf(fid, '%d,%.4f,%.2f', N, opt_temps_K(i), opt_temps_C(i));
                    for k = 1:length(x_opt)
                        fprintf(fid, ',%.6f', x_opt(k));
                    end
                    fprintf(fid, '\n');
                    fclose(fid);
                    
                catch ME
                    fprintf('Optimization failed for N=%d: %s\n', N, ME.message);
                    opt_temps_K(i) = NaN;
                    opt_temps_C(i) = NaN;
                    all_x_opt{i} = [];
                end
            end
            
            obj.plot_results(stages, opt_temps_K, opt_temps_C);
            obj.save_summary(stages, opt_temps_K, opt_temps_C, all_x_opt, actual_var_names);
            fprintf('\n✓ Stage Sweep Complete. Results in %s\n', obj.OutputDir);
        end
        
        function plot_results(obj, stages, temps_K, temps_C)
            % Plot in Kelvin
            h1 = figure('Visible', 'off');
            plot(stages, temps_K, '-s', 'LineWidth', 2, 'MarkerFaceColor', 'b', 'MarkerSize', 8);
            grid on;
            xlabel('Number of Stages');
            ylabel('Optimal T_{max} (K)');
            title('Optimum Temperature vs Number of Stages');
            
            % Find global minimum
            [min_T, idx] = min(temps_K);
            best_N = stages(idx);
            
            hold on;
            plot(best_N, min_T, 'rp', 'MarkerSize', 15, 'MarkerFaceColor', 'r');
            text(best_N, min_T, sprintf('  Best: N=%d, T=%.1f K (%.1f°C)', best_N, min_T, min_T-273.15), ...
                'VerticalAlignment', 'bottom', 'FontSize', 10);
            
            % Add reference lines
            yline(373.15, 'r--', '100°C', 'LineWidth', 1.5);
            yline(358.15, 'g--', '85°C Target', 'LineWidth', 1.5);
            
            saveas(h1, fullfile(obj.OutputDir, 'stage_optimization_curve_K.png'));
            close(h1);
            
            % Plot in Celsius
            h2 = figure('Visible', 'off');
            plot(stages, temps_C, '-s', 'LineWidth', 2, 'MarkerFaceColor', 'b', 'MarkerSize', 8);
            grid on;
            xlabel('Number of Stages');
            ylabel('Optimal T_{max} (°C)');
            title('Optimum Temperature vs Number of Stages');
            
            hold on;
            plot(best_N, temps_C(idx), 'rp', 'MarkerSize', 15, 'MarkerFaceColor', 'r');
            text(best_N, temps_C(idx), sprintf('  Best: N=%d, T=%.1f°C', best_N, temps_C(idx)), ...
                'VerticalAlignment', 'bottom', 'FontSize', 10);
            
            % Add reference lines
            yline(100, 'r--', '100°C', 'LineWidth', 1.5);
            yline(85, 'g--', '85°C Target', 'LineWidth', 1.5);
            
            saveas(h2, fullfile(obj.OutputDir, 'stage_optimization_curve_C.png'));
            close(h2);
            
            fprintf('Stage optimization plots saved.\n');
        end
        
        function save_summary(obj, stages, temps_K, temps_C, all_x_opt, var_names)
            % Save a human-readable summary
            fid = fopen(fullfile(obj.OutputDir, 'summary.txt'), 'w');
            fprintf(fid, 'Stage Sweep Summary\n');
            fprintf(fid, '===================\n\n');
            
            [min_T, idx] = min(temps_K);
            fprintf(fid, 'Best Configuration:\n');
            fprintf(fid, '  Stages: %d\n', stages(idx));
            fprintf(fid, '  T_max: %.2f K (%.1f °C)\n\n', min_T, min_T - 273.15);
            
            fprintf(fid, 'All Results:\n');
            fprintf(fid, '%-8s  %-12s  %-12s\n', 'Stages', 'T_max (K)', 'T_max (°C)');
            fprintf(fid, '----------------------------------------\n');
            for i = 1:length(stages)
                fprintf(fid, '%-8d  %-12.2f  %-12.1f\n', stages(i), temps_K(i), temps_C(i));
            end
            
            fprintf(fid, '\nOptimal Parameters for Best Stage Count (%d):\n', stages(idx));
            fprintf(fid, '----------------------------------------------\n');
            if ~isempty(all_x_opt{idx})
                x_best = all_x_opt{idx};
                for v = 1:min(length(var_names), length(x_best))
                    fprintf(fid, '  %-30s: %.6f\n', var_names{v}, x_best(v));
                end
            end
            
            fclose(fid);
        end
    end
end
