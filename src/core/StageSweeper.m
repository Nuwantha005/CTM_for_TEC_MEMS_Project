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
            opt_temps = zeros(size(stages));
            
            % Prepare summary log
            summary_file = fullfile(obj.OutputDir, 'stage_sweep_summary.csv');
            fid = fopen(summary_file, 'w');
            % We don't know the exact param names order from Optimizer easily without running it, 
            % but we know the order in TECOptimizer.m. 
            % Let's assume: current, k_r, ic_ratio, oc_ratio, ic_ang_ratio, oc_ang_ratio, fill
            fprintf(fid, 'Stages,Min_Temp_K,Opt_Current,Opt_Kr,Opt_IC_Ratio,Opt_OC_Ratio,Opt_IC_Ang,Opt_OC_Ang,Opt_Fill\n');
            fclose(fid);
            
            for i = 1:length(stages)
                N = stages(i);
                fprintf('\n=== Optimizing for N = %d Stages ===\n', N);
                
                % Create a temporary config file or modify the object?
                % TECOptimizer loads from file. Let's modify the file or 
                % better, modify the object after loading.
                % But TECOptimizer constructor loads the file.
                % We can create a temp config file for this stage.
                
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
                    [x_opt, fval, ~] = optimizer.run_optimization();
                    
                    opt_temps(i) = fval;
                    
                    % Log to summary
                    fid = fopen(summary_file, 'a');
                    fprintf(fid, '%d,%.4f', N, fval);
                    for k = 1:length(x_opt)
                        fprintf(fid, ',%.4f', x_opt(k));
                    end
                    fprintf(fid, '\n');
                    fclose(fid);
                    
                catch ME
                    fprintf('Optimization failed for N=%d: %s\n', N, ME.message);
                    opt_temps(i) = NaN;
                end
                
                % Clean up temp config (optional, maybe keep for debug)
                % delete(temp_config_path);
            end
            
            obj.plot_results(stages, opt_temps);
            fprintf('Stage Sweep Complete. Results in %s\n', obj.OutputDir);
        end
        
        function plot_results(obj, stages, temps)
            h = figure('Visible', 'off');
            plot(stages, temps, '-s', 'LineWidth', 2, 'MarkerFaceColor', 'b');
            grid on;
            xlabel('Number of Stages');
            ylabel('Optimal Cold Side Temperature (K)');
            title('Optimum Temperature vs Number of Stages');
            
            % Find global minimum
            [min_T, idx] = min(temps);
            best_N = stages(idx);
            
            hold on;
            plot(best_N, min_T, 'rp', 'MarkerSize', 12, 'MarkerFaceColor', 'r');
            text(best_N, min_T, sprintf('  Best: N=%d, T=%.2fK', best_N, min_T), ...
                'VerticalAlignment', 'bottom');
            
            saveas(h, fullfile(obj.OutputDir, 'stage_optimization_curve.png'));
            close(h);
        end
    end
end
