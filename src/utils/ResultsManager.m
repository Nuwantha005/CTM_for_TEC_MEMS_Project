classdef ResultsManager
    % RESULTSMANAGER Handles I/O, saving results, and plotting.

    properties
        OutputDir
    end

    methods
        function obj = ResultsManager(config)
            % Constructor: Creates output directory with timestamp
            timestamp = datestr(now, 'yyyy-mm-dd_HH-MM-SS');
            base_dir = 'output';
            obj.OutputDir = fullfile(base_dir, timestamp);

            if ~exist(obj.OutputDir, 'dir')
                mkdir(obj.OutputDir);
            end

            % Save config snapshot
            fid = fopen(fullfile(obj.OutputDir, 'config_snapshot.json'), 'w');
            fprintf(fid, '%s', jsonencode(config));
            fclose(fid);
        end

        function save_results(obj, T, Q_out, config)
            % SAVE_RESULTS Saves simulation data to .mat file
            filename = fullfile(obj.OutputDir, 'results.mat');
            save(filename, 'T', 'Q_out', 'config');
            fprintf('Results saved to %s\n', filename);
        end

        function save_results_text(obj, T, Q_out, config)
            % SAVE_RESULTS_TEXT Saves results to a readable text file.
            filename = fullfile(obj.OutputDir, 'results.txt');
            fid = fopen(filename, 'w');

            fprintf(fid, 'Radial TEC Simulation Results\n');
            fprintf(fid, 'Date: %s\n', datestr(now));
            fprintf(fid, '--------------------------------\n');
            fprintf(fid, 'q_flux: %.2f W/m2\n', config.boundary_conditions.q_flux_W_m2);
            fprintf(fid, 'I_current: %.2f A\n', config.operating_conditions.I_current_A);
            fprintf(fid, 'Q_out (Calculated): %.2e W\n', Q_out);
            fprintf(fid, '--------------------------------\n');
            fprintf(fid, 'Node 0 (Center) Temp: %.2f K\n', T(1));
            fprintf(fid, '--------------------------------\n');
            fprintf(fid, 'Stage\tT_Si (K)\tT_c (K)\n');

            N = config.geometry.N_stages;
            for i = 1:N
                idx_Si = i + 1;
                idx_c = N + i + 1;
                fprintf(fid, '%d\t%.2f\t%.2f\n', i, T(idx_Si), T(idx_c));
            end

            fclose(fid);
            fprintf('Text results saved to %s\n', filename);
        end

        function plot_temperature_profile(obj, T, geometry)
            % PLOT_TEMPERATURE_PROFILE Generates and saves a plot

            N = geometry.N_stages;
            % Extract vectors (Node 0 is index 1)
            T_0 = T(1);
            T_Si = T(2:N+1);
            T_c = T(N+2:end);

            % Create radial coordinate vector (simplified)
            % In reality, get actual radii from geometry
            r = 1:N;

            h = figure('Visible', 'off');
            plot(0, T_0, 'ko', 'MarkerFaceColor', 'k'); hold on;
            plot(r, T_Si, 'b-o', 'DisplayName', 'Silicon Layer');
            plot(r, T_c, 'r-s', 'DisplayName', 'TEC Layer');
            xlabel('Stage Index');
            ylabel('Temperature (K)');
            title('Radial Temperature Profile');
            legend('Location', 'best');
            grid on;

            saveas(h, fullfile(obj.OutputDir, 'temp_profile.png'));
            close(h);
        end
    end
end
