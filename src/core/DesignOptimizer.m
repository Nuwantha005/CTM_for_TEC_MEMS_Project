classdef DesignOptimizer < handle
    % DESIGNOPTIMIZER - Robust multi-level optimization for Radial TEC
    % 
    % Strategy:
    %   Level 1: Grid search over discrete params (N_stages, theta, t_TEC)
    %   Level 2: Optimize current I for each configuration (golden section)
    %   Level 3: Filter feasible designs (T_max < T_target)
    %   Level 4: Rank by COP and figure of merit
    
    properties
        BaseConfig          % Base configuration struct
        Results             % Optimization results
        OutputDir           % Output directory for results
        
        % Optimization bounds
        N_stages_range      % [min, max] number of stages
        Theta_range         % [min, max] wedge angle in degrees
        t_TEC_range         % [min, max] TEC thickness in um
        k_r_range           % [min, max] radial expansion factor
        I_range             % [min, max] current in A
        
        % Targets
        T_target_C          % Target max temperature in Celsius
        T_water_K           % Coolant temperature
        
        % Feasibility thresholds
        MinThermalConductance   % Minimum K to avoid ill-conditioning
        MaxTempRise             % Maximum expected temp rise for pre-check
    end
    
    methods
        function obj = DesignOptimizer(config_path)
            % Load base configuration
            text = fileread(config_path);
            obj.BaseConfig = jsondecode(text);
            
            % Set default optimization ranges
            obj.N_stages_range = [2, 8];
            obj.Theta_range = [15, 60];      % degrees
            obj.t_TEC_range = [30, 100];     % um
            obj.k_r_range = [0.8, 1.5];
            obj.I_range = [0.001, 0.5];      % A
            
            % Targets
            obj.T_target_C = 100;            % Target: keep chip below 100°C
            obj.T_water_K = obj.BaseConfig.boundary_conditions.T_water_K;
            
            % Feasibility thresholds - relaxed to allow more configs through
            obj.MinThermalConductance = 1e-9;  % W/K
            obj.MaxTempRise = 50000;           % K above ambient (very relaxed, let solver decide)
            
            % Create output directory
            timestamp = datestr(now, 'yyyy-mm-dd_HH-MM-SS');
            obj.OutputDir = fullfile('output', 'optimizations', timestamp);
            if ~exist(obj.OutputDir, 'dir')
                mkdir(obj.OutputDir);
            end
            
            obj.Results = struct();
        end
        
        function results = runFullOptimization(obj, varargin)
            % RUNFULLOPTIMIZATION - Main optimization routine
            %
            % Optional parameters:
            %   'N_stages_list' - specific values to test [2,3,4,5]
            %   'Theta_list'    - specific angles to test [20,30,45]
            %   'Verbose'       - print progress (default: true)
            
            p = inputParser;
            addParameter(p, 'N_stages_list', 2:6);
            addParameter(p, 'Theta_list', [15, 20, 25, 30, 45]);
            addParameter(p, 't_TEC_list', [40, 50, 60, 80]);
            addParameter(p, 'k_r_list', [0.9, 1.0, 1.1, 1.2]);
            addParameter(p, 'Verbose', true);
            parse(p, varargin{:});
            
            N_list = p.Results.N_stages_list;
            theta_list = p.Results.Theta_list;
            t_list = p.Results.t_TEC_list;
            kr_list = p.Results.k_r_list;
            verbose = p.Results.Verbose;
            
            % Calculate total configurations
            total_configs = length(N_list) * length(theta_list) * length(t_list) * length(kr_list);
            
            fprintf('\n========================================\n');
            fprintf('  RADIAL TEC DESIGN OPTIMIZATION\n');
            fprintf('========================================\n');
            fprintf('Target: T_max < %d°C (%.1f K)\n', obj.T_target_C, obj.T_target_C + 273.15);
            fprintf('Coolant: T_water = %.1f K (%.1f°C)\n', obj.T_water_K, obj.T_water_K - 273.15);
            fprintf('Heat flux: %.0f W/m²\n', obj.BaseConfig.boundary_conditions.q_flux_W_m2);
            fprintf('Total configurations to evaluate: %d\n', total_configs);
            fprintf('----------------------------------------\n\n');
            
            % Pre-allocate results storage
            all_results = [];
            config_idx = 0;
            feasible_count = 0;
            failed_count = 0;
            
            % Progress tracking
            start_time = tic;
            
            for n_idx = 1:length(N_list)
                N = N_list(n_idx);
                
                for theta_idx = 1:length(theta_list)
                    theta_deg = theta_list(theta_idx);
                    
                    for t_idx = 1:length(t_list)
                        t_TEC = t_list(t_idx);
                        
                        for kr_idx = 1:length(kr_list)
                            k_r = kr_list(kr_idx);
                            
                            config_idx = config_idx + 1;
                            
                            % Create configuration
                            config = obj.createConfig(N, theta_deg, t_TEC, k_r);
                            
                            if verbose && mod(config_idx, 10) == 0
                                elapsed = toc(start_time);
                                eta = elapsed / config_idx * (total_configs - config_idx);
                                fprintf('Progress: %d/%d (%.1f%%) - ETA: %.0fs\n', ...
                                    config_idx, total_configs, 100*config_idx/total_configs, eta);
                            end
                            
                            % Check feasibility before solving
                            [is_feasible, reason] = obj.checkFeasibility(config);
                            
                            if ~is_feasible
                                if verbose
                                    fprintf('  [SKIP] N=%d, θ=%.0f°, t=%.0fμm, k_r=%.2f: %s\n', ...
                                        N, theta_deg, t_TEC, k_r, reason);
                                end
                                failed_count = failed_count + 1;
                                continue;
                            end
                            
                            % Optimize current for this configuration
                            try
                                [opt_result, success] = obj.optimizeCurrent(config);
                                
                                if success && opt_result.T_max_C < 1000  % Store results up to 1000°C for analysis
                                    opt_result.N_stages = N;
                                    opt_result.theta_deg = theta_deg;
                                    opt_result.t_TEC_um = t_TEC;
                                    opt_result.k_r = k_r;
                                    opt_result.config_idx = config_idx;
                                    opt_result.meets_target = opt_result.T_max_C < obj.T_target_C;
                                    
                                    all_results = [all_results; opt_result];
                                    
                                    if opt_result.meets_target
                                        feasible_count = feasible_count + 1;
                                        if verbose
                                            fprintf('  [OK] N=%d, θ=%.0f°, t=%.0fμm, k_r=%.2f: T=%.1f°C, I=%.3fA, COP=%.3f\n', ...
                                                N, theta_deg, t_TEC, k_r, opt_result.T_max_C, opt_result.I_opt, opt_result.COP);
                                        end
                                    else
                                        if verbose
                                            fprintf('  [HIGH] N=%d, θ=%.0f°, t=%.0fμm, k_r=%.2f: T=%.1f°C\n', ...
                                                N, theta_deg, t_TEC, k_r, opt_result.T_max_C);
                                        end
                                    end
                                else
                                    failed_count = failed_count + 1;
                                end
                            catch ME
                                if verbose
                                    fprintf('  [FAIL] N=%d, θ=%.0f°, t=%.0fμm, k_r=%.2f: %s\n', ...
                                        N, theta_deg, t_TEC, k_r, ME.message);
                                end
                                failed_count = failed_count + 1;
                            end
                        end
                    end
                end
            end
            
            % Store results
            obj.Results.all_results = all_results;
            obj.Results.summary.total_evaluated = config_idx;
            obj.Results.summary.feasible_count = feasible_count;
            obj.Results.summary.failed_count = failed_count;
            obj.Results.summary.elapsed_time = toc(start_time);
            
            % Print summary
            fprintf('\n========================================\n');
            fprintf('  OPTIMIZATION COMPLETE\n');
            fprintf('========================================\n');
            fprintf('Total configurations: %d\n', config_idx);
            fprintf('Feasible (T < %d°C): %d\n', obj.T_target_C, feasible_count);
            fprintf('Failed/Infeasible: %d\n', failed_count);
            fprintf('Elapsed time: %.1f seconds\n', obj.Results.summary.elapsed_time);
            
            % Rank and display top results
            if ~isempty(all_results)
                obj.rankAndDisplayResults();
                obj.generatePlots();
                obj.saveResults();
            else
                fprintf('\nNo feasible designs found!\n');
                fprintf('Consider: increasing t_TEC, reducing q_flux, or adjusting ranges.\n');
            end
            
            results = obj.Results;
        end
        
        function config = createConfig(obj, N_stages, theta_deg, t_TEC_um, k_r)
            % Create configuration struct with specified parameters
            config = obj.BaseConfig;
            config.geometry.N_tsv_limit = N_stages;
            config.geometry.wedge_angle_deg = theta_deg;
            config.geometry.t_TEC_um = t_TEC_um;
            config.geometry.radial_expansion_factor = k_r;
        end
        
        function [is_feasible, reason] = checkFeasibility(obj, config)
            % Pre-check configuration feasibility before solving
            is_feasible = true;
            reason = '';
            
            try
                % Create geometry to check dimensions
                geo = TECGeometry(config);
                
                % Check 1: Minimum stage length
                min_L = inf;
                for i = 1:geo.N_stages
                    [r_in, L, ~, ~, ~, ~, ~, ~, ~, ~] = geo.get_stage_geometry(i);
                    if L < min_L
                        min_L = L;
                    end
                    
                    % Check for negative or too small dimensions
                    if L < 10e-6  % Less than 10 um
                        is_feasible = false;
                        reason = sprintf('Stage %d length too small: %.2f um', i, L*1e6);
                        return;
                    end
                    
                    if r_in < 0
                        is_feasible = false;
                        reason = sprintf('Stage %d has negative inner radius', i);
                        return;
                    end
                end
                
                % Check 2: Reasonable thermal path
                % Estimate thermal resistance: R_th ~ L / (k * A)
                % If R_th is too high, temperatures will be unrealistic
                theta = geo.WedgeAngle;
                R_base = geo.R_base;
                t_TEC = config.geometry.t_TEC_um * 1e-6;
                k_bi2te3 = 1.2;  % Approximate
                
                % Very rough estimate of total thermal conductance
                A_avg = theta * R_base * t_TEC / 2;
                L_total = R_base - geo.R_cyl;
                K_estimate = k_bi2te3 * A_avg / L_total;
                
                if K_estimate < obj.MinThermalConductance
                    is_feasible = false;
                    reason = sprintf('Thermal conductance too low: %.2e W/K', K_estimate);
                    return;
                end
                
                % Check 3: Estimate max temperature (very rough)
                q_flux = config.boundary_conditions.q_flux_W_m2;
                A_chip = 0.5 * theta * R_base^2;
                Q_in = q_flux * A_chip;
                
                dT_estimate = Q_in / K_estimate;
                if dT_estimate > obj.MaxTempRise
                    is_feasible = false;
                    reason = sprintf('Estimated temp rise too high: %.0f K', dT_estimate);
                    return;
                end
                
            catch ME
                is_feasible = false;
                reason = sprintf('Geometry error: %s', ME.message);
            end
        end
        
        function [result, success] = optimizeCurrent(obj, config)
            % Optimize current for given configuration using golden section search
            result = struct();
            success = false;
            
            % Set up solver components
            try
                mat = MaterialProperties(config);
                geo = TECGeometry(config);
                net = ThermalNetwork(geo, mat, config);
            catch ME
                result.error = ME.message;
                return;
            end
            
            % Golden section search for optimal current
            I_min = obj.I_range(1);
            I_max = obj.I_range(2);
            
            phi = (1 + sqrt(5)) / 2;  % Golden ratio
            tol = 1e-4;               % Current tolerance
            max_iter = 30;
            
            % Objective: minimize T_max
            function T_max = evalCurrent(I)
                config.operating_conditions.I_current_A = I;
                net_temp = ThermalNetwork(geo, mat, config);
                
                N = geo.N_stages;
                dim = 2*N + 1;
                T_init = ones(dim, 1) * 300;
                
                % Iterative solve
                T_current = T_init;
                for iter = 1:50
                    [T_new, ~, ~] = net_temp.solve(T_current);
                    
                    % Check for invalid values
                    if any(isnan(T_new)) || any(T_new < 0) || any(T_new > 5000)
                        T_max = inf;
                        return;
                    end
                    
                    diff = max(abs(T_new - T_current));
                    T_current = 0.5 * T_current + 0.5 * T_new;
                    
                    if diff < 1e-4
                        break;
                    end
                end
                
                T_max = max(T_current);
            end
            
            % Golden section search
            a = I_min;
            b = I_max;
            c = b - (b - a) / phi;
            d = a + (b - a) / phi;
            
            fc = evalCurrent(c);
            fd = evalCurrent(d);
            
            for iter = 1:max_iter
                if abs(b - a) < tol
                    break;
                end
                
                if fc < fd
                    b = d;
                    d = c;
                    fd = fc;
                    c = b - (b - a) / phi;
                    fc = evalCurrent(c);
                else
                    a = c;
                    c = d;
                    fc = fd;
                    d = a + (b - a) / phi;
                    fd = evalCurrent(d);
                end
            end
            
            I_opt = (a + b) / 2;
            
            % Final evaluation at optimal current
            config.operating_conditions.I_current_A = I_opt;
            net_final = ThermalNetwork(geo, mat, config);
            
            N = geo.N_stages;
            dim = 2*N + 1;
            T_init = ones(dim, 1) * 300;
            T_current = T_init;
            
            for iter = 1:100
                [T_new, Q_out, Q_in] = net_final.solve(T_current);
                
                if any(isnan(T_new)) || any(T_new < 0)
                    return;
                end
                
                diff = max(abs(T_new - T_current));
                T_current = 0.5 * T_current + 0.5 * T_new;
                
                if diff < 1e-6
                    break;
                end
            end
            
            % Calculate metrics
            T_max_K = max(T_current);
            T_min_K = min(T_current);
            T_max_C = T_max_K - 273.15;
            
            % Calculate power and COP
            % P = V*I = (S*dT + I*R) * I
            S_eff = 0.0004;  % Approximate Seebeck (p+n)
            dT = T_max_K - obj.T_water_K;
            
            % Rough estimates
            P_estimate = abs(Q_in - Q_out);  % Power input
            if P_estimate > 0 && Q_in > 0
                COP = Q_in / P_estimate;
            else
                COP = 0;
            end
            
            % Figure of merit estimate
            Z_est = S_eff^2 / (1.2 * 1e-5);  % S^2 / (k * rho)
            ZT = Z_est * T_max_K;
            
            % Store results
            result.I_opt = I_opt;
            result.T_max_K = T_max_K;
            result.T_min_K = T_min_K;
            result.T_max_C = T_max_C;
            result.T_distribution = T_current;
            result.Q_in = Q_in;
            result.Q_out = Q_out;
            result.P_estimate = P_estimate;
            result.COP = COP;
            result.ZT = ZT;
            result.dT_cooling = (T_max_K - obj.T_water_K);
            result.config = config;
            
            success = true;
        end
        
        function rankAndDisplayResults(obj)
            % Rank results and display top designs
            all_results = obj.Results.all_results;
            
            if isempty(all_results)
                return;
            end
            
            % Convert to table for easier manipulation
            n = length(all_results);
            
            % Extract key metrics
            T_max_C = [all_results.T_max_C]';
            COP = [all_results.COP]';
            I_opt = [all_results.I_opt]';
            N_stages = [all_results.N_stages]';
            theta_deg = [all_results.theta_deg]';
            t_TEC_um = [all_results.t_TEC_um]';
            k_r = [all_results.k_r]';
            meets_target = [all_results.meets_target]';
            
            % Sort by T_max (ascending)
            [~, sort_idx] = sort(T_max_C);
            
            % Display top 10 results
            fprintf('\n========================================\n');
            fprintf('  TOP 10 DESIGNS (by T_max)\n');
            fprintf('========================================\n');
            fprintf('Rank | N  | θ(°) | t(μm) | k_r  | I(mA) | T_max(°C) | COP   | Target\n');
            fprintf('-----|----|----- |-------|------|-------|-----------|-------|-------\n');
            
            top_n = min(10, n);
            for i = 1:top_n
                idx = sort_idx(i);
                target_str = 'YES';
                if ~meets_target(idx)
                    target_str = 'NO ';
                end
                fprintf('%4d | %2d | %4.0f | %5.0f | %4.2f | %5.1f | %9.1f | %5.3f | %s\n', ...
                    i, N_stages(idx), theta_deg(idx), t_TEC_um(idx), k_r(idx), ...
                    I_opt(idx)*1000, T_max_C(idx), COP(idx), target_str);
            end
            
            % Store sorted results
            obj.Results.sorted_idx = sort_idx;
            obj.Results.top_designs = all_results(sort_idx(1:min(5, n)));
            
            % Display best feasible design details
            feasible_idx = find(meets_target);
            if ~isempty(feasible_idx)
                [~, best_feasible] = min(T_max_C(feasible_idx));
                best_idx = feasible_idx(best_feasible);
                
                fprintf('\n========================================\n');
                fprintf('  BEST FEASIBLE DESIGN\n');
                fprintf('========================================\n');
                best = all_results(best_idx);
                fprintf('Number of stages: %d\n', best.N_stages);
                fprintf('Wedge angle: %.1f°\n', best.theta_deg);
                fprintf('TEC thickness: %.0f μm\n', best.t_TEC_um);
                fprintf('Radial expansion factor: %.2f\n', best.k_r);
                fprintf('Optimal current: %.3f A (%.1f mA)\n', best.I_opt, best.I_opt*1000);
                fprintf('Maximum temperature: %.1f°C (%.1f K)\n', best.T_max_C, best.T_max_K);
                fprintf('Temperature drop: %.1f K\n', best.dT_cooling);
                fprintf('COP: %.3f\n', best.COP);
                fprintf('Heat absorbed: %.4f W\n', best.Q_in);
                fprintf('Heat rejected: %.4f W\n', abs(best.Q_out));
                
                obj.Results.best_feasible = best;
            end
        end
        
        function generatePlots(obj)
            % Generate publication-ready plots
            all_results = obj.Results.all_results;
            
            if isempty(all_results)
                return;
            end
            
            % Extract data
            T_max_C = [all_results.T_max_C]';
            COP = [all_results.COP]';
            N_stages = [all_results.N_stages]';
            theta_deg = [all_results.theta_deg]';
            t_TEC_um = [all_results.t_TEC_um]';
            I_opt = [all_results.I_opt]';
            meets_target = [all_results.meets_target]';
            
            % Figure 1: T_max vs N_stages (grouped by theta)
            fig1 = figure('Position', [100, 100, 800, 600]);
            
            unique_theta = unique(theta_deg);
            colors = lines(length(unique_theta));
            hold on;
            
            for i = 1:length(unique_theta)
                mask = theta_deg == unique_theta(i);
                scatter(N_stages(mask), T_max_C(mask), 60, colors(i,:), 'filled', ...
                    'DisplayName', sprintf('θ = %.0f°', unique_theta(i)));
            end
            
            yline(obj.T_target_C, 'r--', 'LineWidth', 2, 'DisplayName', sprintf('Target (%d°C)', obj.T_target_C));
            
            xlabel('Number of Stages', 'FontSize', 12);
            ylabel('Maximum Temperature (°C)', 'FontSize', 12);
            title('Effect of Stage Count on Cooling Performance', 'FontSize', 14);
            legend('Location', 'best');
            grid on;
            
            saveas(fig1, fullfile(obj.OutputDir, 'Tmax_vs_Nstages.png'));
            saveas(fig1, fullfile(obj.OutputDir, 'Tmax_vs_Nstages.fig'));
            
            % Figure 2: T_max vs TEC thickness
            fig2 = figure('Position', [100, 100, 800, 600]);
            
            unique_N = unique(N_stages);
            colors2 = lines(length(unique_N));
            hold on;
            
            for i = 1:length(unique_N)
                mask = N_stages == unique_N(i);
                scatter(t_TEC_um(mask), T_max_C(mask), 60, colors2(i,:), 'filled', ...
                    'DisplayName', sprintf('N = %d', unique_N(i)));
            end
            
            yline(obj.T_target_C, 'r--', 'LineWidth', 2, 'DisplayName', sprintf('Target (%d°C)', obj.T_target_C));
            
            xlabel('TEC Thickness (μm)', 'FontSize', 12);
            ylabel('Maximum Temperature (°C)', 'FontSize', 12);
            title('Effect of TEC Thickness on Cooling Performance', 'FontSize', 14);
            legend('Location', 'best');
            grid on;
            
            saveas(fig2, fullfile(obj.OutputDir, 'Tmax_vs_thickness.png'));
            saveas(fig2, fullfile(obj.OutputDir, 'Tmax_vs_thickness.fig'));
            
            % Figure 3: Optimal current heatmap
            fig3 = figure('Position', [100, 100, 800, 600]);
            
            unique_t = unique(t_TEC_um);
            [N_grid, T_grid] = meshgrid(unique_N, unique_t);
            I_grid = nan(size(N_grid));
            
            for i = 1:numel(N_grid)
                mask = (N_stages == N_grid(i)) & (t_TEC_um == T_grid(i));
                if any(mask)
                    I_grid(i) = mean(I_opt(mask)) * 1000;  % mA
                end
            end
            
            imagesc(unique_N, unique_t, I_grid);
            colorbar;
            colormap(jet);
            xlabel('Number of Stages', 'FontSize', 12);
            ylabel('TEC Thickness (μm)', 'FontSize', 12);
            title('Optimal Current (mA)', 'FontSize', 14);
            set(gca, 'YDir', 'normal');
            
            saveas(fig3, fullfile(obj.OutputDir, 'optimal_current_heatmap.png'));
            saveas(fig3, fullfile(obj.OutputDir, 'optimal_current_heatmap.fig'));
            
            % Figure 4: Pareto front (T_max vs COP)
            fig4 = figure('Position', [100, 100, 800, 600]);
            
            valid = COP > 0 & T_max_C < 500;
            scatter(T_max_C(valid), COP(valid), 80, N_stages(valid), 'filled');
            colorbar;
            colormap(jet);
            
            xline(obj.T_target_C, 'r--', 'LineWidth', 2);
            
            xlabel('Maximum Temperature (°C)', 'FontSize', 12);
            ylabel('Coefficient of Performance (COP)', 'FontSize', 12);
            title('Pareto Front: Temperature vs Efficiency', 'FontSize', 14);
            c = colorbar;
            c.Label.String = 'Number of Stages';
            grid on;
            
            saveas(fig4, fullfile(obj.OutputDir, 'pareto_front.png'));
            saveas(fig4, fullfile(obj.OutputDir, 'pareto_front.fig'));
            
            % Figure 5: Temperature distribution for best design
            if isfield(obj.Results, 'best_feasible')
                fig5 = figure('Position', [100, 100, 900, 500]);
                
                best = obj.Results.best_feasible;
                T_dist = best.T_distribution;
                N = best.N_stages;
                
                % Create labels
                labels = cell(2*N+1, 1);
                labels{1} = 'Node 0\n(Center)';
                for i = 1:N
                    labels{i+1} = sprintf('Si_%d', i);
                    labels{N+i+1} = sprintf('TEC_%d', i);
                end
                
                % Bar plot
                subplot(1,2,1);
                bar(T_dist - 273.15);
                ylabel('Temperature (°C)', 'FontSize', 12);
                xlabel('Node Index', 'FontSize', 12);
                title(sprintf('Temperature Distribution\n(N=%d, θ=%.0f°, t=%.0fμm)', ...
                    best.N_stages, best.theta_deg, best.t_TEC_um), 'FontSize', 12);
                yline(obj.T_target_C, 'r--', 'LineWidth', 2);
                grid on;
                
                % Schematic view
                subplot(1,2,2);
                
                % Si layer temperatures
                T_Si = T_dist(2:N+1) - 273.15;
                T_TEC = T_dist(N+2:end) - 273.15;
                
                plot(1:N, T_Si, 'b-o', 'LineWidth', 2, 'MarkerSize', 8, 'DisplayName', 'Silicon Layer');
                hold on;
                plot(1:N, T_TEC, 'r-s', 'LineWidth', 2, 'MarkerSize', 8, 'DisplayName', 'TEC Layer');
                yline(obj.T_target_C, 'k--', 'LineWidth', 1.5, 'DisplayName', 'Target');
                yline(obj.T_water_K - 273.15, 'c--', 'LineWidth', 1.5, 'DisplayName', 'Coolant');
                
                xlabel('Stage Number', 'FontSize', 12);
                ylabel('Temperature (°C)', 'FontSize', 12);
                title('Radial Temperature Profile', 'FontSize', 12);
                legend('Location', 'best');
                grid on;
                
                saveas(fig5, fullfile(obj.OutputDir, 'best_design_temperature.png'));
                saveas(fig5, fullfile(obj.OutputDir, 'best_design_temperature.fig'));
            end
            
            % Figure 6: Effect of wedge angle
            fig6 = figure('Position', [100, 100, 800, 600]);
            
            hold on;
            for i = 1:length(unique_N)
                mask = N_stages == unique_N(i);
                if sum(mask) > 1
                    [sorted_theta, sort_idx] = sort(theta_deg(mask));
                    sorted_T = T_max_C(mask);
                    sorted_T = sorted_T(sort_idx);
                    plot(sorted_theta, sorted_T, '-o', 'LineWidth', 2, ...
                        'DisplayName', sprintf('N = %d', unique_N(i)));
                end
            end
            
            yline(obj.T_target_C, 'r--', 'LineWidth', 2, 'DisplayName', 'Target');
            xlabel('Wedge Angle (°)', 'FontSize', 12);
            ylabel('Maximum Temperature (°C)', 'FontSize', 12);
            title('Effect of Wedge Angle on Cooling', 'FontSize', 14);
            legend('Location', 'best');
            grid on;
            
            saveas(fig6, fullfile(obj.OutputDir, 'Tmax_vs_wedge_angle.png'));
            saveas(fig6, fullfile(obj.OutputDir, 'Tmax_vs_wedge_angle.fig'));
            
            fprintf('\nPlots saved to: %s\n', obj.OutputDir);
            
            close all;
        end
        
        function saveResults(obj)
            % Save results to files
            
            % Save MATLAB data
            results = obj.Results;
            save(fullfile(obj.OutputDir, 'optimization_results.mat'), 'results');
            
            % Save summary as text
            fid = fopen(fullfile(obj.OutputDir, 'optimization_summary.txt'), 'w');
            
            fprintf(fid, 'RADIAL TEC DESIGN OPTIMIZATION RESULTS\n');
            fprintf(fid, '=======================================\n\n');
            fprintf(fid, 'Date: %s\n\n', datestr(now));
            
            fprintf(fid, 'OPTIMIZATION PARAMETERS\n');
            fprintf(fid, '-----------------------\n');
            fprintf(fid, 'Target T_max: %d°C\n', obj.T_target_C);
            fprintf(fid, 'Coolant temperature: %.1f K\n', obj.T_water_K);
            fprintf(fid, 'Heat flux: %.0f W/m²\n\n', obj.BaseConfig.boundary_conditions.q_flux_W_m2);
            
            fprintf(fid, 'SUMMARY\n');
            fprintf(fid, '-------\n');
            fprintf(fid, 'Total configurations evaluated: %d\n', obj.Results.summary.total_evaluated);
            fprintf(fid, 'Feasible designs (T < %d°C): %d\n', obj.T_target_C, obj.Results.summary.feasible_count);
            fprintf(fid, 'Failed/Infeasible: %d\n', obj.Results.summary.failed_count);
            fprintf(fid, 'Computation time: %.1f seconds\n\n', obj.Results.summary.elapsed_time);
            
            if isfield(obj.Results, 'best_feasible')
                best = obj.Results.best_feasible;
                fprintf(fid, 'BEST FEASIBLE DESIGN\n');
                fprintf(fid, '--------------------\n');
                fprintf(fid, 'Number of stages (N): %d\n', best.N_stages);
                fprintf(fid, 'Wedge angle (θ): %.1f°\n', best.theta_deg);
                fprintf(fid, 'TEC thickness: %.0f μm\n', best.t_TEC_um);
                fprintf(fid, 'Radial expansion factor (k_r): %.2f\n', best.k_r);
                fprintf(fid, 'Optimal current: %.4f A (%.2f mA)\n', best.I_opt, best.I_opt*1000);
                fprintf(fid, 'Maximum temperature: %.2f°C (%.2f K)\n', best.T_max_C, best.T_max_K);
                fprintf(fid, 'Minimum temperature: %.2f K\n', best.T_min_K);
                fprintf(fid, 'COP: %.4f\n', best.COP);
                fprintf(fid, 'Heat absorbed (Q_in): %.6f W\n', best.Q_in);
                fprintf(fid, 'Heat rejected (Q_out): %.6f W\n\n', abs(best.Q_out));
                
                fprintf(fid, 'TEMPERATURE DISTRIBUTION (K)\n');
                fprintf(fid, '----------------------------\n');
                T = best.T_distribution;
                N = best.N_stages;
                fprintf(fid, 'Node 0 (Center): %.2f K (%.2f°C)\n', T(1), T(1)-273.15);
                for i = 1:N
                    fprintf(fid, 'Stage %d - Si: %.2f K, TEC: %.2f K\n', i, T(i+1), T(N+i+1));
                end
            end
            
            fclose(fid);
            
            % Save best config as JSON for COMSOL
            if isfield(obj.Results, 'best_feasible')
                best_config = obj.Results.best_feasible.config;
                best_config.optimization_results.T_max_C = obj.Results.best_feasible.T_max_C;
                best_config.optimization_results.I_opt_A = obj.Results.best_feasible.I_opt;
                best_config.optimization_results.COP = obj.Results.best_feasible.COP;
                
                json_str = jsonencode(best_config);
                % Pretty print
                json_str = strrep(json_str, ',', sprintf(',\n'));
                json_str = strrep(json_str, '{', sprintf('{\n'));
                json_str = strrep(json_str, '}', sprintf('\n}'));
                
                fid = fopen(fullfile(obj.OutputDir, 'best_design_config.json'), 'w');
                fprintf(fid, '%s', json_str);
                fclose(fid);
            end
            
            % Save all feasible designs as CSV
            if ~isempty(obj.Results.all_results)
                all_r = obj.Results.all_results;
                meets = [all_r.meets_target]';
                feasible = all_r(meets);
                
                if ~isempty(feasible)
                    fid = fopen(fullfile(obj.OutputDir, 'feasible_designs.csv'), 'w');
                    fprintf(fid, 'Rank,N_stages,theta_deg,t_TEC_um,k_r,I_opt_mA,T_max_C,COP,Q_in_W,Q_out_W\n');
                    
                    % Sort by T_max
                    T_all = [feasible.T_max_C];
                    [~, sort_idx] = sort(T_all);
                    
                    for i = 1:length(feasible)
                        r = feasible(sort_idx(i));
                        fprintf(fid, '%d,%d,%.1f,%.0f,%.2f,%.3f,%.2f,%.4f,%.6f,%.6f\n', ...
                            i, r.N_stages, r.theta_deg, r.t_TEC_um, r.k_r, ...
                            r.I_opt*1000, r.T_max_C, r.COP, r.Q_in, abs(r.Q_out));
                    end
                    fclose(fid);
                end
            end
            
            fprintf('Results saved to: %s\n', obj.OutputDir);
        end
    end
end
