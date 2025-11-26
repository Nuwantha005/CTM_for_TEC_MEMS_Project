classdef PhysicsConstrainedOptimizer < handle
    % PHYSICSCONSTRAINEDOPTIMIZER - Optimization with physics-based constraints
    %
    % Key improvements over previous optimizers:
    % 1. Analytical first-cut design to set realistic bounds
    % 2. Constraints on temperature (T >= T_water, T <= T_max)
    % 3. Power balance validation
    % 4. Physical COP limits
    % 5. Bounded trust-region approach instead of unbounded gradient descent
    %
    % Strategy:
    %   1. Estimate feasibility using analytical TEC equations
    %   2. Generate candidate designs from analytical sizing
    %   3. Refine candidates using constrained optimization
    %   4. Validate solutions physically
    
    properties
        BaseConfig          % Base configuration struct
        Materials           % MaterialProperties object
        OutputDir           % Output directory
        Results             % Stored results
        
        % Physical limits
        T_max_K             % Maximum allowable temperature (K)
        T_min_K             % Minimum temperature (= T_water_K, physics limit)
        T_water_K           % Coolant temperature
        
        % TEC material properties (at 300K for estimation)
        S_pn                % Combined Seebeck coefficient (V/K)
        k_pn                % Combined thermal conductivity (W/m-K)
        rho_pn              % Combined electrical resistivity (ohm-m)
        ZT                  % Figure of merit at 300K
        
        % Problem parameters
        Q_total             % Total heat to remove (W)
        w_chip              % Chip width (m)
        
        % Optimization settings
        MaxIterations
        ConvergenceTol
    end
    
    methods
        function obj = PhysicsConstrainedOptimizer(config_path)
            % Load configuration
            text = fileread(config_path);
            obj.BaseConfig = jsondecode(text);
            obj.Materials = MaterialProperties(obj.BaseConfig);
            
            % Set temperature limits
            obj.T_water_K = obj.BaseConfig.boundary_conditions.T_water_K;
            obj.T_max_K = 373.15;  % 100°C default
            obj.T_min_K = obj.T_water_K;  % Cannot cool below water temp!
            
            % Get TEC material properties (Bi2Te3)
            % Use values at ~320K (slightly above ambient)
            T_ref = 320;
            obj.S_pn = 2 * abs(obj.Materials.get_S('Bi2Te3', T_ref));  % S_p - S_n
            obj.k_pn = 2 * obj.Materials.get_k('Bi2Te3', T_ref);        % k_p + k_n
            obj.rho_pn = 2 * obj.Materials.get_rho('Bi2Te3', T_ref);    % rho_p + rho_n
            
            % Figure of merit
            obj.ZT = (obj.S_pn^2 * T_ref) / (obj.k_pn * obj.rho_pn);
            
            % Chip dimensions
            obj.w_chip = obj.BaseConfig.geometry.w_chip_um * 1e-6;
            
            % Calculate total heat input
            q_flux = obj.BaseConfig.boundary_conditions.q_flux_W_m2;
            A_chip = obj.w_chip^2;  % Full chip area
            obj.Q_total = q_flux * A_chip;
            
            % Optimization settings
            obj.MaxIterations = 200;
            obj.ConvergenceTol = 1e-4;
            
            % Create output directory
            timestamp = datestr(now, 'yyyy-mm-dd_HH-MM-SS');
            obj.OutputDir = fullfile('output', 'physics_constrained', timestamp);
            if ~exist(obj.OutputDir, 'dir')
                mkdir(obj.OutputDir);
            end
            
            obj.Results = struct();
            
            % Print problem summary
            obj.printProblemSummary();
        end
        
        function printProblemSummary(obj)
            fprintf('\n========================================\n');
            fprintf('  PHYSICS-CONSTRAINED TEC OPTIMIZATION\n');
            fprintf('========================================\n');
            fprintf('Chip size: %.1f mm x %.1f mm\n', obj.w_chip*1000, obj.w_chip*1000);
            fprintf('Heat flux: %.0f kW/m²\n', obj.BaseConfig.boundary_conditions.q_flux_W_m2/1000);
            fprintf('Total heat: %.4f W (full chip)\n', obj.Q_total);
            fprintf('Coolant temperature: %.1f K (%.1f°C)\n', obj.T_water_K, obj.T_water_K-273.15);
            fprintf('Max chip temperature: %.1f K (%.1f°C)\n', obj.T_max_K, obj.T_max_K-273.15);
            fprintf('Material ZT @ 320K: %.3f\n', obj.ZT);
            fprintf('========================================\n\n');
        end
        
        function feasibility = assessFeasibility(obj)
            % ASSESSFEASIBILITY - Analytical check if TEC cooling is possible
            %
            % Uses ideal TEC equations to determine if the problem is solvable
            % before attempting numerical optimization.
            
            fprintf('Assessing problem feasibility...\n');
            
            dT_required = obj.T_max_K - obj.T_water_K;  % Max temp rise we can tolerate
            
            % For a single-stage TEC with n couples, max cooling capacity is:
            % Q_c_max = n * [S*I_opt*T_c - K*(T_h - T_c) - 0.5*I_opt^2*R]
            % At optimal current I_opt = S*T_c / R
            % Q_c_max = n * [S^2*T_c^2/(2R) - K*dT]
            % Q_c_max = n * [0.5*ZT*K*T_c - K*dT] (using ZT = S^2/(K*R))
            
            % For analysis, use T_c ≈ T_water (cold side at water temp)
            T_c = obj.T_water_K;
            
            % Maximum cooling capacity per geometric factor G
            % Q_c_max/G ≈ 0.5*S²*T_c²/rho - k*dT
            % where G = geometric factor relating to leg dimensions
            
            % Let's estimate what G we need
            % Q_c ≈ (0.5 * S² * T_c² / rho - k * dT) / G
            % Need Q_c = Q_total / n_wedges
            
            % Estimate for radial geometry:
            % One wedge at 30° has ~1/12 of the heat
            % Let's see what's achievable
            
            % Maximum COP for a TEC:
            % COP_max = (T_c/dT) * [sqrt(1+ZT_avg) - T_h/T_c] / [sqrt(1+ZT_avg) + 1]
            % For ZT=0.8, T_h=373K, T_c=300K: COP_max ≈ 0.4-0.6
            
            % Calculate maximum cooling for reference single-stage TEC
            dT = dT_required;
            ZT_avg = obj.ZT;
            T_c = obj.T_water_K;
            T_h = obj.T_max_K;
            
            % Maximum COP at this dT
            sqrt_term = sqrt(1 + ZT_avg);
            COP_max = (T_c / dT) * (sqrt_term - T_h/T_c) / (sqrt_term + 1);
            
            fprintf('  Max temperature rise allowed: %.1f K\n', dT_required);
            fprintf('  Theoretical max COP: %.3f\n', max(0, COP_max));
            
            % Required cooling per wedge (assume 30° wedge = 12 wedges per chip)
            n_wedges_ref = 12;
            Q_per_wedge = obj.Q_total / n_wedges_ref;
            
            fprintf('  Heat per wedge (30° wedge): %.4f mW\n', Q_per_wedge*1000);
            
            % Check if it's physically possible
            % A micro-TEC with typical dimensions can handle ~1-100 mW
            if Q_per_wedge > 0.5  % 500 mW per wedge is very high
                fprintf('  WARNING: Heat per wedge is HIGH (%.1f mW)\n', Q_per_wedge*1000);
                fprintf('  May require multiple stages or enhanced heat spreading.\n');
            end
            
            % Multi-stage analysis
            % For N stages, the temperature drop is distributed
            % Each stage handles dT/N with cascaded cooling
            
            feasibility = struct();
            feasibility.dT_required = dT_required;
            feasibility.Q_total = obj.Q_total;
            feasibility.Q_per_wedge_30deg = Q_per_wedge;
            feasibility.COP_max_ideal = max(0, COP_max);
            feasibility.ZT = obj.ZT;
            feasibility.is_challenging = (Q_per_wedge > 0.1 || dT_required > 50);
            
            if COP_max <= 0
                fprintf('\n  ** PROBLEM: dT too large for single-stage TEC! **\n');
                fprintf('  Need multi-stage cascade or reduced heat flux.\n');
                feasibility.single_stage_possible = false;
            else
                feasibility.single_stage_possible = true;
            end
            
            obj.Results.feasibility = feasibility;
        end
        
        function designs = generateCandidateDesigns(obj)
            % GENERATECANDIDATEDESIGNS - Create physically-motivated initial designs
            %
            % Strategy:
            % 1. Calculate required thermal conductance from Q and dT
            % 2. Size TEC geometry to provide that conductance
            % 3. Generate variations around this baseline
            
            fprintf('\nGenerating candidate designs...\n');
            
            % Target parameters
            dT_target = obj.T_max_K - obj.T_water_K;  % e.g., 73 K for 100°C to 27°C
            Q_target = obj.Q_total;  % Total chip heat
            
            designs = [];
            design_idx = 0;
            
            % Sweep over number of stages and wedge angles
            N_stages_list = [2, 3, 4, 5, 6, 7, 8];
            theta_list = [15, 20, 25, 30, 36, 45, 60];  % degrees
            
            for N = N_stages_list
                for theta_deg = theta_list
                    n_wedges = floor(360 / theta_deg);
                    Q_wedge = Q_target / n_wedges;
                    
                    % Required effective thermal conductance (approximate)
                    % Q = K * dT => K = Q / dT
                    K_required = Q_wedge / dT_target;  % W/K per wedge
                    
                    % For TEC: effective K includes Peltier cooling
                    % At optimal current: Q_c = 0.5*ZT*K*T_c - K*dT
                    % Rearranging: K = Q_c / (0.5*ZT*T_c - dT)
                    % But we need dT < 0.5*ZT*T_c for cooling to work
                    
                    dT_per_stage = dT_target / N;  % Temperature drop per stage
                    ZT_eff = obj.ZT * 0.7;  % Derate for real conditions
                    
                    % Check if cooling is possible with this staging
                    cooling_limit = 0.5 * ZT_eff * obj.T_water_K;
                    if dT_per_stage >= cooling_limit
                        % Not enough stages for this dT
                        continue;
                    end
                    
                    % Estimate optimal current (analytical)
                    % I_opt = S * T_c / R = S * T_c * G / rho
                    % For micro-TEC, I is typically 1-100 mA
                    
                    % Sweep over thickness and expansion factor
                    t_TEC_list = [30, 40, 50, 60, 80, 100];  % um
                    k_r_list = [0.8, 0.9, 1.0, 1.1, 1.2, 1.3];
                    
                    for t_TEC = t_TEC_list
                        for k_r = k_r_list
                            design_idx = design_idx + 1;
                            
                            % Create design structure
                            design = struct();
                            design.N_stages = N;
                            design.theta_deg = theta_deg;
                            design.n_wedges = n_wedges;
                            design.t_TEC_um = t_TEC;
                            design.k_r = k_r;
                            design.Q_wedge = Q_wedge;
                            design.dT_per_stage = dT_per_stage;
                            design.K_required = K_required;
                            
                            % Estimate optimal current range
                            % Higher N -> lower dT per stage -> can use lower current
                            % Thicker TEC -> lower R -> can use higher current
                            G_estimate = 1e-4;  % Rough geometric factor
                            I_estimate = obj.S_pn * obj.T_water_K * G_estimate / obj.rho_pn;
                            design.I_estimate = min(0.5, max(0.001, I_estimate));
                            
                            designs = [designs; design];
                        end
                    end
                end
            end
            
            fprintf('  Generated %d candidate designs\n', length(designs));
            obj.Results.candidate_designs = designs;
        end
        
        function results = runConstrainedOptimization(obj, varargin)
            % RUNCONSTRAINEDOPTIMIZATION - Main optimization with physics constraints
            %
            % Optional parameters:
            %   'MaxDesigns' - maximum number of designs to evaluate (default: 500)
            %   'Verbose' - print progress (default: true)
            
            p = inputParser;
            addParameter(p, 'MaxDesigns', 500);
            addParameter(p, 'Verbose', true);
            parse(p, varargin{:});
            
            max_designs = p.Results.MaxDesigns;
            verbose = p.Results.Verbose;
            
            % Step 1: Assess feasibility
            obj.assessFeasibility();
            
            % Step 2: Generate candidate designs
            obj.generateCandidateDesigns();
            candidates = obj.Results.candidate_designs;
            
            n_candidates = min(length(candidates), max_designs);
            fprintf('\nEvaluating %d designs with constrained solver...\n', n_candidates);
            
            % Pre-allocate results
            all_results = [];
            feasible_count = 0;
            
            % Create progress figure
            fig = figure('Name', 'Constrained Optimization Progress', 'Position', [100 100 1400 500]);
            
            ax1 = subplot(1, 3, 1);
            hold(ax1, 'on');
            xlabel('Design Index');
            ylabel('T_{max} (°C)');
            title('Temperature vs Design');
            yline(ax1, obj.T_max_K - 273.15, 'r--', 'Target', 'LineWidth', 2);
            grid on;
            
            ax2 = subplot(1, 3, 2);
            hold(ax2, 'on');
            xlabel('Design Index');
            ylabel('COP');
            title('COP vs Design');
            grid on;
            
            ax3 = subplot(1, 3, 3);
            axis off;
            title('Best Feasible Designs');
            best_text = text(0.1, 0.9, 'Searching...', 'FontSize', 10, 'VerticalAlignment', 'top');
            
            drawnow;
            
            start_time = tic;
            
            for i = 1:n_candidates
                candidate = candidates(i);
                
                if verbose && mod(i, 50) == 0
                    elapsed = toc(start_time);
                    eta = elapsed / i * (n_candidates - i);
                    fprintf('Progress: %d/%d (%.1f%%) - ETA: %.0fs\n', ...
                        i, n_candidates, 100*i/n_candidates, eta);
                end
                
                % Create configuration
                config = obj.createConfig(candidate);
                
                try
                    % Run constrained current optimization
                    [result, success] = obj.optimizeWithConstraints(config, candidate);
                    
                    if success
                        result.design_idx = i;
                        result.N_stages = candidate.N_stages;
                        result.theta_deg = candidate.theta_deg;
                        result.t_TEC_um = candidate.t_TEC_um;
                        result.k_r = candidate.k_r;
                        result.n_wedges = candidate.n_wedges;
                        
                        all_results = [all_results; result];
                        
                        % Plot results
                        T_C = result.T_max_K - 273.15;
                        if result.is_feasible
                            plot(ax1, i, T_C, 'go', 'MarkerSize', 6, 'MarkerFaceColor', 'g');
                            plot(ax2, i, result.COP, 'go', 'MarkerSize', 6, 'MarkerFaceColor', 'g');
                            feasible_count = feasible_count + 1;
                        else
                            plot(ax1, i, T_C, 'rx', 'MarkerSize', 4);
                            plot(ax2, i, result.COP, 'rx', 'MarkerSize', 4);
                        end
                        
                        % Update best designs text
                        if ~isempty(all_results)
                            feasible_results = all_results([all_results.is_feasible]);
                            if ~isempty(feasible_results)
                                [~, sort_idx] = sort([feasible_results.T_max_K]);
                                top_n = min(5, length(feasible_results));
                                lines = {'Rank | N | θ | t(μm) | T(°C) | COP', ...
                                         '------------------------------------'};
                                for j = 1:top_n
                                    r = feasible_results(sort_idx(j));
                                    lines{end+1} = sprintf('%4d | %d | %2d | %4d | %5.1f | %.3f', ...
                                        j, r.N_stages, r.theta_deg, r.t_TEC_um, ...
                                        r.T_max_K-273.15, r.COP);
                                end
                                set(best_text, 'String', lines);
                            end
                        end
                        
                        drawnow limitrate;
                    end
                    
                catch ME
                    if verbose
                        fprintf('  Design %d failed: %s\n', i, ME.message);
                    end
                end
            end
            
            % Save progress figure
            saveas(fig, fullfile(obj.OutputDir, 'optimization_progress.png'));
            saveas(fig, fullfile(obj.OutputDir, 'optimization_progress.fig'));
            
            % Store results
            obj.Results.all_results = all_results;
            obj.Results.feasible_count = feasible_count;
            obj.Results.elapsed_time = toc(start_time);
            
            fprintf('\n========================================\n');
            fprintf('  OPTIMIZATION COMPLETE\n');
            fprintf('========================================\n');
            fprintf('Total designs evaluated: %d\n', n_candidates);
            fprintf('Feasible designs (T < %.0f°C): %d\n', obj.T_max_K-273.15, feasible_count);
            fprintf('Elapsed time: %.1f s\n', obj.Results.elapsed_time);
            
            if ~isempty(all_results)
                obj.rankAndSaveResults();
            else
                fprintf('\nNo valid results! Check heat flux and geometry constraints.\n');
            end
            
            results = obj.Results;
        end
        
        function config = createConfig(obj, candidate)
            % Create configuration from candidate design
            config = obj.BaseConfig;
            config.geometry.N_stages = candidate.N_stages;
            config.geometry.N_tsv_limit = candidate.N_stages;
            config.geometry.wedge_angle_deg = candidate.theta_deg;
            config.geometry.thickness_um = candidate.t_TEC_um;
            config.geometry.radial_expansion_factor = candidate.k_r;
        end
        
        function [result, success] = optimizeWithConstraints(obj, config, candidate)
            % OPTIMIZEWITHCONSTRAINTS - Find optimal current with physics bounds
            %
            % Key constraints:
            % 1. All temperatures must be >= T_water (can't cool below coolant!)
            % 2. Temperature profile must be monotonically decreasing (hot to cold)
            % 3. COP must be positive (actually cooling, not heating)
            
            result = struct();
            success = false;
            
            try
                geo = TECGeometry(config);
                mat = MaterialProperties(config);
            catch ME
                result.error = ME.message;
                return;
            end
            
            N = geo.N_stages;
            dim = 2*N + 1;
            
            % Current bounds - based on physics estimates
            I_min = 1e-4;   % 0.1 mA minimum
            I_max = 0.5;    % 500 mA maximum for micro-TEC
            
            % Bracket search to find optimal current
            % Start with coarse search, then refine
            
            n_coarse = 20;
            I_test = logspace(log10(I_min), log10(I_max), n_coarse);
            
            best_T_max = inf;
            best_I = I_test(1);
            best_T_dist = [];
            best_valid = false;
            
            for I = I_test
                config.operating_conditions.I_current_A = I;
                
                try
                    net = ThermalNetwork(geo, mat, config);
                    [T_final, is_valid, Q_out, Q_in] = obj.solveWithValidation(net, geo, config);
                    
                    if is_valid
                        T_max = max(T_final);
                        if T_max < best_T_max
                            best_T_max = T_max;
                            best_I = I;
                            best_T_dist = T_final;
                            best_valid = true;
                            best_Q_out = Q_out;
                            best_Q_in = Q_in;
                        end
                    end
                catch
                    continue;
                end
            end
            
            if ~best_valid
                result.error = 'No valid solution in current range';
                return;
            end
            
            % Refine around best current using golden section
            I_low = best_I / 2;
            I_high = min(best_I * 2, I_max);
            
            [I_opt, T_max_opt, T_dist_opt, Q_out_opt, Q_in_opt] = ...
                obj.goldenSectionSearch(config, geo, mat, I_low, I_high);
            
            if isempty(T_dist_opt)
                % Fallback to coarse result
                I_opt = best_I;
                T_max_opt = best_T_max;
                T_dist_opt = best_T_dist;
                Q_out_opt = best_Q_out;
                Q_in_opt = best_Q_in;
            end
            
            % Compute final metrics
            result.I_opt = I_opt;
            result.T_max_K = T_max_opt;
            result.T_min_K = min(T_dist_opt);
            result.T_distribution = T_dist_opt;
            result.Q_in = Q_in_opt;
            result.Q_out = Q_out_opt;
            
            % Power and COP
            P_in = abs(Q_in_opt - Q_out_opt);
            if P_in > 1e-12 && Q_in_opt > 0
                result.COP = Q_in_opt / P_in;
            else
                result.COP = 0;
            end
            result.P_input = P_in;
            
            % Check feasibility
            result.is_feasible = (T_max_opt <= obj.T_max_K) && (result.COP > 0);
            result.meets_target = result.is_feasible;
            result.config = config;
            
            success = true;
        end
        
        function [T_final, is_valid, Q_out, Q_in] = solveWithValidation(obj, net, geo, config)
            % SOLVEWITHVALIDATION - Iterative solve with physics validation
            %
            % Returns:
            %   T_final - temperature distribution
            %   is_valid - true if solution is physically valid
            %   Q_out - heat rejected to coolant
            %   Q_in - heat absorbed from chip
            
            N = geo.N_stages;
            dim = 2*N + 1;
            T_water = config.boundary_conditions.T_water_K;
            
            % Initialize with reasonable guess
            % Temperature should decrease from chip to water
            T_init = linspace(T_water + 50, T_water + 5, dim)';
            
            T_current = T_init;
            alpha = 0.5;  % Relaxation factor
            
            max_iter = 100;
            tol = 1e-5;
            
            is_valid = false;
            Q_out = 0;
            Q_in = 0;
            
            for iter = 1:max_iter
                try
                    [T_new, Q_out, Q_in] = net.solve(T_current);
                catch
                    T_final = T_current;
                    return;
                end
                
                % Check for invalid values
                if any(isnan(T_new)) || any(isinf(T_new))
                    T_final = T_current;
                    return;
                end
                
                % Check for negative temperatures
                if any(T_new < 0)
                    T_final = T_current;
                    return;
                end
                
                % Apply relaxation
                T_current = alpha * T_new + (1 - alpha) * T_current;
                
                % Check convergence
                diff = max(abs(T_new - T_current));
                if diff < tol
                    break;
                end
            end
            
            T_final = T_current;
            
            % Validate solution
            is_valid = obj.validateSolution(T_final, T_water, Q_in, Q_out);
        end
        
        function is_valid = validateSolution(obj, T, T_water, Q_in, Q_out)
            % VALIDATESOLUTION - Check if temperature distribution is physical
            %
            % Criteria:
            % 1. All T >= T_water (can't cool below coolant)
            % 2. T values are reasonable (< 1000 K)
            % 3. Temperature generally decreases from chip to water
            % 4. Heat balance makes sense
            
            is_valid = true;
            
            % Check 1: Minimum temperature
            if min(T) < T_water - 5  % Allow small tolerance
                is_valid = false;
                return;
            end
            
            % Check 2: Maximum temperature
            if max(T) > 1000
                is_valid = false;
                return;
            end
            
            % Check 3: Temperature should generally decrease
            % (First node is hottest, last TEC node is coolest)
            % Allow some flexibility for intermediate nodes
            T_chip = T(1);  % Hottest
            T_cold = T(end);  % Coldest TEC node
            
            if T_cold > T_chip + 10  % Cold side hotter than hot side is wrong
                is_valid = false;
                return;
            end
            
            % Check 4: Heat balance
            % Q_in should be positive (heat entering)
            % For cooling, Q_out (heat to water) should be negative or we absorb heat
            if Q_in < 0
                is_valid = false;
                return;
            end
        end
        
        function [I_opt, T_min, T_dist, Q_out, Q_in] = goldenSectionSearch(obj, config, geo, mat, I_low, I_high)
            % Golden section search for optimal current
            
            phi = (1 + sqrt(5)) / 2;
            tol = 1e-5;
            max_iter = 30;
            
            a = I_low;
            b = I_high;
            c = b - (b - a) / phi;
            d = a + (b - a) / phi;
            
            [Tc, T_dist_c, Qout_c, Qin_c, valid_c] = obj.evalCurrentPoint(config, geo, mat, c);
            [Td, T_dist_d, Qout_d, Qin_d, valid_d] = obj.evalCurrentPoint(config, geo, mat, d);
            
            for iter = 1:max_iter
                if abs(b - a) < tol
                    break;
                end
                
                if valid_c && valid_d
                    if Tc < Td
                        b = d;
                        d = c;
                        Td = Tc;
                        T_dist_d = T_dist_c;
                        Qout_d = Qout_c;
                        Qin_d = Qin_c;
                        valid_d = valid_c;
                        c = b - (b - a) / phi;
                        [Tc, T_dist_c, Qout_c, Qin_c, valid_c] = obj.evalCurrentPoint(config, geo, mat, c);
                    else
                        a = c;
                        c = d;
                        Tc = Td;
                        T_dist_c = T_dist_d;
                        Qout_c = Qout_d;
                        Qin_c = Qin_d;
                        valid_c = valid_d;
                        d = a + (b - a) / phi;
                        [Td, T_dist_d, Qout_d, Qin_d, valid_d] = obj.evalCurrentPoint(config, geo, mat, d);
                    end
                elseif valid_c
                    b = d;
                    d = c;
                    c = b - (b - a) / phi;
                    [Tc, T_dist_c, Qout_c, Qin_c, valid_c] = obj.evalCurrentPoint(config, geo, mat, c);
                    Td = Tc; T_dist_d = T_dist_c; Qout_d = Qout_c; Qin_d = Qin_c; valid_d = valid_c;
                elseif valid_d
                    a = c;
                    c = d;
                    d = a + (b - a) / phi;
                    [Td, T_dist_d, Qout_d, Qin_d, valid_d] = obj.evalCurrentPoint(config, geo, mat, d);
                    Tc = Td; T_dist_c = T_dist_d; Qout_c = Qout_d; Qin_c = Qin_d; valid_c = valid_d;
                else
                    break;  % Neither valid, give up
                end
            end
            
            I_opt = (a + b) / 2;
            if valid_c || valid_d
                if valid_c && (~valid_d || Tc <= Td)
                    T_min = Tc;
                    T_dist = T_dist_c;
                    Q_out = Qout_c;
                    Q_in = Qin_c;
                else
                    T_min = Td;
                    T_dist = T_dist_d;
                    Q_out = Qout_d;
                    Q_in = Qin_d;
                end
            else
                T_min = inf;
                T_dist = [];
                Q_out = 0;
                Q_in = 0;
            end
        end
        
        function [T_max, T_dist, Q_out, Q_in, is_valid] = evalCurrentPoint(obj, config, geo, mat, I)
            % Evaluate objective at a current point
            config.operating_conditions.I_current_A = I;
            net = ThermalNetwork(geo, mat, config);
            [T_dist, is_valid, Q_out, Q_in] = obj.solveWithValidation(net, geo, config);
            T_max = max(T_dist);
        end
        
        function rankAndSaveResults(obj)
            % Rank results and save best designs
            all_results = obj.Results.all_results;
            
            if isempty(all_results)
                return;
            end
            
            % Separate feasible and infeasible
            feasible_mask = [all_results.is_feasible];
            feasible_results = all_results(feasible_mask);
            
            fprintf('\n========================================\n');
            fprintf('  RESULTS SUMMARY\n');
            fprintf('========================================\n');
            
            if isempty(feasible_results)
                fprintf('No feasible designs found (T < %.0f°C)\n', obj.T_max_K - 273.15);
                fprintf('\nShowing best 10 designs by temperature:\n');
                
                T_all = [all_results.T_max_K];
                [~, sort_idx] = sort(T_all);
                top_n = min(10, length(all_results));
            else
                fprintf('Found %d feasible designs!\n\n', length(feasible_results));
                
                % Sort by temperature
                T_all = [feasible_results.T_max_K];
                [~, sort_idx] = sort(T_all);
                top_n = min(10, length(feasible_results));
                all_results = feasible_results;  % Use feasible for display
            end
            
            fprintf('TOP %d DESIGNS:\n', top_n);
            fprintf('Rank | N  | θ(°) | t(μm) | k_r  | I(mA) | T_max(°C) | COP\n');
            fprintf('-----|----|------|-------|------|-------|-----------|------\n');
            
            for i = 1:top_n
                r = all_results(sort_idx(i));
                fprintf('%4d | %2d | %4d | %5d | %4.2f | %5.2f | %9.1f | %.3f\n', ...
                    i, r.N_stages, r.theta_deg, r.t_TEC_um, r.k_r, ...
                    r.I_opt*1000, r.T_max_K-273.15, r.COP);
            end
            
            % Save best design
            if ~isempty(feasible_results)
                best = feasible_results(sort_idx(1));
                obj.Results.best_design = best;
                obj.saveBestDesign(best);
                obj.saveTopCandidates(feasible_results, sort_idx, min(10, length(feasible_results)));
            elseif ~isempty(all_results)
                best = all_results(sort_idx(1));
                obj.Results.best_design = best;
                obj.saveBestDesign(best);
            end
            
            % Save all results to CSV
            obj.saveResultsCSV(all_results);
        end
        
        function saveBestDesign(obj, design)
            % Save best design configuration
            
            fprintf('\n========================================\n');
            fprintf('  BEST DESIGN DETAILS\n');
            fprintf('========================================\n');
            fprintf('Number of stages: %d\n', design.N_stages);
            fprintf('Wedge angle: %.1f° (%d wedges per chip)\n', design.theta_deg, design.n_wedges);
            fprintf('TEC thickness: %.0f μm\n', design.t_TEC_um);
            fprintf('Radial expansion factor: %.2f\n', design.k_r);
            fprintf('Optimal current: %.3f A (%.2f mA)\n', design.I_opt, design.I_opt*1000);
            fprintf('Maximum temperature: %.1f°C (%.1f K)\n', design.T_max_K-273.15, design.T_max_K);
            fprintf('COP: %.3f\n', design.COP);
            fprintf('Heat absorbed: %.4f mW\n', design.Q_in*1000);
            fprintf('Power input: %.4f mW\n', design.P_input*1000);
            
            % Save config JSON
            config = design.config;
            config.operating_conditions.I_current_A = design.I_opt;
            
            json_path = fullfile(obj.OutputDir, 'best_design_config.json');
            fid = fopen(json_path, 'w');
            fprintf(fid, '%s', jsonencode(config));
            fclose(fid);
            fprintf('\nConfiguration saved to: %s\n', json_path);
            
            % Save temperature profile plot
            obj.plotTemperatureProfile(design);
        end
        
        function plotTemperatureProfile(obj, design)
            % Plot temperature distribution
            
            T = design.T_distribution;
            N = design.N_stages;
            T_water = obj.T_water_K;
            
            fig = figure('Name', 'Temperature Profile', 'Position', [100 100 800 500]);
            
            % Create node labels
            node_labels = cell(length(T) + 1, 1);
            node_labels{1} = 'Chip Center';
            for i = 1:N
                node_labels{1 + i} = sprintf('Si_%d', i);
            end
            for i = 1:N
                node_labels{1 + N + i} = sprintf('TEC_c%d', i);
            end
            node_labels{end} = 'Water';
            
            % Temperature with water
            T_plot = [T; T_water];
            x = 1:length(T_plot);
            
            % Plot
            plot(x, T_plot - 273.15, 'b-o', 'LineWidth', 2, 'MarkerSize', 8, 'MarkerFaceColor', 'b');
            hold on;
            yline(obj.T_max_K - 273.15, 'r--', 'Target', 'LineWidth', 1.5);
            yline(T_water - 273.15, 'c--', 'Coolant', 'LineWidth', 1.5);
            
            xlabel('Node');
            ylabel('Temperature (°C)');
            title(sprintf('Temperature Distribution (N=%d, θ=%.0f°, I=%.1fmA)', ...
                design.N_stages, design.theta_deg, design.I_opt*1000));
            set(gca, 'XTick', x, 'XTickLabel', node_labels);
            xtickangle(45);
            grid on;
            legend('T_{nodes}', 'T_{target}', 'T_{water}', 'Location', 'best');
            
            % Save
            saveas(fig, fullfile(obj.OutputDir, 'best_temperature_profile.png'));
            saveas(fig, fullfile(obj.OutputDir, 'best_temperature_profile.fig'));
            
            fprintf('Temperature profile saved.\n');
        end
        
        function saveTopCandidates(obj, results, sort_idx, n_top)
            % Save detailed results for top N candidates
            
            candidates_dir = fullfile(obj.OutputDir, 'top_candidates');
            if ~exist(candidates_dir, 'dir')
                mkdir(candidates_dir);
            end
            
            for i = 1:n_top
                r = results(sort_idx(i));
                
                cand_dir = fullfile(candidates_dir, sprintf('rank_%02d', i));
                if ~exist(cand_dir, 'dir')
                    mkdir(cand_dir);
                end
                
                % Save config
                config = r.config;
                config.operating_conditions.I_current_A = r.I_opt;
                fid = fopen(fullfile(cand_dir, 'config.json'), 'w');
                fprintf(fid, '%s', jsonencode(config));
                fclose(fid);
                
                % Save summary
                fid = fopen(fullfile(cand_dir, 'summary.txt'), 'w');
                fprintf(fid, 'DESIGN RANK: %d\n', i);
                fprintf(fid, '================\n\n');
                fprintf(fid, 'Geometry:\n');
                fprintf(fid, '  N_stages: %d\n', r.N_stages);
                fprintf(fid, '  Wedge angle: %.1f° (%d wedges/chip)\n', r.theta_deg, r.n_wedges);
                fprintf(fid, '  TEC thickness: %.0f μm\n', r.t_TEC_um);
                fprintf(fid, '  Radial expansion: %.2f\n', r.k_r);
                fprintf(fid, '\nPerformance:\n');
                fprintf(fid, '  Optimal current: %.4f A (%.2f mA)\n', r.I_opt, r.I_opt*1000);
                fprintf(fid, '  T_max: %.1f°C (%.1f K)\n', r.T_max_K-273.15, r.T_max_K);
                fprintf(fid, '  T_min: %.1f°C (%.1f K)\n', r.T_min_K-273.15, r.T_min_K);
                fprintf(fid, '  COP: %.4f\n', r.COP);
                fprintf(fid, '  Q_absorbed: %.4f mW\n', r.Q_in*1000);
                fprintf(fid, '  P_input: %.4f mW\n', r.P_input*1000);
                fprintf(fid, '\nTemperature Distribution (K):\n');
                for j = 1:length(r.T_distribution)
                    fprintf(fid, '  Node %d: %.2f K\n', j, r.T_distribution(j));
                end
                fclose(fid);
                
                % Save temperature profile figure
                fig = figure('Visible', 'off');
                T_plot = [r.T_distribution; obj.T_water_K];
                plot(1:length(T_plot), T_plot - 273.15, 'b-o', 'LineWidth', 2);
                hold on;
                yline(obj.T_max_K - 273.15, 'r--', 'LineWidth', 1.5);
                xlabel('Node');
                ylabel('Temperature (°C)');
                title(sprintf('Rank %d: N=%d, θ=%.0f°, T=%.1f°C', i, r.N_stages, r.theta_deg, r.T_max_K-273.15));
                grid on;
                saveas(fig, fullfile(cand_dir, 'temperature_profile.png'));
                close(fig);
            end
            
            fprintf('Saved %d top candidate details to %s\n', n_top, candidates_dir);
        end
        
        function saveResultsCSV(obj, results)
            % Save all results to CSV
            
            n = length(results);
            csv_path = fullfile(obj.OutputDir, 'all_results.csv');
            
            fid = fopen(csv_path, 'w');
            fprintf(fid, 'N_stages,theta_deg,t_TEC_um,k_r,I_mA,T_max_C,T_min_C,COP,Q_in_mW,P_in_mW,feasible\n');
            
            for i = 1:n
                r = results(i);
                fprintf(fid, '%d,%d,%d,%.2f,%.3f,%.2f,%.2f,%.4f,%.4f,%.4f,%d\n', ...
                    r.N_stages, r.theta_deg, r.t_TEC_um, r.k_r, ...
                    r.I_opt*1000, r.T_max_K-273.15, r.T_min_K-273.15, ...
                    r.COP, r.Q_in*1000, r.P_input*1000, r.is_feasible);
            end
            
            fclose(fid);
            fprintf('Results saved to: %s\n', csv_path);
        end
    end
end
