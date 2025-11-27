% run_global_optimization.m
% Global optimization using Genetic Algorithm and Particle Swarm
%
% This finds the globally optimal TEC design by exploring the full
% parameter space, avoiding local minima that gradient-based methods miss.
%
% Optimizes ALL parameters as RATIOS (except N_stages, T_target, N_tsv_limit)
%
% Uses: Global Optimization Toolbox (ga, particleswarm)

clear; clc;
addpath(genpath('../../src'));

fprintf('╔════════════════════════════════════════════════════════════════════╗\n');
fprintf('║         GLOBAL OPTIMIZATION FOR TEC DESIGN (FULL PARAMETERS)       ║\n');
fprintf('╚════════════════════════════════════════════════════════════════════╝\n\n');

%% Configuration - FIXED PARAMETERS (not optimized)
CONFIG = struct();
CONFIG.T_target_C = 85;              % Target max temperature (°C) - FIXED
CONFIG.N_stages = 3;                 % Number of stages - FIXED
CONFIG.N_tsv_limit = 0;              % TSV limit - FIXED

% Fixed boundary conditions (design requirements)
CONFIG.q_flux_W_m2 = 1e6;            % Heat flux (W/m²) - FIXED DESIGN REQUIREMENT
CONFIG.h_conv_W_m2K = 1e6;           % Convection coefficient (W/m²K) - FIXED
CONFIG.T_water_K = 300;              % Coolant temperature (K) - FIXED

% Reference values for ratio-based parameters
CONFIG.ref.w_chip_um = 10000;        % Reference chip width (µm)
CONFIG.ref.thickness_um = 100;       % Reference TEC thickness (µm)

% Output directory
timestamp = datestr(now, 'yyyy-mm-dd_HH-MM-SS');
OUTPUT_DIR = fullfile('../../output', 'global_optimization', timestamp);
if ~exist(OUTPUT_DIR, 'dir')
    mkdir(OUTPUT_DIR);
end

%% Define optimization variables and bounds
% ═══════════════════════════════════════════════════════════════════════════
% OPTIMIZATION VARIABLE CONFIGURATION
% ═══════════════════════════════════════════════════════════════════════════
% Variables are loaded from: src/config/optimization_variables.m
% Edit that file to change bounds, defaults, or enable/disable variables.
% This provides a SINGLE CONTROL POINT for all optimization scripts.
% ═══════════════════════════════════════════════════════════════════════════

[var_names, lb, ub, x0, all_vars] = optimization_variables();
nvars = length(var_names);

% Store all_vars in CONFIG for use in objective function
CONFIG.all_vars = all_vars;

fprintf('Optimization Variables (%d enabled):\n', nvars);
fprintf('═══════════════════════════════════════════════════════════════\n');
fprintf('  %-35s  %10s  %10s  %10s\n', 'Variable', 'Lower', 'Upper', 'Initial');
fprintf('───────────────────────────────────────────────────────────────\n');
for i = 1:nvars
    fprintf('  %-35s  %10.4f  %10.4f  %10.4f\n', var_names{i}, lb(i), ub(i), x0(i));
end
fprintf('═══════════════════════════════════════════════════════════════\n');
fprintf('\nFixed Parameters (Design Requirements):\n');
fprintf('  N_stages:     %d\n', CONFIG.N_stages);
fprintf('  N_tsv_limit:  %d\n', CONFIG.N_tsv_limit);
fprintf('  T_target:     %.0f °C\n', CONFIG.T_target_C);
fprintf('  q_flux:       %.2e W/m² (%.0f kW/m²)\n', CONFIG.q_flux_W_m2, CONFIG.q_flux_W_m2/1e3);
fprintf('  h_conv:       %.2e W/m²K\n', CONFIG.h_conv_W_m2K);
fprintf('  T_water:      %.1f K (%.1f °C)\n', CONFIG.T_water_K, CONFIG.T_water_K - 273.15);
fprintf('───────────────────────────────────────────────────────────────\n\n');

%% Create base configuration
base_config = create_base_config(CONFIG);

%% Define objective function
% Minimize T_max while penalizing infeasible designs

objective = @(x) tec_objective(x, base_config, CONFIG);

%% Genetic Algorithm Options
fprintf('=== PHASE 1: GENETIC ALGORITHM ===\n\n');

ga_options = optimoptions('ga', ...
    'PopulationSize', 100, ...
    'MaxGenerations', 50, ...
    'FunctionTolerance', 1e-4, ...
    'Display', 'iter', ...
    'PlotFcn', {@gaplotbestf, @gaplotdistance}, ...
    'UseParallel', false, ...  % Set true if you have Parallel Computing Toolbox
    'OutputFcn', @(options, state, flag) ga_output(options, state, flag, var_names));

fprintf('Running Genetic Algorithm...\n');
fprintf('Population: 100, Generations: 50\n\n');

tic;
[x_ga, fval_ga, exitflag_ga, output_ga, population_ga, scores_ga] = ...
    ga(objective, nvars, [], [], [], [], lb, ub, [], ga_options);
time_ga = toc;

fprintf('\nGA completed in %.1f seconds\n', time_ga);
fprintf('Best solution: T_max = %.2f °C\n\n', fval_ga);

%% Particle Swarm Optimization
fprintf('=== PHASE 2: PARTICLE SWARM OPTIMIZATION ===\n\n');

pso_options = optimoptions('particleswarm', ...
    'SwarmSize', 100, ...
    'MaxIterations', 100, ...
    'FunctionTolerance', 1e-6, ...
    'Display', 'iter', ...
    'PlotFcn', @pswplotbestf, ...
    'UseParallel', false, ...
    'InitialSwarmMatrix', population_ga);  % Start from GA population

fprintf('Running Particle Swarm...\n');
fprintf('Swarm Size: 100, Max Iterations: 100\n\n');

tic;
[x_pso, fval_pso, exitflag_pso, output_pso] = ...
    particleswarm(objective, nvars, lb, ub, pso_options);
time_pso = toc;

fprintf('\nPSO completed in %.1f seconds\n', time_pso);
fprintf('Best solution: T_max = %.2f °C\n\n', fval_pso);

%% Compare with Local Optimization (fmincon)
fprintf('=== PHASE 3: LOCAL REFINEMENT (fmincon) ===\n\n');

% Use best global solution as starting point
x0 = x_pso;

fmincon_options = optimoptions('fmincon', ...
    'Algorithm', 'sqp', ...
    'Display', 'iter', ...
    'MaxIterations', 100, ...
    'OptimalityTolerance', 1e-8);

tic;
[x_local, fval_local, exitflag_local] = ...
    fmincon(objective, x0, [], [], [], [], lb, ub, [], fmincon_options);
time_local = toc;

fprintf('\nLocal refinement completed in %.1f seconds\n', time_local);
fprintf('Final solution: T_max = %.2f °C\n\n', fval_local);

%% Results Summary
fprintf('╔════════════════════════════════════════════════════════════╗\n');
fprintf('║                    OPTIMIZATION RESULTS                    ║\n');
fprintf('╚════════════════════════════════════════════════════════════╝\n\n');

% Choose best overall
if fval_local <= fval_pso && fval_local <= fval_ga
    x_best = x_local;
    fval_best = fval_local;
    method_best = 'fmincon (local refinement)';
elseif fval_pso <= fval_ga
    x_best = x_pso;
    fval_best = fval_pso;
    method_best = 'Particle Swarm';
else
    x_best = x_ga;
    fval_best = fval_ga;
    method_best = 'Genetic Algorithm';
end

fprintf('Method Comparison:\n');
fprintf('─────────────────────────────────────────────────────────────\n');
fprintf('  %-25s  T_max = %7.2f °C  (%.1f s)\n', 'Genetic Algorithm:', fval_ga, time_ga);
fprintf('  %-25s  T_max = %7.2f °C  (%.1f s)\n', 'Particle Swarm:', fval_pso, time_pso);
fprintf('  %-25s  T_max = %7.2f °C  (%.1f s)\n', 'Local Refinement:', fval_local, time_local);
fprintf('─────────────────────────────────────────────────────────────\n\n');

fprintf('BEST SOLUTION (from %s):\n', method_best);
fprintf('═══════════════════════════════════════════════════════════════\n');
fprintf('Fixed Design Requirements:\n');
fprintf('  Heat Flux:         %.2e W/m² (%.0f kW/m²)\n', CONFIG.q_flux_W_m2, CONFIG.q_flux_W_m2/1e3);
fprintf('  Coolant Temp:      %.1f K (%.1f °C)\n', CONFIG.T_water_K, CONFIG.T_water_K - 273.15);
fprintf('  h_conv:            %.2e W/m²K\n', CONFIG.h_conv_W_m2K);
fprintf('\nOptimized Operating Conditions:\n');
fprintf('  Current:           %.3f A (%.1f mA)\n', x_best(1), x_best(1)*1000);
fprintf('\nOptimized Geometry:\n');
fprintf('  Wedge Angle:       %.1f°\n', x_best(2));
fprintf('  TEC Thickness:     %.1f µm (ratio=%.3f)\n', x_best(3)*CONFIG.ref.thickness_um, x_best(3));
fprintf('  Chip Thickness:    %.1f µm (ratio=%.3f)\n', x_best(4)*100, x_best(4));
fprintf('  R_cyl:             %.1f µm (ratio=%.3f)\n', x_best(5)*CONFIG.ref.w_chip_um, x_best(5));
fprintf('  k_r (radial exp):  %.3f\n', x_best(6));
fprintf('  Fill Factor:       %.3f\n', x_best(7));
fprintf('\nOptimized Interconnects:\n');
fprintf('  Interconnect ratio:          %.3f\n', x_best(8));
fprintf('  Outerconnect ratio:          %.3f\n', x_best(9));
fprintf('  Interconnect angle ratio:    %.3f\n', x_best(10));
fprintf('  Outerconnect angle ratio:    %.3f\n', x_best(11));
fprintf('  Interconnect thickness ratio:%.3f\n', x_best(12));
fprintf('  Outerconnect thickness ratio:%.3f\n', x_best(13));
fprintf('  Insulation width ratio:      %.3f\n', x_best(14));
fprintf('───────────────────────────────────────────────────────────────\n');
fprintf('  T_max:             %.2f °C\n', fval_best);
fprintf('═══════════════════════════════════════════════════════════════\n\n');

%% Validate best solution with physical feasibility checks
fprintf('Validating best solution with physical feasibility checks...\n');

% Validate all candidate solutions
candidates = {
    x_ga, fval_ga, 'Genetic Algorithm';
    x_pso, fval_pso, 'Particle Swarm';
    x_local, fval_local, 'Local Refinement'
};

valid_solutions = {};

% T_water is now fixed from CONFIG
T_water_K_fixed = CONFIG.T_water_K;

for i = 1:size(candidates, 1)
    x_cand = candidates{i, 1};
    fval_cand = candidates{i, 2};
    method_cand = candidates{i, 3};
    
    [is_valid, reason, T_max_cand, T_profile_cand, COP_cand] = ...
        validate_solution(x_cand, base_config, T_water_K_fixed, CONFIG);
    
    if is_valid
        valid_solutions{end+1} = struct(...
            'x', x_cand, ...
            'fval', fval_cand, ...
            'method', method_cand, ...
            'T_max', T_max_cand, ...
            'T_profile', T_profile_cand, ...
            'COP', COP_cand);
        fprintf('  ✓ %s: VALID (T_max = %.2f °C)\n', method_cand, T_max_cand - 273.15);
    else
        fprintf('  ✗ %s: INVALID - %s\n', method_cand, reason);
    end
end

% Select best valid solution
if isempty(valid_solutions)
    fprintf('\n⚠ WARNING: No physically valid solutions found!\n');
    fprintf('  All solutions violate physical constraints.\n');
    fprintf('  Consider adjusting bounds or checking the model.\n\n');
    
    % Fall back to least bad solution
    x_best = x_local;
    fval_best = fval_local;
    method_best = 'Local Refinement (unvalidated)';
    [T_max, T_profile, Q_in, Q_out, COP] = evaluate_design(x_best, base_config, CONFIG);
else
    % Find best among valid solutions
    best_idx = 1;
    best_Tmax = valid_solutions{1}.T_max;
    for i = 2:length(valid_solutions)
        if valid_solutions{i}.T_max < best_Tmax
            best_idx = i;
            best_Tmax = valid_solutions{i}.T_max;
        end
    end
    
    x_best = valid_solutions{best_idx}.x;
    fval_best = valid_solutions{best_idx}.fval;
    method_best = valid_solutions{best_idx}.method;
    T_max = valid_solutions{best_idx}.T_max;
    T_profile = valid_solutions{best_idx}.T_profile;
    COP = valid_solutions{best_idx}.COP;
    
    % Recalculate Q values for the best solution
    [~, ~, Q_in, Q_out, ~] = evaluate_design(x_best, base_config, CONFIG);
    
    fprintf('\n✓ Selected best VALID solution from: %s\n', method_best);
end

fprintf('\nDetailed Results:\n');
fprintf('  T_max:    %.2f °C\n', T_max - 273.15);
fprintf('  T_center: %.2f °C\n', T_profile(1) - 273.15);
fprintf('  T_edge:   %.2f °C\n', T_profile(end) - 273.15);
fprintf('  ΔT:       %.2f K\n', T_profile(1) - T_profile(end));
fprintf('  Q_in:     %.4f W\n', Q_in);
fprintf('  Q_out:    %.4f W\n', Q_out);
fprintf('  COP:      %.3f\n', COP);

% Temperature profile sanity check
fprintf('\nTemperature Profile (should decrease from center to edge):\n');
for i = 1:length(T_profile)
    fprintf('  Node %d: %.2f °C\n', i, T_profile(i) - 273.15);
end

%% Save results
results = struct();
results.timestamp = timestamp;
results.config = CONFIG;
results.variables = var_names;
results.bounds.lb = lb;
results.bounds.ub = ub;

results.ga.x = x_ga;
results.ga.fval = fval_ga;
results.ga.time = time_ga;
results.ga.output = output_ga;

results.pso.x = x_pso;
results.pso.fval = fval_pso;
results.pso.time = time_pso;
results.pso.output = output_pso;

results.local.x = x_local;
results.local.fval = fval_local;
results.local.time = time_local;

results.best.x = x_best;
results.best.fval = fval_best;
results.best.method = method_best;
results.best.T_profile = T_profile;
results.best.COP = COP;

save(fullfile(OUTPUT_DIR, 'global_optimization_results.mat'), 'results');

% Save to CSV for easy viewing - include all parameters
% Fix: var_names is already a column cell array, use it directly
valid_var_names = matlab.lang.makeValidName(var_names);
csv_data = array2table(x_best(:)', 'VariableNames', valid_var_names(:)');
csv_data.T_max_C = fval_best;
csv_data.COP = COP;
writetable(csv_data, fullfile(OUTPUT_DIR, 'optimal_design.csv'));

% Also save a human-readable summary
fid = fopen(fullfile(OUTPUT_DIR, 'optimal_design_summary.txt'), 'w');
fprintf(fid, 'Global Optimization Results\n');
fprintf(fid, '===========================\n');
fprintf(fid, 'Timestamp: %s\n', timestamp);
fprintf(fid, 'Method: %s\n\n', method_best);
fprintf(fid, 'Fixed Parameters:\n');
fprintf(fid, '  N_stages: %d\n', CONFIG.N_stages);
fprintf(fid, '  N_tsv_limit: %d\n', CONFIG.N_tsv_limit);
fprintf(fid, '  T_target: %.0f C\n\n', CONFIG.T_target_C);
fprintf(fid, 'Optimized Parameters:\n');
for i = 1:nvars
    fprintf(fid, '  %s: %.6f\n', var_names{i}, x_best(i));
end
fprintf(fid, '\nResults:\n');
fprintf(fid, '  T_max: %.2f C\n', fval_best);
fprintf(fid, '  COP: %.4f\n', COP);
fclose(fid);

fprintf('\nResults saved to: %s\n', OUTPUT_DIR);

%% Plot convergence comparison
figure('Position', [100, 100, 1000, 400]);

subplot(1,2,1);
bar([fval_ga, fval_pso, fval_local]);
set(gca, 'XTickLabel', {'GA', 'PSO', 'fmincon'});
ylabel('T_{max} (°C)');
title('Optimization Method Comparison');
grid on;

subplot(1,2,2);
bar([time_ga, time_pso, time_local]);
set(gca, 'XTickLabel', {'GA', 'PSO', 'fmincon'});
ylabel('Time (s)');
title('Computation Time');
grid on;

saveas(gcf, fullfile(OUTPUT_DIR, 'optimization_comparison.png'));

%% Plot temperature profile for best solution using ResultsManager
fprintf('\nGenerating temperature profile plot...\n');

% Create a full config for the best solution to get geometry
best_config = build_full_config(x_best, base_config, CONFIG);
best_geometry = TECGeometry(best_config);

% Create a temporary ResultsManager pointing to our output directory
rm = ResultsManager(best_config);
% Override output directory to use our global optimization folder
rm.OutputDir = OUTPUT_DIR;

% T_water is fixed from CONFIG
T_water_best = CONFIG.T_water_K;

% Plot and save the temperature profile
rm.plot_temperature_profile(T_profile, best_geometry, ...
    'best_solution_temp_profile.png', ...
    sprintf('Best Solution (T_{max}=%.1f°C, COP=%.3f)', fval_best, COP), ...
    '', T_water_best);

fprintf('Temperature profile saved to: %s\n', fullfile(OUTPUT_DIR, 'best_solution_temp_profile.png'));

% Also display the figure
h = figure('Name', 'Global Optimization - Best Solution Temperature Profile');
N = best_geometry.N_stages;
T_0 = T_profile(1);
T_Si = T_profile(2:N+1);
T_c = T_profile(N+2:end);

r_chip = 0:N;
T_chip = [T_0; T_Si];
r_tec = 1:(N+1);
T_tec = [T_c; T_water_best];

plot(r_chip, T_chip - 273.15, '-ob', 'LineWidth', 1.8, 'MarkerFaceColor', 'b', 'DisplayName', 'Silicon Layer'); hold on;
plot(r_tec, T_tec - 273.15, '-sr', 'LineWidth', 1.8, 'MarkerFaceColor', 'r', 'DisplayName', 'TEC Cold Side');
plot(N+1, T_water_best - 273.15, 'g^', 'MarkerSize', 10, 'MarkerFaceColor', 'g', 'DisplayName', 'Coolant');

xlabel('Stage Index');
ylabel('Temperature (°C)');
title(sprintf('Best Solution Temperature Profile\nMethod: %s | T_{max}=%.1f°C | COP=%.3f', ...
    method_best, fval_best, COP));
legend('Location', 'best');
grid on;

saveas(h, fullfile(OUTPUT_DIR, 'best_solution_temp_profile_celsius.png'));
fprintf('Temperature profile (Celsius) saved.\n');

fprintf('\n✓ Global optimization complete!\n');

%% ==================== HELPER FUNCTIONS ====================

function config = create_base_config(CONFIG)
    % Create base configuration with FIXED parameters only
    % Optimizable parameters will be set in evaluate_design
    
    config = struct();
    
    % Fixed geometry parameters
    config.geometry.N_stages = CONFIG.N_stages;
    config.geometry.N_tsv_limit = CONFIG.N_tsv_limit;
    config.geometry.w_chip_um = CONFIG.ref.w_chip_um;
    
    % TSV parameters (fixed)
    config.geometry.tsv.R_TSV_um = 10;
    config.geometry.tsv.P_TSV_um = 20;
    config.geometry.tsv.g_rad_um = 10;
    config.geometry.tsv.t_SOI_um = 100;
    
    % Reference values for ratio calculations
    config.ref = CONFIG.ref;
    
    % Material properties (fixed)
    config.materials.Bi2Te3 = struct('k', 1.2, 'rho', 1e-5, 'S', 0.0002);
    config.materials.Cu = struct('k', 400, 'rho', 1.7e-8);
    config.materials.Si = struct('k', 150, 'rho', 0.01);
    config.materials.AlN = struct('k', 170, 'rho', 1e10);
    config.materials.SiO2 = struct('k', 1.4, 'rho', 1e14);
    config.materials.Al2O3 = struct('k', 30, 'rho', 1e12);
end

function config = build_full_config(x, base_config, CONFIG)
    % Build a complete config struct from optimization vector x
    % Uses the same variable names as optimization_variables.m
    
    config = base_config;
    
    % Get variable names from CONFIG
    all_vars = CONFIG.all_vars;
    enabled_mask = [all_vars{:, 5}];
    
    % Create a map of variable values: enabled get x values, disabled get initial values
    var_map = containers.Map();
    x_idx = 1;
    for i = 1:size(all_vars, 1)
        name = all_vars{i, 1};
        if enabled_mask(i)  % enabled
            var_map(name) = x(x_idx);
            x_idx = x_idx + 1;
        else  % disabled - use initial value
            var_map(name) = all_vars{i, 4};
        end
    end
    
    % Helper function to get value (with default)
    get_val = @(name, default) get_var_or_default(var_map, name, default);
    
    % Operating conditions - variable name is 'current' (matching optimization_variables.m)
    config.operating_conditions.I_current_A = get_val('current', 0.025);
    
    % Boundary conditions - FIXED from CONFIG
    config.boundary_conditions.q_flux_W_m2 = CONFIG.q_flux_W_m2;
    config.boundary_conditions.h_conv_W_m2K = CONFIG.h_conv_W_m2K;
    config.boundary_conditions.T_water_K = CONFIG.T_water_K;
    
    % Geometry - main parameters (direct values, not ratios)
    config.geometry.thickness_um = get_val('thickness_um', 200);
    config.geometry.wedge_angle_deg = get_val('wedge_angle_deg', 30);
    config.geometry.R_cyl_um = get_val('R_cyl_um', 1000);
    
    % TSV/SOI
    if isfield(config.geometry, 'tsv')
        config.geometry.tsv.t_SOI_um = get_val('t_SOI_um', 100);
    end
    
    % Expansion and fill
    config.geometry.radial_expansion_factor = get_val('k_r', 1.15);
    config.geometry.fill_factor = get_val('fill_factor', 0.95);
    config.geometry.insulation_width_ratio = get_val('insulation_width_ratio', 0.04);
    
    % Interconnects
    config.geometry.interconnect_ratio = get_val('interconnect_ratio', 0.15);
    config.geometry.outerconnect_ratio = get_val('outerconnect_ratio', 0.15);
    config.geometry.interconnect_angle_ratio = get_val('interconnect_angle_ratio', 0.16);
    config.geometry.outerconnect_angle_ratio = get_val('outerconnect_angle_ratio', 0.16);
    config.geometry.interconnect_thickness_ratio = get_val('interconnect_thickness_ratio', 1.0);
    config.geometry.outerconnect_thickness_ratio = get_val('outerconnect_thickness_ratio', 1.0);
end

function val = get_var_or_default(var_map, name, default)
    if var_map.isKey(name)
        val = var_map(name);
    else
        val = default;
    end
end

function cost = tec_objective(x, base_config, CONFIG)
    % Objective function for optimization
    % Dynamically handles enabled/disabled variables from CONFIG.all_vars
    %
    % FIXED from CONFIG: q_flux, h_conv, T_water, N_stages, N_tsv_limit, T_target
    %
    % Includes guardrails for physically unrealistic results
    
    PENALTY = 1e6;  % Large penalty for infeasible designs
    T_target_C = CONFIG.T_target_C;
    T_water_K = CONFIG.T_water_K;  % Fixed from CONFIG
    
    try
        [T_max_K, T_profile, Q_in, Q_out, COP] = evaluate_design(x, base_config, CONFIG);
        T_max_C = T_max_K - 273.15;
        T_min_K = min(T_profile);
        T_center_K = T_profile(1);  % Center/hot side temperature
        T_edge_K = T_profile(end);  % Edge/cold side temperature
        
        % ============ PHYSICAL FEASIBILITY CHECKS ============
        
        % Check 1: No negative temperatures (absolute zero violation)
        if any(T_profile < 0)
            cost = PENALTY;
            return;
        end
        
        % Check 2: All temperatures must be above absolute zero with margin
        if T_min_K < 200  % Below -73°C is unrealistic for this application
            cost = PENALTY;
            return;
        end
        
        % Check 3: Center should be hotter than edge (heat flows outward)
        % Heat is applied at center, removed at edge
        if T_center_K < T_edge_K
            cost = PENALTY;
            return;
        end
        
        % Check 4: Edge temperature should be near or above coolant temp
        % TEC can cool below ambient, but not excessively
        T_min_allowed = T_water_K - 50;  % Max 50K subcooling is realistic
        if T_edge_K < T_min_allowed
            cost = PENALTY;
            return;
        end
        
        % Check 5: Temperature must not be unrealistically high
        if T_max_C > 500  % Above 500°C is unrealistic
            cost = PENALTY;
            return;
        end
        
        % Check 6: COP sanity check (should be positive for cooling)
        if COP < -1  % Negative COP indicates heating, not cooling
            cost = PENALTY;
            return;
        end
        
        % Check 7: Energy balance - Q_out should be greater than Q_in
        if Q_out < Q_in * 0.5  % Gross energy imbalance
            cost = PENALTY;
            return;
        end
        
        % Check 8: Monotonic temperature profile (heat should flow outward)
        % For radial TEC, temperature should generally decrease from center to edge
        % Allow small violations due to TEC pumping effects
        for i = 2:length(T_profile)
            if T_profile(i) > T_profile(i-1) + 5  % Allow 5K tolerance
                cost = PENALTY;
                return;
            end
        end
        
        % ============ OBJECTIVE CALCULATION ============
        
        % Primary objective: minimize T_max
        cost = T_max_C;
        
        % Soft penalty for exceeding target (gradual, not hard cutoff)
        if T_max_C > T_target_C
            cost = cost + 10 * (T_max_C - T_target_C);
        end
        
        % Bonus for good COP (encourages efficient designs)
        if COP > 0.1
            cost = cost - 0.5 * min(COP, 2);  % Small bonus, capped
        end
        
    catch ME
        % Debug: uncomment to see what's failing
        % fprintf('Objective failed: %s\n', ME.message);
        cost = PENALTY;  % Penalty for failed simulations
    end
end

function [T_max, T_profile, Q_in, Q_out, COP] = evaluate_design(x, base_config, CONFIG)
    % Evaluate a TEC design with dynamically configured optimization parameters
    % Uses CONFIG.all_vars to map x vector to parameter names
    % Uses same variable names as optimization_variables.m
    
    config = base_config;
    
    % Get variable names from CONFIG
    all_vars = CONFIG.all_vars;
    enabled_mask = [all_vars{:, 5}];
    
    % Create a map of variable values: enabled get x values, disabled get initial values
    var_map = containers.Map();
    x_idx = 1;
    for i = 1:size(all_vars, 1)
        name = all_vars{i, 1};
        if enabled_mask(i)  % enabled
            var_map(name) = x(x_idx);
            x_idx = x_idx + 1;
        else  % disabled - use initial value
            var_map(name) = all_vars{i, 4};
        end
    end
    
    % Helper function to get value (with default)
    get_val = @(name, default) get_var_or_default(var_map, name, default);
    
    % Operating conditions - variable name is 'current' (matching optimization_variables.m)
    config.operating_conditions.I_current_A = get_val('current', 0.025);
    
    % Boundary conditions - FIXED from CONFIG
    config.boundary_conditions.q_flux_W_m2 = CONFIG.q_flux_W_m2;
    config.boundary_conditions.h_conv_W_m2K = CONFIG.h_conv_W_m2K;
    config.boundary_conditions.T_water_K = CONFIG.T_water_K;
    
    % Geometry - main parameters (direct values, not ratios)
    config.geometry.thickness_um = get_val('thickness_um', 200);
    config.geometry.wedge_angle_deg = get_val('wedge_angle_deg', 30);
    config.geometry.R_cyl_um = get_val('R_cyl_um', 1000);
    
    % TSV/SOI
    if isfield(config.geometry, 'tsv')
        config.geometry.tsv.t_SOI_um = get_val('t_SOI_um', 100);
    end
    
    % Expansion and fill
    config.geometry.radial_expansion_factor = get_val('k_r', 1.15);
    config.geometry.fill_factor = get_val('fill_factor', 0.95);
    config.geometry.insulation_width_ratio = get_val('insulation_width_ratio', 0.04);
    
    % Interconnects
    config.geometry.interconnect_ratio = get_val('interconnect_ratio', 0.15);
    config.geometry.outerconnect_ratio = get_val('outerconnect_ratio', 0.15);
    config.geometry.interconnect_angle_ratio = get_val('interconnect_angle_ratio', 0.16);
    config.geometry.outerconnect_angle_ratio = get_val('outerconnect_angle_ratio', 0.16);
    config.geometry.interconnect_thickness_ratio = get_val('interconnect_thickness_ratio', 1.0);
    config.geometry.outerconnect_thickness_ratio = get_val('outerconnect_thickness_ratio', 1.0);
    
    % Create solver objects
    materials = MaterialProperties(config);
    geometry = TECGeometry(config);
    network = ThermalNetwork(geometry, materials, config);
    
    N = geometry.N_stages;
    T_water = config.boundary_conditions.T_water_K;
    
    % CRITICAL: Use warm initial guess (like TECOptimizer does)
    % Starting from T_water causes convergence issues
    T = ones(2*N + 1, 1) * (T_water + 50);
    
    % CRITICAL: Use relaxation iteration (like TECOptimizer does)
    % Without relaxation, solver can diverge for extreme parameters
    for iter = 1:100
        T_old = T;
        try
            [T_new, Q_out, Q_in] = network.solve(T);
        catch
            % Solver failed - return penalty values
            T_max = 1e6;
            T_profile = T;
            Q_in = 0;
            Q_out = 0;
            COP = 0;
            return;
        end
        
        % Check for invalid values
        if any(isnan(T_new)) || any(isinf(T_new)) || any(T_new < 0)
            T_max = 1e6;
            T_profile = T;
            Q_in = 0;
            Q_out = 0;
            COP = 0;
            return;
        end
        
        % RELAXATION: Blend old and new (critical for stability)
        T = 0.5 * T_new + 0.5 * T;
        
        if max(abs(T - T_old)) < 1e-6
            break;
        end
    end
    
    T_max = max(T);
    T_profile = T;
    
    % COP = Q_in / (Q_out - Q_in)
    P_elec = Q_out - Q_in;
    if P_elec > 0
        COP = Q_in / P_elec;
    else
        COP = 0;
    end
end

function [state, options, optchanged] = ga_output(options, state, flag, ~)
    % Custom output function for GA
    % var_names argument replaced with ~ as it's unused
    optchanged = false;
    
    if strcmp(flag, 'iter')
        % Could log or save intermediate results here
    end
end

function [is_valid, reason, T_max, T_profile, COP] = validate_solution(x, base_config, T_water_K, CONFIG)
    % Validate a solution for physical feasibility
    % Returns: is_valid (bool), reason (string), and thermal results
    %
    % Physical constraints checked:
    % 1. No negative temperatures
    % 2. Temperatures above reasonable minimum (200 K)
    % 3. Center hotter than edge (correct heat flow direction)
    % 4. Edge temperature reasonable relative to coolant
    % 5. Monotonic temperature profile (no inversions)
    % 6. Reasonable COP
    
    is_valid = false;
    reason = 'Unknown error';  % Default reason
    T_max = NaN;
    T_profile = [];
    COP = NaN;
    
    try
        [T_max, T_profile, Q_in, Q_out, COP] = evaluate_design(x, base_config, CONFIG);
        
        T_min = min(T_profile);
        T_center = T_profile(1);
        T_edge = T_profile(end);
        
        % Check 1: No negative temperatures
        if any(T_profile < 0)
            reason = 'Negative temperature detected';
            return;
        end
        
        % Check 2: Minimum temperature check
        if T_min < 200  % Below -73°C
            reason = sprintf('Temperature too low (%.1f K)', T_min);
            return;
        end
        
        % Check 3: Center should be hotter than edge
        if T_center < T_edge
            reason = sprintf('Inverted profile: center (%.1f K) < edge (%.1f K)', T_center, T_edge);
            return;
        end
        
        % Check 4: Edge temperature relative to coolant
        T_min_allowed = T_water_K - 50;  % Max 50K subcooling
        if T_edge < T_min_allowed
            reason = sprintf('Edge too cold (%.1f K < %.1f K limit)', T_edge, T_min_allowed);
            return;
        end
        
        % Check 5: Maximum temperature check
        if T_max > 773.15  % Above 500°C
            reason = sprintf('Temperature too high (%.1f °C)', T_max - 273.15);
            return;
        end
        
        % Check 6: Monotonic profile check (allow small violations)
        for i = 2:length(T_profile)
            if T_profile(i) > T_profile(i-1) + 10  % 10K tolerance
                reason = sprintf('Non-monotonic profile at node %d', i);
                return;
            end
        end
        
        % Check 7: COP sanity
        if COP < -0.5
            reason = sprintf('Invalid COP (%.3f)', COP);
            return;
        end
        
        % Check 8: Energy balance
        if Q_out < 0 || Q_in < 0
            reason = 'Negative heat flow';
            return;
        end
        
        % All checks passed
        is_valid = true;
        reason = 'Valid';
        
    catch ME
        reason = sprintf('Evaluation failed: %s', ME.message);
    end
end
