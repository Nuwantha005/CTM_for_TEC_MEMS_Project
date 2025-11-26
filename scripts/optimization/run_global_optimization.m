% run_global_optimization.m
% Global optimization using Genetic Algorithm and Particle Swarm
%
% This finds the globally optimal TEC design by exploring the full
% parameter space, avoiding local minima that gradient-based methods miss.
%
% Uses: Global Optimization Toolbox (ga, particleswarm)

clear; clc;
addpath(genpath('../../src'));

fprintf('╔════════════════════════════════════════════════════════════╗\n');
fprintf('║           GLOBAL OPTIMIZATION FOR TEC DESIGN               ║\n');
fprintf('╚════════════════════════════════════════════════════════════╝\n\n');

%% Configuration
CONFIG = struct();
CONFIG.q_flux_W_m2 = 1000;          % Heat flux to optimize for
CONFIG.T_target_C = 85;              % Target max temperature (°C)
CONFIG.N_stages = 3;                 % Fixed for COMSOL template
CONFIG.wedge_angle_deg = 30;         % Fixed for COMSOL template

% Output directory
timestamp = datestr(now, 'yyyy-mm-dd_HH-MM-SS');
OUTPUT_DIR = fullfile('../../output', 'global_optimization', timestamp);
if ~exist(OUTPUT_DIR, 'dir')
    mkdir(OUTPUT_DIR);
end

%% Define optimization variables and bounds
% Variables: [I_current (A), thickness (um), k_r, fill_factor]

var_names = {'I_current (A)', 'thickness (µm)', 'k_r', 'fill_factor'};
nvars = 4;

% Lower bounds
lb = [0.01,   50,  0.8, 0.70];  % 10 mA, 50 µm, k_r=0.8, ff=0.70

% Upper bounds  
ub = [0.20,  400,  1.5, 0.99];  % 200 mA, 400 µm, k_r=1.5, ff=0.99

fprintf('Optimization Variables:\n');
fprintf('─────────────────────────────────────────\n');
for i = 1:nvars
    fprintf('  %s: [%.2f, %.2f]\n', var_names{i}, lb(i), ub(i));
end
fprintf('─────────────────────────────────────────\n\n');

%% Create base configuration
base_config = create_base_config(CONFIG);

%% Define objective function
% Minimize T_max while penalizing infeasible designs

objective = @(x) tec_objective(x, base_config, CONFIG.T_target_C);

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
fprintf('─────────────────────────────────────────────────────────────\n');
fprintf('  Current:        %.1f mA\n', x_best(1) * 1000);
fprintf('  TEC Thickness:  %.0f µm\n', x_best(2));
fprintf('  k_r:            %.3f\n', x_best(3));
fprintf('  Fill Factor:    %.3f\n', x_best(4));
fprintf('  ─────────────────────\n');
fprintf('  T_max:          %.2f °C\n', fval_best);
fprintf('─────────────────────────────────────────────────────────────\n\n');

%% Validate best solution
fprintf('Validating best solution...\n');
[T_max, T_profile, Q_in, Q_out, COP] = evaluate_design(x_best, base_config);

fprintf('\nDetailed Results:\n');
fprintf('  T_max:    %.2f °C\n', T_max - 273.15);
fprintf('  T_center: %.2f °C\n', T_profile(1) - 273.15);
fprintf('  Q_in:     %.4f W\n', Q_in);
fprintf('  Q_out:    %.4f W\n', Q_out);
fprintf('  COP:      %.3f\n', COP);

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

% Save to CSV for easy viewing
T = table(x_best(1)*1000, x_best(2), x_best(3), x_best(4), fval_best, COP, ...
    'VariableNames', {'I_mA', 't_um', 'k_r', 'fill_factor', 'T_max_C', 'COP'});
writetable(T, fullfile(OUTPUT_DIR, 'optimal_design.csv'));

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

fprintf('\n✓ Global optimization complete!\n');

%% ==================== HELPER FUNCTIONS ====================

function config = create_base_config(CONFIG)
    config = struct();
    config.geometry.N_stages = CONFIG.N_stages;
    config.geometry.wedge_angle_deg = CONFIG.wedge_angle_deg;
    config.geometry.w_chip_um = 10000;
    config.geometry.R_cyl_um = 1000;
    config.geometry.t_chip_um = 50;
    config.geometry.interconnect_ratio = 0.15;
    config.geometry.outerconnect_ratio = 0.15;
    config.geometry.insulation_width_ratio = 0.04;
    config.geometry.interconnect_angle_ratio = 0.16;
    config.geometry.outerconnect_angle_ratio = 0.16;
    config.geometry.interconnect_thickness_ratio = 1.0;
    config.geometry.outerconnect_thickness_ratio = 1.0;
    config.geometry.tsv.R_TSV_um = 10;
    config.geometry.tsv.P_TSV_um = 20;
    config.geometry.tsv.g_rad_um = 10;
    config.geometry.tsv.t_SOI_um = 100;
    
    config.boundary_conditions.q_flux_W_m2 = CONFIG.q_flux_W_m2;
    config.boundary_conditions.T_water_K = 300;
    config.boundary_conditions.h_conv_W_m2K = 1e6;
    
    config.materials.Bi2Te3 = struct('k', 1.2, 'rho', 1e-5, 'S', 0.0002);
    config.materials.Cu = struct('k', 400, 'rho', 1.7e-8);
    config.materials.Si = struct('k', 150, 'rho', 0.01);
    config.materials.AlN = struct('k', 170, 'rho', 1e10);
    config.materials.SiO2 = struct('k', 1.4, 'rho', 1e14);
    config.materials.Al2O3 = struct('k', 30, 'rho', 1e12);
end

function cost = tec_objective(x, base_config, T_target_C)
    % Objective function for optimization
    % x = [I_current, thickness, k_r, fill_factor]
    
    try
        [T_max_K, ~, ~, ~, ~] = evaluate_design(x, base_config);
        T_max_C = T_max_K - 273.15;
        
        % Primary objective: minimize T_max
        cost = T_max_C;
        
        % Penalty for exceeding target
        if T_max_C > T_target_C
            cost = cost + 10 * (T_max_C - T_target_C);
        end
        
        % Penalty for negative or unrealistic temperatures
        if T_max_C < 0 || T_max_C > 500
            cost = 1e6;
        end
        
    catch
        cost = 1e6;  % Penalty for failed simulations
    end
end

function [T_max, T_profile, Q_in, Q_out, COP] = evaluate_design(x, base_config)
    % Evaluate a TEC design
    % x = [I_current, thickness, k_r, fill_factor]
    
    config = base_config;
    config.operating_conditions.I_current_A = x(1);
    config.geometry.thickness_um = x(2);
    config.geometry.radial_expansion_factor = x(3);
    config.geometry.fill_factor = x(4);
    
    materials = MaterialProperties(config);
    geometry = TECGeometry(config);
    network = ThermalNetwork(geometry, materials, config);
    
    N = geometry.N_stages;
    T = ones(2*N + 1, 1) * 300;
    
    for iter = 1:100
        T_old = T;
        [T, Q_out, Q_in] = network.solve(T);
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

function [state, options, optchanged] = ga_output(options, state, flag, var_names)
    % Custom output function for GA
    optchanged = false;
    
    if strcmp(flag, 'iter')
        % Could log or save intermediate results here
    end
end
