% run_multiobjective_optimization.m
% Multi-objective optimization: Minimize T_max AND Power consumption
%
% This finds the Pareto-optimal front showing trade-offs between
% cooling performance and electrical power consumption.
%
% Uses: Global Optimization Toolbox (gamultiobj)

clear; clc;
addpath(genpath('../../src'));

fprintf('╔════════════════════════════════════════════════════════════╗\n');
fprintf('║        MULTI-OBJECTIVE OPTIMIZATION FOR TEC DESIGN        ║\n');
fprintf('╚════════════════════════════════════════════════════════════╝\n\n');

%% Configuration
CONFIG = struct();
CONFIG.q_flux_W_m2 = 20000;          % Heat flux
CONFIG.N_stages = 3;                 % Fixed for COMSOL template
CONFIG.wedge_angle_deg = 30;         % Fixed for COMSOL template

% Output directory
timestamp = datestr(now, 'yyyy-mm-dd_HH-MM-SS');
OUTPUT_DIR = fullfile('../../output', 'multiobjective_optimization', timestamp);
if ~exist(OUTPUT_DIR, 'dir')
    mkdir(OUTPUT_DIR);
end

%% Define optimization variables and bounds
% Variables: [I_current (A), thickness (um), k_r, fill_factor]

var_names = {'I_current (A)', 'thickness (µm)', 'k_r', 'fill_factor'};
nvars = 4;

% Lower bounds
lb = [0.01,   50,  0.8, 0.70];

% Upper bounds  
ub = [0.20,  400,  1.5, 0.99];

fprintf('Optimization Variables:\n');
fprintf('─────────────────────────────────────────\n');
for i = 1:nvars
    fprintf('  %s: [%.2f, %.2f]\n', var_names{i}, lb(i), ub(i));
end
fprintf('\nObjectives:\n');
fprintf('  1. Minimize T_max (°C)\n');
fprintf('  2. Minimize Power consumption (W)\n');
fprintf('─────────────────────────────────────────\n\n');

%% Create base configuration
base_config = create_base_config(CONFIG);

%% Define multi-objective function
% Returns [T_max, Power]

multiobjective = @(x) tec_multiobjective(x, base_config);

%% Run Multi-Objective Genetic Algorithm
fprintf('=== RUNNING gamultiobj ===\n\n');
fprintf('Finding Pareto-optimal front...\n');
fprintf('Population: 100, Generations: 100\n\n');

options = optimoptions('gamultiobj', ...
    'PopulationSize', 100, ...
    'MaxGenerations', 100, ...
    'ParetoFraction', 0.35, ...
    'Display', 'iter', ...
    'PlotFcn', {@gaplotpareto}, ...
    'UseParallel', false);

tic;
[x_pareto, fval_pareto, exitflag, output] = ...
    gamultiobj(multiobjective, nvars, [], [], [], [], lb, ub, options);
opt_time = toc;

fprintf('\nOptimization completed in %.1f seconds\n', opt_time);
fprintf('Found %d Pareto-optimal solutions\n\n', size(x_pareto, 1));

%% Analyze Pareto Front
fprintf('╔════════════════════════════════════════════════════════════╗\n');
fprintf('║                    PARETO FRONT ANALYSIS                   ║\n');
fprintf('╚════════════════════════════════════════════════════════════╝\n\n');

T_max_pareto = fval_pareto(:,1);
Power_pareto = fval_pareto(:,2);

% Find extreme points
[~, idx_min_T] = min(T_max_pareto);
[~, idx_min_P] = min(Power_pareto);

% Find knee point (balanced solution) using distance to utopia point
utopia = [min(T_max_pareto), min(Power_pareto)];
distances = sqrt(((T_max_pareto - utopia(1))/range(T_max_pareto)).^2 + ...
                 ((Power_pareto - utopia(2))/range(Power_pareto)).^2);
[~, idx_knee] = min(distances);

fprintf('Key Solutions on Pareto Front:\n');
fprintf('─────────────────────────────────────────────────────────────────────────────────\n');
fprintf('%-15s | %8s | %8s | %8s | %8s | %8s | %8s\n', ...
    'Solution', 'I (mA)', 't (µm)', 'k_r', 'ff', 'T_max(°C)', 'Power(W)');
fprintf('─────────────────────────────────────────────────────────────────────────────────\n');

% Min Temperature
fprintf('%-15s | %8.1f | %8.0f | %8.3f | %8.3f | %8.2f | %8.4f\n', ...
    'Min T_max', x_pareto(idx_min_T,1)*1000, x_pareto(idx_min_T,2), ...
    x_pareto(idx_min_T,3), x_pareto(idx_min_T,4), ...
    T_max_pareto(idx_min_T), Power_pareto(idx_min_T));

% Knee Point (Balanced)
fprintf('%-15s | %8.1f | %8.0f | %8.3f | %8.3f | %8.2f | %8.4f\n', ...
    'Balanced', x_pareto(idx_knee,1)*1000, x_pareto(idx_knee,2), ...
    x_pareto(idx_knee,3), x_pareto(idx_knee,4), ...
    T_max_pareto(idx_knee), Power_pareto(idx_knee));

% Min Power
fprintf('%-15s | %8.1f | %8.0f | %8.3f | %8.3f | %8.2f | %8.4f\n', ...
    'Min Power', x_pareto(idx_min_P,1)*1000, x_pareto(idx_min_P,2), ...
    x_pareto(idx_min_P,3), x_pareto(idx_min_P,4), ...
    T_max_pareto(idx_min_P), Power_pareto(idx_min_P));

fprintf('─────────────────────────────────────────────────────────────────────────────────\n\n');

%% Plot Pareto Front
figure('Position', [100, 100, 1200, 500], 'Name', 'Pareto Front Analysis');

% Pareto front with key points
subplot(1,2,1);
scatter(T_max_pareto, Power_pareto * 1000, 50, 'b', 'filled', 'MarkerFaceAlpha', 0.6);
hold on;
scatter(T_max_pareto(idx_min_T), Power_pareto(idx_min_T)*1000, 200, 'g', 'filled', 'MarkerEdgeColor', 'k', 'LineWidth', 2);
scatter(T_max_pareto(idx_knee), Power_pareto(idx_knee)*1000, 200, 'm', 'filled', 'MarkerEdgeColor', 'k', 'LineWidth', 2);
scatter(T_max_pareto(idx_min_P), Power_pareto(idx_min_P)*1000, 200, 'r', 'filled', 'MarkerEdgeColor', 'k', 'LineWidth', 2);

xlabel('Maximum Temperature (°C)', 'FontSize', 12);
ylabel('Power Consumption (mW)', 'FontSize', 12);
title('Pareto Front: T_{max} vs Power', 'FontSize', 14);
legend('Pareto solutions', 'Min T_{max}', 'Balanced (Knee)', 'Min Power', 'Location', 'northeast');
grid on;

% Temperature constraint line
xline(85, 'r--', 'T_{max} = 85°C', 'LineWidth', 1.5, 'LabelHorizontalAlignment', 'left');

% Parameter distributions
subplot(1,2,2);
param_labels = {'I (mA)', 't (µm)', 'k_r', 'fill'};
param_data = [x_pareto(:,1)*1000, x_pareto(:,2), x_pareto(:,3), x_pareto(:,4)];

% Normalize for visualization
param_norm = (param_data - min(param_data)) ./ (max(param_data) - min(param_data));

boxplot(param_norm, 'Labels', param_labels);
ylabel('Normalized Value', 'FontSize', 12);
title('Parameter Distributions on Pareto Front', 'FontSize', 14);
grid on;

saveas(gcf, fullfile(OUTPUT_DIR, 'pareto_front.png'));
saveas(gcf, fullfile(OUTPUT_DIR, 'pareto_front.fig'));

%% Additional Analysis: Trade-off Curves
figure('Position', [100, 100, 1200, 800], 'Name', 'Trade-off Analysis');

% T_max vs Current
subplot(2,2,1);
scatter(x_pareto(:,1)*1000, T_max_pareto, 50, Power_pareto*1000, 'filled');
xlabel('Current (mA)');
ylabel('T_{max} (°C)');
title('Temperature vs Current');
colorbar; ylabel(colorbar, 'Power (mW)');
grid on;

% T_max vs Thickness
subplot(2,2,2);
scatter(x_pareto(:,2), T_max_pareto, 50, Power_pareto*1000, 'filled');
xlabel('TEC Thickness (µm)');
ylabel('T_{max} (°C)');
title('Temperature vs Thickness');
colorbar; ylabel(colorbar, 'Power (mW)');
grid on;

% Power vs Current
subplot(2,2,3);
scatter(x_pareto(:,1)*1000, Power_pareto*1000, 50, T_max_pareto, 'filled');
xlabel('Current (mA)');
ylabel('Power (mW)');
title('Power vs Current');
colorbar; ylabel(colorbar, 'T_{max} (°C)');
grid on;

% 3D view
subplot(2,2,4);
scatter3(x_pareto(:,1)*1000, x_pareto(:,2), T_max_pareto, 50, Power_pareto*1000, 'filled');
xlabel('Current (mA)');
ylabel('Thickness (µm)');
zlabel('T_{max} (°C)');
title('3D Parameter Space');
colorbar; ylabel(colorbar, 'Power (mW)');
grid on;
view(45, 30);

saveas(gcf, fullfile(OUTPUT_DIR, 'tradeoff_analysis.png'));
saveas(gcf, fullfile(OUTPUT_DIR, 'tradeoff_analysis.fig'));

%% Save Results
results = struct();
results.timestamp = timestamp;
results.config = CONFIG;
results.variables = var_names;
results.bounds.lb = lb;
results.bounds.ub = ub;
results.optimization_time = opt_time;

results.pareto.x = x_pareto;
results.pareto.fval = fval_pareto;
results.pareto.T_max = T_max_pareto;
results.pareto.Power = Power_pareto;

results.key_solutions.min_temp.x = x_pareto(idx_min_T,:);
results.key_solutions.min_temp.T_max = T_max_pareto(idx_min_T);
results.key_solutions.min_temp.Power = Power_pareto(idx_min_T);

results.key_solutions.balanced.x = x_pareto(idx_knee,:);
results.key_solutions.balanced.T_max = T_max_pareto(idx_knee);
results.key_solutions.balanced.Power = Power_pareto(idx_knee);

results.key_solutions.min_power.x = x_pareto(idx_min_P,:);
results.key_solutions.min_power.T_max = T_max_pareto(idx_min_P);
results.key_solutions.min_power.Power = Power_pareto(idx_min_P);

save(fullfile(OUTPUT_DIR, 'multiobjective_results.mat'), 'results');

% Save Pareto front to CSV
T = table(x_pareto(:,1)*1000, x_pareto(:,2), x_pareto(:,3), x_pareto(:,4), ...
    T_max_pareto, Power_pareto*1000, ...
    'VariableNames', {'I_mA', 't_um', 'k_r', 'fill_factor', 'T_max_C', 'Power_mW'});
writetable(T, fullfile(OUTPUT_DIR, 'pareto_front.csv'));

fprintf('\nResults saved to: %s\n', OUTPUT_DIR);
fprintf('\n✓ Multi-objective optimization complete!\n');

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

function f = tec_multiobjective(x, base_config)
    % Multi-objective function
    % Returns: [T_max (°C), Power (W)]
    
    try
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
        
        T_max_C = max(T) - 273.15;
        Power_W = Q_out - Q_in;  % Electrical power input
        
        % Penalize invalid solutions
        if T_max_C < 0 || T_max_C > 500 || Power_W < 0
            f = [1e6, 1e6];
        else
            f = [T_max_C, Power_W];
        end
        
    catch
        f = [1e6, 1e6];  % Penalty for failed simulations
    end
end
