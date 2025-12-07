%% Compare MATLAB CTM vs COMSOL FEM
% Using actual COMSOL parametric sweep data from comsol_param_sweep.txt
% Single 30° wedge comparison

clear all; close all; clc;

% Get script directory for relative paths
script_dir = fileparts(mfilename('fullpath'));
project_root = fullfile(script_dir, '..', '..');
addpath(genpath(fullfile(project_root, 'src')));

fprintf('╔═══════════════════════════════════════════════════════════════╗\n');
fprintf('║     MATLAB vs COMSOL COMPARISON (Single Wedge)               ║\n');
fprintf('╚═══════════════════════════════════════════════════════════════╝\n\n');

%% Load COMSOL parametric study results from file
% Data from: data/comsol_param_sweep.txt
% Columns: I_0 (A), T_chip_avg (K), T_max (K), Q_in (W), Q_out (W)
comsol_data = struct();
comsol_data.I = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10];  % Current (A)
comsol_data.T_chip_avg = [524.55, 509.85, 504.63, 507.58, 517.75, 534.20, 556.44, 583.70, 615.56, 651.64, 691.57];  % K
comsol_data.T_max = [701.65, 657.18, 631.90, 623.74, 632.15, 653.98, 685.80, 725.61, 772.29, 824.87, 882.49];  % K
comsol_data.Q_in = [2.6179, 2.6180, 2.6179, 2.6179, 2.6180, 2.6179, 2.6181, 2.6184, 2.6188, 2.6191, 2.6174];  % W
comsol_data.Q_out = [2.6166, 2.5991, 2.6809, 2.8543, 3.1169, 3.4629, 3.8898, 4.3942, 4.9748, 5.6299, 6.3587];  % W

fprintf('COMSOL Data Loaded:\n');
fprintf('  At I=0 A: T_chip = %.2f K, T_max = %.2f K, Q_in = %.4f W\n', ...
        comsol_data.T_chip_avg(1), comsol_data.T_max(1), comsol_data.Q_in(1));
fprintf('  T_water = 300 K (assumed)\n\n');

%% Setup MATLAB model with base configuration
[~, ~, ~, ~, ~, CONFIG] = optimization_variables();

% Create base configuration matching COMSOL model
% Parameters from data/COMSOL_parameters.txt
config = struct();
config.geometry.N_stages = 3;
config.geometry.wedge_angle_deg = 30;           % LL_theta = 30
config.geometry.w_chip_um = 10000;              % 10 mm chip (derived from geometry)
config.geometry.t_chip_um = 50;                 % LL_t_chip = 50
config.geometry.R_cyl_um = 1000;                % LL_R_cyl = 1000
config.geometry.thickness_um = 1000;            % LL_t_TEC = 1000 (NOT 500!)
config.geometry.t_ins_um = 10;                  % Vertical insulator
config.geometry.radial_expansion_factor = 1.15; % LL_k_r = 1.15
config.geometry.azimuthal_gap_um = 20;          % W_az
config.geometry.insulation_width_ratio = 0.05;  % W_is/L (LL_w_is = 50, L~1000 → ratio ~ 0.05)
config.geometry.fill_factor = 0.9;              % LL_fill_factor = 0.9
config.geometry.interconnect_ratio = 0.1;       % LL_ic_w_r = 0.1
config.geometry.interconnect_thickness_ratio = 0.6;  % LL_ic_t_r = 0.6
config.geometry.interconnect_angle_ratio = 0.5;      % LL_ic_angle_r = 0.5
config.geometry.outerconnect_ratio = 0.1;            % LL_oc_w_r = 0.1
config.geometry.outerconnect_thickness_ratio = 0.6;  % LL_oc_t_r = 0.6
config.geometry.outerconnect_angle_ratio = 0.5;      % LL_oc_angle_r = 0.5

config.boundary_conditions.T_water_K = 300;
config.boundary_conditions.h_conv_W_m2K = 1e6;
config.boundary_conditions.q_flux_W_m2 = 200e3;
config.operating_conditions.I_current_A = 0;
config.materials = CONFIG.materials;

% Initialize components
mat_props = MaterialProperties(config);
geom = TECGeometry(config);

fprintf('MATLAB Geometry:\n');
fprintf('  N_stages = %d, θ = %.1f°\n', geom.N_stages, config.geometry.wedge_angle_deg);
fprintf('  R_base = %.2f mm, R_cyl = %.2f mm\n', geom.R_base*1e3, geom.R_cyl*1e3);
fprintf('  t_TEC = %.0f µm, L_1 = %.1f µm\n', geom.Thickness*1e6, geom.L_1*1e6);
fprintf('\n');

%% Loop over current sweep
n_current = length(comsol_data.I);
T_center_matlab = zeros(1, n_current);
T_max_matlab = zeros(1, n_current);
Q_out_matlab = zeros(1, n_current);
Q_in_matlab = zeros(1, n_current);

fprintf('Running MATLAB thermal network solver...\n');

for idx = 1:n_current
    I_current = comsol_data.I(idx);
    
    % Update current in config
    config.operating_conditions.I_current_A = I_current;
    
    % Rebuild thermal network with new current
    thermal_net = ThermalNetwork(geom, mat_props, config);
    
    N = geom.N_stages;
    T_water = config.boundary_conditions.T_water_K;
    
    % Warm initial guess
    T = ones(2*N + 1, 1) * (T_water + 100);
    
    % Solve with relaxation
    converged = false;
    for iter = 1:200
        T_old = T;
        try
            [T_new, Q_out_calc, Q_in_calc] = thermal_net.solve(T);
        catch ME
            fprintf('  I=%.0f A: Solver error: %s\n', I_current, ME.message);
            break;
        end
        
        if any(isnan(T_new)) || any(isinf(T_new)) || any(T_new < 0)
            break;
        end
        
        % Relaxation
        T = 0.3 * T_new + 0.7 * T;
        
        if max(abs(T - T_old)) < 1e-4
            converged = true;
            break;
        end
    end
    
    % Extract results
    T_center_matlab(idx) = T(1);  % Center temperature (T_h_0)
    T_max_matlab(idx) = max(T);
    Q_out_matlab(idx) = Q_out_calc;
    Q_in_matlab(idx) = Q_in_calc;
    
    status = '';
    if ~converged, status = ' (not converged)'; end
    fprintf('  I=%2.0f A: T_chip=%.1f K, T_max=%.1f K, Q_out=%.3f W%s\n', ...
            I_current, T(1), max(T), Q_out_calc, status);
end

%% Print Comparison Table
fprintf('\n');
fprintf('╔═══════════════════════════════════════════════════════════════════════════════╗\n');
fprintf('║                    TEMPERATURE COMPARISON TABLE                              ║\n');
fprintf('╠═══════════════════════════════════════════════════════════════════════════════╣\n');
fprintf('║ I (A) │ MATLAB T_chip │ COMSOL T_chip │ Error (K) │ Error (%%)  │ Q_out Err  ║\n');
fprintf('╠═══════════════════════════════════════════════════════════════════════════════╣\n');

for idx = 1:n_current
    I_val = comsol_data.I(idx);
    T_matlab = T_center_matlab(idx);
    T_comsol = comsol_data.T_chip_avg(idx);
    error_K = T_matlab - T_comsol;
    error_pct = 100 * error_K / T_comsol;
    Q_err_pct = 100 * (Q_out_matlab(idx) - comsol_data.Q_out(idx)) / comsol_data.Q_out(idx);
    
    fprintf('║ %5.0f │ %13.1f │ %13.1f │ %9.1f │ %10.1f │ %10.1f ║\n', ...
            I_val, T_matlab, T_comsol, error_K, error_pct, Q_err_pct);
end
fprintf('╚═══════════════════════════════════════════════════════════════════════════════╝\n');

%% Summary Statistics
fprintf('\n=== SUMMARY ===\n');

% Thermal resistance at I=0
T_water = 300;
Q_in_comsol = comsol_data.Q_in(1);
dT_comsol = comsol_data.T_chip_avg(1) - T_water;
dT_matlab = T_center_matlab(1) - T_water;
R_th_comsol = dT_comsol / Q_in_comsol;
R_th_matlab = dT_matlab / Q_in_comsol;

fprintf('\nThermal Resistance (at I=0):\n');
fprintf('  COMSOL: dT = %.1f K, R_th = %.1f K/W\n', dT_comsol, R_th_comsol);
fprintf('  MATLAB: dT = %.1f K, R_th = %.1f K/W\n', dT_matlab, R_th_matlab);
fprintf('  Ratio (MATLAB/COMSOL): %.2f×\n', R_th_matlab / R_th_comsol);

% Average errors
avg_T_error = mean(T_center_matlab - comsol_data.T_chip_avg);
avg_T_error_pct = 100 * mean((T_center_matlab - comsol_data.T_chip_avg) ./ comsol_data.T_chip_avg);
fprintf('\nAverage T_chip error: %.1f K (%.1f%%)\n', avg_T_error, avg_T_error_pct);

%% Plot Comparison
figure('Position', [100 100 1400 400]);

% T_chip vs Current
subplot(1, 3, 1);
plot(comsol_data.I, comsol_data.T_chip_avg, 'bo-', 'LineWidth', 2, 'MarkerSize', 8); hold on;
plot(comsol_data.I, T_center_matlab, 'rs--', 'LineWidth', 2, 'MarkerSize', 8);
xlabel('Current (A)', 'FontSize', 12);
ylabel('T_{chip} (K)', 'FontSize', 12);
title('Chip Surface Temperature', 'FontSize', 14);
legend('COMSOL FEM', 'MATLAB CTM', 'Location', 'best', 'FontSize', 11);
grid on;
xlim([0 10]);

% T_max vs Current
subplot(1, 3, 2);
plot(comsol_data.I, comsol_data.T_max, 'bo-', 'LineWidth', 2, 'MarkerSize', 8); hold on;
plot(comsol_data.I, T_max_matlab, 'rs--', 'LineWidth', 2, 'MarkerSize', 8);
xlabel('Current (A)', 'FontSize', 12);
ylabel('T_{max} (K)', 'FontSize', 12);
title('Maximum Temperature', 'FontSize', 14);
legend('COMSOL FEM', 'MATLAB CTM', 'Location', 'best', 'FontSize', 11);
grid on;
xlim([0 10]);

% Q_out vs Current
subplot(1, 3, 3);
plot(comsol_data.I, comsol_data.Q_out, 'bo-', 'LineWidth', 2, 'MarkerSize', 8); hold on;
plot(comsol_data.I, Q_out_matlab, 'rs--', 'LineWidth', 2, 'MarkerSize', 8);
xlabel('Current (A)', 'FontSize', 12);
ylabel('Q_{out} (W)', 'FontSize', 12);
title('Heat Output (Cold Side)', 'FontSize', 14);
legend('COMSOL FEM', 'MATLAB CTM', 'Location', 'best', 'FontSize', 11);
grid on;
xlim([0 10]);

sgtitle('MATLAB CTM vs COMSOL FEM Comparison (Single 30° Wedge)', 'FontSize', 16, 'FontWeight', 'bold');

% Save outputs
output_dir = fullfile(project_root, 'output', 'diagnosis');
if ~exist(output_dir, 'dir'), mkdir(output_dir); end
savefig(fullfile(output_dir, 'matlab_vs_comsol_comparison.fig'));
print(fullfile(output_dir, 'matlab_vs_comsol_comparison.png'), '-dpng', '-r150');

fprintf('\n✓ Plots saved to: %s\n', output_dir);
