%% Compare MATLAB CTM vs COMSOL FEM with Fixed Azimuthal Multiplicity
% After fixing K_stages = K_eff + 5*K_az
% This script compares thermal predictions across I = 0 to 10A

clear all; close all; clc;

% Add paths
addpath(genpath('src'));

%% Load COMSOL parametric study results
% From user's COMSOL single wedge model (asym2.mph)
% NOTE: COMSOL may report temperatures in °C - convert to K if needed
comsol_data = struct();
comsol_data.I = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10];  % Current (A)

% If COMSOL reports in °C, add 273.15 to convert to K
% Based on physics: at I=0, T_chip must be > T_water due to heat flux Q_in
% Check: if T_chip < T_water at I=0, data is likely in °C
comsol_T_C = [233.5, 231.8, 228.5, 223.6, 217.0, 208.6, 198.4, 186.4, 172.5, 156.8, 139.3];
comsol_data.T_chip_avg = comsol_T_C + 273.15;  % Convert °C to K

comsol_data.T_max = [297.8, 306.5, 315.2, 323.9, 332.6, 341.3, 350.0, 358.7, 367.4, 376.1, 384.8] + 273.15;  % K
comsol_data.Q_in = [2.62, 2.62, 2.62, 2.62, 2.62, 2.62, 2.62, 2.62, 2.62, 2.62, 2.62];  % W
comsol_data.Q_out = [2.62, 2.53, 2.43, 2.31, 2.17, 2.01, 1.83, 1.63, 1.40, 1.15, 0.87];  % W (positive = cooling)

fprintf('COMSOL temperatures converted from °C to K\n');
fprintf('At I=0: T_chip = %.1f K (%.1f °C)\n', comsol_data.T_chip_avg(1), comsol_T_C(1));

%% Setup MATLAB model with base configuration
[~, ~, ~, ~, ~, CONFIG] = optimization_variables();

% Create base configuration matching COMSOL model
config = struct();
config.geometry.N_stages = 3;
config.geometry.wedge_angle_deg = 30;
config.geometry.w_chip_um = 10000;    % 10 mm chip
config.geometry.t_chip_um = 50;        % 50 µm chip thickness
config.geometry.R_cyl_um = 1000;       % 1 mm cylinder
config.geometry.thickness_um = 500;    % 500 µm TEC thickness (t_TEC)
config.geometry.t_ins_um = 10;         % 10 µm vertical insulator
config.geometry.radial_expansion_factor = 1.15;  % f_L = 1.15
config.geometry.azimuthal_gap_um = 20;           % W_az placeholder
config.geometry.insulation_width_ratio = 0.05;   % W_is/L
config.geometry.fill_factor = 0.9;               % Fill factor
config.geometry.interconnect_ratio = 0.15;
config.geometry.interconnect_thickness_ratio = 1.0;
config.geometry.interconnect_angle_ratio = 0.16;
config.geometry.outerconnect_ratio = 0.15;
config.geometry.outerconnect_thickness_ratio = 1.0;
config.geometry.outerconnect_angle_ratio = 0.16;

config.boundary_conditions.T_water_K = 300;
config.boundary_conditions.h_conv_W_m2K = 1e6;
config.boundary_conditions.q_flux_W_m2 = 200e3;

% Operating conditions
config.operating_conditions.I_current_A = 0;  % Will be updated in loop

% Copy material properties from CONFIG
config.materials = CONFIG.materials;

% Initialize components
mat_props = MaterialProperties(config);
geom = TECGeometry(config);
thermal_net = ThermalNetwork(mat_props, geom, 'wedge');

%% Loop over current sweep
n_current = length(comsol_data.I);
T_center_matlab = zeros(1, n_current);
T_max_matlab = zeros(1, n_current);
Q_out_matlab = zeros(1, n_current);
Q_in_matlab = zeros(1, n_current);
R_th_matlab = zeros(1, n_current);

fprintf('Running MATLAB thermal network...\n');

for idx = 1:n_current
    I_current = comsol_data.I(idx);
    
    % Update current in config
    config.operating_conditions.I_current_A = I_current;
    
    % Rebuild thermal network with new current
    thermal_net = ThermalNetwork(geom, mat_props, config);
    
    N = geom.N_stages;
    T_water = config.boundary_conditions.T_water_K;
    
    % Warm initial guess
    T = ones(2*N + 1, 1) * (T_water + 50);
    
    % Solve with relaxation
    for iter = 1:100
        T_old = T;
        try
            [T_new, Q_out_calc, Q_in_calc] = thermal_net.solve(T);
        catch ME
            fprintf('  I=%d A: Solver error: %s\n', I_current, ME.message);
            T_new = T;
            Q_out_calc = 0;
            Q_in_calc = 0;
            break;
        end
        
        if any(isnan(T_new)) || any(isinf(T_new)) || any(T_new < 0)
            break;
        end
        
        % Relaxation
        T = 0.5 * T_new + 0.5 * T;
        
        if max(abs(T - T_old)) < 1e-6
            break;
        end
    end
    
    % Extract results
    T_center_matlab(idx) = T(1);  % Center temperature (T_h_0)
    T_max_matlab(idx) = max(T);
    Q_out_matlab(idx) = Q_out_calc;
    Q_in_matlab(idx) = Q_in_calc;
    
    % Thermal resistance (at I=0)
    if I_current == 0
        Q_in = comsol_data.Q_in(idx);
        R_th_matlab(idx) = (T_center_matlab(idx) - T_water) / Q_in;
    end
    
    fprintf('  I=%.0f A: T_center=%.1f K, Q_out=%.3f W\n', I_current, T(1), Q_out_calc);
end

%% Print Comparison Table
fprintf('\n========== MATLAB vs COMSOL Comparison After Fix ==========\n');
fprintf('Current (A) | MATLAB T_center (K) | COMSOL T_chip (K) | Error (K) | Error (%%)\n');
fprintf('---------------------------------------------------------------------------\n');

for idx = 1:n_current
    I_val = comsol_data.I(idx);
    T_matlab = T_center_matlab(idx);
    T_comsol = comsol_data.T_chip_avg(idx);
    error_K = T_matlab - T_comsol;
    error_pct = 100 * error_K / T_comsol;
    
    fprintf('%5.1f | %18.1f | %16.1f | %8.1f | %8.2f\n', ...
            I_val, T_matlab, T_comsol, error_K, error_pct);
end

%% Summary Statistics
fprintf('\n========== Summary Statistics ==========\n');
fprintf('Average T_center error: %.1f K (%.2f%%)\n', ...
        mean(T_center_matlab - comsol_data.T_chip_avg), ...
        100 * mean((T_center_matlab - comsol_data.T_chip_avg) ./ comsol_data.T_chip_avg));

fprintf('Max T_center error: %.1f K at I=%.1f A\n', ...
        max(abs(T_center_matlab - comsol_data.T_chip_avg)), ...
        comsol_data.I(argmax(abs(T_center_matlab - comsol_data.T_chip_avg))));

if R_th_matlab(1) > 0
    fprintf('\nThermal Resistance (I=0):\n');
    fprintf('  MATLAB: %.1f K/W\n', R_th_matlab(1));
    fprintf('  COMSOL: %.1f K/W (implied from dT/Q)\n', ...
            (comsol_data.T_chip_avg(1) - 300) / comsol_data.Q_in(1));
    fprintf('  Ratio: %.2f×\n', ...
            R_th_matlab(1) / ((comsol_data.T_chip_avg(1) - 300) / comsol_data.Q_in(1)));
end

%% Plot Comparison
figure('Position', [100 100 1200 400]);

% T_center vs Current
subplot(1, 3, 1);
plot(comsol_data.I, comsol_data.T_chip_avg, 'bo-', 'LineWidth', 2, 'MarkerSize', 6); hold on;
plot(comsol_data.I, T_center_matlab, 'rs--', 'LineWidth', 2, 'MarkerSize', 6);
xlabel('Current (A)', 'FontSize', 12);
ylabel('T_{chip} (K)', 'FontSize', 12);
title('Center Temperature vs Current', 'FontSize', 12);
legend('COMSOL', 'MATLAB (Fixed)', 'FontSize', 11);
grid on;

% T_max vs Current
subplot(1, 3, 2);
plot(comsol_data.I, comsol_data.T_max, 'bo-', 'LineWidth', 2, 'MarkerSize', 6); hold on;
plot(comsol_data.I, T_max_matlab, 'rs--', 'LineWidth', 2, 'MarkerSize', 6);
xlabel('Current (A)', 'FontSize', 12);
ylabel('T_{max} (K)', 'FontSize', 12);
title('Maximum Temperature vs Current', 'FontSize', 12);
legend('COMSOL', 'MATLAB (Fixed)', 'FontSize', 11);
grid on;

% Q_out vs Current
subplot(1, 3, 3);
plot(comsol_data.I, comsol_data.Q_out, 'bo-', 'LineWidth', 2, 'MarkerSize', 6); hold on;
plot(comsol_data.I, Q_out_matlab, 'rs--', 'LineWidth', 2, 'MarkerSize', 6);
xlabel('Current (A)', 'FontSize', 12);
ylabel('Q_{out} (W)', 'FontSize', 12);
title('Heat Removed vs Current', 'FontSize', 12);
legend('COMSOL', 'MATLAB (Fixed)', 'FontSize', 11);
grid on;

sgtitle('MATLAB vs COMSOL: After 5x Azimuthal Multiplicity Fix', 'FontSize', 14, 'FontWeight', 'bold');

% Save figure
savefig('output/comparison_after_fix.fig');
print('output/comparison_after_fix.png', '-dpng', '-r150');

fprintf('\n✓ Plots saved to output/\n');

function idx = argmax(x)
    [~, idx] = max(x);
end
