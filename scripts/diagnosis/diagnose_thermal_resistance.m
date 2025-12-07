%% Diagnostic Script: Analyze Thermal Resistance Breakdown
% This script computes the individual resistances to understand why 
% MATLAB predicts 15× higher thermal resistance than COMSOL

clear all; close all; clc;
addpath(genpath('src'));

fprintf('╔═══════════════════════════════════════════════════════════════╗\n');
fprintf('║     THERMAL RESISTANCE DIAGNOSTIC                            ║\n');
fprintf('╚═══════════════════════════════════════════════════════════════╝\n\n');

%% Load configuration
[~, ~, ~, ~, ~, CONFIG] = optimization_variables();

config = struct();
config.geometry.N_stages = 3;
config.geometry.wedge_angle_deg = 30;
config.geometry.w_chip_um = 10000;
config.geometry.t_chip_um = 50;
config.geometry.R_cyl_um = 1000;
config.geometry.thickness_um = 500;  % t_TEC = 500 µm
config.geometry.t_ins_um = 10;
config.geometry.radial_expansion_factor = 1.15;
config.geometry.azimuthal_gap_um = 20;
config.geometry.insulation_width_ratio = 0.05;
config.geometry.fill_factor = 0.9;
config.geometry.interconnect_ratio = 0.15;
config.geometry.interconnect_thickness_ratio = 1.0;
config.geometry.interconnect_angle_ratio = 0.16;
config.geometry.outerconnect_ratio = 0.15;
config.geometry.outerconnect_thickness_ratio = 1.0;
config.geometry.outerconnect_angle_ratio = 0.16;
config.boundary_conditions.T_water_K = 300;
config.boundary_conditions.h_conv_W_m2K = 1e6;
config.boundary_conditions.q_flux_W_m2 = 200e3;
config.operating_conditions.I_current_A = 0;
config.materials = CONFIG.materials;

mat = MaterialProperties(config);
geom = TECGeometry(config);

%% Print geometry parameters
fprintf('=== GEOMETRY PARAMETERS ===\n');
fprintf('N_stages:        %d\n', geom.N_stages);
fprintf('Wedge angle:     %.4f rad (%.1f°)\n', geom.WedgeAngle, geom.WedgeAngle*180/pi);
fprintf('R_base:          %.6f m (%.1f µm)\n', geom.R_base, geom.R_base*1e6);
fprintf('R_cyl:           %.6f m (%.1f µm)\n', geom.R_cyl, geom.R_cyl*1e6);
fprintf('t_TEC:           %.6f m (%.1f µm)\n', geom.Thickness, geom.Thickness*1e6);
fprintf('L_1:             %.6f m (%.1f µm)\n', geom.L_1, geom.L_1*1e6);
fprintf('\n');

%% Material properties at 350 K (typical operating point)
T_ref = 350;  % K
k_Bi2Te3 = mat.get_k('Bi2Te3', T_ref);
k_AlN = mat.get_k('AlN', T_ref);
k_SiO2 = mat.get_k('SiO2', T_ref);
k_Si = mat.get_k('Si', T_ref);

fprintf('=== MATERIAL PROPERTIES (at %.0f K) ===\n', T_ref);
fprintf('k_Bi2Te3:        %.2f W/mK\n', k_Bi2Te3);
fprintf('k_AlN:           %.1f W/mK\n', k_AlN);
fprintf('k_SiO2:          %.2f W/mK\n', k_SiO2);
fprintf('k_Si:            %.1f W/mK\n', k_Si);
fprintf('\n');

%% Calculate resistances for each stage
fprintf('=== STAGE-BY-STAGE THERMAL RESISTANCE ===\n');
fprintf('---------------------------------------------------------------------------\n');
fprintf('Stage | r_in (µm) | L (µm)  | G (1/m)    | K_legs | R_is   | K_az   | K_total\n');
fprintf('---------------------------------------------------------------------------\n');

theta = geom.WedgeAngle;
t_tec = geom.Thickness;

K_stages = zeros(geom.N_stages, 1);
R_stages = zeros(geom.N_stages, 1);

for i = 1:geom.N_stages
    [r_in, L, w_ic, t_ic, beta_ic, w_oc, t_oc, beta_oc, w_az, w_is] = geom.get_stage_geometry(i);
    r_out = r_in + L;
    
    % Calculate G factor
    G = geom.calculate_G(r_in, L, w_ic, t_ic, beta_ic, w_oc, t_oc, beta_oc, w_az, w_is);
    
    % Thermal conductance of TE legs
    k_p = k_Bi2Te3;
    k_n = k_Bi2Te3;
    K_legs = (k_p + k_n) / G;
    
    % Insulator resistance
    R_is = geom.calculate_R_thermal_insulator(r_out, w_is, t_tec, theta, k_AlN);
    
    % Azimuthal conductance
    K_az = geom.calculate_K_azimuthal(r_in, L, w_az, t_tec, k_SiO2, theta);
    
    % Combined stage conductance (with 5x azimuthal)
    R_eff = R_is + 1/K_legs;
    K_eff = 1/R_eff;
    K_stage = K_eff + 5*K_az;
    
    K_stages(i) = K_stage;
    R_stages(i) = 1/K_stage;
    
    fprintf('%5d | %9.2f | %7.2f | %10.2e | %6.4f | %6.2f | %6.5f | %7.5f\n', ...
            i, r_in*1e6, L*1e6, G, K_legs, R_is, K_az, K_stage);
end

fprintf('---------------------------------------------------------------------------\n');

%% Total thermal resistance through TEC stages (in series)
R_total_TEC = sum(R_stages);
K_total_TEC = 1/R_total_TEC;

fprintf('\n=== SUMMARY ===\n');
fprintf('Individual stage resistances: R_1 = %.1f K/W, R_2 = %.1f K/W, R_3 = %.1f K/W\n', ...
        R_stages(1), R_stages(2), R_stages(3));
fprintf('Total TEC resistance (3 stages in series): %.1f K/W\n', R_total_TEC);
fprintf('Total TEC conductance: %.5f W/K\n', K_total_TEC);

%% Compare to simple estimate
% Simple conduction: R = L / (k*A)
L_total = 0;
A_avg = 0;
for i = 1:geom.N_stages
    [r_in, L, ~, ~, ~, ~, ~, ~, ~, ~] = geom.get_stage_geometry(i);
    r_out = r_in + L;
    A_i = (theta/2) * (r_out^2 - r_in^2) * t_tec;  % Area
    L_total = L_total + L;
    A_avg = A_avg + A_i;
end
A_avg = A_avg / geom.N_stages;

R_simple = L_total / (2 * k_Bi2Te3 * A_avg);  % 2 legs in parallel

fprintf('\n=== SIMPLE CONDUCTION ESTIMATE ===\n');
fprintf('Total path length: %.1f µm\n', L_total*1e6);
fprintf('Average cross-sectional area: %.3e m²\n', A_avg);
fprintf('Simple R estimate (2 legs parallel): %.1f K/W\n', R_simple);

%% Expected from COMSOL
Q_in_comsol = 2.62;  % W
T_chip_comsol = 506.6;  % K (233.5°C converted to K)
T_water = 300;
dT_comsol = T_chip_comsol - T_water;
R_apparent_comsol = dT_comsol / Q_in_comsol;

fprintf('\n=== COMSOL REFERENCE (I=0) ===\n');
fprintf('Q_in = %.2f W\n', Q_in_comsol);
fprintf('T_chip = %.1f K (%.1f °C)\n', T_chip_comsol, T_chip_comsol - 273.15);
fprintf('T_water = %.1f K\n', T_water);
fprintf('dT = %.1f K\n', dT_comsol);
fprintf('Apparent R_th = dT/Q = %.1f K/W\n', R_apparent_comsol);

fprintf('\n=== DISCREPANCY ===\n');
fprintf('MATLAB R_total = %.1f K/W\n', R_total_TEC);
fprintf('COMSOL R_apparent = %.1f K/W\n', R_apparent_comsol);
fprintf('Ratio = %.1f×\n', R_total_TEC / R_apparent_comsol);

%% Investigate G factor in detail
fprintf('\n=== G FACTOR DETAILED BREAKDOWN ===\n');
[r_in, L, w_ic, t_ic, beta_ic, w_oc, t_oc, beta_oc, w_az, w_is] = geom.get_stage_geometry(1);
r_out = r_in + L;

fprintf('Stage 1 geometry:\n');
fprintf('  r_in = %.2f µm, r_out = %.2f µm, L = %.2f µm\n', r_in*1e6, r_out*1e6, L*1e6);
fprintf('  w_ic = %.2f µm, t_ic = %.2f µm, beta_ic = %.4f rad\n', w_ic*1e6, t_ic*1e6, beta_ic);
fprintf('  w_oc = %.2f µm, t_oc = %.2f µm, beta_oc = %.4f rad\n', w_oc*1e6, t_oc*1e6, beta_oc);
fprintf('  w_az = %.2f µm\n', w_az*1e6);
fprintf('  w_is = %.2f µm\n', w_is*1e6);

% Effective cross-sectional area
% For a wedge: A_eff = (theta/2) * (r_out^2 - r_in^2) * t_tec (without interconnect losses)
A_wedge_full = (theta/2) * (r_out^2 - r_in^2) * t_tec;
fprintf('\n  Full wedge area: %.3e m²\n', A_wedge_full);

% Simple G = L / A for comparison
G_simple = L / A_wedge_full;
G_actual = geom.calculate_G(r_in, L, w_ic, t_ic, beta_ic, w_oc, t_oc, beta_oc, w_az, w_is);
fprintf('  Simple G (L/A): %.3e 1/m\n', G_simple);
fprintf('  Actual G (with interconnects): %.3e 1/m\n', G_actual);
fprintf('  Ratio (actual/simple): %.2f\n', G_actual / G_simple);
