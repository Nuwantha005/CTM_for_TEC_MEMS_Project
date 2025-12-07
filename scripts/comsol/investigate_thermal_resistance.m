%% INVESTIGATE THERMAL RESISTANCE DISCREPANCY
% This script investigates the 5x thermal resistance discrepancy between
% MATLAB CTM and COMSOL FEM models.
%
% Using t_TEC = 500 um to match COMSOL parametric study data
% Sweeps current from 0 to 10 A to compare with COMSOL results
%
% COMSOL Variables used:
%   - dom1: Global average temperature
%   - bnd1: Average temperature on heat flux chip surface  
%   - bnd2: Total heat flux in from chip surface
%   - bnd3: Total heat flux out from constant temperature boundary
%
% Usage:
%   1. Run: clear all
%   2. Start COMSOL server: comsolmphserver -port 2036
%   3. Run this script

clear; clc;

%% ============ CONFIGURATION ============
script_path = mfilename('fullpath');
if isempty(script_path)
    PROJECT_ROOT = pwd;
else
    [script_dir, ~, ~] = fileparts(script_path);
    PROJECT_ROOT = fileparts(fileparts(script_dir));
end

% Add paths
addpath(genpath(fullfile(PROJECT_ROOT, 'src')));
addpath('F:\EngineeringSoftware\COMSOL\COMSOL63\Multiphysics\mli');

COMSOL_PORT = 2036;
COMSOL_MODEL_PATH = 'E:\Semester 7\ME4311 - MicroNano Electro Mechanical Systems and Nanotechnology\Project\COMSOL\Test_Wedge\asym2.mph';

% Output to main project output folder
OUTPUT_DIR = fullfile(PROJECT_ROOT, 'output', 'thermal_resistance_study');
if ~exist(OUTPUT_DIR, 'dir')
    mkdir(OUTPUT_DIR);
end

%% ============ TEST PARAMETERS (Matching COMSOL Study) ============
% Fixed parameters
params.M = 12;                    % Number of wedges
params.theta_deg = 30;            % Wedge angle
params.t_TEC_um = 500;            % TEC thickness (matching COMSOL study)
params.t_chip_um = 50;            % Chip thickness
params.t_SOI_um = 20;             % SOI/insulator thickness
params.R_cyl_um = 1000;           % Cylinder radius
params.f_L = 1.15;                % Length ratio
params.fill_factor = 0.9;         % Fill factor
params.w_is_um = 50;              % Radial insulator width

% Interconnect parameters
params.ic_w_r = 0.1;
params.ic_t_r = 0.6;
params.ic_angle_r = 0.5;
params.oc_w_r = 0.1;
params.oc_t_r = 0.6;
params.oc_angle_r = 0.5;

% Heat flux (200 kW/m² to match COMSOL Q_in ~ 2.62 W)
params.q_Wm2 = 2e5;

% Current sweep (matching COMSOL parametric study)
I_sweep = 0:1:10;  % 0 to 10 A

% Calculate derived parameters
params.theta_rad = deg2rad(params.theta_deg);
w_chip = 10000e-6;  % 10 mm chip
R_cyl = params.R_cyl_um * 1e-6;
R_base = w_chip / sqrt(2);
N_stages = 3;

% Calculate L_1 using geometric series
w_is = params.w_is_um * 1e-6;
L_total = (R_base - R_cyl) - (N_stages + 1) * w_is;
f_L = params.f_L;
if abs(f_L - 1) < 1e-6
    L_1 = L_total / N_stages;
else
    L_1 = L_total * (1 - f_L) / (1 - f_L^N_stages);
end
params.L_1_um = L_1 * 1e6;

%% ============ COMSOL REFERENCE DATA ============
% From your parametric study (t_TEC = 500 um)
comsol_ref.I_A = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10];
comsol_ref.T_chip_avg_K = [524.55, 509.85, 504.63, 507.58, 517.75, 534.20, 556.44, 583.70, 615.56, 651.64, 691.57];
comsol_ref.T_max_K = [701.65, 657.18, 631.90, 623.74, 632.15, 653.98, 685.80, 725.61, 772.29, 824.87, 882.49];
comsol_ref.Q_in_W = [2.6179, 2.6180, 2.6179, 2.6179, 2.6180, 2.6179, 2.6181, 2.6184, 2.6188, 2.6191, 2.6174];
comsol_ref.Q_out_W = [2.6166, 2.5991, 2.6809, 2.8543, 3.1169, 3.4629, 3.8898, 4.3942, 4.9748, 5.6299, 6.3587];

fprintf('========================================================================\n');
fprintf('     THERMAL RESISTANCE INVESTIGATION\n');
fprintf('========================================================================\n');
fprintf('     Comparing MATLAB CTM vs COMSOL FEM\n');
fprintf('     t_TEC = %d um, q = %.0f kW/m²\n', params.t_TEC_um, params.q_Wm2/1000);
fprintf('========================================================================\n\n');

%% ============ RUN MATLAB MODEL FOR CURRENT SWEEP ============
fprintf('Running MATLAB CTM for current sweep...\n\n');

% Build config structure
config = struct();
config.geometry.N_stages = N_stages;
config.geometry.M_wedges = params.M;
config.geometry.wedge_angle_deg = params.theta_deg;
config.geometry.thickness_um = params.t_TEC_um;
config.geometry.radial_expansion_factor = params.f_L;
config.geometry.w_chip_um = 10000;
config.geometry.R_cyl_um = params.R_cyl_um;
config.geometry.t_chip_um = params.t_chip_um;
config.geometry.t_ins_um = params.t_SOI_um;
config.geometry.interconnect_ratio = params.ic_w_r;
config.geometry.outerconnect_ratio = params.oc_w_r;
config.geometry.insulation_width_um = params.w_is_um;
config.geometry.interconnect_angle_ratio = params.ic_angle_r;
config.geometry.outerconnect_angle_ratio = params.oc_angle_r;
config.geometry.interconnect_thickness_ratio = params.ic_t_r;
config.geometry.outerconnect_thickness_ratio = params.oc_t_r;
config.geometry.fill_factor = params.fill_factor;

config.boundary_conditions.q_flux_W_m2 = params.q_Wm2;
config.boundary_conditions.T_water_K = 300;
config.boundary_conditions.h_conv_W_m2K = 1e6;

config.materials.Bi2Te3 = struct('k', 1.2, 'rho', 1e-5, 'S', 0.0002);
config.materials.Cu = struct('k', 400, 'rho', 1.7e-8);
config.materials.Si = struct('k', 150, 'rho', 0.01);
config.materials.AlN = struct('k', 170, 'rho', 1e10);
config.materials.SiO2 = struct('k', 1.4, 'rho', 1e14);
config.materials.Al2O3 = struct('k', 30, 'rho', 1e12);

% Initialize materials and geometry once
materials = MaterialProperties(config);
geometry = TECGeometry(config);

% Print geometry info
fprintf('Geometry (single wedge):\n');
fprintf('  R_cyl = %.0f um\n', params.R_cyl_um);
fprintf('  R_base = %.0f um\n', R_base*1e6);
fprintf('  L_1 = %.2f um\n', params.L_1_um);
fprintf('  t_TEC = %d um\n', params.t_TEC_um);
fprintf('  theta = %.1f deg\n', params.theta_deg);
fprintf('  fill_factor = %.2f\n\n', params.fill_factor);

for i = 1:N_stages
    [r_in, L, ~, ~, ~, ~, ~, ~, w_az, ~] = geometry.get_stage_geometry(i);
    G = geometry.calculate_G(r_in, L, params.ic_w_r*L, params.ic_t_r*params.t_TEC_um*1e-6, ...
        params.ic_angle_r*params.theta_rad, params.oc_w_r*L, params.oc_t_r*params.t_TEC_um*1e-6, ...
        params.oc_angle_r*params.theta_rad, w_az, params.w_is_um*1e-6);
    K_legs = 2*1.2/G;
    fprintf('  Stage %d: r_in=%.0f um, L=%.0f um, w_az=%.1f um, R_legs=%.0f K/W\n', ...
        i, r_in*1e6, L*1e6, w_az*1e6, 1/K_legs);
end
fprintf('\n');

% Storage for results
matlab_results.I_A = I_sweep;
matlab_results.T_center_K = zeros(size(I_sweep));
matlab_results.T_max_K = zeros(size(I_sweep));
matlab_results.Q_out_W = zeros(size(I_sweep));
matlab_results.Q_in_W = zeros(size(I_sweep));

% Run for each current value
for idx = 1:length(I_sweep)
    I_current = I_sweep(idx);
    config.operating_conditions.I_current_A = I_current;
    
    % Rebuild network with new current
    network = ThermalNetwork(geometry, materials, config);
    
    % Solve iteratively
    dim = 2*N_stages + 1;
    T = ones(dim, 1) * 300;
    
    for iter = 1:100
        T_old = T;
        [T, Q_out, Q_in] = network.solve(T);
        if max(abs(T - T_old)) < 1e-6
            break;
        end
    end
    
    matlab_results.T_center_K(idx) = T(1);
    matlab_results.T_max_K(idx) = max(T);
    matlab_results.Q_out_W(idx) = Q_out;
    matlab_results.Q_in_W(idx) = Q_in;
    
    fprintf('  I = %2.0f A: T_center = %7.2f K, T_max = %7.2f K, Q_out = %.4f W\n', ...
        I_current, T(1), max(T), Q_out);
end

%% ============ COMPARISON ANALYSIS ============
fprintf('\n');
fprintf('========================================================================\n');
fprintf('                    COMPARISON: MATLAB vs COMSOL\n');
fprintf('========================================================================\n');
fprintf('  I(A)   MATLAB_Tavg    COMSOL_Tavg    Diff(K)    Diff(%%)\n');
fprintf('------------------------------------------------------------------------\n');

for idx = 1:length(comsol_ref.I_A)
    I = comsol_ref.I_A(idx);
    T_matlab = matlab_results.T_center_K(idx);
    T_comsol = comsol_ref.T_chip_avg_K(idx);
    diff_K = T_matlab - T_comsol;
    diff_pct = 100 * diff_K / T_comsol;
    fprintf('  %2.0f     %7.2f K      %7.2f K    %+7.1f    %+6.1f%%\n', ...
        I, T_matlab, T_comsol, diff_K, diff_pct);
end

fprintf('------------------------------------------------------------------------\n\n');

fprintf('T_max Comparison:\n');
fprintf('  I(A)   MATLAB_Tmax    COMSOL_Tmax    Diff(K)    Diff(%%)\n');
fprintf('------------------------------------------------------------------------\n');

for idx = 1:length(comsol_ref.I_A)
    I = comsol_ref.I_A(idx);
    T_matlab = matlab_results.T_max_K(idx);
    T_comsol = comsol_ref.T_max_K(idx);
    diff_K = T_matlab - T_comsol;
    diff_pct = 100 * diff_K / T_comsol;
    fprintf('  %2.0f     %7.2f K      %7.2f K    %+7.1f    %+6.1f%%\n', ...
        I, T_matlab, T_comsol, diff_K, diff_pct);
end

fprintf('------------------------------------------------------------------------\n\n');

fprintf('Q_out Comparison:\n');
fprintf('  I(A)   MATLAB_Qout    COMSOL_Qout    Diff(W)    Diff(%%)\n');
fprintf('------------------------------------------------------------------------\n');

for idx = 1:length(comsol_ref.I_A)
    I = comsol_ref.I_A(idx);
    Q_matlab = matlab_results.Q_out_W(idx);
    Q_comsol = comsol_ref.Q_out_W(idx);
    diff_W = Q_matlab - Q_comsol;
    diff_pct = 100 * diff_W / Q_comsol;
    fprintf('  %2.0f     %7.4f W      %7.4f W    %+7.4f    %+6.1f%%\n', ...
        I, Q_matlab, Q_comsol, diff_W, diff_pct);
end

fprintf('------------------------------------------------------------------------\n\n');

%% ============ THERMAL RESISTANCE ANALYSIS ============
fprintf('========================================================================\n');
fprintf('                    THERMAL RESISTANCE ANALYSIS\n');
fprintf('========================================================================\n');

% At I = 0 A (no Peltier effect), thermal resistance is simply dT/Q
T_water = 300;  % K

idx_I0 = 1;  % I = 0 A
dT_matlab_I0 = matlab_results.T_center_K(idx_I0) - T_water;
dT_comsol_I0 = comsol_ref.T_chip_avg_K(idx_I0) - T_water;
Q_in_matlab = matlab_results.Q_in_W(idx_I0);
Q_in_comsol = comsol_ref.Q_in_W(idx_I0);

R_th_matlab = dT_matlab_I0 / Q_in_matlab;
R_th_comsol = dT_comsol_I0 / Q_in_comsol;

fprintf('\nAt I = 0 A (pure conduction, no Peltier):\n');
fprintf('  MATLAB: dT = %.2f K, Q_in = %.4f W, R_th = %.1f K/W\n', ...
    dT_matlab_I0, Q_in_matlab, R_th_matlab);
fprintf('  COMSOL: dT = %.2f K, Q_in = %.4f W, R_th = %.1f K/W\n', ...
    dT_comsol_I0, Q_in_comsol, R_th_comsol);
fprintf('  Ratio R_th_matlab / R_th_comsol = %.2f\n\n', R_th_matlab / R_th_comsol);

% Check heat balance
fprintf('Heat Balance at I = 0 A:\n');
fprintf('  MATLAB: Q_in = %.4f W, Q_out = %.4f W, balance = %.4f W\n', ...
    Q_in_matlab, matlab_results.Q_out_W(idx_I0), Q_in_matlab - matlab_results.Q_out_W(idx_I0));
fprintf('  COMSOL: Q_in = %.4f W, Q_out = %.4f W, balance = %.4f W\n', ...
    Q_in_comsol, comsol_ref.Q_out_W(idx_I0), Q_in_comsol - comsol_ref.Q_out_W(idx_I0));

%% ============ DETAILED BREAKDOWN ============
fprintf('\n========================================================================\n');
fprintf('                    DETAILED THERMAL BREAKDOWN\n');
fprintf('========================================================================\n');

% Recalculate thermal resistances in MATLAB model
fprintf('\nMATLAB CTM thermal resistances (t_TEC = %d um):\n', params.t_TEC_um);

theta = params.theta_rad;
t_TEC = params.t_TEC_um * 1e-6;
t_chip = params.t_chip_um * 1e-6;
t_SOI = params.t_SOI_um * 1e-6;

k_BiTe = 1.2;
k_Si = 150;
k_AlN = 170;
k_SiO2 = 1.4;

R_total_TEC = 0;
for i = 1:N_stages
    [r_in, L, w_ic, t_ic, beta_ic, w_oc, t_oc, beta_oc, w_az, w_is_stage] = geometry.get_stage_geometry(i);
    G = geometry.calculate_G(r_in, L, w_ic, t_ic, beta_ic, w_oc, t_oc, beta_oc, w_az, w_is_stage);
    K_legs = 2*k_BiTe/G;
    R_legs = 1/K_legs;
    
    r_out = r_in + L;
    R_is = (1/(k_AlN * t_TEC * theta)) * log((r_out + w_is_stage)/r_out);
    
    K_az = k_SiO2 * (w_az * t_TEC) / L;
    
    R_stage = 1/(1/(R_legs + R_is) + K_az);
    R_total_TEC = R_total_TEC + R_stage;
    
    fprintf('  Stage %d:\n', i);
    fprintf('    R_TE_legs = %.1f K/W (G = %.2e)\n', R_legs, G);
    fprintf('    R_radial_insulator = %.2f K/W\n', R_is);
    fprintf('    K_azimuthal = %.2e W/K (parallel path)\n', K_az);
    fprintf('    R_stage_effective = %.1f K/W\n', R_stage);
end

fprintf('\n  Total TEC resistance (3 stages series) ~ %.1f K/W\n', R_total_TEC);
fprintf('  (This is approximate - actual network is more complex)\n');

%% ============ SAVE RESULTS ============
fprintf('\n========================================================================\n');
fprintf('Saving results...\n');

results = struct();
results.params = params;
results.matlab = matlab_results;
results.comsol_ref = comsol_ref;
results.analysis.R_th_matlab = R_th_matlab;
results.analysis.R_th_comsol = R_th_comsol;
results.analysis.R_th_ratio = R_th_matlab / R_th_comsol;
results.timestamp = datestr(now);

filename = sprintf('thermal_resistance_study_%s.mat', datestr(now, 'yyyy-mm-dd_HH-MM-SS'));
save(fullfile(OUTPUT_DIR, filename), 'results');
fprintf('  Saved: %s\n', fullfile(OUTPUT_DIR, filename));

% Save comparison CSV
csv_file = fullfile(OUTPUT_DIR, 'comparison_results.csv');
fid = fopen(csv_file, 'w');
fprintf(fid, 'I_A,MATLAB_Tcenter_K,COMSOL_Tavg_K,MATLAB_Tmax_K,COMSOL_Tmax_K,MATLAB_Qout_W,COMSOL_Qout_W\n');
for idx = 1:length(I_sweep)
    fprintf(fid, '%.1f,%.2f,%.2f,%.2f,%.2f,%.4f,%.4f\n', ...
        I_sweep(idx), matlab_results.T_center_K(idx), comsol_ref.T_chip_avg_K(idx), ...
        matlab_results.T_max_K(idx), comsol_ref.T_max_K(idx), ...
        matlab_results.Q_out_W(idx), comsol_ref.Q_out_W(idx));
end
fclose(fid);
fprintf('  Saved: %s\n', csv_file);

fprintf('\n========================================================================\n');
fprintf('                         INVESTIGATION COMPLETE\n');
fprintf('========================================================================\n');
