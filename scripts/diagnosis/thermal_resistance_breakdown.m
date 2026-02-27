%% Thermal Resistance Breakdown Analysis
% Analyze individual resistance components to find the source of 14x discrepancy
% between MATLAB CTM and COMSOL FEM

clear all; close all; clc;

script_dir = fileparts(mfilename('fullpath'));
project_root = fullfile(script_dir, '..', '..');
addpath(genpath(fullfile(project_root, 'src')));

fprintf('╔═══════════════════════════════════════════════════════════════╗\n');
fprintf('║     THERMAL RESISTANCE BREAKDOWN ANALYSIS                    ║\n');
fprintf('╚═══════════════════════════════════════════════════════════════╝\n\n');

%% Load configuration
[~, ~, ~, ~, ~, CONFIG] = optimization_variables();

config = struct();
config.geometry.N_stages = 3;
config.geometry.wedge_angle_deg = 30;
config.geometry.w_chip_um = 10000;
config.geometry.t_chip_um = 50;
config.geometry.R_cyl_um = 1000;
config.geometry.thickness_um = 1000;  % t_TEC = 1000 µm (from COMSOL_parameters.txt)
config.geometry.t_ins_um = 10;
config.geometry.radial_expansion_factor = 1.15;
config.geometry.azimuthal_gap_um = 20;
config.geometry.insulation_width_ratio = 0.05;
config.geometry.fill_factor = 0.9;
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

mat = MaterialProperties(config);
geom = TECGeometry(config);

%% Print geometry
fprintf('=== GEOMETRY ===\n');
fprintf('N_stages = %d\n', geom.N_stages);
fprintf('θ = %.4f rad (%.1f°)\n', geom.WedgeAngle, geom.WedgeAngle*180/pi);
fprintf('R_base = %.4f m (%.2f mm) - chip outer edge (inscribed circle)\n', geom.R_base, geom.R_base*1e3);
fprintf('R_cyl = %.4f m (%.2f mm) - inner cylinder\n', geom.R_cyl, geom.R_cyl*1e3);
fprintf('t_TEC = %.0f µm\n', geom.Thickness*1e6);
fprintf('t_chip = %.0f µm\n', config.geometry.t_chip_um);
fprintf('L_1 = %.1f µm (first stage length)\n', geom.L_1*1e6);
fprintf('\n');

%% Material properties at reference temperature
T_ref = 400;  % K - approximate operating point
k_TE = mat.get_k('Bi2Te3', T_ref);
k_AlN = mat.get_k('AlN', T_ref);
k_SiO2 = mat.get_k('SiO2', T_ref);
k_Si = mat.get_k('Si', T_ref);

fprintf('=== MATERIAL PROPERTIES (at %.0f K) ===\n', T_ref);
fprintf('k_Bi2Te3 = %.2f W/mK\n', k_TE);
fprintf('k_AlN = %.1f W/mK (radial insulator)\n', k_AlN);
fprintf('k_SiO2 = %.2f W/mK (azimuthal insulator)\n', k_SiO2);
fprintf('k_Si = %.1f W/mK (chip)\n', k_Si);
fprintf('\n');

%% Calculate resistances for each stage
fprintf('=== STAGE-BY-STAGE ANALYSIS ===\n\n');

theta = geom.WedgeAngle;
t_tec = geom.Thickness;
t_chip = config.geometry.t_chip_um * 1e-6;

total_R_TE = 0;
total_R_is = 0;

for i = 1:geom.N_stages
    [r_in, L, w_ic, t_ic, beta_ic, w_oc, t_oc, beta_oc, w_az, w_is] = geom.get_stage_geometry(i);
    r_out = r_in + L;
    
    fprintf('--- Stage %d ---\n', i);
    fprintf('r_in = %.2f µm, r_out = %.2f µm, L = %.2f µm\n', r_in*1e6, r_out*1e6, L*1e6);
    
    % G factor (current MATLAB implementation)
    G = geom.calculate_G(r_in, L, w_ic, t_ic, beta_ic, w_oc, t_oc, beta_oc, w_az, w_is);
    
    % Thermal resistance of TE legs: R = G / (k_p + k_n)
    R_TE = G / (2 * k_TE);  % Two legs (P and N) with same k
    K_TE = 1 / R_TE;
    
    % Simple radial conduction estimate for comparison
    % R = (1/kθt) * ln(r_out/r_in) for radial flow through annular sector
    R_TE_simple = (1 / (k_TE * theta * t_tec)) * log(r_out / r_in);
    
    % For HALF wedge (one leg), then 2 legs in parallel
    R_TE_half = (1 / (k_TE * (theta/2) * t_tec)) * log(r_out / r_in);
    R_TE_parallel = R_TE_half / 2;  % Two legs in parallel
    
    % Insulator resistance
    R_is = geom.calculate_R_thermal_insulator(r_out, w_is, t_tec, theta, k_AlN);
    
    % Azimuthal conductance
    K_az = geom.calculate_K_azimuthal(r_in, L, w_az, t_tec, k_SiO2, theta);
    R_az = 1 / K_az;
    
    fprintf('  G factor = %.2e 1/m\n', G);
    fprintf('  R_TE (MATLAB, from G) = %.1f K/W\n', R_TE);
    fprintf('  R_TE (simple radial, full θ) = %.1f K/W\n', R_TE_simple);
    fprintf('  R_TE (half θ, 2 parallel) = %.1f K/W\n', R_TE_parallel);
    fprintf('  R_insulator = %.2f K/W\n', R_is);
    fprintf('  R_azimuthal = %.1f K/W (K_az = %.2e W/K)\n', R_az, K_az);
    
    % Combined stage resistance (current implementation with 5x azimuthal)
    R_eff = R_is + R_TE;
    K_eff = 1 / R_eff;
    K_stage = K_eff + 5 * K_az;  % With 5x multiplier
    R_stage = 1 / K_stage;
    
    fprintf('  R_stage (R_is + R_TE) = %.1f K/W\n', R_eff);
    fprintf('  R_stage (with 5x K_az parallel) = %.1f K/W\n', R_stage);
    fprintf('\n');
    
    total_R_TE = total_R_TE + R_TE;
    total_R_is = total_R_is + R_is;
end

%% Total resistance (stages in series)
fprintf('=== TOTAL THERMAL RESISTANCE (3 stages in series) ===\n');
fprintf('Sum of R_TE (all stages): %.1f K/W\n', total_R_TE);
fprintf('Sum of R_insulator (all stages): %.2f K/W\n', total_R_is);
fprintf('\n');

%% COMSOL reference
fprintf('=== COMSOL REFERENCE (I=0) ===\n');
Q_in = 2.6179;  % W
T_chip_comsol = 524.55;  % K
T_water = 300;  % K
dT = T_chip_comsol - T_water;
R_apparent = dT / Q_in;

fprintf('Q_in = %.4f W\n', Q_in);
fprintf('T_chip = %.2f K\n', T_chip_comsol);
fprintf('T_water = %.0f K\n', T_water);
fprintf('dT = %.1f K\n', dT);
fprintf('R_apparent = dT/Q_in = %.1f K/W\n', R_apparent);
fprintf('\n');

%% What resistance SHOULD give COMSOL's result?
fprintf('=== REQUIRED RESISTANCE TO MATCH COMSOL ===\n');
fprintf('To get T_chip = %.1f K with Q_in = %.3f W:\n', T_chip_comsol, Q_in);
fprintf('  Required R_th = %.1f K/W\n', R_apparent);
fprintf('\n');
fprintf('MATLAB predicts R_th ≈ %.0f K/W (from sum of stage resistances)\n', total_R_TE + total_R_is);
fprintf('Ratio = %.1fx\n', (total_R_TE + total_R_is) / R_apparent);
fprintf('\n');

%% Possible explanations
fprintf('=== INVESTIGATION: WHAT COULD CAUSE 14x LOWER RESISTANCE? ===\n\n');

% 1. Check if MATLAB is using correct area
fprintf('1. Cross-sectional area check:\n');
r1 = geom.R_cyl + geom.L_1 * 0.05;  % Start of stage 1
r2 = r1 + geom.L_1;
A_annular = (theta/2) * (r2^2 - r1^2);  % Annular sector area (one leg)
fprintf('   Annular sector area (one leg, stage 1): %.3e m²\n', A_annular);
fprintf('   Cross-section at mid-radius: θ/2 * r_mid * t = %.3e m²\n', (theta/2) * ((r1+r2)/2) * t_tec);
L_stage1 = geom.L_1;
G_simple = L_stage1 / (A_annular * t_tec);
fprintf('   Simple G = L/A = %.2e 1/m (if using annular area × t)\n', G_simple);
fprintf('   This doesn''t match either...\n\n');

% 2. Check units
fprintf('2. Units check:\n');
fprintf('   G is in 1/m (correct for R = rho*G or R = G/k)\n');
fprintf('   k is in W/mK (correct)\n');
fprintf('   R = G/k has units [1/m] / [W/mK] = K/W (correct)\n\n');

% 3. Compare to simple estimate
fprintf('3. Simple thermal resistance estimate:\n');
% For a cylinder shell: R = ln(r_out/r_in) / (2πkL)
% For a wedge (angle θ): R = ln(r_out/r_in) / (k*θ*t)
r_in_total = geom.R_cyl + geom.L_1 * 0.05;
r_out_total = geom.R_base;
R_simple_total = (1 / (k_TE * theta * t_tec)) * log(r_out_total / r_in_total);
fprintf('   Simple radial R (full wedge, r_cyl to R_base): %.1f K/W\n', R_simple_total);
fprintf('   This is also very high compared to COMSOL (%.1f K/W)\n\n', R_apparent);

% 4. What if effective cross-section is much larger?
fprintf('4. What effective cross-section would give COMSOL''s R?\n');
L_total = r_out_total - r_in_total;
A_required = L_total / (k_TE * R_apparent);
fprintf('   Required A_eff = L_total / (k * R_apparent) = %.3e m²\n', A_required);
fprintf('   Actual annular area = %.3e m²\n', (theta/2) * (r_out_total^2 - r_in_total^2) * t_tec);
fprintf('   Ratio = %.1fx\n', A_required / ((theta/2) * (r_out_total^2 - r_in_total^2) * t_tec));
fprintf('\n');

% 5. Heat spreading in chip
fprintf('5. Chip layer contribution:\n');
R_chip_spreading = (1 / (k_Si * theta * t_chip)) * log(r_out_total / r_in_total);
fprintf('   Chip spreading resistance: %.1f K/W\n', R_chip_spreading);
fprintf('   This is much lower than TEC resistance due to high k_Si\n\n');

% 6. What if we have geometry parameter mismatch?
fprintf('6. Checking geometry parameters vs COMSOL:\n');
fprintf('   Does COMSOL use same R_cyl = %.0f µm?\n', geom.R_cyl*1e6);
fprintf('   Does COMSOL use same t_TEC = %.0f µm?\n', t_tec*1e6);
fprintf('   Does COMSOL use same θ = %.1f°?\n', theta*180/pi);
fprintf('   Does COMSOL use same k_Bi2Te3 = %.2f W/mK?\n', k_TE);
fprintf('\n');

fprintf('=== CONCLUSION ===\n');
fprintf('The 14x discrepancy cannot be explained by:\n');
fprintf('  - Azimuthal conduction (already added 5x multiplier, minimal effect)\n');
fprintf('  - Vertical conduction (insulated in COMSOL)\n');
fprintf('  - Simple formula errors (both methods give similar high R)\n\n');
fprintf('Most likely causes:\n');
fprintf('  1. Geometry parameter mismatch between MATLAB and COMSOL\n');
fprintf('  2. Different thermal conductivity values\n');
fprintf('  3. Different definition of where heat enters (chip area vs TEC area)\n');
fprintf('  4. COMSOL includes conduction paths not in MATLAB model\n');
