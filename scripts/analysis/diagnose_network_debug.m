% diagnose_network_debug.m
% Deep dive into thermal network equations

clear; clc;
addpath(genpath('src'));

%% Load config
config = jsondecode(fileread('src/config/default_params.json'));
config.boundary_conditions.q_flux_W_m2 = 100;  % Very low flux
config.operating_conditions.I_current_A = 0.05;  % Low current

fprintf('=== THERMAL NETWORK DEBUG ===\n\n');

%% Setup
geometry = TECGeometry(config);
materials = MaterialProperties(config);
N = geometry.N_stages;

fprintf('--- Configuration ---\n');
fprintf('N_stages: %d\n', N);
fprintf('q_flux: %.0f W/m²\n', config.boundary_conditions.q_flux_W_m2);
fprintf('I: %.3f A\n', config.operating_conditions.I_current_A);
fprintf('T_water: %.0f K\n', config.boundary_conditions.T_water_K);

%% Build thermal network manually step by step
T_water = config.boundary_conditions.T_water_K;
q_flux = config.boundary_conditions.q_flux_W_m2;
I = config.operating_conditions.I_current_A;
theta = geometry.WedgeAngle;
R_base = geometry.R_base;
R_cyl = geometry.R_cyl;
t_chip = config.geometry.t_chip_um * 1e-6;
t_tec = geometry.Thickness;

%% Calculate heat input
A_total = 0.5 * theta * R_base^2;
A_cyl = 0.5 * theta * R_cyl^2;
Q_total = q_flux * A_total;
Q_cyl = q_flux * A_cyl;

fprintf('\n--- Heat Input ---\n');
fprintf('Total wedge area: %.4f mm²\n', A_total*1e6);
fprintf('Center cylinder area: %.4f mm²\n', A_cyl*1e6);
fprintf('Total heat input: %.4f mW\n', Q_total*1000);
fprintf('Heat at center: %.4f mW\n', Q_cyl*1000);

%% Collect stage properties
fprintf('\n--- Stage Properties ---\n');
K_stages = zeros(N, 1);
Re_stages = zeros(N, 1);
R_TSV = zeros(N, 1);
R_lat_Si = zeros(N, 1);
Q_gen = zeros(N, 1);
S_stages = zeros(N, 1);

for i = 1:N
    [r_in, L, w_ic, t_ic, beta_ic, w_oc, t_oc, beta_oc, w_az, w_is] = geometry.get_stage_geometry(i);
    r_out = r_in + L;
    
    T_avg = 320;
    k_te = materials.get_k('Bi2Te3', T_avg);
    rho_te = materials.get_rho('Bi2Te3', T_avg);
    S_te = materials.get_S('Bi2Te3', T_avg);
    k_is = materials.get_k('AlN', T_avg);
    k_az = materials.get_k('SiO2', T_avg);
    k_Si = materials.get_k('Si', T_avg);
    k_Cu = materials.get_k('Cu', T_avg);
    rho_Cu = materials.get_rho('Cu', T_avg);
    
    % G factor and thermal conductance
    G = geometry.calculate_G(r_in, L, w_ic, t_ic, beta_ic, w_oc, t_oc, beta_oc, w_az, w_is);
    K_legs = 2 * k_te / G;
    
    % Thermal insulator
    r_end_leg = r_in + L - w_is;
    R_is = geometry.calculate_R_thermal_insulator(r_end_leg, w_is, t_tec, theta, k_is);
    K_eff = 1 / (R_is + 1/K_legs);
    
    % Azimuthal
    K_az = geometry.calculate_K_azimuthal(r_in, L, w_az, t_tec, k_az, theta);
    
    K_stages(i) = K_eff + K_az;
    
    % Electrical resistance
    Re_stages(i) = 2 * rho_te * G;
    
    % Seebeck
    S_stages(i) = 2 * S_te;
    
    % TSV vertical resistance
    [~, R_TSV(i)] = geometry.calculate_TSV_vertical_resistance(r_in, w_ic, beta_ic, k_Cu, i);
    
    % Lateral Si resistance
    if i < N
        [r_in_next, L_next] = geometry.get_stage_geometry(i+1);
        r_mid_i = r_in + L/2;
        r_mid_next = r_in_next + L_next/2;
        R_lat_Si(i) = log(r_mid_next/r_mid_i) / (k_Si * theta * t_chip);
    else
        R_lat_Si(i) = inf;
    end
    
    % Heat generation at this stage
    A_stage = 0.5 * theta * (r_out^2 - r_in^2);
    Q_gen(i) = q_flux * A_stage;
    
    fprintf('Stage %d:\n', i);
    fprintf('  K_stage = %.4e W/K (thermal)\n', K_stages(i));
    fprintf('  Re = %.4f Ohm (electrical)\n', Re_stages(i));
    fprintf('  S = %.1f uV/K\n', S_stages(i)*1e6);
    fprintf('  R_TSV = %.4e K/W\n', R_TSV(i));
    fprintf('  R_lat_Si = %.4e K/W\n', R_lat_Si(i));
    fprintf('  Q_gen = %.4f mW\n', Q_gen(i)*1000);
end

%% Analyze heat paths
fprintf('\n--- Heat Flow Analysis ---\n');

% If all heat went through TEC stages in series:
R_TEC_total = sum(1./K_stages);
K_TEC_total = 1/R_TEC_total;
dT_through_TEC = Q_total / K_TEC_total;
fprintf('If heat flowed only through TEC:\n');
fprintf('  Total TEC thermal conductance: %.4e W/K\n', K_TEC_total);
fprintf('  dT = %.1f K\n', dT_through_TEC);

% TSV path (parallel to TEC)
R_TSV_total = sum(R_TSV(isfinite(R_TSV)));
fprintf('TSV path (stages with TSVs): R = %.4e K/W\n', R_TSV_total);

%% The actual matrix system
fprintf('\n--- Matrix Analysis ---\n');

network = ThermalNetwork(geometry, materials, config);
dim = 2*N + 1;
T_init = ones(dim, 1) * 300;

[M, B] = network.assemble_system(T_init);

fprintf('Matrix M (%dx%d):\n', dim, dim);
disp(full(M));
fprintf('\nVector B:\n');
disp(B);

% Check matrix conditioning
cond_M = condest(M);
fprintf('Matrix condition number: %.2e\n', cond_M);

% Solve
T_sol = M \ B;
fprintf('\nSolution:\n');
fprintf('T_0 = %.1f K (%.1f°C)\n', T_sol(1), T_sol(1)-273.15);
for i = 1:N
    fprintf('T_Si,%d = %.1f K, T_c,%d = %.1f K\n', ...
        i, T_sol(1+i), i, T_sol(N+1+i));
end

%% Check energy balance
fprintf('\n--- Energy Balance Check ---\n');
fprintf('Total heat in (q*A): %.4f mW\n', Q_total*1000);

% Heat removed by final stage TEC
T_c_N = T_sol(N+1+N);  % Cold side of outer stage
Q_peltier = S_stages(N) * I * T_c_N;
Q_joule = I^2 * Re_stages(N);
Q_backconduction = K_stages(N) * (T_water - T_c_N);
Q_out = Q_peltier + Q_backconduction - Q_joule/2;

fprintf('Heat out (TEC N):\n');
fprintf('  Peltier: %.4f mW\n', Q_peltier*1000);
fprintf('  Back-cond: %.4f mW\n', Q_backconduction*1000);
fprintf('  Joule/2: %.4f mW\n', Q_joule/2*1000);
fprintf('  Net out: %.4f mW\n', Q_out*1000);

%% Check if there's a sign error in the Peltier term
fprintf('\n--- Checking Peltier Direction ---\n');
fprintf('Current direction: flowing from center to edge (cold to hot side)\n');
fprintf('For TEC cooling, Peltier should ABSORB heat at cold side.\n');
fprintf('Q_peltier = S*I*T_c = %.4f mW\n', Q_peltier*1000);
if Q_peltier > 0
    status = 'POSITIVE';
else
    status = 'NEGATIVE';
end
fprintf('This is %s (should be positive for heat absorption)\n', status);

%% The issue: check the matrix coefficient for Peltier
fprintf('\n--- Matrix Coefficient Check ---\n');
fprintf('Looking at equation for cold side node...\n');
fprintf('The S*I term appears in M matrix as coefficient of T.\n');
fprintf('If S*I*T_c absorbs heat at cold side:\n');
fprintf('  - It should ADD to the node (positive in energy balance)\n');
fprintf('  - In matrix form: this means NEGATIVE coefficient in M\n');

% Check actual matrix structure
fprintf('\nActual M matrix coefficients for T_c nodes:\n');
for i = 1:N
    idx_c = N+1+i;
    fprintf('Row %d (T_c,%d): ', idx_c, i);
    fprintf('diag=%.4f, ', M(idx_c, idx_c));
    if i > 1
        fprintf('prev=%.4f, ', M(idx_c, idx_c-1));
    end
    if i < N
        fprintf('next=%.4f', M(idx_c, idx_c+1));
    end
    fprintf('\n');
end

%% Try without TEC (I=0) to check passive thermal resistance
fprintf('\n--- Without TEC (I=0) ---\n');
config2 = config;
config2.operating_conditions.I_current_A = 0;

network2 = ThermalNetwork(geometry, materials, config2);
[M2, B2] = network2.assemble_system(T_init);
T_passive = M2 \ B2;

fprintf('Passive thermal solution:\n');
fprintf('T_0 = %.1f K (%.1f°C)\n', T_passive(1), T_passive(1)-273.15);
fprintf('dT from water: %.1f K\n', T_passive(1) - T_water);

% Expected dT from simple thermal resistance
fprintf('\nExpected dT = Q/K = %.1f K\n', Q_total / K_TEC_total);
