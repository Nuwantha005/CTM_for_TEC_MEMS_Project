% diagnose_optimal_current.m
% Find optimal current for maximum cooling

clear; clc;
addpath(genpath('src'));

%% Load config
config = jsondecode(fileread('src/config/default_params.json'));
config.boundary_conditions.q_flux_W_m2 = 1000;  % 1 kW/m²

fprintf('=== FINDING OPTIMAL CURRENT ===\n\n');

%% Setup
geometry = TECGeometry(config);
materials = MaterialProperties(config);

N = geometry.N_stages;
theta = geometry.WedgeAngle;
R_base = geometry.R_base;
A_wedge = 0.5 * theta * R_base^2;

T_avg = 350;  % Estimate

%% Get stage 1 properties (outer stage, where cooling happens)
[r_in, L, w_ic, t_ic, beta_ic, w_oc, t_oc, beta_oc, w_az, w_is] = geometry.get_stage_geometry(N);
t = geometry.Thickness;

k_te = materials.get_k('Bi2Te3', T_avg);
S = materials.get_S('Bi2Te3', T_avg);
rho_te = materials.get_rho('Bi2Te3', T_avg);
k_is = materials.get_k('AlN', T_avg);

G = geometry.calculate_G(r_in, L, w_ic, t_ic, beta_ic, w_oc, t_oc, beta_oc, w_az, w_is);
K_legs = 2 * k_te / G;
Re_legs = 2 * rho_te * G;

[R_ic, R_oc] = geometry.calculate_R_electrical_interconnects(r_in, L, w_ic, t_ic, beta_ic, w_oc, t_oc, beta_oc, materials.get_rho('Cu', T_avg));
R_total = Re_legs + R_ic + R_oc;

r_end_leg = r_in + L - w_is;
R_is = geometry.calculate_R_thermal_insulator(r_end_leg, w_is, t, theta, k_is);
K_eff = 1 / (R_is + 1/K_legs);

fprintf('Stage %d (outer stage) properties:\n', N);
fprintf('K_eff = %.4e W/K\n', K_eff);
fprintf('R_elec = %.4f Ohm\n', R_total);
fprintf('2S = %.1f uV/K\n', 2*S*1e6);

%% Optimal current analysis
% For a simple TEC: Q_c = S*I*T_c - 0.5*I²R - K*(T_h - T_c)
% dQ_c/dI = S*T_c - I*R = 0
% I_opt = S*T_c/R

T_c_estimate = 350;  % Cold side temperature estimate
I_opt = 2*S * T_c_estimate / R_total;
fprintf('\n--- Optimal Current ---\n');
fprintf('For T_c = %.0f K, I_opt = %.4f A\n', T_c_estimate, I_opt);
fprintf('Current in config: %.4f A\n', config.operating_conditions.I_current_A);

%% Scan currents
fprintf('\n--- Current Scan ---\n');
fprintf('I (mA) | Q_peltier | Q_joule/2 | Q_back | Q_net (mW)\n');
fprintf('-------|-----------|-----------|--------|----------\n');

T_c = 373;  % 100C
T_h = 300;  % Water temp

currents = linspace(0.001, 0.2, 20);
Q_nets = zeros(size(currents));

for j = 1:length(currents)
    I = currents(j);
    Q_peltier = 2*S * I * T_c;
    Q_joule = I^2 * R_total;
    Q_back = K_eff * (T_h - T_c);  % Negative since T_h < T_c
    Q_net = Q_peltier - Q_joule/2 + Q_back;
    Q_nets(j) = Q_net;
    
    fprintf('%6.1f | %9.3f | %9.3f | %6.3f | %9.3f\n', ...
        I*1000, Q_peltier*1000, Q_joule/2*1000, Q_back*1000, Q_net*1000);
end

[Q_max, idx] = max(Q_nets);
I_max_Q = currents(idx);
fprintf('\nMax Q_net = %.3f mW at I = %.1f mA\n', Q_max*1000, I_max_Q*1000);

%% Check ZT factor
fprintf('\n--- ZT Analysis ---\n');
ZT = S^2 * T_avg / (k_te * rho_te);
fprintf('ZT (single leg) = %.3f\n', ZT);
fprintf('ZT (TE pair, 2S²) = %.3f\n', 4*S^2 * T_avg / (2*k_te * 2*rho_te));

% Maximum COP for single-stage TEC
dT = T_c - T_h;  % Note: hot side is colder (water)!
% Actually for cooling, cold side should be at chip, hot at water
% But our geometry is inverted - chip is hot, water is cold heat sink

fprintf('\nNOTE: In this geometry:\n');
fprintf('  - Chip (center) = HOT side with heat load\n');
fprintf('  - Water (outside) = COLD side (heat sink at 300K)\n');
fprintf('  - TEC must pump heat OUTWARD (from center to edge)\n');

% For cooling chip, we need T_chip < some target
% TEC absorbs heat at T_h (chip) and rejects at T_c (water)
% Q_h = Q_c + electrical power

%% The real issue: thermal resistance
fprintf('\n--- THE REAL PROBLEM ---\n');
Q_load = config.boundary_conditions.q_flux_W_m2 * A_wedge;
fprintf('Heat load per wedge: %.3f mW\n', Q_load*1000);
fprintf('Max TEC cooling: %.3f mW\n', Q_max*1000);

if Q_max < Q_load
    fprintf('\n*** FUNDAMENTAL LIMIT ***\n');
    fprintf('TEC cannot pump enough heat!\n');
    fprintf('Deficit: %.3f mW\n', (Q_load - Q_max)*1000);
    fprintf('This heat must be conducted, causing dT = %.1f K\n', (Q_load - Q_max)/K_eff);
else
    fprintf('\nTEC CAN handle this heat load.\n');
end

%% What heat flux CAN be handled?
q_feasible = Q_max / A_wedge;
fprintf('\n--- Feasible Heat Flux ---\n');
fprintf('Max cooling capacity: %.3f mW per wedge\n', Q_max*1000);
fprintf('This corresponds to q = %.0f W/m²\n', q_feasible);
fprintf('Or total chip power = %.2f mW\n', q_feasible * A_wedge * 12 * 1000);

%% Try with optimal current in full model
fprintf('\n--- Running with Optimal Current ---\n');
config.operating_conditions.I_current_A = I_max_Q;
config.boundary_conditions.q_flux_W_m2 = 100;  % Very low flux

geometry2 = TECGeometry(config);
network = ThermalNetwork(geometry2, materials, config);

dim = 2*N + 1;
T_init = ones(dim, 1) * 300;
T_current = T_init;

for iter = 1:100
    T_old = T_current;
    [T_new, ~, ~] = network.solve(T_current);
    T_current = T_new;
    
    if max(abs(T_new - T_old)) < 1e-4
        break;
    end
end

fprintf('At q=100 W/m², I=%.1f mA:\n', I_max_Q*1000);
fprintf('T_chip = %.1f K (%.1f°C)\n', T_current(1), T_current(1)-273.15);
