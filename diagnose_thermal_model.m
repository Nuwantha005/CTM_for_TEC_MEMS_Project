% diagnose_thermal_model.m
% Debug script to understand why temperatures are so high

clear; clc;
addpath(genpath('src'));

%% Load config
config = jsondecode(fileread('src/config/default_params.json'));

% Use low heat flux for testing
config.boundary_conditions.q_flux_W_m2 = 1000;  % 1 kW/m²

fprintf('=== THERMAL MODEL DIAGNOSTICS ===\n\n');

%% Geometry calculations
geometry = TECGeometry(config);
fprintf('--- Geometry ---\n');
fprintf('Chip width: %.2f mm\n', config.geometry.w_chip_um/1000);
fprintf('TEC thickness: %.1f um\n', config.geometry.thickness_um);
fprintf('Number of stages: %d\n', geometry.N_stages);
fprintf('Wedge angle: %.1f deg\n', rad2deg(geometry.WedgeAngle));
fprintf('R_cyl: %.2f mm\n', geometry.R_cyl * 1000);
fprintf('R_base: %.2f mm\n', geometry.R_base * 1000);
fprintf('L_1: %.2f um\n', geometry.L_1 * 1e6);

%% Calculate areas and heat loads
theta = geometry.WedgeAngle;
R_base = geometry.R_base;
A_wedge = 0.5 * theta * R_base^2;
n_wedges = round(2*pi / theta);

fprintf('\n--- Heat Load ---\n');
fprintf('Area per wedge: %.4f mm²\n', A_wedge * 1e6);
fprintf('Number of wedges: %d\n', n_wedges);
fprintf('Total chip area: %.4f mm²\n', A_wedge * n_wedges * 1e6);
fprintf('Heat flux: %.0f W/m²\n', config.boundary_conditions.q_flux_W_m2);
fprintf('Heat per wedge: %.4f mW\n', config.boundary_conditions.q_flux_W_m2 * A_wedge * 1000);
fprintf('Total heat: %.2f mW\n', config.boundary_conditions.q_flux_W_m2 * A_wedge * n_wedges * 1000);

%% Stage geometry
fprintf('\n--- Stage Geometry ---\n');
for i = 1:geometry.N_stages
    [r_in, L, w_ic, t_ic, beta_ic, w_oc, t_oc, beta_oc, w_az, w_is] = geometry.get_stage_geometry(i);
    fprintf('Stage %d: r_in=%.2f um, L=%.2f um, w_ic=%.2f um\n', ...
        i, r_in*1e6, L*1e6, w_ic*1e6);
end

%% Materials
materials = MaterialProperties(config);
T_avg = 320;  % K

fprintf('\n--- Material Properties @ %.0fK ---\n', T_avg);
fprintf('Bi2Te3: k=%.2f W/mK, S=%.1f uV/K, rho=%.2e Ohm-m\n', ...
    materials.get_k('Bi2Te3', T_avg), ...
    materials.get_S('Bi2Te3', T_avg)*1e6, ...
    materials.get_rho('Bi2Te3', T_avg));
fprintf('Cu: k=%.0f W/mK, rho=%.2e Ohm-m\n', ...
    materials.get_k('Cu', T_avg), materials.get_rho('Cu', T_avg));
fprintf('AlN: k=%.0f W/mK\n', materials.get_k('AlN', T_avg));

%% Calculate thermal resistances
fprintf('\n--- Thermal Resistances ---\n');
[r_in, L, w_ic, t_ic, beta_ic, w_oc, t_oc, beta_oc, w_az, w_is] = geometry.get_stage_geometry(1);
k_te = materials.get_k('Bi2Te3', T_avg);
k_is = materials.get_k('AlN', T_avg);
t = geometry.Thickness;

% G factor
G = geometry.calculate_G(r_in, L, w_ic, t_ic, beta_ic, w_oc, t_oc, beta_oc, w_az, w_is);
fprintf('G factor (stage 1): %.4e\n', G);

K_legs = 2 * k_te / G;
fprintf('K_legs (2 legs): %.4e W/K\n', K_legs);

% Thermal insulator resistance
r_end_leg = r_in + L - w_is;
R_is = geometry.calculate_R_thermal_insulator(r_end_leg, w_is, t, theta, k_is);
fprintf('R_insulator: %.4e K/W\n', R_is);

% Effective thermal conductance
R_eff = R_is + 1/K_legs;
K_eff = 1/R_eff;
fprintf('K_effective (stage 1): %.4e W/K\n', K_eff);

%% Expected temperature rise (simple estimate)
Q_wedge = config.boundary_conditions.q_flux_W_m2 * A_wedge;
dT_estimate = Q_wedge / K_eff;
fprintf('\n--- Simple Estimate ---\n');
fprintf('Heat to remove per wedge: %.4f mW\n', Q_wedge * 1000);
fprintf('If all heat went through K_eff, dT ~ %.1f K\n', dT_estimate);

%% TEC cooling capacity
S = materials.get_S('Bi2Te3', T_avg);
rho_te = materials.get_rho('Bi2Te3', T_avg);
I = config.operating_conditions.I_current_A;
T_c = 373;  % Target 100C
T_h = config.boundary_conditions.T_water_K;  % 300K

% Get electrical resistance
Re = 2 * rho_te * G;  % Two legs in series
[R_ic, R_oc] = geometry.calculate_R_electrical_interconnects(r_in, L, w_ic, t_ic, beta_ic, w_oc, t_oc, beta_oc, materials.get_rho('Cu', T_avg));
R_total = Re + R_ic + R_oc;

fprintf('\n--- TEC Electrical ---\n');
fprintf('Current: %.3f A\n', I);
fprintf('Seebeck coefficient (2S): %.1f uV/K\n', 2*S*1e6);
fprintf('Electrical resistance (legs): %.4e Ohm\n', Re);
fprintf('Interconnect resistance: %.4e Ohm\n', R_ic + R_oc);
fprintf('Total resistance: %.4e Ohm\n', R_total);

% Peltier cooling
Q_peltier = 2 * S * I * T_c;
Q_joule = I^2 * R_total;
Q_backconduction = K_eff * (T_h - T_c);  % Note: T_h < T_c so this is negative

fprintf('\n--- TEC Heat Terms @ T_c=%.0fK, I=%.3fA ---\n', T_c, I);
fprintf('Peltier cooling: %.4f mW (removes heat)\n', Q_peltier * 1000);
fprintf('Joule heating (half at cold): %.4f mW (adds heat)\n', Q_joule/2 * 1000);
fprintf('Back-conduction (T_h < T_c): %.4f mW\n', Q_backconduction * 1000);

Q_net_cold = Q_peltier - Q_joule/2 + Q_backconduction;
fprintf('Net heat pumped at cold: %.4f mW\n', Q_net_cold * 1000);

%% Run actual solve
fprintf('\n--- Running Thermal Network Solve ---\n');
network = ThermalNetwork(geometry, materials, config);

% Initial temperature guess
N = geometry.N_stages;
dim = 2*N + 1;
T_init = ones(dim, 1) * config.simulation.T_initial_guess;

% Iterate
T_current = T_init;
for iter = 1:50
    T_old = T_current;
    [T_new, ~, ~] = network.solve(T_current);
    T_current = T_new;
    
    err = max(abs(T_new - T_old));
    if iter <= 5 || mod(iter, 10) == 0
        fprintf('Iter %d: max_err=%.2e, T_chip=%.1fK (%.1f°C)\n', ...
            iter, err, T_new(1), T_new(1) - 273.15);
    end
    
    if err < 1e-4
        fprintf('Converged at iteration %d\n', iter);
        break;
    end
end

fprintf('\n--- Final Temperature Distribution ---\n');
fprintf('T_0 (chip center): %.1f K (%.1f°C)\n', T_current(1), T_current(1)-273.15);
for i = 1:N
    fprintf('T_Si,%d: %.1f K, T_c,%d: %.1f K\n', ...
        i, T_current(1+i), i, T_current(N+1+i));
end
fprintf('T_water: %.1f K\n', config.boundary_conditions.T_water_K);

%% Check TSV resistance
fprintf('\n--- TSV Analysis ---\n');
k_TSV = materials.get_k('Cu', T_avg);
for i = 1:min(3, N)
    [r_in, ~, w_ic, ~, beta_ic, ~, ~, ~, ~, ~] = geometry.get_stage_geometry(i);
    [N_TSV, R_TSV] = geometry.calculate_TSV_vertical_resistance(r_in, w_ic, beta_ic, k_TSV, i);
    fprintf('Stage %d: N_TSV=%d, R_TSV=%.4e K/W\n', i, N_TSV, R_TSV);
end

%% Check if TSV limit is causing issues
fprintf('\n--- Effect of N_tsv_limit ---\n');
fprintf('Current N_tsv_limit: %d\n', config.geometry.N_tsv_limit);
if config.geometry.N_tsv_limit < N
    fprintf('WARNING: TSV limit < N_stages means some stages have R_TSV = 1e9!\n');
    fprintf('This creates extremely high thermal resistance.\n');
end
