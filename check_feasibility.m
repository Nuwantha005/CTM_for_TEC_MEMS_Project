%% QUICK FEASIBILITY CHECK FOR TEC DESIGN
% This script quickly checks if a TEC design is feasible for your heat flux
% and provides analytical estimates before running full optimization.
%
% Author: Auto-generated
% Date: 2025-11-26

clear; clc;

%% Add paths
addpath(genpath('src'));

%% User Inputs
fprintf('\n========================================\n');
fprintf('  TEC FEASIBILITY QUICK CHECK\n');
fprintf('========================================\n');

% Chip parameters
w_chip_mm = 10;         % Chip width (mm)
q_flux_kW_m2 = 100;     % Heat flux (kW/m²)

% Cooling parameters
T_water_C = 27;         % Water temperature (°C)
T_max_C = 100;          % Max allowable chip temperature (°C)

% TEC material (Bi2Te3 typical values)
S_pn = 0.4e-3;          % Seebeck coefficient (V/K) - p+n combined
k_pn = 2.4;             % Thermal conductivity (W/m-K) - p+n
rho_pn = 2e-5;          % Electrical resistivity (ohm-m) - p+n

%% Convert units
w_chip = w_chip_mm * 1e-3;  % m
q_flux = q_flux_kW_m2 * 1e3; % W/m²
T_water_K = T_water_C + 273.15;
T_max_K = T_max_C + 273.15;

%% Calculate problem parameters
A_chip = w_chip^2;              % Chip area (m²)
Q_total = q_flux * A_chip;      % Total heat (W)
dT_max = T_max_K - T_water_K;   % Maximum temperature rise (K)

fprintf('\nProblem Parameters:\n');
fprintf('  Chip size: %.1f mm x %.1f mm = %.2f cm²\n', w_chip_mm, w_chip_mm, A_chip*1e4);
fprintf('  Heat flux: %.0f kW/m² = %.0f W/cm²\n', q_flux_kW_m2, q_flux/1e4);
fprintf('  Total chip heat: %.2f W\n', Q_total);
fprintf('  Required cooling: keep chip below %.0f°C\n', T_max_C);
fprintf('  Max temperature rise: %.1f K\n', dT_max);

%% TEC Material Figure of Merit
ZT = (S_pn^2 * T_water_K) / (k_pn * rho_pn);
fprintf('\nMaterial Properties (Bi2Te3):\n');
fprintf('  ZT @ %.0f K: %.2f\n', T_water_K, ZT);

%% Theoretical TEC Limits
% Maximum cooling capacity for ideal single-stage TEC
% Q_c_max per unit G = 0.5 * S² * T_c² / ρ - k * dT
% where G is the geometric factor (length/area)

% At optimal current: Q_c = 0.5 * ZT * K * T_c - K * dT
% For net cooling: 0.5 * ZT * T_c > dT
cooling_limit_dT = 0.5 * ZT * T_water_K;

fprintf('\nTheoretical Limits:\n');
fprintf('  Max dT for single-stage cooling: %.1f K\n', cooling_limit_dT);
fprintf('  Your required dT: %.1f K\n', dT_max);

if dT_max > cooling_limit_dT
    N_min = ceil(dT_max / cooling_limit_dT);
    fprintf('  ** Single stage INSUFFICIENT! **\n');
    fprintf('  Minimum stages required: %d (estimated)\n', N_min);
else
    fprintf('  Single stage MAY be sufficient.\n');
end

%% Maximum COP estimate
% COP = (T_c / dT) * [sqrt(1+ZT) - T_h/T_c] / [sqrt(1+ZT) + 1]
sqrt_term = sqrt(1 + ZT);
if dT_max > 0
    COP_max = (T_water_K / dT_max) * (sqrt_term - T_max_K/T_water_K) / (sqrt_term + 1);
else
    COP_max = inf;
end

fprintf('  Theoretical max COP: %.3f\n', max(0, COP_max));

if COP_max <= 0
    fprintf('  ** WARNING: COP <= 0 means TEC cannot provide net cooling! **\n');
    fprintf('  The temperature gradient is too steep for effective TEC cooling.\n');
end

%% Analyze wedge configurations
fprintf('\n========================================\n');
fprintf('  WEDGE CONFIGURATION ANALYSIS\n');
fprintf('========================================\n');
fprintf('θ (°) | Wedges | Q/wedge (mW) | Feasibility\n');
fprintf('------|--------|--------------|------------\n');

theta_list = [15, 20, 25, 30, 36, 45, 60, 72, 90];
for theta = theta_list
    n_wedges = floor(360 / theta);
    Q_wedge = Q_total / n_wedges * 1000;  % mW
    
    % Rough feasibility: micro-TEC can handle ~10-100 mW typically
    if Q_wedge < 50
        feasibility = 'Good';
    elseif Q_wedge < 200
        feasibility = 'Marginal';
    elseif Q_wedge < 500
        feasibility = 'Challenging';
    else
        feasibility = 'Difficult';
    end
    
    fprintf('%5.0f | %6d | %12.1f | %s\n', theta, n_wedges, Q_wedge, feasibility);
end

%% Multi-stage analysis
fprintf('\n========================================\n');
fprintf('  MULTI-STAGE ANALYSIS\n');
fprintf('========================================\n');
fprintf('N | dT/stage (K) | Within limit? | Estimated feasibility\n');
fprintf('--|--------------|---------------|----------------------\n');

for N = 1:10
    dT_stage = dT_max / N;
    within = dT_stage < cooling_limit_dT;
    
    if within
        margin = (cooling_limit_dT - dT_stage) / cooling_limit_dT * 100;
        if margin > 30
            feas = 'Likely feasible';
        else
            feas = 'Possible';
        end
    else
        feas = 'Not feasible';
    end
    
    within_str = 'Yes';
    if ~within
        within_str = 'No';
    end
    
    fprintf('%d | %12.1f | %13s | %s\n', N, dT_stage, within_str, feas);
    
    if within && margin > 50
        break;  % Found good solution
    end
end

%% Recommendations
fprintf('\n========================================\n');
fprintf('  RECOMMENDATIONS\n');
fprintf('========================================\n');

if Q_total > 1  % More than 1W total
    fprintf('1. Heat load is HIGH (%.2f W). Consider:\n', Q_total);
    fprintf('   - Reducing heat flux at the source\n');
    fprintf('   - Using larger chip area to spread heat\n');
    fprintf('   - Combining TEC with other cooling methods\n');
end

if dT_max > cooling_limit_dT
    fprintf('2. Temperature gradient requires MULTI-STAGE design:\n');
    fprintf('   - Minimum %d stages recommended\n', ceil(dT_max / cooling_limit_dT));
    fprintf('   - More stages = lower dT per stage = easier cooling\n');
end

if COP_max < 0.1
    fprintf('3. COP is very low (%.3f). TEC will consume significant power.\n', max(0, COP_max));
    fprintf('   - Consider if active TEC cooling is the right approach\n');
    fprintf('   - Passive cooling or hybrid approaches may be better\n');
end

fprintf('\n4. For COMSOL simulation, try these starting parameters:\n');
if dT_max > cooling_limit_dT
    N_rec = min(8, max(2, ceil(dT_max / cooling_limit_dT * 1.5)));
else
    N_rec = 2;
end
theta_rec = 30;
t_rec = 50;

fprintf('   - N_stages: %d\n', N_rec);
fprintf('   - Wedge angle: %d°\n', theta_rec);
fprintf('   - TEC thickness: %d μm\n', t_rec);
fprintf('   - Start with I = 0.01-0.1 A and sweep\n');

fprintf('\n========================================\n');
fprintf('  Run run_physics_constrained.m for full optimization\n');
fprintf('========================================\n');
