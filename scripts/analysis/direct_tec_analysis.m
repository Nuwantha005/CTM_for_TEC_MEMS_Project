% direct_tec_analysis.m
% Direct analytical calculation of TEC cooling capacity
% Bypasses the thermal network matrix issues
%
% Uses simplified 1D TEC model with radial geometry corrections

clear; clc;

fprintf('╔════════════════════════════════════════════════════════════════╗\n');
fprintf('║     DIRECT TEC COOLING CAPACITY ANALYSIS                      ║\n');
fprintf('║     Analytical model bypassing matrix solver issues           ║\n');
fprintf('╚════════════════════════════════════════════════════════════════╝\n\n');

%% TEC Material Properties (Bi2Te3 at 300K)
S_p = 200e-6;    % Seebeck p-type, V/K
S_n = -200e-6;   % Seebeck n-type, V/K (negative)
S = S_p - S_n;   % Total Seebeck = 400 µV/K
k_te = 1.2;      % Thermal conductivity, W/m-K
rho_te = 1e-5;   % Electrical resistivity, Ohm-m

% Figure of merit
ZT = S^2 / (rho_te * k_te) * 300;
fprintf('Bi2Te3 Material:\n');
fprintf('  Seebeck (total): %.0f µV/K\n', S*1e6);
fprintf('  Thermal k: %.2f W/m-K\n', k_te);
fprintf('  Electrical ρ: %.1e Ω-m\n', rho_te);
fprintf('  ZT @ 300K: %.3f\n\n', ZT);

%% Geometry parameters
fprintf('Geometry Parameters:\n');
fprintf('─────────────────────────────────────────────────────────────────\n');

% Base geometry
R_chip = 1e-3;      % Chip radius, m (1 mm)
theta = deg2rad(30); % Wedge angle
n_wedges = round(2*pi / theta);

fprintf('  Chip radius: %.1f mm\n', R_chip*1000);
fprintf('  Wedge angle: %.0f° (%d wedges)\n', rad2deg(theta), n_wedges);

%% TEC Cooling Analysis
% For a single-stage TEC, the maximum cooling power is:
% Q_c_max = 0.5 * S² * T_c² / (ρ * L / A) - K * ΔT
%
% Where:
% - L = TEC leg length (thickness)
% - A = TEC leg cross-section area
% - K = thermal conductance = k * A / L
% - ΔT = T_h - T_c (temperature difference across TEC)
%
% Optimal current: I_opt = S * T_c / R
% Maximum cooling: Q_c_max = 0.5 * Z * T_c² / (L/A) - K * ΔT
%                         = (S² * T_c²) / (2 * ρ * L / A) - (k * A / L) * ΔT
%
% Per unit area:
% q_c_max = (S² * T_c²) / (2 * ρ * L) - (k / L) * ΔT

fprintf('\n═══════════════════════════════════════════════════════════════════\n');
fprintf('        MAXIMUM COOLING CAPACITY vs TEC THICKNESS\n');
fprintf('═══════════════════════════════════════════════════════════════════\n\n');

T_c = 353;  % Cold side (chip) at 80°C = 353K
T_h = 300;  % Hot side (coolant) at 27°C = 300K
DeltaT = T_c - T_h;  % = 53K (we're cooling FROM 80°C TO 27°C)

% Wait - if T_c > T_h, we're not cooling, we're heating!
% For cooling: T_c should be < T_h (chip colder than coolant)
% Rethink: We want the chip to be at max 80°C, with coolant at 27°C
% TEC pumps heat FROM chip (cold side) TO coolant (hot side)
% So T_c = T_chip, T_h = T_coolant + TEC temp rise
% If chip is at 80°C and coolant at 27°C, the TEC must pump heat 
% AGAINST a 53K gradient, which is very challenging for Bi2Te3!

% Let me reconsider:
% - Coolant is at 300K (27°C)
% - Chip should be at max 353K (80°C)
% - The TEC can only cool IF T_c < T_h
% - But the PASSIVE case: heat flows chip→coolant naturally if chip is hotter
% - TEC ACTIVE cooling: can make chip COLDER than coolant

fprintf('Scenario: Passive cooling (no TEC)\n');
fprintf('  Heat naturally flows from chip (80°C) to coolant (27°C)\n');
fprintf('  The TEC in passive mode just adds thermal resistance!\n\n');

% In ACTIVE TEC mode, we want: T_chip < T_coolant (chip colder than coolant)
% Let's find maximum heat flux when T_chip = T_coolant (ΔT = 0)
% and when T_chip = T_coolant - 10K (actual cooling)

T_coolant = 300;  % K

% For TEC legs, the cooling capacity per unit area at optimal current:
% q_c = (S²*T_c²)/(2*ρ*L) - (k/L)*ΔT
% At ΔT = 0: q_c_max = (S²*T_c²)/(2*ρ*L)

fprintf('Maximum cooling flux at ΔT=0 (T_chip = T_coolant):\n');
fprintf('─────────────────────────────────────────────────────────────────\n');
fprintf('  t_TEC (µm)    q_c_max (W/cm²)    I_opt (mA)    Power (W/cm²)\n');
fprintf('─────────────────────────────────────────────────────────────────\n');

thicknesses = [50, 100, 200, 500, 1000, 2000, 5000] * 1e-6;

for L = thicknesses
    % Maximum cooling at ΔT = 0
    T_c = T_coolant;  % At ΔT = 0
    q_c_max = (S^2 * T_c^2) / (2 * rho_te * L);  % W/m²
    
    % Optimal current density
    % I = S*T_c / R, where R = ρ*L/A, so I/A = S*T_c / (ρ*L)
    j_opt = S * T_c / (rho_te * L);  % A/m² (current density)
    
    % Voltage and power
    V = S * (T_h - T_c) + j_opt * rho_te * L;  % V (per leg pair)
    P_elec = j_opt * V;  % W/m² electrical power input
    
    % For a typical TEC leg cross-section (e.g., 100µm x 100µm)
    A_leg = (100e-6)^2;  % m²
    I_opt = j_opt * A_leg * 1000;  % mA
    
    fprintf('  %6.0f        %8.1f          %6.1f         %6.1f\n', ...
        L*1e6, q_c_max/10000, I_opt, P_elec/10000);
end

fprintf('\n');
fprintf('Note: These are THEORETICAL MAXIMUMS assuming perfect heat extraction\n');
fprintf('      from the TEC hot side to the coolant.\n\n');

%% Realistic cooling with temperature difference
fprintf('═══════════════════════════════════════════════════════════════════\n');
fprintf('        COOLING WITH ΔT (Chip hotter than coolant)\n');
fprintf('═══════════════════════════════════════════════════════════════════\n\n');

% If chip is hotter than coolant, we don't need TEC for cooling!
% The TEC would add thermal resistance and REDUCE passive cooling.
%
% TEC makes sense when:
% 1. We want chip COLDER than coolant (active cooling)
% 2. We want to boost cooling beyond passive conduction
%
% For 3D IC with hot top die, the TEC should:
% - Pump heat FROM bottom die TO coolant
% - Maintain bottom die below some threshold
%
% Let's calculate: Given a heat flux from top die onto bottom die,
% what's the achievable bottom die temperature?

fprintf('Scenario: Heat flux from top die impinging on bottom die\n');
fprintf('Goal: Keep bottom die temperature < 80°C\n');
fprintf('Coolant: 27°C (300K)\n\n');

T_max = 353;  % 80°C max
T_coolant = 300;  % 27°C

% For given heat flux q_in (W/m²), find if TEC can maintain T_chip < T_max
% Energy balance at chip:
% q_in = q_TEC + q_passive
% q_TEC = S*I*T_c - 0.5*I²*R - K*(T_c - T_h)
%
% But this is complex. Let's use COP approach instead.

fprintf('Maximum heat that TEC can remove at different ΔT:\n');
fprintf('─────────────────────────────────────────────────────────────────\n');

L = 200e-6;  % 200 µm TEC thickness
fprintf('TEC thickness: %.0f µm\n\n', L*1e6);

fprintf('  ΔT (K)    T_chip(°C)    q_c (W/cm²)    COP\n');
fprintf('─────────────────────────────────────────────────────────────────\n');

DeltaT_values = [-20, -10, -5, 0, 5, 10, 20, 30, 40, 50];

for DeltaT = DeltaT_values
    T_c = T_coolant + DeltaT;  % Chip temp
    T_h = T_coolant;            % Hot side (at coolant)
    
    if DeltaT < 0
        % Active cooling: chip colder than coolant
        % q_c = S²*T_c²/(2*ρ*L) - k/L * (T_h - T_c)
        q_c = (S^2 * T_c^2) / (2 * rho_te * L) - (k_te / L) * (T_h - T_c);
        
        % COP = q_c / P_elec (rough estimate)
        % At optimal current: COP ≈ T_c / (2*(T_h - T_c)) for small ZT
        if T_h ~= T_c
            COP = T_c / (2 * abs(T_h - T_c));
        else
            COP = inf;
        end
    else
        % Passive or negative case: chip hotter than coolant
        % Heat flows naturally from chip to coolant
        % TEC can augment this, but passive conduction dominates
        % q_c = k/L * (T_c - T_h) + TEC boost
        q_passive = (k_te / L) * (T_c - T_h);
        
        % TEC boost at optimal current
        q_tec_boost = (S^2 * T_c^2) / (2 * rho_te * L);
        q_c = q_passive + q_tec_boost;
        COP = inf;  % Not applicable for natural gradient
    end
    
    fprintf('  %+5.0f      %6.1f         %8.1f      %6.2f\n', ...
        DeltaT, T_c - 273, q_c/10000, COP);
end

fprintf('\n');
fprintf('Key insight:\n');
fprintf('  When ΔT > 0 (chip hotter than coolant), heat flows naturally!\n');
fprintf('  TEC ADDS to this by pumping more heat.\n');
fprintf('  At ΔT = +53K (80°C chip, 27°C coolant):\n');

DeltaT = 53;
T_c = T_coolant + DeltaT;
q_passive = (k_te / L) * DeltaT;
q_tec_boost = (S^2 * T_c^2) / (2 * rho_te * L);
q_total = q_passive + q_tec_boost;

fprintf('    Passive conduction: %.1f W/cm²\n', q_passive/10000);
fprintf('    TEC boost:          %.1f W/cm²\n', q_tec_boost/10000);
fprintf('    TOTAL:              %.1f W/cm²\n', q_total/10000);

fprintf('\n');
fprintf('═══════════════════════════════════════════════════════════════════\n');
fprintf('                    CONCLUSION\n');
fprintf('═══════════════════════════════════════════════════════════════════\n\n');

fprintf('For T_chip = 80°C, T_coolant = 27°C, t_TEC = 200µm:\n');
fprintf('  Maximum removable heat flux: %.1f W/cm² (%.0f kW/m²)\n', q_total/10000, q_total/1000);
fprintf('\n');
fprintf('This is in the range needed for:\n');
fprintf('  ✓ DRAM (1-5 W/cm²)\n');
fprintf('  ✓ SRAM Cache (10 W/cm²)\n');
fprintf('  ✓ IO Controllers (20 W/cm²)\n');
fprintf('  ✓ Media Engines (20 W/cm²)\n');
fprintf('  ✓ NPU (30 W/cm²)\n');

if q_total/10000 >= 50
    fprintf('  ✓ NPU MAC hotspots (50 W/cm²)\n');
else
    fprintf('  ✗ NPU MAC hotspots (50 W/cm²) - need thinner TEC or better materials\n');
end

%% Optimize TEC thickness for maximum heat removal at 80°C chip
fprintf('\n═══════════════════════════════════════════════════════════════════\n');
fprintf('        OPTIMIZATION: Best TEC thickness for 80°C chip\n');
fprintf('═══════════════════════════════════════════════════════════════════\n\n');

DeltaT = 53;  % Fixed: 80°C chip, 27°C coolant
T_c = T_coolant + DeltaT;

thicknesses_fine = logspace(-5, -2, 100);  % 10µm to 10mm

q_total_arr = zeros(size(thicknesses_fine));
for idx = 1:length(thicknesses_fine)
    L = thicknesses_fine(idx);
    q_passive = (k_te / L) * DeltaT;
    q_tec = (S^2 * T_c^2) / (2 * rho_te * L);
    q_total_arr(idx) = q_passive + q_tec;
end

[q_max, idx_max] = max(q_total_arr);
L_opt = thicknesses_fine(idx_max);

fprintf('Optimal TEC thickness: %.0f µm\n', L_opt*1e6);
fprintf('Maximum heat flux: %.1f W/cm² (%.0f kW/m²)\n', q_max/10000, q_max/1000);

% But wait - this is wrong! As L→0, both terms go to infinity!
% This means the model breaks down for very thin TEC.
%
% The real constraint is:
% 1. Minimum manufacturable thickness (~50µm)
% 2. Contact resistance dominates at thin films
% 3. Current carrying capacity limits

fprintf('\nNote: Thinner is better in this simple model, but:\n');
fprintf('  - Minimum practical thickness: ~50 µm\n');
fprintf('  - Contact resistance becomes significant at thin films\n');
fprintf('  - Current density limits for reliability\n');

L = 50e-6;  % Minimum practical
q_passive = (k_te / L) * DeltaT;
q_tec = (S^2 * T_c^2) / (2 * rho_te * L);
q_total = q_passive + q_tec;

fprintf('\nAt minimum practical thickness (50 µm):\n');
fprintf('  Maximum heat flux: %.1f W/cm² (%.0f kW/m²)\n', q_total/10000, q_total/1000);
