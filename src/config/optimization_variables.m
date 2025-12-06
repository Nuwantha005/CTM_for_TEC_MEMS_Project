function [var_names, lb, ub, x0, all_vars, CONFIG] = optimization_variables()
% OPTIMIZATION_VARIABLES - Central configuration for all optimization scripts
%
% Uses notation from: Thermal_Network_For_Radial_TEC.tex
%
% This function provides a SINGLE CONTROL POINT for:
%   1. Optimization variables (bounds, defaults, enable/disable)
%   2. Fixed/Integer parameters (N, M, T_target)
%   3. Boundary conditions (q_flux, h_conv, T_water)
%   4. Base geometry defaults
%
% IMPORTANT: M (number of wedges) is an INTEGER parameter, not continuous.
%   θ = 360°/M must divide evenly, so M must be an integer.
%   Examples: M=4 → θ=90°, M=6 → θ=60°, M=12 → θ=30°, M=36 → θ=10°
%
% Used by: TECOptimizer, run_global_optimization, run_multiobjective_optimization,
%          StageSweeper, run_optimization
%
% Usage:
%   [var_names, lb, ub, x0, all_vars, CONFIG] = optimization_variables();
%
% Returns:
%   var_names - Cell array of enabled variable names
%   lb        - Lower bounds for enabled variables
%   ub        - Upper bounds for enabled variables
%   x0        - Initial values for enabled variables
%   all_vars  - Full configuration table (for reference)
%   CONFIG    - Struct with fixed parameters and boundary conditions
%
% ═══════════════════════════════════════════════════════════════════════════
% FIXED / INTEGER PARAMETERS (Design Requirements)
% These are NOT continuous optimization variables
% ═══════════════════════════════════════════════════════════════════════════

CONFIG = struct();

% --- Target & Constraints ---
CONFIG.T_target_C = 85;              % Target max temperature (°C)

% --- INTEGER PARAMETERS (must be swept, not optimized continuously) ---
CONFIG.N = 3;                        % N: Number of TEC stages (integer)
CONFIG.M = 12;                       % M: Number of wedges (integer)
% θ = 360°/M = 30° for M=12

% Derived from M: wedge angle in degrees and radians
CONFIG.theta_deg = 360 / CONFIG.M;   % θ = 360°/M
CONFIG.theta_rad = 2*pi / CONFIG.M;  % θ = 2π/M

% Valid M values for sweeps (divisors of 360 for nice angles)
CONFIG.M_valid = [4, 6, 8, 9, 10, 12, 15, 18, 20, 24, 30, 36, 40, 45, 60, 72, 90, 120, 180, 360];

% --- Boundary Conditions ---
CONFIG.q_flux_W_m2 = 50000;          % Heat flux at chip (W/m²)
CONFIG.h_conv_W_m2K = 1e6;           % Convection coefficient (W/m²K)
CONFIG.T_water_K = 300;              % T_water: Coolant temperature (K)

% --- Chip Geometry (Fixed) ---
CONFIG.L_chip_um = 10000;            % L_chip: Chip length (µm)
CONFIG.W_chip_um = 10000;            % W_chip: Chip width (µm)
CONFIG.t_chip_um = 50;               % t_chip: Chip thickness (µm)

% --- Material Properties ---
% Thermoelectric: Bi2Te3
CONFIG.materials.Bi2Te3 = struct('k', 1.2, 'rho', 1e-5, 'S', 0.0002);
% Interconnects: Cu
CONFIG.materials.Cu = struct('k', 400, 'rho', 1.7e-8);
% Chip/Substrate: Si
CONFIG.materials.Si = struct('k', 150, 'rho', 0.01);
% Radial Insulator: AlN
CONFIG.materials.AlN = struct('k', 170, 'rho', 1e10);
% Azimuthal Insulator: SiO2
CONFIG.materials.SiO2 = struct('k', 1.4, 'rho', 1e14);
% Vertical Insulator: Al2O3 or SiO2
CONFIG.materials.Al2O3 = struct('k', 30, 'rho', 1e12);

% ═══════════════════════════════════════════════════════════════════════════
% CONTINUOUS OPTIMIZATION VARIABLES
% NOTE: M (number of wedges) is NOT included here - it's an integer parameter
% ═══════════════════════════════════════════════════════════════════════════
% Set 'enabled' (5th column) to true/false to include/exclude each variable.
% Disabled variables use their initial value as a constant.
%
% Format: {name, lower_bound, upper_bound, initial_value, enabled, description}
% ═══════════════════════════════════════════════════════════════════════════

all_vars = {
    % ═══════════════════════════════════════════════════════════════════════
    % ELECTRICAL (Operating Current)
    % ═══════════════════════════════════════════════════════════════════════
    'I',                            0.0005, 1.0,   0.025,  true,   'I: Operating current [A]';

    % ═══════════════════════════════════════════════════════════════════════
    % GEOMETRY - PRIMARY DIMENSIONS
    % NOTE: theta_deg is REMOVED - use M (integer) instead
    % ═══════════════════════════════════════════════════════════════════════
    't_TEC_um',                     50,     500,  200,    true,   't_TEC: TEC layer thickness [µm]';
    'r_cyl_um',                     500,    5000,  1000,   true,   'r_cyl: Central cylinder radius [µm]';
    't_ins_um',                     5,      100,   10,     true,   't_ins: Insulator layer thickness [µm]';

    % ═══════════════════════════════════════════════════════════════════════
    % GEOMETRY - DIMENSIONAL REDUCTION FACTORS
    % ═══════════════════════════════════════════════════════════════════════
    'f_L',                          0.2,    3.0,   1.15,   true,   'f_L: Length ratio L_{i+1}/L_i';
    'W_az_um',                      5,      100,   20,     true,   'W_az: Azimuthal insulator width [µm]';
    'W_is_ratio',                   0.02,   0.15,  0.05,   true,   'W_is/L: Radial insulator width ratio';

    % ═══════════════════════════════════════════════════════════════════════
    % INTERCONNECTOR (ic) - Cold side electrical contact
    % W_ic = f_{ic,W} · L, t_ic = f_{ic,t} · t_TEC, β_ic = f_{ic,β} · θ
    % ═══════════════════════════════════════════════════════════════════════
    'f_ic_W',                       0.1,    0.35,  0.15,   true,   'f_{ic,W}: Interconnector width ratio W_ic/L';
    'f_ic_t',                       0.5,    2.0,   1.0,    true,   'f_{ic,t}: Interconnector thickness ratio t_ic/t_TEC';
    'f_ic_beta',                    0.1,    0.4,   0.16,   true,   'f_{ic,β}: Interconnector angle ratio β_ic/θ';

    % ═══════════════════════════════════════════════════════════════════════
    % OUTERCONNECTOR (oc) - Hot side electrical contact
    % W_oc = f_{oc,W} · L, t_oc = f_{oc,t} · t_TEC, β_oc = f_{oc,β} · θ
    % ═══════════════════════════════════════════════════════════════════════
    'f_oc_W',                       0.1,    0.35,  0.15,   true,   'f_{oc,W}: Outerconnector width ratio W_oc/L';
    'f_oc_t',                       0.5,    2.0,   1.0,    true,   'f_{oc,t}: Outerconnector thickness ratio t_oc/t_TEC';
    'f_oc_beta',                    0.1,    0.4,   0.16,   true,   'f_{oc,β}: Outerconnector angle ratio β_oc/θ';
    };

% ═══════════════════════════════════════════════════════════════════════════
% Filter to only enabled variables
% ═══════════════════════════════════════════════════════════════════════════
enabled_mask = [all_vars{:, 5}];
vars = all_vars(enabled_mask, 1:4);

var_names = vars(:, 1);
lb = [vars{:, 2}]';
ub = [vars{:, 3}]';
x0 = [vars{:, 4}]';

% Store all_vars in CONFIG for objective functions
CONFIG.all_vars = all_vars;
end
