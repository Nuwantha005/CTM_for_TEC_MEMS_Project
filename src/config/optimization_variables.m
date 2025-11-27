function [var_names, lb, ub, x0, all_vars, CONFIG] = optimization_variables()
    % OPTIMIZATION_VARIABLES - Central configuration for all optimization scripts
    %
    % This function provides a SINGLE CONTROL POINT for:
    %   1. Optimization variables (bounds, defaults, enable/disable)
    %   2. Fixed parameters (N_stages, N_tsv_limit, T_target)
    %   3. Boundary conditions (q_flux, h_conv, T_water)
    %   4. Base geometry defaults
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
    % FIXED PARAMETERS (Design Requirements)
    % ═══════════════════════════════════════════════════════════════════════════
    
    CONFIG = struct();
    
    % --- Target & Constraints ---
    CONFIG.T_target_C = 85;              % Target max temperature (°C)
    CONFIG.N_stages = 3;                 % Number of TEC stages
    CONFIG.N_tsv_limit = 0;              % TSV limit (0 = no TSVs)
    
    % --- Boundary Conditions ---
    CONFIG.q_flux_W_m2 = 2000;            % Heat flux at chip (W/m²)
    CONFIG.h_conv_W_m2K = 1e6;           % Convection coefficient (W/m²K)
    CONFIG.T_water_K = 300;              % Coolant temperature (K)
    
    % --- Base Geometry (non-optimized defaults) ---
    CONFIG.w_chip_um = 10000;            % Chip width (µm)
    CONFIG.t_chip_um = 50;               % Chip thickness (µm)
    
    % --- TSV Parameters ---
    CONFIG.tsv.R_TSV_um = 10;            % TSV radius (µm)
    CONFIG.tsv.P_TSV_um = 20;            % TSV pitch (µm)
    CONFIG.tsv.g_rad_um = 10;            % Radial gap (µm)
    
    % --- Material Properties ---
    CONFIG.materials.Bi2Te3 = struct('k', 1.2, 'rho', 1e-5, 'S', 0.0002);
    CONFIG.materials.Cu = struct('k', 400, 'rho', 1.7e-8);
    CONFIG.materials.Si = struct('k', 150, 'rho', 0.01);
    CONFIG.materials.AlN = struct('k', 170, 'rho', 1e10);
    CONFIG.materials.SiO2 = struct('k', 1.4, 'rho', 1e14);
    CONFIG.materials.Al2O3 = struct('k', 30, 'rho', 1e12);
    
    % ═══════════════════════════════════════════════════════════════════════════
    % OPTIMIZATION VARIABLE CONFIGURATION
    % ═══════════════════════════════════════════════════════════════════════════
    % Set 'enabled' (5th column) to true/false to include/exclude each variable.
    % Disabled variables use their initial value as a constant.
    %
    % Format: {name, lower_bound, upper_bound, initial_value, enabled}
    % ═══════════════════════════════════════════════════════════════════════════
    
    all_vars = {
        % === ELECTRICAL ===
        'current',                      0.0005, 1.0,   0.025,  true;   % TEC current [A]
        
        % === GEOMETRY - MAIN ===
        'thickness_um',                 50,     2000,  200,    true;   % TEC element thickness [µm]
        'wedge_angle_deg',              10,     90,    30,     true;   % Wedge angle [deg]
        'R_cyl_um',                     500,    5000,  1000,   true;   % Cylinder radius [µm]
        't_SOI_um',                     20,     500,   100,    true;   % SOI thickness [µm]
        
        % === GEOMETRY - RATIOS ===
        'k_r',                          0.2,    3.0,   1.15,   true;   % Radial expansion factor
        'fill_factor',                  0.8,    0.99,  0.95,   true;   % TE material fill factor
        'insulation_width_ratio',       0.02,   0.1,   0.04,   true;   % Insulation width ratio
        
        % === INTERCONNECTS ===
        'interconnect_ratio',           0.1,    0.35,  0.15,   true;   % Inner interconnect width ratio
        'outerconnect_ratio',           0.1,    0.35,  0.15,   true;   % Outer interconnect width ratio
        'interconnect_angle_ratio',     0.1,    0.4,   0.16,   true;   % Inner interconnect angle ratio
        'outerconnect_angle_ratio',     0.1,    0.4,   0.16,   true;   % Outer interconnect angle ratio
        'interconnect_thickness_ratio', 0.5,    2.0,   1.0,    true;   % Inner interconnect thickness ratio
        'outerconnect_thickness_ratio', 0.5,    2.0,   1.0,    true;   % Outer interconnect thickness ratio
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
