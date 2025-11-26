% explore_design_improvements.m
% Explore how design changes affect maximum heat flux capacity

clear; clc;
addpath(genpath('src'));

fprintf('=== EXPLORING DESIGN IMPROVEMENTS ===\n\n');

%% Base config
config_base = jsondecode(fileread('src/config/default_params.json'));
materials = MaterialProperties(config_base);

%% Helper function to find max feasible heat flux
function [q_max, T_at_max] = find_max_q(config, materials, target_T_C)
    q_low = 100;
    q_high = 100000;
    
    for iter = 1:20
        q_test = (q_low + q_high) / 2;
        config.boundary_conditions.q_flux_W_m2 = q_test;
        
        try
            geometry = TECGeometry(config);
            network = ThermalNetwork(geometry, materials, config);
            
            N = geometry.N_stages;
            dim = 2*N + 1;
            T_current = ones(dim, 1) * 300;
            
            for j = 1:50
                T_old = T_current;
                [T_new, ~, ~] = network.solve(T_current);
                T_current = T_new;
                if max(abs(T_new - T_old)) < 1e-4
                    break;
                end
            end
            
            T_chip = T_current(1) - 273.15;
            
            if T_chip > 0 && T_chip < target_T_C
                q_low = q_test;
            else
                q_high = q_test;
            end
        catch
            q_high = q_test;
        end
    end
    
    q_max = q_low;
    T_at_max = T_chip;
end

%% 1. Effect of Number of Stages
fprintf('=== 1. EFFECT OF NUMBER OF STAGES ===\n');
fprintf('N_stages | Max q (W/m²) | Total Power (W)\n');
fprintf('---------|--------------|----------------\n');

for N = 1:6
    config = config_base;
    config.geometry.N_stages = N;
    config.operating_conditions.I_current_A = 0.02;
    
    [q_max, ~] = find_max_q(config, materials, 100);
    P_total = q_max * 100e-6;  % 10mm x 10mm = 100 mm²
    
    fprintf('%8d | %12.0f | %14.4f\n', N, q_max, P_total);
end

%% 2. Effect of Wedge Angle (more wedges = more cooling area)
fprintf('\n=== 2. EFFECT OF WEDGE ANGLE ===\n');
fprintf('Angle (deg) | N_wedges | Max q (W/m²) | Total Power (W)\n');
fprintf('------------|----------|--------------|----------------\n');

for angle = [45, 30, 20, 15, 10, 5]
    config = config_base;
    config.geometry.wedge_angle_deg = angle;
    config.operating_conditions.I_current_A = 0.02;
    
    n_wedges = round(360 / angle);
    [q_max, ~] = find_max_q(config, materials, 100);
    P_total = q_max * 100e-6;
    
    fprintf('%11d | %8d | %12.0f | %14.4f\n', angle, n_wedges, q_max, P_total);
end

%% 3. Effect of TEC Thickness
fprintf('\n=== 3. EFFECT OF TEC THICKNESS ===\n');
fprintf('Thickness (um) | Max q (W/m²) | Total Power (W)\n');
fprintf('---------------|--------------|----------------\n');

for t = [10, 20, 50, 100, 200, 500]
    config = config_base;
    config.geometry.thickness_um = t;
    config.operating_conditions.I_current_A = 0.02;
    
    [q_max, ~] = find_max_q(config, materials, 100);
    P_total = q_max * 100e-6;
    
    fprintf('%14d | %12.0f | %14.4f\n', t, q_max, P_total);
end

%% 4. Effect of Fill Factor
fprintf('\n=== 4. EFFECT OF FILL FACTOR ===\n');
fprintf('Fill Factor | Max q (W/m²) | Total Power (W)\n');
fprintf('------------|--------------|----------------\n');

for ff = [0.5, 0.7, 0.8, 0.9, 0.95, 0.99]
    config = config_base;
    config.geometry.fill_factor = ff;
    config.operating_conditions.I_current_A = 0.02;
    
    [q_max, ~] = find_max_q(config, materials, 100);
    P_total = q_max * 100e-6;
    
    fprintf('%11.2f | %12.0f | %14.4f\n', ff, q_max, P_total);
end

%% 5. Combined optimizations
fprintf('\n=== 5. COMBINED OPTIMIZATIONS ===\n');

% Baseline
config = config_base;
config.operating_conditions.I_current_A = 0.02;
[q_baseline, ~] = find_max_q(config, materials, 100);
fprintf('Baseline: q_max = %.0f W/m², P = %.4f W\n', q_baseline, q_baseline * 100e-6);

% Best single modifications
config = config_base;
config.geometry.wedge_angle_deg = 10;  % More wedges
config.operating_conditions.I_current_A = 0.02;
[q_angle, ~] = find_max_q(config, materials, 100);
fprintf('Small angle (10°): q_max = %.0f W/m², P = %.4f W\n', q_angle, q_angle * 100e-6);

config = config_base;
config.geometry.thickness_um = 20;  % Thinner TEC
config.operating_conditions.I_current_A = 0.02;
[q_thin, ~] = find_max_q(config, materials, 100);
fprintf('Thin TEC (20um): q_max = %.0f W/m², P = %.4f W\n', q_thin, q_thin * 100e-6);

% Combined
config = config_base;
config.geometry.wedge_angle_deg = 10;
config.geometry.thickness_um = 20;
config.geometry.fill_factor = 0.95;
config.geometry.N_stages = 2;
config.operating_conditions.I_current_A = 0.02;
[q_combined, ~] = find_max_q(config, materials, 100);
fprintf('Combined optimizations: q_max = %.0f W/m², P = %.4f W\n', q_combined, q_combined * 100e-6);

%% 6. Aggressive optimization for high heat flux
fprintf('\n=== 6. PUSHING FOR HIGHER HEAT FLUX ===\n');

config = config_base;
config.geometry.wedge_angle_deg = 5;   % 72 wedges
config.geometry.thickness_um = 10;     % Very thin
config.geometry.fill_factor = 0.98;
config.geometry.N_stages = 2;
config.geometry.R_cyl_um = 500;        % Smaller center
config.operating_conditions.I_current_A = 0.015;

geometry = TECGeometry(config);
fprintf('Configuration:\n');
fprintf('  Wedge angle: %.0f° (%d wedges)\n', config.geometry.wedge_angle_deg, round(360/config.geometry.wedge_angle_deg));
fprintf('  TEC thickness: %.0f um\n', config.geometry.thickness_um);
fprintf('  Fill factor: %.2f\n', config.geometry.fill_factor);
fprintf('  N_stages: %d\n', config.geometry.N_stages);
fprintf('  R_cyl: %.0f um\n', config.geometry.R_cyl_um);

[q_aggressive, T_chip] = find_max_q(config, materials, 100);
fprintf('Result: q_max = %.0f W/m², P = %.4f W, T_chip = %.1f°C\n', ...
    q_aggressive, q_aggressive * 100e-6, T_chip);

%% Summary
fprintf('\n=== SUMMARY ===\n');
fprintf('To achieve 100 kW/m² (10W chip power), need ~%.0fx improvement\n', 100000/q_aggressive);
fprintf('Current best: %.0f W/m² (%.3f W)\n', q_aggressive, q_aggressive * 100e-6);
fprintf('\nRecommendations:\n');
fprintf('1. Use smaller wedge angles (more parallel cooling paths)\n');
fprintf('2. Use thinner TEC layers (lower thermal resistance)\n');
fprintf('3. Consider heat spreading to increase effective chip area\n');
fprintf('4. For 100 kW/m², explore alternative cooling (forced convection, microchannels)\n');
