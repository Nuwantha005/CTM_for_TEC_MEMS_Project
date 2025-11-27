%% FIND FEASIBLE OPERATING POINT
% This script finds what heat flux levels ARE feasible with current TEC design
% and what temperature targets can be achieved.
%
% Key insight: At 100 kW/m² with 10mm chip, total heat = 10W
%              Divided among 12 wedges (30°) = 833 mW per wedge
%              This is VERY high for a micro-TEC!

clear; clc; close all;

%% Add paths
% Resolve config_path and add src to path robustly
config_path = get_config_path();
src_dir = fileparts(fileparts(config_path));
addpath(genpath(src_dir));

%% Parameters
config = jsondecode(fileread(config_path));

fprintf('\n========================================\n');
fprintf('  FINDING FEASIBLE OPERATING POINT\n');
fprintf('========================================\n');

% Fixed parameters
w_chip_m = config.geometry.w_chip_um * 1e-6;
T_water_K = config.boundary_conditions.T_water_K;

fprintf('Chip: %.1f mm x %.1f mm\n', w_chip_m*1000, w_chip_m*1000);
fprintf('Water temp: %.1f K (%.1f°C)\n', T_water_K, T_water_K - 273.15);

%% Test different heat flux levels
fprintf('\n--- Testing Different Heat Flux Levels ---\n');
fprintf('q (kW/m²) | Q_total (mW) | Q/wedge (mW) | T_max (°C) | Feasible?\n');
fprintf('----------|--------------|--------------|------------|----------\n');

q_flux_test = [1000, 5000, 10000, 20000, 50000, 100000, 200000];  % W/m²

for q = q_flux_test
    % Update config
    config.boundary_conditions.q_flux_W_m2 = q;
    
    % Save and create solver
    fid = fopen(config_path, 'w');
    fprintf(fid, '%s', jsonencode(config));
    fclose(fid);
    
    % Calculate heat
    A_chip = w_chip_m^2;
    Q_total = q * A_chip;
    n_wedges = floor(360 / config.geometry.wedge_angle_deg);
    Q_wedge = Q_total / n_wedges * 1000;  % mW
    
    % Try to solve
    try
        mat = MaterialProperties(config);
        geo = TECGeometry(config);
        net = ThermalNetwork(geo, mat, config);
        
        N = geo.N_stages;
        dim = 2*N + 1;
        T_init = ones(dim, 1) * (T_water_K + 30);
        T_current = T_init;
        
        for iter = 1:100
            [T_new, ~, ~] = net.solve(T_current);
            if any(isnan(T_new)) || any(T_new < 0)
                break;
            end
            T_current = 0.5 * T_new + 0.5 * T_current;
            if max(abs(T_new - T_current)) < 1e-5
                break;
            end
        end
        
        T_max_C = max(T_current) - 273.15;
        T_min_C = min(T_current) - 273.15;
        
        if T_min_C < T_water_K - 273.15 - 5
            status = 'INVALID';
        elseif T_max_C < 100
            status = 'YES';
        elseif T_max_C < 150
            status = 'MARGINAL';
        else
            status = 'NO';
        end
        
        fprintf('%9.0f | %12.1f | %12.1f | %10.1f | %s\n', ...
            q/1000, Q_total*1000, Q_wedge, T_max_C, status);
    catch ME
        fprintf('%9.0f | %12.1f | %12.1f | %10s | ERROR\n', ...
            q/1000, Q_total*1000, Q_wedge, '-');
    end
end

%% Now find the maximum feasible heat flux for T < 100°C
fprintf('\n--- Binary Search for Max Feasible Heat Flux ---\n');

q_low = 1000;
q_high = 100000;
T_target = 100;  % °C

for iter = 1:20
    q_mid = (q_low + q_high) / 2;
    
    config.boundary_conditions.q_flux_W_m2 = q_mid;
    
    try
        mat = MaterialProperties(config);
        geo = TECGeometry(config);
        net = ThermalNetwork(geo, mat, config);
        
        N = geo.N_stages;
        dim = 2*N + 1;
        T_current = ones(dim, 1) * (T_water_K + 30);
        
        for solve_iter = 1:100
            [T_new, ~, ~] = net.solve(T_current);
            if any(isnan(T_new)) || any(T_new < 0)
                T_current = inf * ones(dim, 1);
                break;
            end
            T_current = 0.5 * T_new + 0.5 * T_current;
            if max(abs(T_new - T_current)) < 1e-5
                break;
            end
        end
        
        T_max_C = max(T_current) - 273.15;
        
        if T_max_C < T_target
            q_low = q_mid;
        else
            q_high = q_mid;
        end
        
        fprintf('  q = %.1f kW/m² -> T_max = %.1f°C\n', q_mid/1000, T_max_C);
        
    catch
        q_high = q_mid;
    end
    
    if abs(q_high - q_low) < 100
        break;
    end
end

q_max_feasible = q_low;
fprintf('\n** Maximum feasible heat flux for T < 100°C: %.1f kW/m² **\n', q_max_feasible/1000);

%% Test with more stages
fprintf('\n--- Effect of Number of Stages ---\n');
fprintf('At q = 50 kW/m²:\n');
fprintf('N_stages | T_max (°C)\n');
fprintf('---------|----------\n');

config.boundary_conditions.q_flux_W_m2 = 50000;

for N_test = 2:8
    config.geometry.N_stages = N_test;
    config.geometry.N_tsv_limit = N_test;
    
    try
        mat = MaterialProperties(config);
        geo = TECGeometry(config);
        net = ThermalNetwork(geo, mat, config);
        
        dim = 2*N_test + 1;
        T_current = ones(dim, 1) * (T_water_K + 30);
        
        for solve_iter = 1:100
            [T_new, ~, ~] = net.solve(T_current);
            if any(isnan(T_new)) || any(T_new < 0)
                break;
            end
            T_current = 0.5 * T_new + 0.5 * T_current;
            if max(abs(T_new - T_current)) < 1e-5
                break;
            end
        end
        
        T_max_C = max(T_current) - 273.15;
        fprintf('%8d | %9.1f\n', N_test, T_max_C);
    catch
        fprintf('%8d | ERROR\n', N_test);
    end
end

%% Recommendations
fprintf('\n========================================\n');
fprintf('  RECOMMENDATIONS\n');
fprintf('========================================\n');
fprintf('\nFor your 100 kW/m² requirement:\n');
fprintf('  - Current TEC design cannot achieve T < 100°C\n');
fprintf('  - Max feasible with current design: ~%.0f kW/m²\n', q_max_feasible/1000);
fprintf('\nOptions:\n');
fprintf('  1. Reduce heat flux to %.0f kW/m² or lower\n', q_max_feasible/1000);
fprintf('  2. Accept higher chip temperature (e.g., 125-150°C)\n');
fprintf('  3. Use smaller wedge angles (more wedges = more cooling area)\n');
fprintf('  4. Combine TEC with better heat spreading (larger base)\n');
fprintf('  5. Use higher ZT materials (not Bi2Te3)\n');

% Save config back to reasonable value
config.boundary_conditions.q_flux_W_m2 = 10000;
config.geometry.N_stages = 3;
fid = fopen(config_path, 'w');
fprintf(fid, '%s', jsonencode(config));
fclose(fid);
