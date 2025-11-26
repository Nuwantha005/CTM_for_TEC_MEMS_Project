% test_various_conditions.m
% Test the thermal model at various conditions to find feasible operating points

clear; clc;
addpath(genpath('src'));

fprintf('=== TESTING VARIOUS CONDITIONS ===\n\n');

%% Base config
config = jsondecode(fileread('src/config/default_params.json'));
materials = MaterialProperties(config);

T_target = 373.15;  % 100°C in K
T_water = config.boundary_conditions.T_water_K;

%% Test matrix of conditions
heat_fluxes = [100, 500, 1000, 2000, 5000, 10000, 50000, 100000];  % W/m²
currents = [0, 0.01, 0.02, 0.05, 0.1, 0.15, 0.2];  % A

fprintf('Heat Flux (W/m²) | Current (mA) | T_chip (°C) | Status\n');
fprintf('-----------------|--------------|-------------|-------\n');

results = [];

for q = heat_fluxes
    for I = currents
        config.boundary_conditions.q_flux_W_m2 = q;
        config.operating_conditions.I_current_A = I;
        
        try
            geometry = TECGeometry(config);
            network = ThermalNetwork(geometry, materials, config);
            
            N = geometry.N_stages;
            dim = 2*N + 1;
            T_init = ones(dim, 1) * 300;
            
            T_current = T_init;
            converged = false;
            
            for iter = 1:50
                T_old = T_current;
                [T_new, ~, ~] = network.solve(T_current);
                T_current = T_new;
                
                if max(abs(T_new - T_old)) < 1e-4
                    converged = true;
                    break;
                end
            end
            
            T_chip = T_current(1);
            T_chip_C = T_chip - 273.15;
            
            if T_chip_C < 100
                status = 'GOOD';
            elseif T_chip_C < 150
                status = 'OK';
            else
                status = 'HIGH';
            end
            
            % Check for unphysical solutions
            if T_chip < T_water
                status = 'UNPHYS';
                T_chip_C = -999;
            end
            
            fprintf('%16.0f | %12.1f | %11.1f | %s\n', q, I*1000, T_chip_C, status);
            
            results(end+1, :) = [q, I*1000, T_chip_C];
            
        catch ME
            fprintf('%16.0f | %12.1f | %11s | ERROR\n', q, I*1000, 'ERR');
        end
    end
    fprintf('\n');
end

%% Find best operating points
fprintf('\n=== BEST OPERATING POINTS ===\n');

% For each heat flux, find best current
unique_q = unique(results(:,1));
for q = unique_q'
    idx = results(:,1) == q;
    subset = results(idx, :);
    
    % Find minimum temperature
    [T_min, min_idx] = min(subset(:,3));
    I_best = subset(min_idx, 2);
    
    if T_min < 100 && T_min > 0
        status = 'FEASIBLE';
    elseif T_min < 150 && T_min > 0
        status = 'MARGINAL';
    else
        status = 'INFEASIBLE';
    end
    
    fprintf('q = %6.0f W/m²: Best I = %5.1f mA, T_chip = %6.1f°C [%s]\n', ...
        q, I_best, T_min, status);
end

%% Detail on feasible points
fprintf('\n=== DETAILED FEASIBLE DESIGNS ===\n');

feasible_designs = results(results(:,3) > 0 & results(:,3) < 100, :);
if isempty(feasible_designs)
    fprintf('No fully feasible designs found (T < 100°C)!\n\n');
    
    % Show marginal designs
    marginal = results(results(:,3) > 0 & results(:,3) < 150, :);
    if ~isempty(marginal)
        fprintf('Marginal designs (T < 150°C):\n');
        for i = 1:size(marginal, 1)
            fprintf('  q = %.0f W/m², I = %.1f mA -> T = %.1f°C\n', ...
                marginal(i,1), marginal(i,2), marginal(i,3));
        end
    end
else
    fprintf('Found %d feasible designs:\n', size(feasible_designs,1));
    for i = 1:size(feasible_designs, 1)
        fprintf('  q = %.0f W/m², I = %.1f mA -> T = %.1f°C\n', ...
            feasible_designs(i,1), feasible_designs(i,2), feasible_designs(i,3));
    end
end

%% Calculate max feasible heat flux for specific current
fprintf('\n=== MAX FEASIBLE HEAT FLUX BY CURRENT ===\n');

for I = [0, 0.02, 0.05, 0.1]
    % Binary search for max heat flux that gives T < 100C
    q_low = 100;
    q_high = 100000;
    
    for iter = 1:20
        q_test = (q_low + q_high) / 2;
        
        config.boundary_conditions.q_flux_W_m2 = q_test;
        config.operating_conditions.I_current_A = I;
        
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
        
        if T_chip < 100 && T_chip > T_water - 273.15
            q_low = q_test;
        else
            q_high = q_test;
        end
    end
    
    fprintf('I = %5.1f mA: Max q = %.0f W/m² (%.2f W total chip)\n', ...
        I*1000, q_low, q_low * 100e-6);  % 100 mm² = 100e-6 m²
end
