% generate_comsol_design.m
% Generate a feasible design for COMSOL validation
% Based on comprehensive parameter exploration

clear; clc;
addpath(genpath('src'));

fprintf('=== GENERATING COMSOL-READY DESIGN ===\n\n');

%% Best configuration from exploration
% Key insights:
% - Thicker TEC (100-500um) reduces electrical resistance, helps cooling
% - 3-5 stages optimal
% - 30-45° wedge angle optimal (not too many or few wedges)
% - High fill factor helps slightly

config = jsondecode(fileread('src/config/default_params.json'));

% Optimized parameters
config.geometry.N_stages = 4;
config.geometry.wedge_angle_deg = 30;
config.geometry.thickness_um = 200;        % Thicker = lower electrical R
config.geometry.fill_factor = 0.95;
config.geometry.interconnect_ratio = 0.15;
config.geometry.outerconnect_ratio = 0.15;
config.geometry.insulation_width_ratio = 0.04;
config.geometry.radial_expansion_factor = 1.15;

% Optimal current (found from search)
config.operating_conditions.I_current_A = 0.025;

% Water cooling
config.boundary_conditions.h_conv_W_m2K = 1e6;  % Effective for water cooling
config.boundary_conditions.T_water_K = 300;

fprintf('=== OPTIMIZED CONFIGURATION ===\n');
fprintf('Geometry:\n');
fprintf('  Chip: %.0f mm x %.0f mm\n', config.geometry.w_chip_um/1000, config.geometry.w_chip_um/1000);
fprintf('  N_stages: %d\n', config.geometry.N_stages);
fprintf('  Wedge angle: %.0f° (%d wedges)\n', config.geometry.wedge_angle_deg, round(360/config.geometry.wedge_angle_deg));
fprintf('  TEC thickness: %.0f um\n', config.geometry.thickness_um);
fprintf('  Fill factor: %.2f\n', config.geometry.fill_factor);
fprintf('\nOperating conditions:\n');
fprintf('  Current: %.1f mA\n', config.operating_conditions.I_current_A * 1000);
fprintf('  Water temp: %.0f K (%.0f°C)\n', config.boundary_conditions.T_water_K, config.boundary_conditions.T_water_K - 273.15);

%% Find maximum feasible heat flux for T < 100°C
materials = MaterialProperties(config);

q_low = 100;
q_high = 10000;
I_optimal = 0.025;

fprintf('\n=== SEARCHING FOR MAXIMUM FEASIBLE HEAT FLUX ===\n');

% First pass: find approximate range
for I = [0.01, 0.02, 0.025, 0.03, 0.05, 0.1]
    config.operating_conditions.I_current_A = I;
    
    for iter = 1:25
        q_test = (q_low + q_high) / 2;
        config.boundary_conditions.q_flux_W_m2 = q_test;
        
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
        
        T_chip_C = T_current(1) - 273.15;
        
        if T_chip_C > 0 && T_chip_C < 100
            q_low = q_test;
        else
            q_high = q_test;
        end
    end
    
    if q_low > I_optimal * 1000  % Track best result
        I_optimal = I;
    end
    
    fprintf('I = %5.1f mA: Max q = %.0f W/m² -> T = %.1f°C\n', I*1000, q_low, T_chip_C);
    q_high = 10000;  % Reset for next current
end

%% Final optimization with best current
config.operating_conditions.I_current_A = 0.025;  % Best from above
q_low = 100;
q_high = 5000;

for iter = 1:30
    q_test = (q_low + q_high) / 2;
    config.boundary_conditions.q_flux_W_m2 = q_test;
    
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
    
    T_chip_C = T_current(1) - 273.15;
    
    if T_chip_C > 0 && T_chip_C < 100
        q_low = q_test;
    else
        q_high = q_test;
    end
end

q_max_feasible = q_low;
config.boundary_conditions.q_flux_W_m2 = q_max_feasible;

%% Generate final solution at feasible point
geometry = TECGeometry(config);
network = ThermalNetwork(geometry, materials, config);

N = geometry.N_stages;
dim = 2*N + 1;
T_current = ones(dim, 1) * 300;

for iter = 1:100
    T_old = T_current;
    [T_new, Q_out, Q_in] = network.solve(T_current);
    T_current = T_new;
    if max(abs(T_new - T_old)) < 1e-6
        break;
    end
end

fprintf('\n=== FINAL FEASIBLE DESIGN ===\n');
fprintf('Heat flux: %.0f W/m²\n', q_max_feasible);
fprintf('Total chip power: %.2f mW (%.4f W)\n', q_max_feasible * 100e-6 * 1000, q_max_feasible * 100e-6);
fprintf('\nTemperature distribution:\n');
fprintf('  T_chip (center): %.1f K (%.1f°C)\n', T_current(1), T_current(1) - 273.15);
for i = 1:N
    fprintf('  Stage %d: T_Si = %.1f K, T_cold = %.1f K\n', ...
        i, T_current(1+i), T_current(N+1+i));
end
fprintf('  T_water: %.1f K\n', config.boundary_conditions.T_water_K);

fprintf('\nHeat flow:\n');
fprintf('  Q_in (per wedge): %.4f mW\n', Q_in * 1000);
fprintf('  Q_out (per wedge): %.4f mW\n', Q_out * 1000);
n_wedges = round(360 / config.geometry.wedge_angle_deg);
fprintf('  Total Q_in: %.4f mW\n', Q_in * n_wedges * 1000);
fprintf('  Total Q_out: %.4f mW\n', Q_out * n_wedges * 1000);

%% Export COMSOL parameters
fprintf('\n=== COMSOL PARAMETERS ===\n');
fprintf('Save this to a parameter file for COMSOL import:\n\n');

% Print in COMSOL-compatible format
fprintf('// Geometry Parameters\n');
fprintf('w_chip %.6e[m]\n', config.geometry.w_chip_um * 1e-6);
fprintf('N_stages %d\n', config.geometry.N_stages);
fprintf('theta %.6f[deg]\n', config.geometry.wedge_angle_deg);
fprintf('t_tec %.6e[m]\n', config.geometry.thickness_um * 1e-6);
fprintf('t_chip %.6e[m]\n', config.geometry.t_chip_um * 1e-6);
fprintf('R_cyl %.6e[m]\n', config.geometry.R_cyl_um * 1e-6);
fprintf('fill_factor %.4f\n', config.geometry.fill_factor);
fprintf('k_r %.4f\n', config.geometry.radial_expansion_factor);
fprintf('interconnect_ratio %.4f\n', config.geometry.interconnect_ratio);
fprintf('outerconnect_ratio %.4f\n', config.geometry.outerconnect_ratio);
fprintf('insulation_width_ratio %.4f\n', config.geometry.insulation_width_ratio);

fprintf('\n// Operating Conditions\n');
fprintf('I_current %.6e[A]\n', config.operating_conditions.I_current_A);
fprintf('q_flux %.6e[W/m^2]\n', config.boundary_conditions.q_flux_W_m2);
fprintf('T_water %.6f[K]\n', config.boundary_conditions.T_water_K);
fprintf('h_conv %.6e[W/(m^2*K)]\n', config.boundary_conditions.h_conv_W_m2K);

fprintf('\n// Stage Geometry (computed)\n');
for i = 1:N
    [r_in, L, w_ic, t_ic, beta_ic, w_oc, t_oc, beta_oc, w_az, w_is] = geometry.get_stage_geometry(i);
    fprintf('// Stage %d\n', i);
    fprintf('r_in_%d %.6e[m]\n', i, r_in);
    fprintf('L_%d %.6e[m]\n', i, L);
    fprintf('w_ic_%d %.6e[m]\n', i, w_ic);
    fprintf('w_oc_%d %.6e[m]\n', i, w_oc);
    fprintf('w_is_%d %.6e[m]\n', i, w_is);
end

%% Save to file
output_dir = 'output/comsol_designs';
if ~exist(output_dir, 'dir')
    mkdir(output_dir);
end

timestamp = datestr(now, 'yyyy-mm-dd_HH-MM-SS');
filename = fullfile(output_dir, sprintf('comsol_design_%s.txt', timestamp));

fid = fopen(filename, 'w');
fprintf(fid, '%% Radial TEC Design for COMSOL\n');
fprintf(fid, '%% Generated: %s\n', timestamp);
fprintf(fid, '%% Max heat flux for T < 100°C: %.0f W/m²\n\n', q_max_feasible);

fprintf(fid, '%% Geometry Parameters\n');
fprintf(fid, 'w_chip %.6e[m]\n', config.geometry.w_chip_um * 1e-6);
fprintf(fid, 'N_stages %d\n', config.geometry.N_stages);
fprintf(fid, 'theta %.6f[deg]\n', config.geometry.wedge_angle_deg);
fprintf(fid, 't_tec %.6e[m]\n', config.geometry.thickness_um * 1e-6);
fprintf(fid, 't_chip %.6e[m]\n', config.geometry.t_chip_um * 1e-6);
fprintf(fid, 'R_cyl %.6e[m]\n', config.geometry.R_cyl_um * 1e-6);
fprintf(fid, 'fill_factor %.4f\n', config.geometry.fill_factor);
fprintf(fid, 'k_r %.4f\n', config.geometry.radial_expansion_factor);
fprintf(fid, 'interconnect_ratio %.4f\n', config.geometry.interconnect_ratio);
fprintf(fid, 'outerconnect_ratio %.4f\n', config.geometry.outerconnect_ratio);
fprintf(fid, 'insulation_width_ratio %.4f\n', config.geometry.insulation_width_ratio);

fprintf(fid, '\n%% Operating Conditions\n');
fprintf(fid, 'I_current %.6e[A]\n', config.operating_conditions.I_current_A);
fprintf(fid, 'q_flux %.6e[W/m^2]\n', config.boundary_conditions.q_flux_W_m2);
fprintf(fid, 'T_water %.6f[K]\n', config.boundary_conditions.T_water_K);
fprintf(fid, 'h_conv %.6e[W/(m^2*K)]\n', config.boundary_conditions.h_conv_W_m2K);

fprintf(fid, '\n%% Stage Geometry (computed)\n');
for i = 1:N
    [r_in, L, w_ic, t_ic, beta_ic, w_oc, t_oc, beta_oc, w_az, w_is] = geometry.get_stage_geometry(i);
    fprintf(fid, '%% Stage %d\n', i);
    fprintf(fid, 'r_in_%d %.6e[m]\n', i, r_in);
    fprintf(fid, 'L_%d %.6e[m]\n', i, L);
    fprintf(fid, 'w_ic_%d %.6e[m]\n', i, w_ic);
    fprintf(fid, 'w_oc_%d %.6e[m]\n', i, w_oc);
    fprintf(fid, 'w_is_%d %.6e[m]\n', i, w_is);
end

fprintf(fid, '\n%% Expected Results\n');
fprintf(fid, '%% T_chip_max = %.1f °C\n', T_current(1) - 273.15);
fprintf(fid, '%% Total power = %.4f W\n', q_max_feasible * 100e-6);

fclose(fid);
fprintf('\n\nDesign saved to: %s\n', filename);

%% Also save JSON config
json_filename = fullfile(output_dir, sprintf('comsol_config_%s.json', timestamp));
json_str = jsonencode(config);
% Pretty print
json_str = strrep(json_str, ',', sprintf(',\n  '));
json_str = strrep(json_str, '{', sprintf('{\n  '));
json_str = strrep(json_str, '}', sprintf('\n}'));

fid = fopen(json_filename, 'w');
fprintf(fid, '%s', json_str);
fclose(fid);
fprintf('Config saved to: %s\n', json_filename);
