% find_true_max_heat_flux.m
% Find the TRUE maximum heat flux capability using binary search
% This properly accounts for the thermal physics
%
% Uses binary search to find the exact q where T_chip = 80°C

clear; clc;
addpath(genpath('../../src'));

fprintf('╔════════════════════════════════════════════════════════════════╗\n');
fprintf('║     TRUE MAXIMUM HEAT FLUX CAPABILITY ANALYSIS                ║\n');
fprintf('║     Using binary search to find exact thermal limits          ║\n');
fprintf('╚════════════════════════════════════════════════════════════════════╝\n\n');

%% Configuration
T_CHIP_MAX_C = 80;      % Maximum allowed chip temperature (°C)
T_CHIP_MAX_K = T_CHIP_MAX_C + 273.15;

% Output directory
timestamp = datestr(now, 'yyyy-mm-dd_HH-MM-SS');
OUTPUT_DIR = fullfile('../../output', 'max_heat_flux', timestamp);
if ~exist(OUTPUT_DIR, 'dir')
    mkdir(OUTPUT_DIR);
end

%% Define search space
% Discrete parameters
N_stages_options = [1, 2, 3, 4, 5];
wedge_angle_options = [15, 20, 30, 45, 60];

% Continuous parameter ranges for optimization
I_range = [0.05, 0.40];      % Amps (50-400 mA)
t_range = [100, 500];        % µm
k_r_range = [0.8, 1.5];
ff_range = [0.80, 0.99];

fprintf('Search Parameters:\n');
fprintf('─────────────────────────────────────────────────────────────\n');
fprintf('  N_stages:      1 - 5\n');
fprintf('  wedge_angle:   15°, 20°, 30°, 45°, 60°\n');
fprintf('  I_current:     %.0f - %.0f mA\n', I_range(1)*1000, I_range(2)*1000);
fprintf('  thickness:     %.0f - %.0f µm\n', t_range(1), t_range(2));
fprintf('  k_r:           %.2f - %.2f\n', k_r_range(1), k_r_range(2));
fprintf('  fill_factor:   %.2f - %.2f\n', ff_range(1), ff_range(2));
fprintf('─────────────────────────────────────────────────────────────\n');
fprintf('Target: T_chip ≤ %.0f°C\n\n', T_CHIP_MAX_C);

%% Create base config
base_config = create_base_config(2, 30);

%% Find optimal parameters for each stage/wedge combination
all_results = [];

fprintf('=== SEARCHING FOR MAXIMUM HEAT FLUX ===\n\n');
total_combos = length(N_stages_options) * length(wedge_angle_options);
combo_count = 0;

for N_stages = N_stages_options
    for wedge_angle = wedge_angle_options
        combo_count = combo_count + 1;
        
        fprintf('Configuration %d/%d: N_stages = %d, θ = %d°\n', ...
            combo_count, total_combos, N_stages, wedge_angle);
        
        % Find optimal continuous parameters for max q
        base_config = create_base_config(N_stages, wedge_angle);
        
        % Multi-start optimization
        best_q = 0;
        best_params = [];
        best_T = 0;
        
        % Grid search over starting points
        I_starts = linspace(I_range(1), I_range(2), 4);
        t_starts = linspace(t_range(1), t_range(2), 3);
        
        for I0 = I_starts
            for t0 = t_starts
                x0 = [I0, t0, 1.0, 0.90];
                
                % Optimize for maximum q
                [q_opt, x_opt, T_opt, success] = optimize_for_max_q(...
                    base_config, x0, I_range, t_range, k_r_range, ff_range, T_CHIP_MAX_K);
                
                if success && q_opt > best_q
                    best_q = q_opt;
                    best_params = x_opt;
                    best_T = T_opt;
                end
            end
        end
        
        if best_q > 100
            result = struct();
            result.N_stages = N_stages;
            result.wedge_angle_deg = wedge_angle;
            result.I_mA = best_params(1) * 1000;
            result.thickness_um = best_params(2);
            result.k_r = best_params(3);
            result.fill_factor = best_params(4);
            result.q_max_Wm2 = best_q;
            result.T_chip_C = best_T - 273.15;
            
            % Calculate derived metrics
            R_cyl_um = base_config.geometry.R_cyl_um;
            result.Q_total_mW = best_q * pi * (R_cyl_um * 1e-6)^2 * 1000;
            result.N_wedges = 360 / wedge_angle;
            
            all_results = [all_results; result];
            
            fprintf('  ✓ q_max = %.0f W/m² @ I=%.0fmA, t=%.0fµm, T=%.1f°C\n\n', ...
                best_q, result.I_mA, result.thickness_um, result.T_chip_C);
        else
            fprintf('  ✗ No valid solution found\n\n');
        end
    end
end

%% Sort by maximum heat flux
if ~isempty(all_results)
    [~, sort_idx] = sort([all_results.q_max_Wm2], 'descend');
    all_results = all_results(sort_idx);
end

%% Display Results
fprintf('\n');
fprintf('╔════════════════════════════════════════════════════════════════╗\n');
fprintf('║              TOP DESIGN CANDIDATES                            ║\n');
fprintf('║         (Ranked by Maximum Heat Flux Capability)              ║\n');
fprintf('╚════════════════════════════════════════════════════════════════════╝\n\n');

fprintf('┌──────┬────────┬───────┬────────┬────────┬───────┬───────┬───────────┬──────────┬─────────┐\n');
fprintf('│ Rank │ Stages │ θ(°)  │ I(mA)  │ t(µm)  │  k_r  │  ff   │ q_max     │ Q_total  │ T_chip  │\n');
fprintf('│      │        │       │        │        │       │       │ (W/m²)    │ (mW)     │ (°C)    │\n');
fprintf('├──────┼────────┼───────┼────────┼────────┼───────┼───────┼───────────┼──────────┼─────────┤\n');

n_display = min(20, length(all_results));
for i = 1:n_display
    r = all_results(i);
    fprintf('│ %4d │ %6d │ %5.0f │ %6.0f │ %6.0f │ %5.2f │ %5.2f │ %9.0f │ %8.2f │ %7.1f │\n', ...
        i, r.N_stages, r.wedge_angle_deg, r.I_mA, r.thickness_um, ...
        r.k_r, r.fill_factor, r.q_max_Wm2, r.Q_total_mW, r.T_chip_C);
end
fprintf('└──────┴────────┴───────┴────────┴────────┴───────┴───────┴───────────┴──────────┴─────────┘\n\n');

%% Best design summary
if ~isempty(all_results)
    best = all_results(1);
    
    fprintf('═══════════════════════════════════════════════════════════════════\n');
    fprintf('                    BEST DESIGN (Maximum q)\n');
    fprintf('═══════════════════════════════════════════════════════════════════\n');
    fprintf('  Configuration:\n');
    fprintf('    Stages:        %d\n', best.N_stages);
    fprintf('    Wedge Angle:   %d° (%d wedges)\n', best.wedge_angle_deg, best.N_wedges);
    fprintf('    Current:       %.0f mA\n', best.I_mA);
    fprintf('    TEC Thickness: %.0f µm\n', best.thickness_um);
    fprintf('    k_r:           %.3f\n', best.k_r);
    fprintf('    Fill Factor:   %.3f\n', best.fill_factor);
    fprintf('\n');
    fprintf('  Performance:\n');
    fprintf('    Max Heat Flux: %.0f W/m² (%.2f W/cm²)\n', best.q_max_Wm2, best.q_max_Wm2/10000);
    fprintf('    Total Power:   %.2f mW\n', best.Q_total_mW);
    fprintf('    Chip Temp:     %.1f°C (limit: %d°C)\n', best.T_chip_C, T_CHIP_MAX_C);
    fprintf('═══════════════════════════════════════════════════════════════════\n\n');
end

%% Component Matching
fprintf('╔════════════════════════════════════════════════════════════════╗\n');
fprintf('║           COMPONENT MATCHING ANALYSIS                         ║\n');
fprintf('╚════════════════════════════════════════════════════════════════════╝\n\n');

components = {
    'SRAM Cache (256KB)',       500,   0.5;
    'SRAM Cache (1MB)',         1200,  2.0;
    'SRAM Cache (4MB)',         2500,  8.0;
    'L3 Cache Slice',           2000,  3.0;
    'Low-Power MCU Core',       1000,  0.3;
    'DSP Block',                1500,  0.5;
    'Media Engine (H.264)',     2500,  1.5;
    'Media Engine (H.265)',     4000,  3.0;
    'NPU/AI Accelerator',       3000,  2.0;
    'ARM Cortex-A55',           5000,  1.5;
    'ARM Cortex-A78',           15000, 4.0;
    'GPU Shader Core',          10000, 3.0;
    'WiFi 6 Baseband',          1500,  0.3;
    '5G Modem Baseband',        12000, 4.0;
    'DDR5 PHY',                 8000,  2.0;
    'HBM PHY',                  15000, 5.0;
    'Voltage Regulator',        20000, 3.0;
    'RF Power Amplifier',       25000, 2.0;
    'VCSEL Array',              40000, 2.0;
};

if ~isempty(all_results)
    q_max_available = best.q_max_Wm2;
    
    fprintf('Based on best design (q_max = %.0f W/m²):\n\n', q_max_available);
    
    fprintf('✓ CAN COOL:\n');
    for i = 1:size(components, 1)
        q_typ = components{i, 2};
        if q_typ <= q_max_available * 0.9
            fprintf('  %-25s  %5.0f W/m²\n', components{i, 1}, q_typ);
        end
    end
    
    fprintf('\n⚠ MARGINAL:\n');
    for i = 1:size(components, 1)
        q_typ = components{i, 2};
        if q_typ > q_max_available * 0.9 && q_typ <= q_max_available
            fprintf('  %-25s  %5.0f W/m²\n', components{i, 1}, q_typ);
        end
    end
    
    fprintf('\n✗ CANNOT COOL:\n');
    for i = 1:size(components, 1)
        q_typ = components{i, 2};
        if q_typ > q_max_available
            fprintf('  %-25s  %5.0f W/m²\n', components{i, 1}, q_typ);
        end
    end
end

%% Save results
results_struct = struct();
results_struct.timestamp = timestamp;
results_struct.T_chip_max_C = T_CHIP_MAX_C;
results_struct.all_candidates = all_results;
if ~isempty(all_results)
    results_struct.best = best;
end

save(fullfile(OUTPUT_DIR, 'true_max_heat_flux_results.mat'), 'results_struct');

fprintf('\nResults saved to: %s\n', OUTPUT_DIR);
fprintf('\n✓ Analysis complete!\n');

%% ==================== HELPER FUNCTIONS ====================

function config = create_base_config(N_stages, wedge_angle)
    config = struct();
    config.geometry.N_stages = N_stages;
    config.geometry.wedge_angle_deg = wedge_angle;
    config.geometry.w_chip_um = 10000;
    config.geometry.R_cyl_um = 1000;
    config.geometry.t_chip_um = 50;
    config.geometry.interconnect_ratio = 0.15;
    config.geometry.outerconnect_ratio = 0.15;
    config.geometry.insulation_width_ratio = 0.04;
    config.geometry.interconnect_angle_ratio = 0.16;
    config.geometry.outerconnect_angle_ratio = 0.16;
    config.geometry.interconnect_thickness_ratio = 1.0;
    config.geometry.outerconnect_thickness_ratio = 1.0;
    config.geometry.tsv.R_TSV_um = 10;
    config.geometry.tsv.P_TSV_um = 20;
    config.geometry.tsv.g_rad_um = 10;
    config.geometry.tsv.t_SOI_um = 100;
    
    config.boundary_conditions.T_water_K = 300;
    config.boundary_conditions.h_conv_W_m2K = 1e6;
    
    config.materials.Bi2Te3 = struct('k', 1.2, 'rho', 1e-5, 'S', 0.0002);
    config.materials.Cu = struct('k', 400, 'rho', 1.7e-8);
    config.materials.Si = struct('k', 150, 'rho', 0.01);
    config.materials.AlN = struct('k', 170, 'rho', 1e10);
    config.materials.SiO2 = struct('k', 1.4, 'rho', 1e14);
    config.materials.Al2O3 = struct('k', 30, 'rho', 1e12);
end

function [q_max, x_opt, T_at_max, success] = optimize_for_max_q(base_config, x0, I_range, t_range, k_r_range, ff_range, T_max_K)
    % Find maximum q by searching for parameters that maximize cooling capacity
    
    lb = [I_range(1), t_range(1), k_r_range(1), ff_range(1)];
    ub = [I_range(2), t_range(2), k_r_range(2), ff_range(2)];
    
    % Objective: find params that allow highest q while T_chip <= T_max_K
    % We use binary search within the objective to find q_max for each param set
    
    best_q = 0;
    best_x = x0;
    best_T = T_max_K;
    
    % For each parameter combination, find max q using binary search
    options = optimoptions('fmincon', ...
        'Algorithm', 'sqp', ...
        'Display', 'off', ...
        'MaxIterations', 100, ...
        'MaxFunctionEvaluations', 1000);
    
    % We want to maximize q, which is returned by binary search
    objective = @(x) -find_max_q_binary(x, base_config, T_max_K);
    
    try
        [x_opt, neg_q, exitflag] = fmincon(objective, x0, [], [], [], [], lb, ub, [], options);
        
        if exitflag > 0
            q_max = -neg_q;
            T_at_max = get_chip_temp_at_q(x_opt, base_config, q_max);
            success = true;
        else
            q_max = 0;
            x_opt = x0;
            T_at_max = T_max_K;
            success = false;
        end
    catch
        q_max = 0;
        x_opt = x0;
        T_at_max = T_max_K;
        success = false;
    end
end

function q_max = find_max_q_binary(x, base_config, T_max_K)
    % Binary search to find max q where T_chip = T_max_K
    
    config = base_config;
    config.operating_conditions.I_current_A = x(1);
    config.geometry.thickness_um = x(2);
    config.geometry.radial_expansion_factor = x(3);
    config.geometry.fill_factor = x(4);
    
    q_low = 100;
    q_high = 100000;  % 100 kW/m² upper limit
    
    % Check if q_low already exceeds temperature
    T_low = simulate_temp(config, q_low);
    if T_low > T_max_K
        q_max = 0;
        return;
    end
    
    % Check if we can handle q_high
    T_high = simulate_temp(config, q_high);
    if T_high <= T_max_K
        q_max = q_high;
        return;
    end
    
    % Binary search
    for iter = 1:30
        q_mid = (q_low + q_high) / 2;
        T_mid = simulate_temp(config, q_mid);
        
        if abs(T_mid - T_max_K) < 0.5  % Within 0.5K
            q_max = q_mid;
            return;
        end
        
        if T_mid < T_max_K
            q_low = q_mid;
        else
            q_high = q_mid;
        end
    end
    
    q_max = q_low;  % Conservative estimate
end

function T_chip = simulate_temp(config, q_flux)
    config.boundary_conditions.q_flux_W_m2 = q_flux;
    
    try
        materials = MaterialProperties(config);
        geometry = TECGeometry(config);
        network = ThermalNetwork(geometry, materials, config);
        
        N = geometry.N_stages;
        T = ones(2*N + 1, 1) * 300;
        
        for iter = 1:100
            T_old = T;
            [T, ~, ~] = network.solve(T);
            if max(abs(T - T_old)) < 1e-6
                break;
            end
        end
        
        T_chip = max(T);
    catch
        T_chip = 1e6;
    end
end

function T_chip = get_chip_temp_at_q(x, base_config, q_flux)
    config = base_config;
    config.operating_conditions.I_current_A = x(1);
    config.geometry.thickness_um = x(2);
    config.geometry.radial_expansion_factor = x(3);
    config.geometry.fill_factor = x(4);
    
    T_chip = simulate_temp(config, q_flux);
end
