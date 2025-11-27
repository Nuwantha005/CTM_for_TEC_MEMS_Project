% find_enhanced_max_heat_flux.m
% Enhanced optimization with:
%   1. Extended TEC thickness (up to 5mm)
%   2. TSV stages as optimization parameter
%   3. Larger chip radius options
%   4. Improved thermal coupling analysis
%
% Target: Find designs that can handle 10-50 W/cm² heat flux

clear; clc;
addpath(genpath('../../src'));

fprintf('╔════════════════════════════════════════════════════════════════╗\n');
fprintf('║     ENHANCED MAXIMUM HEAT FLUX OPTIMIZATION                   ║\n');
fprintf('║     Targeting 10-50 W/cm² for real component cooling          ║\n');
fprintf('╚════════════════════════════════════════════════════════════════╝\n\n');

%% Real component requirements (from research)
fprintf('REAL COMPONENT HEAT FLUX REQUIREMENTS:\n');
fprintf('─────────────────────────────────────────────────────────────────\n');
fprintf('  DRAM (base):           1-5 W/cm²    (10,000-50,000 W/m²)\n');
fprintf('  SRAM Cache:            10 W/cm²     (100,000 W/m²)\n');
fprintf('  IO Controllers:        20 W/cm²     (200,000 W/m²)\n');
fprintf('  Media Engines:         20 W/cm²     (200,000 W/m²)\n');
fprintf('  NPU (average):         30 W/cm²     (300,000 W/m²)\n');
fprintf('  NPU (MAC hotspots):    50 W/cm²     (500,000 W/m²)\n');
fprintf('─────────────────────────────────────────────────────────────────\n\n');

%% Configuration
T_CHIP_MAX_C = 80;
T_CHIP_MAX_K = T_CHIP_MAX_C + 273.15;

timestamp = datestr(now, 'yyyy-mm-dd_HH-MM-SS');
OUTPUT_DIR = fullfile('../../output', 'enhanced_max_heat_flux', timestamp);
if ~exist(OUTPUT_DIR, 'dir')
    mkdir(OUTPUT_DIR);
end

%% Extended parameter space
fprintf('EXTENDED PARAMETER SPACE:\n');
fprintf('─────────────────────────────────────────────────────────────────\n');

% Discrete parameters
N_stages_options = [2, 3, 4, 5];
N_TSV_stages_options = [1, 2, 3, 4, 5];  % NEW: TSV penetration depth
wedge_angle_options = [15, 30, 45];
R_cyl_um_options = [1000, 2000, 3000, 5000];  % NEW: Larger chip radii

% Continuous parameter ranges - EXTENDED
I_range = [0.01, 1.0];           % Amps (10mA - 1A)
t_range = [100, 5000];           % µm (up to 5mm!)
k_r_range = [0.8, 2.0];
ff_range = [0.5, 0.99];

fprintf('  N_stages:       2 - 5\n');
fprintf('  N_TSV_stages:   1 - 5 (TSV penetration depth)\n');
fprintf('  wedge_angle:    15°, 30°, 45°\n');
fprintf('  R_cyl:          1, 2, 3, 5 mm (chip radius)\n');
fprintf('  I_current:      10 mA - 1 A\n');
fprintf('  thickness:      100 µm - 5 mm\n');
fprintf('  k_r:            0.8 - 2.0\n');
fprintf('  fill_factor:    0.5 - 0.99\n');
fprintf('─────────────────────────────────────────────────────────────────\n');
fprintf('Target: T_chip ≤ %.0f°C\n\n', T_CHIP_MAX_C);

%% Results storage
all_results = [];

%% Exhaustive search over discrete parameters
total_discrete = length(N_stages_options) * length(N_TSV_stages_options) * ...
                 length(wedge_angle_options) * length(R_cyl_um_options);
fprintf('Total discrete combinations to test: %d\n\n', total_discrete);

combo_count = 0;
valid_count = 0;

for N_stages = N_stages_options
    for N_TSV = N_TSV_stages_options
        % TSV stages can't exceed total stages
        if N_TSV > N_stages
            continue;
        end
        
        for wedge_angle = wedge_angle_options
            for R_cyl_um = R_cyl_um_options
                combo_count = combo_count + 1;
                
                fprintf('Testing %d/%d: N=%d, N_TSV=%d, θ=%d°, R=%.0fmm ... ', ...
                    combo_count, total_discrete, N_stages, N_TSV, wedge_angle, R_cyl_um/1000);
                
                % Create config for this discrete combination
                base_config = create_enhanced_config(N_stages, N_TSV, wedge_angle, R_cyl_um);
                
                % Optimize continuous parameters
                [q_max, x_opt, T_opt, success] = optimize_enhanced(...
                    base_config, I_range, t_range, k_r_range, ff_range, T_CHIP_MAX_K);
                
                if success && q_max > 1000
                    valid_count = valid_count + 1;
                    
                    result = struct();
                    result.N_stages = N_stages;
                    result.N_TSV_stages = N_TSV;
                    result.wedge_angle_deg = wedge_angle;
                    result.R_cyl_mm = R_cyl_um / 1000;
                    result.I_mA = x_opt(1) * 1000;
                    result.thickness_um = x_opt(2);
                    result.k_r = x_opt(3);
                    result.fill_factor = x_opt(4);
                    result.q_max_Wm2 = q_max;
                    result.q_max_Wcm2 = q_max / 10000;
                    result.T_chip_C = T_opt - 273.15;
                    
                    % Derived metrics
                    A_chip = pi * (R_cyl_um * 1e-6)^2;
                    result.Q_total_mW = q_max * A_chip * 1000;
                    result.chip_area_mm2 = A_chip * 1e6;
                    
                    all_results = [all_results; result];
                    
                    fprintf('✓ q_max = %.0f W/m² (%.1f W/cm²)\n', q_max, q_max/10000);
                else
                    fprintf('✗\n');
                end
            end
        end
    end
end

fprintf('\n');
fprintf('Valid configurations found: %d / %d\n\n', valid_count, combo_count);

%% Sort by maximum heat flux
if ~isempty(all_results)
    [~, sort_idx] = sort([all_results.q_max_Wm2], 'descend');
    all_results = all_results(sort_idx);
end

%% Display top results
fprintf('╔════════════════════════════════════════════════════════════════╗\n');
fprintf('║              TOP DESIGN CANDIDATES                            ║\n');
fprintf('║         (Ranked by Maximum Heat Flux Capability)              ║\n');
fprintf('╚════════════════════════════════════════════════════════════════╝\n\n');

fprintf('┌──────┬────────┬───────┬───────┬────────┬────────┬────────┬────────────┬────────────┬─────────┐\n');
fprintf('│ Rank │ Stages │ N_TSV │ θ(°)  │ R(mm)  │ I(mA)  │ t(µm)  │ q_max      │ q_max      │ T_chip  │\n');
fprintf('│      │        │       │       │        │        │        │ (W/m²)     │ (W/cm²)    │ (°C)    │\n');
fprintf('├──────┼────────┼───────┼───────┼────────┼────────┼────────┼────────────┼────────────┼─────────┤\n');

n_display = min(25, length(all_results));
for i = 1:n_display
    r = all_results(i);
    fprintf('│ %4d │ %6d │ %5d │ %5.0f │ %6.1f │ %6.0f │ %6.0f │ %10.0f │ %10.1f │ %7.1f │\n', ...
        i, r.N_stages, r.N_TSV_stages, r.wedge_angle_deg, r.R_cyl_mm, ...
        r.I_mA, r.thickness_um, r.q_max_Wm2, r.q_max_Wcm2, r.T_chip_C);
end
fprintf('└──────┴────────┴───────┴───────┴────────┴────────┴────────┴────────────┴────────────┴─────────┘\n\n');

%% Feasibility analysis
fprintf('╔════════════════════════════════════════════════════════════════╗\n');
fprintf('║              COMPONENT FEASIBILITY ANALYSIS                   ║\n');
fprintf('╚════════════════════════════════════════════════════════════════╝\n\n');

if ~isempty(all_results)
    best = all_results(1);
    q_max = best.q_max_Wcm2;
    
    fprintf('Best design achieves: %.1f W/cm² (%.0f W/m²)\n\n', q_max, best.q_max_Wm2);
    
    components = {
        'DRAM (base)',              1;
        'DRAM (stacked)',           5;
        'SRAM Cache',               10;
        'IO Controllers',           20;
        'Media Engines',            20;
        'NPU (average)',            30;
        'NPU (MAC hotspots)',       50;
    };
    
    fprintf('Component Cooling Feasibility:\n');
    for i = 1:size(components, 1)
        name = components{i, 1};
        q_req = components{i, 2};
        
        if q_max >= q_req * 1.1  % 10% margin
            status = '✓ YES';
            margin = (q_max / q_req - 1) * 100;
            note = sprintf('(%.0f%% margin)', margin);
        elseif q_max >= q_req
            status = '⚠ MARGINAL';
            note = '';
        else
            status = '✗ NO';
            gap = (q_req / q_max - 1) * 100;
            note = sprintf('(need %.0f%% more)', gap);
        end
        
        fprintf('  %-25s %3.0f W/cm² → %s %s\n', name, q_req, status, note);
    end
end

%% Best design summary
if ~isempty(all_results)
    best = all_results(1);
    
    fprintf('\n═══════════════════════════════════════════════════════════════════\n');
    fprintf('                    BEST DESIGN CONFIGURATION\n');
    fprintf('═══════════════════════════════════════════════════════════════════\n');
    fprintf('  Geometry:\n');
    fprintf('    Stages:         %d\n', best.N_stages);
    fprintf('    TSV Stages:     %d (penetration depth)\n', best.N_TSV_stages);
    fprintf('    Wedge Angle:    %d°\n', best.wedge_angle_deg);
    fprintf('    Chip Radius:    %.1f mm (Area = %.2f mm²)\n', best.R_cyl_mm, best.chip_area_mm2);
    fprintf('    TEC Thickness:  %.0f µm (%.2f mm)\n', best.thickness_um, best.thickness_um/1000);
    fprintf('\n');
    fprintf('  Operating Parameters:\n');
    fprintf('    Current:        %.0f mA\n', best.I_mA);
    fprintf('    k_r:            %.3f\n', best.k_r);
    fprintf('    Fill Factor:    %.3f\n', best.fill_factor);
    fprintf('\n');
    fprintf('  Performance:\n');
    fprintf('    Max Heat Flux:  %.0f W/m² (%.1f W/cm²)\n', best.q_max_Wm2, best.q_max_Wcm2);
    fprintf('    Total Power:    %.2f mW (%.3f W)\n', best.Q_total_mW, best.Q_total_mW/1000);
    fprintf('    Chip Temp:      %.1f°C (limit: %d°C)\n', best.T_chip_C, T_CHIP_MAX_C);
    fprintf('═══════════════════════════════════════════════════════════════════\n');
end

%% Save results
results_struct = struct();
results_struct.timestamp = timestamp;
results_struct.T_chip_max_C = T_CHIP_MAX_C;
results_struct.all_candidates = all_results;
if ~isempty(all_results)
    results_struct.best = best;
end

save(fullfile(OUTPUT_DIR, 'enhanced_results.mat'), 'results_struct');

% Save CSV
if ~isempty(all_results)
    T = table([all_results.N_stages]', [all_results.N_TSV_stages]', ...
        [all_results.wedge_angle_deg]', [all_results.R_cyl_mm]', ...
        [all_results.I_mA]', [all_results.thickness_um]', ...
        [all_results.k_r]', [all_results.fill_factor]', ...
        [all_results.q_max_Wm2]', [all_results.q_max_Wcm2]', ...
        [all_results.T_chip_C]', ...
        'VariableNames', {'N_stages', 'N_TSV', 'wedge_deg', 'R_mm', ...
                          'I_mA', 't_um', 'k_r', 'ff', ...
                          'q_max_Wm2', 'q_max_Wcm2', 'T_chip_C'});
    writetable(T, fullfile(OUTPUT_DIR, 'enhanced_candidates.csv'));
end

fprintf('\nResults saved to: %s\n', OUTPUT_DIR);
fprintf('\n✓ Enhanced optimization complete!\n');

%% ==================== HELPER FUNCTIONS ====================

function config = create_enhanced_config(N_stages, N_TSV_stages, wedge_angle, R_cyl_um)
    config = struct();
    config.geometry.N_stages = N_stages;
    config.geometry.N_TSV_stages = N_TSV_stages;  % NEW parameter
    config.geometry.wedge_angle_deg = wedge_angle;
    config.geometry.w_chip_um = 10000;
    config.geometry.R_cyl_um = R_cyl_um;
    config.geometry.t_chip_um = 50;
    config.geometry.interconnect_ratio = 0.15;
    config.geometry.outerconnect_ratio = 0.15;
    config.geometry.insulation_width_ratio = 0.04;
    config.geometry.interconnect_angle_ratio = 0.16;
    config.geometry.outerconnect_angle_ratio = 0.16;
    config.geometry.interconnect_thickness_ratio = 1.0;
    config.geometry.outerconnect_thickness_ratio = 1.0;
    
    % TSV parameters - increased density for better thermal coupling
    config.geometry.tsv.R_TSV_um = 15;     % Slightly larger TSVs
    config.geometry.tsv.P_TSV_um = 30;     % Pitch
    config.geometry.tsv.g_rad_um = 15;
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

function [q_max, x_opt, T_at_max, success] = optimize_enhanced(base_config, I_range, t_range, k_r_range, ff_range, T_max_K)
    % Optimize for maximum heat flux
    
    lb = [I_range(1), t_range(1), k_r_range(1), ff_range(1)];
    ub = [I_range(2), t_range(2), k_r_range(2), ff_range(2)];
    
    best_q = 0;
    best_x = [];
    best_T = T_max_K;
    
    % Multi-start optimization
    I_starts = linspace(I_range(1), I_range(2), 4);
    t_starts = [200, 500, 1000, 2000, 4000];
    
    for I0 = I_starts
        for t0 = t_starts
            if t0 > t_range(2)
                continue;
            end
            
            x0 = [I0, t0, 1.0, 0.85];
            
            try
                q_found = find_max_q_binary(x0, base_config, T_max_K);
                
                if q_found > best_q
                    best_q = q_found;
                    best_x = x0;
                    best_T = get_chip_temp(x0, base_config, q_found);
                end
            catch
                % Skip failed attempts
            end
        end
    end
    
    q_max = best_q;
    x_opt = best_x;
    T_at_max = best_T;
    success = best_q > 1000;
end

function q_max = find_max_q_binary(x, base_config, T_max_K)
    % Binary search to find max q where T_chip = T_max_K
    
    q_low = 1000;
    q_high = 1000000;  % 100 W/cm² upper limit
    
    % Check if q_low already exceeds temperature
    T_low = get_chip_temp(x, base_config, q_low);
    if T_low > T_max_K
        q_max = 0;
        return;
    end
    
    % Check if we can handle q_high
    T_high = get_chip_temp(x, base_config, q_high);
    if T_high <= T_max_K
        q_max = q_high;
        return;
    end
    
    % Binary search
    for iter = 1:30
        q_mid = (q_low + q_high) / 2;
        T_mid = get_chip_temp(x, base_config, q_mid);
        
        if abs(T_mid - T_max_K) < 0.5
            q_max = q_mid;
            return;
        end
        
        if T_mid < T_max_K
            q_low = q_mid;
        else
            q_high = q_mid;
        end
    end
    
    q_max = q_low;
end

function T_chip = get_chip_temp(x, base_config, q_flux)
    config = base_config;
    config.operating_conditions.I_current_A = x(1);
    config.geometry.thickness_um = x(2);
    config.geometry.radial_expansion_factor = x(3);
    config.geometry.fill_factor = x(4);
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
