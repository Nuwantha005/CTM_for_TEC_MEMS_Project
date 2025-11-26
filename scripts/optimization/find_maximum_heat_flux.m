% find_maximum_heat_flux.m
% Find the maximum heat flux the TEC can handle while keeping T_chip < 80°C
%
% This optimizes ALL parameters including number of stages to find
% the design envelope - what the TEC technology CAN do.
%
% Output: List of candidate designs with their maximum heat flux capability

clear; clc;
addpath(genpath('../../src'));

fprintf('╔════════════════════════════════════════════════════════════════╗\n');
fprintf('║     MAXIMUM HEAT FLUX CAPABILITY ANALYSIS                     ║\n');
fprintf('║     Finding what this TEC technology CAN cool                 ║\n');
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

%% Define optimization space
% Now including N_stages as a variable!

fprintf('Optimization Variables:\n');
fprintf('─────────────────────────────────────────────────────────────\n');
fprintf('  N_stages:      1 - 5 (discrete)\n');
fprintf('  I_current:     10 - 250 mA\n');
fprintf('  thickness:     50 - 500 µm\n');
fprintf('  k_r:           0.8 - 1.5\n');
fprintf('  fill_factor:   0.7 - 0.99\n');
fprintf('  wedge_angle:   15 - 60° (discrete: 15, 20, 30, 45, 60)\n');
fprintf('─────────────────────────────────────────────────────────────\n');
fprintf('Constraint: T_chip < %.0f°C\n\n', T_CHIP_MAX_C);

%% Exhaustive search over discrete parameters + optimization for continuous
% This is more reliable than trying to optimize discrete + continuous together

% Discrete parameter options
N_stages_options = [1, 2, 3, 4, 5];
wedge_angle_options = [15, 20, 30, 45, 60];  % degrees

% Results storage
all_results = [];

fprintf('=== SEARCHING DESIGN SPACE ===\n\n');
fprintf('Testing %d stage configurations × %d wedge angles = %d combinations\n\n', ...
    length(N_stages_options), length(wedge_angle_options), ...
    length(N_stages_options) * length(wedge_angle_options));

total_combos = length(N_stages_options) * length(wedge_angle_options);
combo_count = 0;

for N_stages = N_stages_options
    for wedge_angle = wedge_angle_options
        combo_count = combo_count + 1;
        
        fprintf('─────────────────────────────────────────────────────────────\n');
        fprintf('Configuration %d/%d: N_stages = %d, θ = %d°\n', ...
            combo_count, total_combos, N_stages, wedge_angle);
        fprintf('─────────────────────────────────────────────────────────────\n');
        
        % Create base config for this N_stages and wedge_angle
        base_config = create_base_config(N_stages, wedge_angle);
        
        % Find maximum heat flux for this configuration
        [q_max, x_opt, T_achieved, exitflag] = find_max_q_for_config(base_config, T_CHIP_MAX_K);
        
        if exitflag > 0 && q_max > 100  % Valid solution with q > 100 W/m²
            result = struct();
            result.N_stages = N_stages;
            result.wedge_angle_deg = wedge_angle;
            result.I_mA = x_opt(1) * 1000;
            result.thickness_um = x_opt(2);
            result.k_r = x_opt(3);
            result.fill_factor = x_opt(4);
            result.q_max_Wm2 = q_max;
            result.T_chip_C = T_achieved - 273.15;
            result.Q_total_mW = q_max * pi * (base_config.geometry.R_cyl_um * 1e-6)^2 * 1000;
            result.exitflag = exitflag;
            
            % Calculate derived metrics
            result.N_wedges = 360 / wedge_angle;
            result.chip_area_mm2 = pi * (base_config.geometry.R_cyl_um / 1000)^2;
            
            all_results = [all_results; result];
            
            fprintf('  ✓ q_max = %.0f W/m² (%.2f mW total)\n', q_max, result.Q_total_mW);
            fprintf('    I = %.0f mA, t = %.0f µm, k_r = %.2f, ff = %.2f\n', ...
                result.I_mA, result.thickness_um, result.k_r, result.fill_factor);
            fprintf('    T_chip = %.1f°C\n\n', result.T_chip_C);
        else
            fprintf('  ✗ No valid solution found\n\n');
        end
    end
end

%% Sort results by maximum heat flux
if ~isempty(all_results)
    [~, sort_idx] = sort([all_results.q_max_Wm2], 'descend');
    all_results = all_results(sort_idx);
end

%% Display Top Candidates
fprintf('\n');
fprintf('╔════════════════════════════════════════════════════════════════╗\n');
fprintf('║              TOP DESIGN CANDIDATES                            ║\n');
fprintf('║         (Ranked by Maximum Heat Flux Capability)              ║\n');
fprintf('╚════════════════════════════════════════════════════════════════════╝\n\n');

fprintf('┌──────┬────────┬───────┬────────┬────────┬───────┬───────┬───────────┬──────────┬─────────┐\n');
fprintf('│ Rank │ Stages │ θ(°)  │ I(mA)  │ t(µm)  │  k_r  │  ff   │ q_max     │ Q_total  │ T_chip  │\n');
fprintf('│      │        │       │        │        │       │       │ (W/m²)    │ (mW)     │ (°C)    │\n');
fprintf('├──────┼────────┼───────┼────────┼────────┼───────┼───────┼───────────┼──────────┼─────────┤\n');

n_display = min(15, length(all_results));
for i = 1:n_display
    r = all_results(i);
    fprintf('│ %4d │ %6d │ %5.0f │ %6.0f │ %6.0f │ %5.2f │ %5.2f │ %9.0f │ %8.2f │ %7.1f │\n', ...
        i, r.N_stages, r.wedge_angle_deg, r.I_mA, r.thickness_um, ...
        r.k_r, r.fill_factor, r.q_max_Wm2, r.Q_total_mW, r.T_chip_C);
end
fprintf('└──────┴────────┴───────┴────────┴────────┴───────┴───────┴───────────┴──────────┴─────────┘\n\n');

%% Best overall design
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

%% Component Matching - What can we cool?
fprintf('╔════════════════════════════════════════════════════════════════╗\n');
fprintf('║           COMPONENT MATCHING ANALYSIS                         ║\n');
fprintf('║       What can we put on the bottom die?                      ║\n');
fprintf('╚════════════════════════════════════════════════════════════════════╝\n\n');

% Component database (typical heat flux values)
components = {
    'SRAM Cache (256KB)',       500,   0.5,   'Memory';
    'SRAM Cache (512KB)',       800,   1.0,   'Memory';
    'SRAM Cache (1MB)',         1200,  2.0,   'Memory';
    'SRAM Cache (2MB)',         1800,  4.0,   'Memory';
    'L3 Cache Slice',           2000,  3.0,   'Memory';
    'Low-Power MCU Core',       1000,  0.3,   'Logic';
    'DSP Block',                1500,  0.5,   'Logic';
    'Media Engine (H.264)',     2500,  1.5,   'Logic';
    'NPU/AI Accelerator',       3000,  2.0,   'Logic';
    'Bluetooth LE',             300,   0.01,  'RF';
    'WiFi 6 Baseband',          1500,  0.3,   'RF';
    'Audio Codec',              200,   0.05,  'Analog';
    'ADC/DAC Block',            500,   0.1,   'Analog';
    'Power Management IC',      800,   0.2,   'Power';
    'Sensor Hub',               400,   0.05,  'Sensor';
    'Security Engine',          600,   0.1,   'Security';
    'USB Controller',           700,   0.15,  'I/O';
    'PCIe Controller',          1200,  0.25,  'I/O';
};

if ~isempty(all_results)
    q_max_available = best.q_max_Wm2;
    
    fprintf('Based on best design (q_max = %.0f W/m²):\n\n', q_max_available);
    
    fprintf('✓ CAN COOL (with margin):\n');
    fprintf('─────────────────────────────────────────────────────────────\n');
    can_cool = {};
    for i = 1:size(components, 1)
        name = components{i, 1};
        q_typ = components{i, 2};
        power = components{i, 3};
        category = components{i, 4};
        
        margin = (q_max_available - q_typ) / q_max_available * 100;
        
        if q_typ <= q_max_available * 0.9  % 10% safety margin
            fprintf('  %-25s  %5.0f W/m²  (%.1f W)  [%.0f%% margin]\n', ...
                name, q_typ, power, margin);
            can_cool{end+1} = name;
        end
    end
    
    fprintf('\n⚠ MARGINAL (needs validation):\n');
    fprintf('─────────────────────────────────────────────────────────────\n');
    for i = 1:size(components, 1)
        name = components{i, 1};
        q_typ = components{i, 2};
        power = components{i, 3};
        
        if q_typ > q_max_available * 0.9 && q_typ <= q_max_available
            fprintf('  %-25s  %5.0f W/m²  (%.1f W)\n', name, q_typ, power);
        end
    end
    
    fprintf('\n✗ CANNOT COOL (exceeds capability):\n');
    fprintf('─────────────────────────────────────────────────────────────\n');
    for i = 1:size(components, 1)
        name = components{i, 1};
        q_typ = components{i, 2};
        power = components{i, 3};
        
        if q_typ > q_max_available
            fprintf('  %-25s  %5.0f W/m²  (%.1f W)  [%.0f%% over]\n', ...
                name, q_typ, power, (q_typ - q_max_available)/q_max_available * 100);
        end
    end
end

%% Save Results
fprintf('\n\n');

% Save all results
results_struct = struct();
results_struct.timestamp = timestamp;
results_struct.T_chip_max_C = T_CHIP_MAX_C;
results_struct.all_candidates = all_results;
if ~isempty(all_results)
    results_struct.best = best;
end

save(fullfile(OUTPUT_DIR, 'max_heat_flux_results.mat'), 'results_struct');

% Save to CSV
if ~isempty(all_results)
    T = table([all_results.N_stages]', [all_results.wedge_angle_deg]', ...
        [all_results.I_mA]', [all_results.thickness_um]', ...
        [all_results.k_r]', [all_results.fill_factor]', ...
        [all_results.q_max_Wm2]', [all_results.Q_total_mW]', ...
        [all_results.T_chip_C]', ...
        'VariableNames', {'N_stages', 'wedge_deg', 'I_mA', 't_um', ...
                          'k_r', 'fill_factor', 'q_max_Wm2', 'Q_mW', 'T_chip_C'});
    writetable(T, fullfile(OUTPUT_DIR, 'design_candidates.csv'));
end

fprintf('Results saved to: %s\n', OUTPUT_DIR);
fprintf('\n✓ Maximum heat flux analysis complete!\n');

%% Plot results
if length(all_results) > 1
    figure('Position', [100, 100, 1400, 600], 'Name', 'Maximum Heat Flux Analysis');
    
    % q_max by configuration
    subplot(1,3,1);
    stages = [all_results.N_stages];
    q_vals = [all_results.q_max_Wm2];
    scatter(stages, q_vals, 100, [all_results.wedge_angle_deg], 'filled');
    xlabel('Number of Stages');
    ylabel('Maximum Heat Flux (W/m²)');
    title('Max q by Stage Count');
    colorbar; ylabel(colorbar, 'Wedge Angle (°)');
    grid on;
    
    % q_max vs current
    subplot(1,3,2);
    scatter([all_results.I_mA], [all_results.q_max_Wm2], 100, stages, 'filled');
    xlabel('Optimal Current (mA)');
    ylabel('Maximum Heat Flux (W/m²)');
    title('Max q vs Current');
    colorbar; ylabel(colorbar, 'Stages');
    grid on;
    
    % q_max vs thickness
    subplot(1,3,3);
    scatter([all_results.thickness_um], [all_results.q_max_Wm2], 100, stages, 'filled');
    xlabel('TEC Thickness (µm)');
    ylabel('Maximum Heat Flux (W/m²)');
    title('Max q vs Thickness');
    colorbar; ylabel(colorbar, 'Stages');
    grid on;
    
    saveas(gcf, fullfile(OUTPUT_DIR, 'max_heat_flux_analysis.png'));
    saveas(gcf, fullfile(OUTPUT_DIR, 'max_heat_flux_analysis.fig'));
end

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

function [q_max, x_opt, T_achieved, exitflag] = find_max_q_for_config(base_config, T_max_K)
    % Find maximum heat flux for given configuration
    % Optimization variables: [I, thickness, k_r, fill_factor, q]
    
    % Variable bounds
    %       I(A)    t(um)   k_r    ff      q(W/m²)
    lb = [0.01,    50,    0.8,   0.70,   100];
    ub = [0.25,   500,    1.5,   0.99,   50000];
    
    % Initial guess
    x0 = [0.1, 200, 1.0, 0.9, 2000];
    
    % Objective: MAXIMIZE q (minimize -q)
    objective = @(x) -x(5);  % Negative because we minimize
    
    % Constraint: T_chip <= T_max_K
    nonlcon = @(x) tec_temperature_constraint(x, base_config, T_max_K);
    
    % Optimization options
    options = optimoptions('fmincon', ...
        'Algorithm', 'sqp', ...
        'Display', 'off', ...
        'MaxIterations', 200, ...
        'MaxFunctionEvaluations', 2000, ...
        'OptimalityTolerance', 1e-6);
    
    % Try multiple starting points
    best_q = 0;
    best_x = x0;
    best_T = T_max_K;
    best_exitflag = 0;
    
    % Starting points
    I_starts = [0.05, 0.1, 0.15, 0.2];
    t_starts = [100, 200, 300];
    
    for I0 = I_starts
        for t0 = t_starts
            x0_try = [I0, t0, 1.0, 0.9, 2000];
            
            try
                [x, fval, exitflag, ~] = fmincon(objective, x0_try, ...
                    [], [], [], [], lb, ub, nonlcon, options);
                
                if exitflag > 0 && -fval > best_q
                    % Verify the solution
                    [c, ~] = nonlcon(x);
                    if c <= 0
                        best_q = -fval;
                        best_x = x;
                        best_exitflag = exitflag;
                        
                        % Get actual temperature
                        [T_chip, ~] = simulate_tec(x(1:4), x(5), base_config);
                        best_T = T_chip;
                    end
                end
            catch
                % Skip failed attempts
            end
        end
    end
    
    q_max = best_q;
    x_opt = best_x(1:4);
    T_achieved = best_T;
    exitflag = best_exitflag;
end

function [c, ceq] = tec_temperature_constraint(x, base_config, T_max_K)
    % Nonlinear constraint: T_chip <= T_max_K
    
    I = x(1);
    thickness = x(2);
    k_r = x(3);
    fill_factor = x(4);
    q_flux = x(5);
    
    try
        [T_chip, valid] = simulate_tec([I, thickness, k_r, fill_factor], q_flux, base_config);
        
        if valid
            c = T_chip - T_max_K;  % Must be <= 0
        else
            c = 1e6;  % Infeasible
        end
    catch
        c = 1e6;  % Infeasible
    end
    
    ceq = [];  % No equality constraints
end

function [T_chip, valid] = simulate_tec(x_cont, q_flux, base_config)
    % Simulate TEC and return chip temperature
    
    config = base_config;
    config.operating_conditions.I_current_A = x_cont(1);
    config.geometry.thickness_um = x_cont(2);
    config.geometry.radial_expansion_factor = x_cont(3);
    config.geometry.fill_factor = x_cont(4);
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
        valid = (T_chip > 250) && (T_chip < 500);  % Sanity check
        
    catch
        T_chip = 1e6;
        valid = false;
    end
end
