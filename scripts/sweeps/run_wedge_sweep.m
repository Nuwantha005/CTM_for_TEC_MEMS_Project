% run_wedge_sweep.m
% Sweep over M (number of wedges) to find optimal packing
% Notation follows: Thermal_Network_For_Radial_TEC.tex
%
% INTEGER PARAMETERS:
%   N: Number of stages (fixed)
%   M: Number of wedges (swept in this script)
%
% θ = 360°/M, so M determines the wedge angle

format longG;
clc;
clear all;

% Get the directory where this script is located
script_dir = fileparts(mfilename('fullpath'));
project_root = fullfile(script_dir, '..', '..');

addpath(genpath(fullfile(project_root, 'src')));

% Ensure output directory exists
timestamp = datestr(now, 'yyyy-mm-dd_HH-MM-SS');
output_dir = fullfile(project_root, 'output', 'wedge_sweep', timestamp);
if ~exist(output_dir, 'dir')
    mkdir(output_dir);
end

% Load centralized configuration
[var_names, lb, ub, x0, all_vars, CONFIG] = optimization_variables();

fprintf('══════════════════════════════════════════════════════════════\n');
fprintf('         WEDGE SWEEP (M) - Paper Notation Reference           \n');
fprintf('══════════════════════════════════════════════════════════════\n');
fprintf('=== Settings from optimization_variables.m ===\n');
fprintf('Heat flux (q_flux):     %.0e W/m²\n', CONFIG.q_flux_W_m2);
fprintf('Coolant temp (T_water): %.1f K\n', CONFIG.T_water_K);
fprintf('Target max temp:        %d°C\n', CONFIG.T_target_C);
fprintf('\nInteger Parameters:\n');
fprintf('  N (stages):   %d  (fixed)\n', CONFIG.N);
fprintf('  M (wedges):   will sweep over valid divisors of 360\n');
fprintf('\nValid M values (θ = 360°/M):\n');
fprintf('  ');
for m = CONFIG.M_valid
    fprintf('%d (%.1f°), ', m, 360/m);
end
fprintf('\n══════════════════════════════════════════════════════════════\n\n');

% Wedge values to sweep (subset of M_valid for practical angles)
M_values = [6, 8, 10, 12, 15, 18, 20, 24, 30, 36];

fprintf('Sweeping M = ');
fprintf('%d ', M_values);
fprintf('\n\n');

% Store results
results = struct();
results.M_values = M_values;
results.theta_values = 360 ./ M_values;
results.T_max = zeros(size(M_values));
results.COP = zeros(size(M_values));
results.converged = false(size(M_values));

% Run sweep
for idx = 1:length(M_values)
    M = M_values(idx);
    theta_deg = 360 / M;

    fprintf('═══════════════════════════════════════════════════════════════\n');
    fprintf('M = %d  →  θ = %.1f°\n', M, theta_deg);
    fprintf('═══════════════════════════════════════════════════════════════\n');

    try
        % Create config with this M value
        config = struct();
        config.geometry.N_stages = CONFIG.N;
        config.geometry.M_wedges = M;
        config.geometry.wedge_angle_deg = theta_deg;
        config.geometry.w_chip_um = CONFIG.W_chip_um;
        config.geometry.t_chip_um = CONFIG.t_chip_um;

        % Use initial values from optimization variables
        config.geometry.R_cyl_um = x0(strcmp(var_names, 'r_cyl_um'));
        config.geometry.thickness_um = x0(strcmp(var_names, 't_TEC_um'));
        config.geometry.t_ins_um = x0(strcmp(var_names, 't_ins_um'));
        config.geometry.radial_expansion_factor = x0(strcmp(var_names, 'f_L'));
        config.geometry.azimuthal_gap_um = x0(strcmp(var_names, 'W_az_um'));
        config.geometry.insulation_width_ratio = x0(strcmp(var_names, 'W_is_ratio'));
        config.geometry.interconnect_ratio = x0(strcmp(var_names, 'f_ic_W'));
        config.geometry.interconnect_thickness_ratio = x0(strcmp(var_names, 'f_ic_t'));
        config.geometry.interconnect_angle_ratio = x0(strcmp(var_names, 'f_ic_beta'));
        config.geometry.outerconnect_ratio = x0(strcmp(var_names, 'f_oc_W'));
        config.geometry.outerconnect_thickness_ratio = x0(strcmp(var_names, 'f_oc_t'));
        config.geometry.outerconnect_angle_ratio = x0(strcmp(var_names, 'f_oc_beta'));

        config.boundary_conditions.q_flux_W_m2 = CONFIG.q_flux_W_m2;
        config.boundary_conditions.T_water_K = CONFIG.T_water_K;
        config.boundary_conditions.h_conv_W_m2K = CONFIG.h_conv_W_m2K;

        config.operating_conditions.I_current_A = x0(strcmp(var_names, 'I'));

        config.materials = CONFIG.materials;

        % Run solver
        materials = MaterialProperties(config);
        geometry = TECGeometry(config);
        network = ThermalNetwork(geometry, materials, config);

        N = geometry.N_stages;
        T_water = config.boundary_conditions.T_water_K;
        T = ones(2*N + 1, 1) * (T_water + 50);

        for iter = 1:100
            T_old = T;
            [T_new, Q_out, Q_in] = network.solve(T);
            T = 0.5 * T_new + 0.5 * T;
            if max(abs(T - T_old)) < 1e-6
                break;
            end
        end

        T_max = max(T);
        P_elec = Q_out - Q_in;
        if P_elec > 0
            COP = Q_in / P_elec;
        else
            COP = 0;
        end

        results.T_max(idx) = T_max;
        results.COP(idx) = COP;
        results.converged(idx) = true;

        fprintf('  T_max = %.2f °C, COP = %.3f\n', T_max - 273.15, COP);

    catch ME
        fprintf('  FAILED: %s\n', ME.message);
        results.T_max(idx) = NaN;
        results.COP(idx) = NaN;
        results.converged(idx) = false;
    end
end

%% Plot results
fprintf('\n═══════════════════════════════════════════════════════════════\n');
fprintf('                    WEDGE SWEEP RESULTS                        \n');
fprintf('═══════════════════════════════════════════════════════════════\n');

figure('Position', [100, 100, 1200, 500], 'Name', 'Wedge Sweep Results');

subplot(1, 2, 1);
yyaxis left;
plot(M_values, results.T_max - 273.15, '-ob', 'LineWidth', 2, 'MarkerFaceColor', 'b');
ylabel('T_{max} (°C)');

yyaxis right;
plot(M_values, results.COP, '-sr', 'LineWidth', 2, 'MarkerFaceColor', 'r');
ylabel('COP');

xlabel('M (Number of Wedges)');
title('T_{max} and COP vs M');
grid on;
legend('T_{max}', 'COP', 'Location', 'best');

% Add theta labels on top x-axis
ax1 = gca;
ax1_pos = ax1.Position;
ax2 = axes('Position', ax1_pos, 'XAxisLocation', 'top', 'Color', 'none');
ax2.XLim = ax1.XLim;
ax2.XTick = M_values;
ax2.XTickLabel = arrayfun(@(m) sprintf('%.0f°', 360/m), M_values, 'UniformOutput', false);
ax2.XLabel.String = 'θ (Wedge Angle)';
ax2.YTick = [];

subplot(1, 2, 2);
bar(results.theta_values, results.T_max - 273.15);
xlabel('θ (Wedge Angle, °)');
ylabel('T_{max} (°C)');
title('T_{max} vs Wedge Angle');
grid on;

saveas(gcf, fullfile(output_dir, 'wedge_sweep_results.png'));
saveas(gcf, fullfile(output_dir, 'wedge_sweep_results.fig'));

%% Save results
save(fullfile(output_dir, 'wedge_sweep_results.mat'), 'results');

% Save to CSV
results_table = table(M_values', results.theta_values', results.T_max' - 273.15, results.COP', results.converged', ...
    'VariableNames', {'M', 'theta_deg', 'T_max_C', 'COP', 'Converged'});
writetable(results_table, fullfile(output_dir, 'wedge_sweep_results.csv'));

fprintf('\nResults saved to: %s\n', output_dir);
fprintf('\n✓ Wedge sweep complete!\n');
