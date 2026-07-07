%% REPLAY 1D VALIDATION SWEEPS USING CACHED COMSOL RESULTS
% Recomputes CTM with the current implementation and reuses COMSOL_Tmax from
% existing CSV files. This avoids re-running COMSOL.

clear; clc; close all;

%% ============ PATH SETUP ============
root_dir = fullfile(fileparts(mfilename('fullpath')), '..', '..');
cd(root_dir);
addpath(genpath('src'));
addpath(genpath('scripts'));

%% ============ INPUT/OUTPUT ============
% Folder containing old 1D sweep CSVs.
INPUT_DIR = fullfile('output', 'validity_sweeps', '2026-03-28_14-17-50');

timestamp = datestr(now, 'yyyy-mm-dd_HH-MM-SS');
OUTPUT_DIR = fullfile(INPUT_DIR, ['replay_new_ctm_', timestamp]);
if ~exist(OUTPUT_DIR, 'dir'), mkdir(OUTPUT_DIR); end

% Baseline constants used by original sweep script.
base = struct();
base.f_L = 1.15;
base.t_chip_um = 50;
base.t_SOI_um = 20;
base.R_cyl_um = 1000;
base.w_is_um = 50;
base.I_A = 0.10;
base.T_water = 293.15;
base.fill_factor = 0.9;
base.ic_w_r = 0.1;
base.ic_t_r = 0.6;
base.ic_angle_r = 0.5;
base.oc_w_r = 0.1;
base.oc_t_r = 0.6;
base.oc_angle_r = 0.5;

files = dir(fullfile(INPUT_DIR, '*_sweep_results.csv'));
if isempty(files)
    error('No 1D sweep result CSV files found in: %s', INPUT_DIR);
end

all_results = table();

for f = 1:numel(files)
    in_csv = fullfile(files(f).folder, files(f).name);
    T = readtable(in_csv);

    fprintf('\nProcessing 1D replay for %s\n', files(f).name);

    new_ctm = nan(height(T), 1);
    new_abs = nan(height(T), 1);
    new_pct = nan(height(T), 1);

    for i = 1:height(T)
        q_now = get_col_value(T, i, 'q_Wm2', NaN);
        M_now = get_col_value(T, i, 'M', NaN);
        t_now = get_col_value(T, i, 't_TEC_um', NaN);

        if isnan(q_now) || isnan(M_now) || isnan(t_now)
            fprintf('  Row %d skipped (missing q/M/t).\n', i);
            continue;
        end

        try
            ctm_tmax = run_ctm_case(q_now, M_now, t_now, base);
            new_ctm(i) = ctm_tmax;

            if ismember('COMSOL_Tmax', T.Properties.VariableNames) && ~isnan(T.COMSOL_Tmax(i))
                new_abs(i) = T.COMSOL_Tmax(i) - ctm_tmax;
                new_pct(i) = 100 * abs(new_abs(i)) / max(abs(T.COMSOL_Tmax(i)), eps);
            end
        catch ME
            fprintf('  Row %d CTM error: %s\n', i, ME.message);
        end
    end

    T.CTM_New_Tmax = new_ctm;
    T.Error_Abs_New = new_abs;
    T.Error_Pct_New = new_pct;

    out_csv = fullfile(OUTPUT_DIR, strrep(files(f).name, '.csv', '_replay.csv'));
    writetable(T, out_csv);

    plot_1d_replay(T, files(f).name, OUTPUT_DIR);

    all_results = [all_results; T];
end

master_csv = fullfile(OUTPUT_DIR, 'compiled_all_sweeps_replay.csv');
writetable(all_results, master_csv);

fprintf('\n=== 1D replay complete ===\n');
fprintf('Saved results: %s\n', master_csv);


function v = get_col_value(T, row_idx, col_name, default_v)
    if nargin < 4
        default_v = NaN;
    end

    if ismember(col_name, T.Properties.VariableNames)
        v = T.(col_name)(row_idx);
    else
        v = default_v;
    end
end


function ctm_tmax = run_ctm_case(q_now, M_now, t_now, base)
    theta_deg = 360 / M_now;

    config = struct();
    config.geometry.N_stages = 3;
    config.geometry.wedge_angle_deg = theta_deg;
    config.geometry.thickness_um = t_now;
    config.geometry.radial_expansion_factor = base.f_L;
    config.geometry.w_chip_um = 10000;
    config.geometry.R_cyl_um = base.R_cyl_um;
    config.geometry.t_chip_um = base.t_chip_um;
    config.geometry.t_ins_um = base.t_SOI_um;
    config.geometry.interconnect_ratio = base.ic_w_r;
    config.geometry.outerconnect_ratio = base.oc_w_r;
    config.geometry.insulation_width_um = base.w_is_um;
    config.geometry.interconnect_angle_ratio = base.ic_angle_r;
    config.geometry.outerconnect_angle_ratio = base.oc_angle_r;
    config.geometry.interconnect_thickness_ratio = base.ic_t_r;
    config.geometry.outerconnect_thickness_ratio = base.oc_t_r;
    config.geometry.fill_factor = base.fill_factor;

    config.operating_conditions.I_current_A = base.I_A;
    config.boundary_conditions.q_flux_W_m2 = q_now;
    config.boundary_conditions.T_water_K = base.T_water;
    config.boundary_conditions.h_conv_W_m2K = 1e6;

    config.materials.Bi2Te3 = struct('k', 1.6, 'rho', 1.15e-5, 'S', 210e-6);
    config.materials.Cu    = struct('k', 400, 'rho', 1.667e-8);
    config.materials.Si    = struct('k', 130, 'rho', 0.01);
    config.materials.AlN   = struct('k', 170, 'rho', 1e10);
    config.materials.SiO2  = struct('k', 1.4, 'rho', 1e14);
    config.materials.Al2O3 = struct('k', 35, 'rho', 1e12);

    materials = MaterialProperties(config);
    geometry = TECGeometry(config);
    network = ThermalNetwork(geometry, materials, config);

    T_state = ones(2 * config.geometry.N_stages + 1, 1) * 300;
    for iter = 1:100
        T_old = T_state;
        [T_state, ~, ~] = network.solve(T_state);
        if max(abs(T_state - T_old)) < 1e-6
            break;
        end
    end

    ctm_tmax = max(T_state) - 273.15;
end


function plot_1d_replay(T, source_name, out_dir)
    if ~ismember('Value', T.Properties.VariableNames)
        return;
    end

    x = T.Value;
    [x_sorted, idx] = sort(x);

    have_comsol = ismember('COMSOL_Tmax', T.Properties.VariableNames);
    have_new = ismember('CTM_New_Tmax', T.Properties.VariableNames);

    label_name = source_name;
    label_name = strrep(label_name, '_sweep_results.csv', '');
    label_name = strrep(label_name, '_', ' ');

    if have_comsol && have_new
        fig1 = figure('Visible', 'off');
        plot(x_sorted, T.COMSOL_Tmax(idx), 'k-o', 'LineWidth', 1.2, 'MarkerSize', 5); hold on;
        plot(x_sorted, T.CTM_New_Tmax(idx), 'r-s', 'LineWidth', 1.2, 'MarkerSize', 5);
        grid on;
        xlabel('Sweep value');
        ylabel('T_{max} (C)');
        title(sprintf('Replay: COMSOL vs New CTM (%s)', label_name));
        legend({'COMSOL', 'New CTM'}, 'Location', 'best');
        saveas(fig1, fullfile(out_dir, [strrep(source_name, '.csv', ''), '_compare.png']));
        savefig(fig1, fullfile(out_dir, [strrep(source_name, '.csv', ''), '_compare.fig']));
        close(fig1);

        fig2 = figure('Visible', 'off');
        plot(x_sorted, T.Error_Pct_New(idx), 'b-d', 'LineWidth', 1.2, 'MarkerSize', 5);
        grid on;
        xlabel('\textsf{Number of Wedges} ($N$)', 'Interpreter', 'latex');
        ylabel('Relative Error (%)');
        title(sprintf('Relative Error (%%):  CTM vs COMSOL'));
        saveas(fig2, fullfile(out_dir, [strrep(source_name, '.csv', ''), '_error.png']));
        savefig(fig2, fullfile(out_dir, [strrep(source_name, '.csv', ''), '_error.fig']));
        close(fig2);
    end
end
