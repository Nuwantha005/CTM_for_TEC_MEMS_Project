%% REPLAY 2D VALIDATION SWEEP USING CACHED COMSOL RESULTS
% Recomputes CTM with the current implementation and reuses COMSOL_Tmax from
% existing 2D sweep CSV. Generates side-by-side contour plots.

clear; clc; close all;

%% ============ PATH SETUP ============
root_dir = fullfile(fileparts(mfilename('fullpath')), '..', '..');
cd(root_dir);
addpath(genpath('src'));
addpath(genpath('scripts'));

%% ============ INPUT/OUTPUT ============
INPUT_CSV = fullfile('output', 'validity_sweeps', '2026-03-28_14-59-39', 'HeatFlux_vs_Thickness_2D_results.csv');

if ~isfile(INPUT_CSV)
    error('Input CSV not found: %s', INPUT_CSV);
end

timestamp = datestr(now, 'yyyy-mm-dd_HH-MM-SS');
output_parent = fileparts(INPUT_CSV);
OUTPUT_DIR = fullfile(output_parent, ['replay_new_ctm_', timestamp]);
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

T = readtable(INPUT_CSV);

v1_col = find_var_col(T.Properties.VariableNames, 'Var1_');
v2_col = find_var_col(T.Properties.VariableNames, 'Var2_');
if isempty(v1_col) || isempty(v2_col)
    error('Could not detect Var1_/Var2_ columns in %s', INPUT_CSV);
end

var1_name = erase(v1_col, 'Var1_');
var2_name = erase(v2_col, 'Var2_');

fprintf('Replaying 2D sweep from cached COMSOL data:\n');
fprintf('  Var1: %s\n', var1_name);
fprintf('  Var2: %s\n', var2_name);

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

out_csv = fullfile(OUTPUT_DIR, 'HeatFlux_vs_Thickness_2D_results_replay.csv');
writetable(T, out_csv);

plot_2d_replay(T, v1_col, v2_col, var1_name, var2_name, OUTPUT_DIR);

fprintf('\n=== 2D replay complete ===\n');
fprintf('Saved results: %s\n', out_csv);


function name = find_var_col(var_names, prefix)
    idx = find(startsWith(var_names, prefix), 1, 'first');
    if isempty(idx)
        name = '';
    else
        name = var_names{idx};
    end
end


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


function plot_2d_replay(T, v1_col, v2_col, var1_name, var2_name, out_dir)
    v1_vals = unique(T.(v1_col));
    v2_vals = unique(T.(v2_col));

    [V1, V2] = meshgrid(v1_vals, v2_vals);
    ctm_grid = nan(size(V1));
    comsol_grid = nan(size(V1));
    err_grid = nan(size(V1));

    for i = 1:height(T)
        x = T.(v1_col)(i);
        y = T.(v2_col)(i);

        ix = find(abs(v1_vals - x) < 1e-12, 1, 'first');
        iy = find(abs(v2_vals - y) < 1e-12, 1, 'first');
        if isempty(ix) || isempty(iy)
            continue;
        end

        ctm_grid(iy, ix) = T.CTM_New_Tmax(i);
        comsol_grid(iy, ix) = T.COMSOL_Tmax(i);
        err_grid(iy, ix) = T.Error_Pct_New(i);
    end

    % Use a shared color scale so differences are visually comparable.
    all_vals = [comsol_grid(:); ctm_grid(:)];
    all_vals = all_vals(~isnan(all_vals));
    if isempty(all_vals)
        clim_shared = [0, 1];
    else
        clim_shared = [min(all_vals), max(all_vals)];
        if abs(clim_shared(2) - clim_shared(1)) < eps
            clim_shared = clim_shared + [-0.5, 0.5];
        end
    end

    fig = figure('Visible', 'off', 'Position', [80, 80, 1400, 520]);

    ax1 = subplot(1, 2, 1);
    contourf(V1, V2, comsol_grid, 20, 'LineStyle', 'none');
    colormap(ax1, hot);
    clim(ax1, clim_shared);
    xlabel(pretty_label(var1_name));
    ylabel(pretty_label(var2_name));
    title('COMSOL T_{max} (C)');

    ax2 = subplot(1, 2, 2);
    contourf(V1, V2, ctm_grid, 20, 'LineStyle', 'none');
    colormap(ax2, hot);
    clim(ax2, clim_shared);
    xlabel(pretty_label(var1_name));
    ylabel(pretty_label(var2_name));
    title('New CTM T_{max} (C)');

    cb = colorbar(ax2, 'eastoutside');
    cb.Label.String = 'T_{max} (C)';
    linkaxes([ax1, ax2], 'xy');

    saveas(fig, fullfile(out_dir, 'Replay_2D_COMSOL_vs_CTM_side_by_side.png'));
    savefig(fig, fullfile(out_dir, 'Replay_2D_COMSOL_vs_CTM_side_by_side.fig'));
    close(fig);

    fig_err = figure('Visible', 'off');
    contourf(V1, V2, err_grid, 20, 'LineStyle', 'none');
    colorbar;
    colormap(gca, jet);
    xlabel(pretty_label(var1_name));
    ylabel(pretty_label(var2_name));
    title('New CTM Relative Error vs COMSOL (%)');
    saveas(fig_err, fullfile(out_dir, 'Replay_2D_Error_percent.png'));
    savefig(fig_err, fullfile(out_dir, 'Replay_2D_Error_percent.fig'));
    close(fig_err);
end


function s = pretty_label(var_name)
    s = var_name;
    s = strrep(s, 'q_Wm2', 'Heat Flux (W/m^2)');
    s = strrep(s, 't_TEC_um', 'TEC Thickness (um)');
    s = strrep(s, 'M', 'M (Number of Wedges)');
    s = strrep(s, '_', ' ');
end
