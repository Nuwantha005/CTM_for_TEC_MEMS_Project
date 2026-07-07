%% ONE-TIME REGENERATION OF PUBLICATION PLOTS FROM REPLAY CSV DATA
% Rebuilds side-by-side CTM vs COMSOL and error contour plots from existing
% replay CSV(s) without rerunning COMSOL or CTM solves.

clear; clc; close all;
fprintf('Running replot_publication_from_replay_oneoff\n');

%% ============ USER INPUT (ONE-TIME) ============
TARGET_DIR = 'E:\Semester 7\ME4311 - MicroNano Electro Mechanical Systems and Nanotechnology\Project\Preliminary Optimization\Algorithm\output\validity_sweeps\2026-03-28_14-59-39\replay_new_ctm_2026-04-14_09-43-52';
N_CONTOURS = 20;

if ~isfolder(TARGET_DIR)
    error('Target folder not found: %s', TARGET_DIR);
end

csv_files = dir(fullfile(TARGET_DIR, '*_results_replay.csv'));
if isempty(csv_files)
    csv_files = dir(fullfile(TARGET_DIR, '*.csv'));
end

if isempty(csv_files)
    error('No CSV files found in: %s', TARGET_DIR);
end

fprintf('Found %d CSV file(s) in %s\n', numel(csv_files), TARGET_DIR);

for k = 1:numel(csv_files)
    csv_path = fullfile(csv_files(k).folder, csv_files(k).name);
    fprintf('\nProcessing: %s\n', csv_files(k).name);

    T = readtable(csv_path);

    v1_col = find_var_col(T.Properties.VariableNames, 'Var1_');
    v2_col = find_var_col(T.Properties.VariableNames, 'Var2_');
    if isempty(v1_col) || isempty(v2_col)
        warning('Skipping %s: Var1_/Var2_ columns not found.', csv_files(k).name);
        continue;
    end

    var1_name = erase(v1_col, 'Var1_');
    var2_name = erase(v2_col, 'Var2_');

    ctm_col = pick_first_existing(T.Properties.VariableNames, {'CTM_New_Tmax', 'MATLAB_Tmax'});
    comsol_col = pick_first_existing(T.Properties.VariableNames, {'COMSOL_Tmax'});
    err_col = pick_first_existing(T.Properties.VariableNames, {'Error_Pct_New', 'Error_Pct'});

    if isempty(ctm_col) || isempty(comsol_col)
        warning('Skipping %s: required CTM/COMSOL columns missing.', csv_files(k).name);
        continue;
    end

    [V1, V2, ctm_grid, comsol_grid, err_grid] = build_grids(T, v1_col, v2_col, ctm_col, comsol_col, err_col);

    pair_name = derive_pair_name(T, csv_files(k).name);
    make_publication_plots(TARGET_DIR, pair_name, V1, V2, ctm_grid, comsol_grid, err_grid, var1_name, var2_name, N_CONTOURS);
end

fprintf('\nDone. Publication-ready replay plots regenerated in:\n%s\n', TARGET_DIR);


function name = find_var_col(var_names, prefix)
    idx = find(startsWith(var_names, prefix), 1, 'first');
    if isempty(idx)
        name = '';
    else
        name = var_names{idx};
    end
end


function col_name = pick_first_existing(var_names, candidates)
    col_name = '';
    for i = 1:numel(candidates)
        if ismember(candidates{i}, var_names)
            col_name = candidates{i};
            return;
        end
    end
end


function [V1, V2, ctm_grid, comsol_grid, err_grid] = build_grids(T, v1_col, v2_col, ctm_col, comsol_col, err_col)
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

        ctm_grid(iy, ix) = T.(ctm_col)(i);
        comsol_grid(iy, ix) = T.(comsol_col)(i);

        if ~isempty(err_col)
            err_grid(iy, ix) = T.(err_col)(i);
        else
            if ~isnan(ctm_grid(iy, ix)) && ~isnan(comsol_grid(iy, ix))
                err_grid(iy, ix) = 100 * abs(comsol_grid(iy, ix) - ctm_grid(iy, ix)) / max(abs(comsol_grid(iy, ix)), eps);
            end
        end
    end
end


function pair_name = derive_pair_name(T, csv_name)
    pair_name = erase(csv_name, '.csv');
    pair_name = erase(pair_name, '_2D_results_replay');
    pair_name = erase(pair_name, '_results_replay');

    if ismember('SweepType', T.Properties.VariableNames) && ~isempty(T.SweepType)
        sweep_val = T.SweepType(1);
        if iscategorical(sweep_val)
            pair_name = char(string(sweep_val));
        else
            pair_name = char(sweep_val);
        end
    end

    if isempty(pair_name)
        pair_name = 'Replay_2D';
    end
end


function make_publication_plots(out_dir, pair_name, V1, V2, ctm_grid, comsol_grid, err_grid, var1_name, var2_name, n_contours)
    xlab = get_latex_param_label(var1_name);
    ylab = get_latex_param_label(var2_name);
    pair_title = prettify_pair_name(pair_name);

    fig_temp = figure('Visible', 'off', 'Color', 'w', 'Position', [80, 80, 1300, 640]);
    tl = tiledlayout(1, 2, 'TileSpacing', 'compact', 'Padding', 'compact');

    ax1 = nexttile(tl, 1);
    contourf(ax1, V1, V2, ctm_grid, n_contours, 'LineStyle', 'none');
    colormap(ax1, hot);
    cb1 = colorbar(ax1);
    cb1.Label.String = 'CTM $T_{max}$ ($^\circ$C)';
    cb1.Label.Interpreter = 'latex';
    xlabel(ax1, xlab, 'Interpreter', 'latex');
    ylabel(ax1, ylab, 'Interpreter', 'latex');
    title(ax1, 'CTM $T_{max}$', 'Interpreter', 'latex');
    axis(ax1, 'square');
    grid(ax1, 'on');

    ax2 = nexttile(tl, 2);
    contourf(ax2, V1, V2, comsol_grid, n_contours, 'LineStyle', 'none');
    colormap(ax2, hot);
    cb2 = colorbar(ax2);
    cb2.Label.String = 'COMSOL $T_{max}$ ($^\circ$C)';
    cb2.Label.Interpreter = 'latex';
    xlabel(ax2, xlab, 'Interpreter', 'latex');
    ylabel(ax2, ylab, 'Interpreter', 'latex');
    title(ax2, 'COMSOL $T_{max}$', 'Interpreter', 'latex');
    axis(ax2, 'square');
    grid(ax2, 'on');

    sgtitle(tl, sprintf('2D Temperature Comparison: %s', pair_title), 'Interpreter', 'none', 'FontWeight', 'bold');

    out_base = sprintf('%s_publication_replot', pair_name);
    saveas(fig_temp, fullfile(out_dir, sprintf('%s_TemperatureComparison_plot.png', out_base)));
    savefig(fig_temp, fullfile(out_dir, sprintf('%s_TemperatureComparison_plot.fig', out_base)));
    close(fig_temp);

    fig_err = figure('Visible', 'off', 'Color', 'w', 'Position', [80, 80, 760, 700]);
    ax = axes(fig_err);
    contourf(ax, V1, V2, err_grid, n_contours, 'LineStyle', 'none');
    colormap(ax, jet);
    cb = colorbar(ax);
    cb.Label.String = 'Relative Error (%)';
    cb.Label.Interpreter = 'tex';
    xlabel(ax, xlab, 'Interpreter', 'latex');
    ylabel(ax, ylab, 'Interpreter', 'latex');
    title(ax, sprintf('Relative Error (%%): %s', pair_title), 'Interpreter', 'none', 'FontWeight', 'bold');
    axis(ax, 'square');
    grid(ax, 'on');

    saveas(fig_err, fullfile(out_dir, sprintf('%s_Error_plot.png', out_base)));
    savefig(fig_err, fullfile(out_dir, sprintf('%s_Error_plot.fig', out_base)));
    close(fig_err);
end


function label = get_latex_param_label(var_name)
    switch var_name
        case 'q_Wm2'
            label = '$q_{\mathrm{i}}\ [\mathrm{W/m^2}]$';
        case 'M'
            label = '$N$';
        case 't_TEC_um'
            label = '$t_{\mathrm{TEC}}\ [\mu\mathrm{m}]$';
        case 'f_L'
            label = '$f_{L}$';
        case 'fill_factor'
            label = '$f_{\mathrm{f}}$';
        case 'w_is_um'
            label = '$W_{\mathrm{is}}\ [\mu\mathrm{m}]$';
        case 't_chip_um'
            label = '$t_{\mathrm{gen}}\ [\mu\mathrm{m}]$';
        otherwise
            label = strrep(var_name, '_', '\_');
    end
end


function s = prettify_pair_name(pair_name)
    s = strrep(pair_name, '_vs_', ' vs ');
    s = strrep(s, '_', ' ');
    s = strrep(s, 'LengthRatio', 'Length Ratio');
    s = strrep(s, 'FillFactor', 'Fill Factor');
    s = strrep(s, 'RadialInsWidth', 'Radial Ins. Width');
    s = strrep(s, 'HeatFlux', 'Heat Flux');
    s = strrep(s, 'Thickness', 'Thickness');
end
