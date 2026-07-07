%% ANALYZE 2D SWEEPS (RANK METRICS + PUBLICATION PLOTS)
% Reads existing *_2D_results.csv files, computes Spearman/Kendall rank
% correlation between CTM and COMSOL Tmax, and regenerates publication-ready
% square plots without modifying existing result files.

clear; clc; close all;
fprintf('Running analyze_2d_sweeps_rank_and_replot v4\n');

%% ============ PATH SETUP ============
root_dir = fullfile(fileparts(mfilename('fullpath')), '..', '..');
cd(root_dir);
addpath(genpath('src'));
addpath(genpath('scripts'));

%% ============ USER CONFIG ============
% Point these to the batch folders produced by overnight runs.
BATCH_DIRS = {
    fullfile('output', 'validity_sweeps', '2026-04-17_23-17-47')
    fullfile('output', 'validity_sweeps', '2026-04-18_04-33-37')
};

REGENERATE_PLOTS = true;
N_CONTOURS = 20;

%% ============ PROCESS EACH BATCH ============
for b = 1:numel(BATCH_DIRS)
    batch_dir = BATCH_DIRS{b};
    if ~isfolder(batch_dir)
        fprintf('Skipping missing folder: %s\n', batch_dir);
        continue;
    end

    csv_files = dir(fullfile(batch_dir, '*', '*_2D_results.csv'));
    if isempty(csv_files)
        fprintf('No 2D result CSV files found in: %s\n', batch_dir);
        continue;
    end

    fprintf('\n============================================================\n');
    fprintf('Batch: %s\n', batch_dir);
    fprintf('Found %d pair CSV files.\n', numel(csv_files));

    summary = table();

    for k = 1:numel(csv_files)
        csv_path = fullfile(csv_files(k).folder, csv_files(k).name);
        pair_folder = csv_files(k).folder;
        pair_name = erase(csv_files(k).name, '_2D_results.csv');

        try
            T = readtable(csv_path);
        catch ME
            fprintf('  Failed to read %s: %s\n', csv_files(k).name, ME.message);
            continue;
        end

        [var1_col, var2_col] = find_var_columns(T.Properties.VariableNames);
        if isempty(var1_col) || isempty(var2_col)
            fprintf('  Skipping %s (Var1_/Var2_ columns not found).\n', csv_files(k).name);
            continue;
        end

        var1_name = erase(var1_col, 'Var1_');
        var2_name = erase(var2_col, 'Var2_');

        if ~ismember('MATLAB_Tmax', T.Properties.VariableNames) || ~ismember('COMSOL_Tmax', T.Properties.VariableNames)
            fprintf('  Skipping %s (MATLAB_Tmax/COMSOL_Tmax missing).\n', csv_files(k).name);
            continue;
        end

        valid = ~isnan(T.MATLAB_Tmax) & ~isnan(T.COMSOL_Tmax);
        n_valid = nnz(valid);
        if n_valid < 3
            fprintf('  Skipping %s (not enough valid points).\n', csv_files(k).name);
            continue;
        end

        x = T.MATLAB_Tmax(valid);
        y = T.COMSOL_Tmax(valid);

        rho = corr(x, y, 'Type', 'Spearman', 'Rows', 'complete');
        % Manual Spearman p-value (t-approximation) to inspect tiny-tail values.
        % Guard edge cases where approximation is not reliable.
        nu = n_valid - 2;
        if nu > 0 && abs(rho) < 1 - 1e-12
            t_stat = rho * sqrt(nu / max(1 - rho^2, eps));
            p_s_manual = 2 * tcdf(abs(t_stat), nu, 'upper');
        else
            p_s_manual = NaN;
        end
        [tau, p_k] = corr(x, y, 'Type', 'Kendall',  'Rows', 'complete');
        trend_label = classify_trend_quality(rho, tau);

        summary = [summary; table(categorical({prettify_pair_name(pair_name)}), ...
                                  rho, tau, p_s_manual, p_k, n_valid, categorical({trend_label}), ...
                                  'VariableNames', {'DesignSpace_2D_Sweep', 'Spearman_rho', 'Kendall_tau', ...
                                                    'Spearman_p', 'Kendall_p', 'N_valid', 'Trend_Quality'})];

        fprintf('  %-40s rho=% .4f  tau=% .4f  p_s=% .3e  p_k=% .3e  N=%d\n', ...
            pair_name, rho, tau, p_s_manual, p_k, n_valid);

        if REGENERATE_PLOTS
            try
                [V1, V2, ctm_grid, comsol_grid, err_grid] = build_grids(T, var1_col, var2_col);
                make_publication_plots(pair_folder, pair_name, V1, V2, ctm_grid, comsol_grid, err_grid, ...
                                       var1_name, var2_name, N_CONTOURS);
            catch ME
                fprintf('    Plot regen failed for %s: %s\n', pair_name, ME.message);
            end
        end
    end

    if isempty(summary)
        fprintf('No valid pair analyses for this batch.\n');
        continue;
    end

    avg_rho = mean(summary.Spearman_rho, 'omitnan');
    avg_tau = mean(summary.Kendall_tau, 'omitnan');
    global_quality = classify_global_quality(avg_rho, avg_tau);

    summary = [summary; table(categorical({'Global Average'}), avg_rho, avg_tau, NaN, NaN, ...
                              sum(summary.N_valid, 'omitnan'), categorical({global_quality}), ...
                              'VariableNames', summary.Properties.VariableNames)];

    out_csv = fullfile(batch_dir, 'rank_correlation_summary.csv');
    writetable(summary, out_csv);

    fprintf('\nSaved rank summary: %s\n', out_csv);
    disp(summary(:, {'DesignSpace_2D_Sweep', 'Spearman_rho', 'Kendall_tau', 'Spearman_p', 'Kendall_p', 'Trend_Quality'}));

    try
        make_rank_barplot(summary, batch_dir);
    catch ME
        fprintf('Could not create rank summary plot: %s\n', ME.message);
    end
end

fprintf('\nDone. Existing result files were not modified.\n');


function [var1_col, var2_col] = find_var_columns(var_names)
    v1_idx = find(startsWith(var_names, 'Var1_'), 1, 'first');
    v2_idx = find(startsWith(var_names, 'Var2_'), 1, 'first');

    if isempty(v1_idx), var1_col = ''; else, var1_col = var_names{v1_idx}; end
    if isempty(v2_idx), var2_col = ''; else, var2_col = var_names{v2_idx}; end
end


function [V1, V2, ctm_grid, comsol_grid, err_grid] = build_grids(T, var1_col, var2_col)
    v1_vals = unique(T.(var1_col));
    v2_vals = unique(T.(var2_col));

    [V1, V2] = meshgrid(v1_vals, v2_vals);
    ctm_grid = nan(size(V1));
    comsol_grid = nan(size(V1));
    err_grid = nan(size(V1));

    for i = 1:height(T)
        ix = find(v1_vals == T.(var1_col)(i), 1, 'first');
        iy = find(v2_vals == T.(var2_col)(i), 1, 'first');
        if isempty(ix) || isempty(iy)
            continue;
        end

        if ismember('MATLAB_Tmax', T.Properties.VariableNames)
            ctm_grid(iy, ix) = T.MATLAB_Tmax(i);
        end
        if ismember('COMSOL_Tmax', T.Properties.VariableNames)
            comsol_grid(iy, ix) = T.COMSOL_Tmax(i);
        end

        if ismember('Error_Pct', T.Properties.VariableNames)
            err_grid(iy, ix) = T.Error_Pct(i);
        else
            if ~isnan(ctm_grid(iy, ix)) && ~isnan(comsol_grid(iy, ix))
                err_grid(iy, ix) = 100 * abs(comsol_grid(iy, ix) - ctm_grid(iy, ix)) / max(abs(comsol_grid(iy, ix)), eps);
            end
        end
    end
end


function make_publication_plots(out_dir, pair_name, V1, V2, ctm_grid, comsol_grid, err_grid, var1_name, var2_name, n_contours)
    xlab = get_latex_param_label(var1_name);
    ylab = get_latex_param_label(var2_name);
    pair_title = prettify_pair_name(pair_name);

    % Side-by-side CTM vs COMSOL (both square with separate colorbars)
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

    saveas(fig_temp, fullfile(out_dir, sprintf('%s_TemperatureComparison_plot_pub.png', pair_name)));
    savefig(fig_temp, fullfile(out_dir, sprintf('%s_TemperatureComparison_plot_pub.fig', pair_name)));
    close(fig_temp);

    % Error plot (square + LaTeX colorbar label)
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

    saveas(fig_err, fullfile(out_dir, sprintf('%s_Error_plot_pub.png', pair_name)));
    savefig(fig_err, fullfile(out_dir, sprintf('%s_Error_plot_pub.fig', pair_name)));
    close(fig_err);
end


function make_rank_barplot(summary_table, out_dir)
    is_global = strcmp(string(summary_table.DesignSpace_2D_Sweep), "Global Average");
    T = summary_table(~is_global, :);
    if isempty(T)
        return;
    end

    fig = figure('Visible', 'off', 'Color', 'w', 'Position', [80, 80, 1200, 520]);
    vals = [T.Spearman_rho, T.Kendall_tau];
    bh = bar(vals, 'grouped');
    bh(1).FaceColor = [0.20 0.45 0.70];
    bh(2).FaceColor = [0.85 0.45 0.10];
    ylim([0, 1]);
    grid on;

    ax = gca;
    ax.XTick = 1:height(T);
    ax.XTickLabel = cellstr(T.DesignSpace_2D_Sweep);
    ax.XTickLabelRotation = 30;
    ylabel('Rank Correlation Coefficient');
    title('Per-Pair Rank Correlation Summary');
    legend({'Spearman \rho', 'Kendall \tau'}, 'Location', 'best');

    saveas(fig, fullfile(out_dir, 'rank_correlation_summary_plot.png'));
    savefig(fig, fullfile(out_dir, 'rank_correlation_summary_plot.fig'));
    close(fig);
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
    s = strrep(s, 'LengthRatio', 'Length Ratio');
    s = strrep(s, 'FillFactor', 'Fill Factor');
    s = strrep(s, 'RadialInsWidth', 'Radial Ins. Width');
    s = strrep(s, 'HeatFlux', 'Heat Flux');
    s = strrep(s, 'Thickness', 'Thickness');
end


function q = classify_trend_quality(rho, tau)
    score = min(abs([rho, tau]));
    if score >= 0.90
        q = 'Excellent';
    elseif score >= 0.75
        q = 'Strong';
    elseif score >= 0.55
        q = 'Moderate';
    else
        q = 'Weak';
    end
end


function q = classify_global_quality(avg_rho, avg_tau)
    score = min(abs([avg_rho, avg_tau]));
    if score >= 0.80
        q = 'High Fidelity';
    elseif score >= 0.65
        q = 'Moderate Fidelity';
    else
        q = 'Low Fidelity';
    end
end
