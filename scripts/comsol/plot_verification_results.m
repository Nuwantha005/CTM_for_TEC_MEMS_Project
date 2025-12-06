%% PLOT VERIFICATION RESULTS
% This script generates comparison plots from the verification results
% Run this after completing all 5 test cases with run_verification_case.m
%
% Usage:
%   1. Complete all 5 test cases
%   2. Run: plot_verification_results

clear; clc; close all;

%% ============ LOAD RESULTS ============
OUTPUT_DIR = fullfile('output', 'model_verification');
csv_file = fullfile(OUTPUT_DIR, 'verification_results.csv');

if ~exist(csv_file, 'file')
    error('Results file not found: %s\nRun all test cases first using run_verification_case.m', csv_file);
end

% Read CSV
data = readtable(csv_file);
fprintf('Loaded %d test case results from:\n  %s\n\n', height(data), csv_file);

%% ============ DISPLAY SUMMARY TABLE ============
fprintf('VERIFICATION RESULTS SUMMARY\n');
fprintf('================================================================================\n');
fprintf('%-5s %-25s %12s %12s %12s %10s\n', 'ID', 'Test Case', 'MATLAB [C]', 'COMSOL [C]', 'Error [C]', 'Error [%]');
fprintf('--------------------------------------------------------------------------------\n');

for i = 1:height(data)
    fprintf('%-5d %-25s %12.2f %12.2f %12.2f %10.1f%%\n', ...
        data.CaseID(i), char(data.Name(i)), data.MATLAB_Tmax(i), ...
        data.COMSOL_Tmax(i), data.Error_Tmax(i), data.Error_pct(i));
end
fprintf('================================================================================\n\n');

% Statistics
valid_idx = ~isnan(data.Error_pct);
if any(valid_idx)
    errors = abs(data.Error_pct(valid_idx));
    fprintf('Statistics:\n');
    fprintf('  Mean Absolute Error: %.2f%%\n', mean(errors));
    fprintf('  Max Absolute Error:  %.2f%%\n', max(errors));
    fprintf('  Min Absolute Error:  %.2f%%\n', min(errors));
    fprintf('  Std Deviation:       %.2f%%\n\n', std(errors));
end

%% ============ CREATE FIGURE ============
fig = figure('Position', [100, 100, 1400, 600], 'Color', 'w');

%% ============ PLOT 1: Temperature Comparison ============
subplot(2, 3, 1);
x = 1:height(data);
bar_data = [data.MATLAB_Tmax, data.COMSOL_Tmax];
b = bar(x, bar_data, 'grouped');
b(1).FaceColor = [0.2 0.4 0.8];  % Blue for MATLAB
b(2).FaceColor = [0.8 0.3 0.2];  % Red for COMSOL

set(gca, 'XTick', x);
labels = cellfun(@(s) strrep(s, '_', ' '), data.Name, 'UniformOutput', false);
set(gca, 'XTickLabel', labels);
xtickangle(30);
ylabel('Maximum Temperature [°C]');
title('Temperature Comparison: MATLAB vs COMSOL');
legend('MATLAB CTM', 'COMSOL FEM', 'Location', 'northwest');
grid on;
box on;

%% ============ PLOT 2: Error Analysis ============
subplot(2, 3, 2);
bar(x, data.Error_Tmax, 'FaceColor', [0.4 0.7 0.4]);
set(gca, 'XTick', x);
set(gca, 'XTickLabel', labels);
xtickangle(30);
ylabel('Error (COMSOL - MATLAB) [°C]');
title('Temperature Prediction Error');
grid on;
box on;

% Add reference lines
hold on;
yline(0, 'k-', 'LineWidth', 1.5);
yline(mean(data.Error_Tmax), 'r--', 'LineWidth', 1.5);
text(0.5, mean(data.Error_Tmax)+0.5, sprintf('Mean: %.2f°C', mean(data.Error_Tmax)), ...
    'Color', 'r', 'FontSize', 9);
hold off;

%% ============ PLOT 3: Percentage Error ============
subplot(2, 3, 3);
bar(x, abs(data.Error_pct), 'FaceColor', [0.9 0.6 0.2]);
set(gca, 'XTick', x);
set(gca, 'XTickLabel', labels);
xtickangle(30);
ylabel('Absolute Error [%]');
title('Percentage Error');
grid on;
box on;

% Add threshold lines
hold on;
yline(5, 'g--', '5% threshold', 'LineWidth', 1.5, 'LabelHorizontalAlignment', 'left');
yline(10, 'r--', '10% threshold', 'LineWidth', 1.5, 'LabelHorizontalAlignment', 'left');
hold off;

%% ============ PLOT 4: Correlation Plot ============
subplot(2, 3, 4);
scatter(data.MATLAB_Tmax, data.COMSOL_Tmax, 100, 'filled', 'MarkerFaceColor', [0.2 0.6 0.4]);
hold on;

% Perfect correlation line
min_T = min([data.MATLAB_Tmax; data.COMSOL_Tmax]) - 2;
max_T = max([data.MATLAB_Tmax; data.COMSOL_Tmax]) + 2;
plot([min_T max_T], [min_T max_T], 'k--', 'LineWidth', 1.5);

% Error bounds (±5%)
plot([min_T max_T], [min_T*1.05 max_T*1.05], 'r:', 'LineWidth', 1);
plot([min_T max_T], [min_T*0.95 max_T*0.95], 'r:', 'LineWidth', 1);

xlabel('MATLAB T_{max} [°C]');
ylabel('COMSOL T_{max} [°C]');
title('Correlation: MATLAB vs COMSOL');
legend('Data Points', 'Perfect Correlation', '±5% Bounds', 'Location', 'southeast');
grid on;
box on;
axis equal;
xlim([min_T max_T]);
ylim([min_T max_T]);
hold off;

%% ============ PLOT 5: Parameter Sensitivity ============
subplot(2, 3, 5);
% Plot error vs heat flux
scatter(data.q_Wm2, abs(data.Error_pct), 100, data.I_mA, 'filled');
colorbar;
xlabel('Heat Flux [W/m²]');
ylabel('Absolute Error [%]');
title('Error vs Heat Flux (color = Current)');
grid on;
box on;

%% ============ PLOT 6: Error Distribution ============
subplot(2, 3, 6);
if height(data) >= 3
    histogram(data.Error_pct, 'FaceColor', [0.5 0.5 0.8], 'EdgeColor', 'w');
    xlabel('Error [%]');
    ylabel('Count');
    title('Error Distribution');
    grid on;
    box on;
    
    % Add statistics
    hold on;
    xline(mean(data.Error_pct), 'r-', 'LineWidth', 2);
    xline(mean(data.Error_pct) + std(data.Error_pct), 'r--', 'LineWidth', 1);
    xline(mean(data.Error_pct) - std(data.Error_pct), 'r--', 'LineWidth', 1);
    legend('Error Distribution', 'Mean', '±1 Std Dev', 'Location', 'best');
    hold off;
else
    % Not enough data for histogram, show text summary
    text(0.5, 0.5, sprintf('Mean Error: %.2f%%\nStd Dev: %.2f%%', ...
        mean(data.Error_pct), std(data.Error_pct)), ...
        'HorizontalAlignment', 'center', 'FontSize', 14);
    axis off;
    title('Error Statistics');
end

%% ============ ADD OVERALL TITLE ============
sgtitle('Mathematical Model Verification Results', 'FontSize', 14, 'FontWeight', 'bold');

%% ============ SAVE FIGURE ============
savefig(fullfile(OUTPUT_DIR, 'verification_plots.fig'));
saveas(fig, fullfile(OUTPUT_DIR, 'verification_plots.png'));
fprintf('Plots saved to:\n');
fprintf('  %s\n', fullfile(OUTPUT_DIR, 'verification_plots.fig'));
fprintf('  %s\n', fullfile(OUTPUT_DIR, 'verification_plots.png'));

%% ============ DETAILED ANALYSIS ============
fprintf('\n\nDETAILED ANALYSIS\n');
fprintf('================================================================================\n');

fprintf('\nTest Cases by Error Magnitude:\n');
[~, sort_idx] = sort(abs(data.Error_pct));
for i = 1:height(data)
    idx = sort_idx(i);
    if abs(data.Error_pct(idx)) < 5
        status = 'EXCELLENT';
    elseif abs(data.Error_pct(idx)) < 10
        status = 'GOOD';
    else
        status = 'NEEDS REVIEW';
    end
    fprintf('  %d. Case %d (%-25s): %.1f%% - %s\n', i, data.CaseID(idx), ...
        char(data.Name(idx)), abs(data.Error_pct(idx)), status);
end

fprintf('\n================================================================================\n');
fprintf('CONCLUSION:\n');
if mean(abs(data.Error_pct)) < 5
    fprintf('  The compact thermal model (CTM) shows EXCELLENT agreement with COMSOL FEM.\n');
    fprintf('  Mean absolute error: %.2f%% (< 5%% threshold)\n', mean(abs(data.Error_pct)));
    fprintf('  The model is well-suited for preliminary design optimization.\n');
elseif mean(abs(data.Error_pct)) < 10
    fprintf('  The compact thermal model (CTM) shows GOOD agreement with COMSOL FEM.\n');
    fprintf('  Mean absolute error: %.2f%% (< 10%% threshold)\n', mean(abs(data.Error_pct)));
    fprintf('  The model is acceptable for preliminary design, with some limitations.\n');
else
    fprintf('  The compact thermal model (CTM) shows MODERATE agreement with COMSOL FEM.\n');
    fprintf('  Mean absolute error: %.2f%% (> 10%% threshold)\n', mean(abs(data.Error_pct)));
    fprintf('  Consider model refinements for improved accuracy.\n');
end
fprintf('================================================================================\n\n');
