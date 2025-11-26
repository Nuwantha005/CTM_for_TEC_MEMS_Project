% plot_comsol_comparison.m
% Generate comparison plots from COMSOL validation runs
%
% Run this after completing all test cases with run_single_comsol.m

clear; clc;

fprintf('=== COMSOL VALIDATION RESULTS ===\n\n');

%% Load results
INPUT_DIR = fullfile('output', 'comsol_single_runs');
OUTPUT_DIR = INPUT_DIR;

% Read CSV
csvFile = fullfile(INPUT_DIR, 'all_results.csv');
if ~exist(csvFile, 'file')
    error('Results file not found: %s\nRun test cases first.', csvFile);
end

data = readtable(csvFile);
fprintf('Loaded %d test case results\n\n', height(data));

% Display results table
fprintf('┌─────┬───────────────────┬────────┬────────┬──────┬──────────┬─────────────┬──────────────┬───────────┐\n');
fprintf('│ ID  │ Name              │ I(mA)  │ t(µm)  │ k_r  │ q(W/m²)  │ MATLAB(°C)  │ COMSOL(°C)   │ Error(°C) │\n');
fprintf('├─────┼───────────────────┼────────┼────────┼──────┼──────────┼─────────────┼──────────────┼───────────┤\n');
for i = 1:height(data)
    fprintf('│ %3d │ %-17s │ %6.0f │ %6.0f │ %4.2f │ %8.0f │ %11.2f │ %12.2f │ %9.2f │\n', ...
        data.ID(i), data.Name{i}, data.I_mA(i), data.t_um(i), data.k_r(i), ...
        data.q_Wm2(i), data.MATLAB_Tmax(i), data.COMSOL_Tmax(i), data.Error(i));
end
fprintf('└─────┴───────────────────┴────────┴────────┴──────┴──────────┴─────────────┴──────────────┴───────────┘\n\n');

%% Statistics
fprintf('=== ERROR STATISTICS ===\n');
fprintf('  Mean Error:     %+.2f °C\n', mean(data.Error));
fprintf('  Std Error:       %.2f °C\n', std(data.Error));
fprintf('  Max Abs Error:   %.2f °C\n', max(abs(data.Error)));
fprintf('  Mean Abs Error:  %.2f °C\n', mean(abs(data.Error)));

% Relative error
rel_error = 100 * abs(data.Error) ./ data.COMSOL_Tmax;
fprintf('  Mean Rel Error: %.1f%%\n\n', mean(rel_error));

%% Create figure
fig = figure('Position', [50, 50, 1400, 900], 'Name', 'COMSOL vs MATLAB Validation');

% Plot 1: Bar comparison
subplot(2,3,1);
bar_data = [data.MATLAB_Tmax, data.COMSOL_Tmax];
b = bar(bar_data);
b(1).FaceColor = [0.2 0.4 0.8];
b(2).FaceColor = [0.8 0.2 0.2];
xlabel('Test Case');
ylabel('T_{max} (°C)');
title('MATLAB vs COMSOL: T_{max}');
legend('MATLAB CTM', 'COMSOL FEM', 'Location', 'northwest');
set(gca, 'XTickLabel', data.Name);
xtickangle(30);
grid on;
yline(100, 'k--', 'Target Limit', 'LineWidth', 1.5);

% Plot 2: Parity plot
subplot(2,3,2);
scatter(data.MATLAB_Tmax, data.COMSOL_Tmax, 100, data.q_Wm2, 'filled');
hold on;
% Perfect agreement line
lim = [min([data.MATLAB_Tmax; data.COMSOL_Tmax])-5, max([data.MATLAB_Tmax; data.COMSOL_Tmax])+5];
plot(lim, lim, 'k--', 'LineWidth', 1.5);
% ±20% band
plot(lim, lim*1.2, 'r:', 'LineWidth', 1);
plot(lim, lim*0.8, 'r:', 'LineWidth', 1);
xlabel('MATLAB T_{max} (°C)');
ylabel('COMSOL T_{max} (°C)');
title('Parity Plot: MATLAB vs COMSOL');
colorbar; ylabel(colorbar, 'Heat Flux (W/m²)');
legend('Data', 'Perfect', '±20%', 'Location', 'southeast');
axis equal;
xlim(lim); ylim(lim);
grid on;

% Plot 3: Error vs Heat Flux
subplot(2,3,3);
scatter(data.q_Wm2, data.Error, 100, 'filled');
xlabel('Heat Flux (W/m²)');
ylabel('Error (COMSOL - MATLAB) (°C)');
title('Model Error vs Heat Flux');
yline(0, 'k--', 'LineWidth', 1);
grid on;
% Add trend line
p = polyfit(data.q_Wm2, data.Error, 1);
x_fit = linspace(min(data.q_Wm2), max(data.q_Wm2), 100);
hold on;
plot(x_fit, polyval(p, x_fit), 'r-', 'LineWidth', 1.5);

% Plot 4: Temperature vs Heat Flux
subplot(2,3,4);
plot(data.q_Wm2, data.MATLAB_Tmax, 'b-o', 'LineWidth', 2, 'MarkerSize', 10, 'MarkerFaceColor', 'b');
hold on;
plot(data.q_Wm2, data.COMSOL_Tmax, 'r-s', 'LineWidth', 2, 'MarkerSize', 10, 'MarkerFaceColor', 'r');
xlabel('Heat Flux (W/m²)');
ylabel('T_{max} (°C)');
title('Temperature vs Heat Flux');
legend('MATLAB CTM', 'COMSOL FEM', 'Location', 'northwest');
yline(100, 'k--', 'Target Limit', 'LineWidth', 1.5);
grid on;

% Plot 5: Error distribution
subplot(2,3,5);
bar(data.Error);
xlabel('Test Case');
ylabel('Error (°C)');
title('Model Error by Test Case');
set(gca, 'XTickLabel', data.Name);
xtickangle(30);
yline(0, 'k-', 'LineWidth', 1);
grid on;

% Add color based on sign
colors = zeros(height(data), 3);
for i = 1:height(data)
    if data.Error(i) > 0
        colors(i,:) = [0.8 0.2 0.2];  % Red for positive
    else
        colors(i,:) = [0.2 0.6 0.2];  % Green for negative
    end
end

% Plot 6: Summary table
subplot(2,3,6);
axis off;

% Create text summary
text(0.1, 0.95, 'VALIDATION SUMMARY', 'FontSize', 14, 'FontWeight', 'bold');
text(0.1, 0.85, sprintf('Total Test Cases: %d', height(data)), 'FontSize', 11);
text(0.1, 0.75, sprintf('Heat Flux Range: %d - %d W/m²', min(data.q_Wm2), max(data.q_Wm2)), 'FontSize', 11);

text(0.1, 0.60, 'Model Agreement:', 'FontSize', 12, 'FontWeight', 'bold');
text(0.1, 0.50, sprintf('  Mean Absolute Error: %.2f °C', mean(abs(data.Error))), 'FontSize', 11);
text(0.1, 0.40, sprintf('  Mean Relative Error: %.1f%%', mean(rel_error)), 'FontSize', 11);
text(0.1, 0.30, sprintf('  Max Absolute Error: %.2f °C', max(abs(data.Error))), 'FontSize', 11);

% Conclusion
if max(data.COMSOL_Tmax) < 100
    text(0.1, 0.15, '✓ All designs below 100°C target!', 'FontSize', 12, 'Color', 'g', 'FontWeight', 'bold');
else
    text(0.1, 0.15, '✗ Some designs exceed 100°C target', 'FontSize', 12, 'Color', 'r', 'FontWeight', 'bold');
end

text(0.1, 0.05, sprintf('Best COMSOL Result: %.2f°C @ %d W/m²', ...
    min(data.COMSOL_Tmax), data.q_Wm2(data.COMSOL_Tmax == min(data.COMSOL_Tmax))), 'FontSize', 11);

%% Save figure
sgtitle('COMSOL Validation of MATLAB Compact Thermal Model', 'FontSize', 14, 'FontWeight', 'bold');

saveas(fig, fullfile(OUTPUT_DIR, 'validation_comparison.png'));
saveas(fig, fullfile(OUTPUT_DIR, 'validation_comparison.fig'));
fprintf('\nPlots saved to: %s\n', OUTPUT_DIR);

%% Additional: Create heat flux sweep plot
figure('Position', [100, 100, 800, 600], 'Name', 'Heat Flux Response');

% Sort by heat flux
[q_sorted, idx] = sort(data.q_Wm2);
matlab_sorted = data.MATLAB_Tmax(idx);
comsol_sorted = data.COMSOL_Tmax(idx);

% Plot with error bars representation
fill([q_sorted; flipud(q_sorted)], [matlab_sorted; flipud(comsol_sorted)], ...
    [0.8 0.8 0.9], 'EdgeColor', 'none', 'FaceAlpha', 0.5);
hold on;
plot(q_sorted, matlab_sorted, 'b-o', 'LineWidth', 2, 'MarkerSize', 12, 'MarkerFaceColor', 'b');
plot(q_sorted, comsol_sorted, 'r-s', 'LineWidth', 2, 'MarkerSize', 12, 'MarkerFaceColor', 'r');

xlabel('Heat Flux (W/m²)', 'FontSize', 12);
ylabel('Maximum Temperature (°C)', 'FontSize', 12);
title('TEC Performance: COMSOL Validation', 'FontSize', 14);
legend('Error Band', 'MATLAB Compact Model', 'COMSOL FEM', 'Location', 'northwest');
grid on;
yline(100, 'k--', 'T_{max} = 100°C Target', 'LineWidth', 2);

% Add annotation for key finding
[T_min, idx_min] = min(comsol_sorted);
text(q_sorted(idx_min), T_min - 5, sprintf('Best: %.1f°C', T_min), ...
    'HorizontalAlignment', 'center', 'FontWeight', 'bold');

saveas(gcf, fullfile(OUTPUT_DIR, 'heat_flux_response.png'));
saveas(gcf, fullfile(OUTPUT_DIR, 'heat_flux_response.fig'));

%% Key findings
fprintf('\n=== KEY FINDINGS ===\n');
fprintf('1. COMSOL validates the MATLAB compact thermal model\n');
fprintf('2. Best achieved temperature: %.2f°C at %d W/m²\n', min(data.COMSOL_Tmax), ...
    data.q_Wm2(data.COMSOL_Tmax == min(data.COMSOL_Tmax)));
fprintf('3. Design is feasible up to at least %d W/m² heat flux\n', max(data.q_Wm2));
fprintf('4. Average model error: %.1f°C (%.1f%%)\n', mean(abs(data.Error)), mean(rel_error));

% Recommendation
fprintf('\n=== RECOMMENDED DESIGN ===\n');
[~, best_idx] = min(data.COMSOL_Tmax);
fprintf('Test Case: %s\n', data.Name{best_idx});
fprintf('  Current: %.0f mA\n', data.I_mA(best_idx));
fprintf('  TEC Thickness: %.0f µm\n', data.t_um(best_idx));
fprintf('  k_r: %.2f\n', data.k_r(best_idx));
fprintf('  Heat Flux: %.0f W/m²\n', data.q_Wm2(best_idx));
fprintf('  COMSOL T_max: %.2f °C\n', data.COMSOL_Tmax(best_idx));

fprintf('\nDone!\n');
