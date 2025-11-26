% extract_comsol_data.m
% Extract data from COMSOL and generate analysis plots
%
% This script connects to COMSOL, runs simulations with the optimal parameters
% found from MATLAB optimization, and extracts the results for comparison.
%
% Prerequisites:
%   1. Start COMSOL server: 
%      "F:\EngineeringSoftware\COMSOL\COMSOL63\Multiphysics\bin\win64\comsolmphserver -port 2036"
%   2. Have the COMSOL template model (.mph) ready
%   3. Model should have "Node 0 Temp" probe for center temperature

clear; clc;
addpath(genpath('src'));

fprintf('=== COMSOL DATA EXTRACTION ===\n\n');

%% Configuration - EDIT THESE
COMSOL_PORT = 2036;
COMSOL_MODEL_PATH = 'E:\Semester 7\ME4311 - MicroNano Electro Mechanical Systems and Nanotechnology\Project\COMSOL\MATLAB_API\template_3_stage.mph';

% Add COMSOL LiveLink path
COMSOL_MLI_PATH = 'F:\EngineeringSoftware\COMSOL\COMSOL63\Multiphysics\mli';
if exist(COMSOL_MLI_PATH, 'dir')
    addpath(COMSOL_MLI_PATH);
    fprintf('Added COMSOL LiveLink to path\n');
end

%% Output setup
timestamp = datestr(now, 'yyyy-mm-dd_HH-MM-SS');
OUTPUT_DIR = fullfile('output', 'comsol_data', timestamp);
if ~exist(OUTPUT_DIR, 'dir')
    mkdir(OUTPUT_DIR);
end
fprintf('Output: %s\n\n', OUTPUT_DIR);

%% Get model path
if isempty(COMSOL_MODEL_PATH)
    fprintf('Please select the COMSOL template model (.mph file)...\n');
    [file, path] = uigetfile('*.mph', 'Select COMSOL Template Model');
    if file == 0
        error('No model selected');
    end
    COMSOL_MODEL_PATH = fullfile(path, file);
end
fprintf('Model: %s\n\n', COMSOL_MODEL_PATH);

%% Define test cases based on MATLAB optimization results
% These are the best configurations from our optimization
% Constrained to 3 stages with TSV in stage 1 only

test_cases = struct([]);

% Test Case 1: Best from optimization (low flux)
test_cases(1).name = 'Optimal_LowFlux';
test_cases(1).I_current_A = 0.025;
test_cases(1).thickness_um = 200;
test_cases(1).k_r = 1.15;
test_cases(1).q_flux_W_m2 = 500;
test_cases(1).wedge_angle_deg = 30;

% Test Case 2: Higher current
test_cases(2).name = 'HighCurrent';
test_cases(2).I_current_A = 0.1;
test_cases(2).thickness_um = 200;
test_cases(2).k_r = 1.15;
test_cases(2).q_flux_W_m2 = 500;
test_cases(2).wedge_angle_deg = 30;

% Test Case 3: Thicker TEC
test_cases(3).name = 'ThickTEC';
test_cases(3).I_current_A = 0.05;
test_cases(3).thickness_um = 300;
test_cases(3).k_r = 1.15;
test_cases(3).q_flux_W_m2 = 500;
test_cases(3).wedge_angle_deg = 30;

% Test Case 4: Higher heat flux
test_cases(4).name = 'HighFlux';
test_cases(4).I_current_A = 0.05;
test_cases(4).thickness_um = 200;
test_cases(4).k_r = 1.15;
test_cases(4).q_flux_W_m2 = 1000;
test_cases(4).wedge_angle_deg = 30;

% Test Case 5: Maximum feasible flux
test_cases(5).name = 'MaxFlux';
test_cases(5).I_current_A = 0.05;
test_cases(5).thickness_um = 300;
test_cases(5).k_r = 1.0;
test_cases(5).q_flux_W_m2 = 2000;
test_cases(5).wedge_angle_deg = 30;

%% Connect to COMSOL
fprintf('Connecting to COMSOL server on port %d...\n', COMSOL_PORT);

try
    mphstart(COMSOL_PORT);
    fprintf('Connected successfully!\n\n');
catch ME
    fprintf('*** COMSOL CONNECTION FAILED ***\n');
    fprintf('Error: %s\n\n', ME.message);
    fprintf('Please ensure COMSOL server is running:\n');
    fprintf('  F:\\EngineeringSoftware\\COMSOL\\COMSOL63\\Multiphysics\\bin\\win64\\comsolmphserver -port %d\n\n', COMSOL_PORT);
    return;
end

%% Load model
fprintf('Loading model...\n');
try
    model = mphload(COMSOL_MODEL_PATH);
    fprintf('Model loaded successfully!\n\n');
catch ME
    fprintf('Failed to load model: %s\n', ME.message);
    return;
end

%% Display current parameters
fprintf('Current model parameters:\n');
try
    params = mphgetexpressions(model.param);
    param_names = fieldnames(params);
    for i = 1:min(20, length(param_names))  % Show first 20
        fprintf('  %s = %s\n', param_names{i}, params.(param_names{i}));
    end
    fprintf('\n');
catch
    fprintf('  Could not retrieve parameters\n\n');
end

%% Run test cases
fprintf('=== RUNNING TEST CASES ===\n\n');

results = struct([]);

for i = 1:length(test_cases)
    tc = test_cases(i);
    fprintf('--- Test Case %d: %s ---\n', i, tc.name);
    
    % Set parameters
    try
        model.param.set('I0', sprintf('%g[A]', tc.I_current_A));
        model.param.set('LL_t_TEC', sprintf('%g[um]', tc.thickness_um));
        model.param.set('LL_k_r', sprintf('%g', tc.k_r));
        model.param.set('q', sprintf('%g[W/m^2]', tc.q_flux_W_m2));
        model.param.set('LL_theta', sprintf('%g[deg]', tc.wedge_angle_deg));
        
        fprintf('  Parameters set:\n');
        fprintf('    I = %.3f A, t = %d um, k_r = %.2f, q = %d W/m², θ = %d°\n', ...
            tc.I_current_A, tc.thickness_um, tc.k_r, tc.q_flux_W_m2, tc.wedge_angle_deg);
    catch ME
        fprintf('  Warning: Could not set some parameters: %s\n', ME.message);
    end
    
    % Run study
    fprintf('  Running simulation...\n');
    tic;
    try
        model.study('std1').run();
        elapsed = toc;
        fprintf('  Completed in %.1f seconds\n', elapsed);
        
        % Extract Node 0 temperature
        try
            % Try probe first
            T_node0 = mphglobal(model, 'comp1.ppb1', 'dataset', 'dset1');
            fprintf('  Node 0 Temp (probe): %.2f K (%.2f°C)\n', T_node0, T_node0 - 273.15);
        catch
            % Try max temperature
            try
                T_data = mpheval(model, 'T', 'dataset', 'dset1');
                T_node0 = max(T_data.d1(:));
                fprintf('  T_max: %.2f K (%.2f°C)\n', T_node0, T_node0 - 273.15);
            catch ME2
                T_node0 = NaN;
                fprintf('  Could not extract temperature: %s\n', ME2.message);
            end
        end
        
        % Store results
        results(i).name = tc.name;
        results(i).params = tc;
        results(i).T_K = T_node0;
        results(i).T_C = T_node0 - 273.15;
        results(i).time = elapsed;
        results(i).success = true;
        
    catch ME
        fprintf('  FAILED: %s\n', ME.message);
        results(i).name = tc.name;
        results(i).params = tc;
        results(i).T_K = NaN;
        results(i).T_C = NaN;
        results(i).time = toc;
        results(i).success = false;
    end
    
    fprintf('\n');
end

%% Parametric sweeps
fprintf('=== PARAMETRIC SWEEPS ===\n\n');

sweeps = struct();

% Reset to baseline
model.param.set('I0', '0.05[A]');
model.param.set('LL_t_TEC', '200[um]');
model.param.set('LL_k_r', '1.15');
model.param.set('q', '1000[W/m^2]');

% Current sweep
fprintf('Current sweep...\n');
I_values = [0.01, 0.02, 0.03, 0.05, 0.075, 0.1, 0.15, 0.2];
T_vs_I = zeros(size(I_values));

for j = 1:length(I_values)
    model.param.set('I0', sprintf('%g[A]', I_values(j)));
    try
        model.study('std1').run();
        T = mphglobal(model, 'comp1.ppb1', 'dataset', 'dset1');
        T_vs_I(j) = T - 273.15;
        fprintf('  I = %.3f A: T = %.1f°C\n', I_values(j), T_vs_I(j));
    catch
        T_vs_I(j) = NaN;
    end
end
sweeps.current.I = I_values;
sweeps.current.T = T_vs_I;

% Reset current
model.param.set('I0', '0.05[A]');

% Heat flux sweep
fprintf('\nHeat flux sweep...\n');
q_values = [100, 250, 500, 750, 1000, 1500, 2000, 2500, 3000];
T_vs_q = zeros(size(q_values));

for j = 1:length(q_values)
    model.param.set('q', sprintf('%g[W/m^2]', q_values(j)));
    try
        model.study('std1').run();
        T = mphglobal(model, 'comp1.ppb1', 'dataset', 'dset1');
        T_vs_q(j) = T - 273.15;
        fprintf('  q = %d W/m²: T = %.1f°C\n', q_values(j), T_vs_q(j));
    catch
        T_vs_q(j) = NaN;
    end
end
sweeps.flux.q = q_values;
sweeps.flux.T = T_vs_q;

% Reset flux
model.param.set('q', '1000[W/m^2]');

% Thickness sweep
fprintf('\nThickness sweep...\n');
t_values = [50, 100, 150, 200, 250, 300, 400, 500];
T_vs_t = zeros(size(t_values));

for j = 1:length(t_values)
    model.param.set('LL_t_TEC', sprintf('%g[um]', t_values(j)));
    try
        model.study('std1').run();
        T = mphglobal(model, 'comp1.ppb1', 'dataset', 'dset1');
        T_vs_t(j) = T - 273.15;
        fprintf('  t = %d um: T = %.1f°C\n', t_values(j), T_vs_t(j));
    catch
        T_vs_t(j) = NaN;
    end
end
sweeps.thickness.t = t_values;
sweeps.thickness.T = T_vs_t;

%% Save all data
fprintf('\n=== SAVING DATA ===\n\n');

save(fullfile(OUTPUT_DIR, 'comsol_results.mat'), 'results', 'sweeps', 'test_cases');
fprintf('Saved MATLAB data to: comsol_results.mat\n');

% Export to CSV
% Test case results
T = table({results.name}', [results.T_C]', ...
    [test_cases.I_current_A]'*1000, [test_cases.thickness_um]', ...
    [test_cases.q_flux_W_m2]', [results.time]', ...
    'VariableNames', {'TestCase', 'T_C', 'I_mA', 'Thickness_um', 'HeatFlux_Wm2', 'Time_s'});
writetable(T, fullfile(OUTPUT_DIR, 'test_results.csv'));

% Sweep data
writetable(table(sweeps.current.I'*1000, sweeps.current.T', 'VariableNames', {'I_mA', 'T_C'}), ...
    fullfile(OUTPUT_DIR, 'sweep_current.csv'));
writetable(table(sweeps.flux.q', sweeps.flux.T', 'VariableNames', {'q_Wm2', 'T_C'}), ...
    fullfile(OUTPUT_DIR, 'sweep_flux.csv'));
writetable(table(sweeps.thickness.t', sweeps.thickness.T', 'VariableNames', {'t_um', 'T_C'}), ...
    fullfile(OUTPUT_DIR, 'sweep_thickness.csv'));

fprintf('Saved CSV files\n');

%% Generate plots
fprintf('\nGenerating plots...\n');

figure('Position', [100, 100, 1400, 900], 'Name', 'COMSOL Results');

% Test case comparison
subplot(2,3,1);
bar([results.T_C]);
set(gca, 'XTickLabel', {results.name}, 'XTickLabelRotation', 45);
ylabel('Temperature (°C)');
title('Test Case Results');
yline(100, 'r--', 'Target');
grid on;

% Current sweep
subplot(2,3,2);
plot(sweeps.current.I*1000, sweeps.current.T, '-o', 'LineWidth', 2, 'MarkerSize', 8);
xlabel('Current (mA)');
ylabel('Temperature (°C)');
title('Effect of TEC Current');
grid on;
[~, idx] = min(sweeps.current.T);
hold on;
plot(sweeps.current.I(idx)*1000, sweeps.current.T(idx), 'r*', 'MarkerSize', 15, 'LineWidth', 2);
text(sweeps.current.I(idx)*1000 + 5, sweeps.current.T(idx), sprintf('Optimal: %.0f mA', sweeps.current.I(idx)*1000));

% Heat flux sweep
subplot(2,3,3);
plot(sweeps.flux.q, sweeps.flux.T, '-o', 'LineWidth', 2, 'MarkerSize', 8);
xlabel('Heat Flux (W/m²)');
ylabel('Temperature (°C)');
title('Effect of Heat Flux');
yline(100, 'r--', 'T = 100°C');
grid on;
% Find max feasible
idx = find(sweeps.flux.T < 100, 1, 'last');
if ~isempty(idx)
    hold on;
    plot(sweeps.flux.q(idx), sweeps.flux.T(idx), 'g*', 'MarkerSize', 15, 'LineWidth', 2);
    text(sweeps.flux.q(idx), sweeps.flux.T(idx) - 5, sprintf('Max: %.0f W/m²', sweeps.flux.q(idx)));
end

% Thickness sweep
subplot(2,3,4);
plot(sweeps.thickness.t, sweeps.thickness.T, '-o', 'LineWidth', 2, 'MarkerSize', 8);
xlabel('TEC Thickness (µm)');
ylabel('Temperature (°C)');
title('Effect of TEC Thickness');
grid on;
[~, idx] = min(sweeps.thickness.T);
hold on;
plot(sweeps.thickness.t(idx), sweeps.thickness.T(idx), 'r*', 'MarkerSize', 15, 'LineWidth', 2);
text(sweeps.thickness.t(idx) + 20, sweeps.thickness.T(idx), sprintf('Optimal: %d µm', sweeps.thickness.t(idx)));

% Performance map (Current vs Flux)
subplot(2,3,5);
% Create a simple contour (would need 2D sweep for real data)
imagesc(sweeps.flux.q, sweeps.current.I*1000, repmat(sweeps.flux.T, length(sweeps.current.I), 1));
xlabel('Heat Flux (W/m²)');
ylabel('Current (mA)');
title('Operating Space (illustrative)');
colorbar;
set(gca, 'YDir', 'normal');

% Summary stats
subplot(2,3,6);
axis off;
text_str = {
    'COMSOL VALIDATION SUMMARY', 
    '', 
    sprintf('Test cases run: %d', length(results)),
    sprintf('Successful: %d', sum([results.success])),
    '',
    sprintf('Best result: %.1f°C', min([results.T_C])),
    sprintf('At q = %d W/m², I = %.0f mA', ...
        results(find([results.T_C] == min([results.T_C]), 1)).params.q_flux_W_m2, ...
        results(find([results.T_C] == min([results.T_C]), 1)).params.I_current_A * 1000),
    '',
    'Optimal parameters:',
    sprintf('  Current: %.0f mA', sweeps.current.I(find(sweeps.current.T == min(sweeps.current.T), 1)) * 1000),
    sprintf('  Thickness: %d µm', sweeps.thickness.t(find(sweeps.thickness.T == min(sweeps.thickness.T), 1))),
    '',
    sprintf('Max flux for T<100°C: %.0f W/m²', ...
        sweeps.flux.q(find(sweeps.flux.T < 100, 1, 'last')))
};
text(0.1, 0.9, text_str, 'VerticalAlignment', 'top', 'FontSize', 11);

saveas(gcf, fullfile(OUTPUT_DIR, 'comsol_analysis.png'));
saveas(gcf, fullfile(OUTPUT_DIR, 'comsol_analysis.fig'));

fprintf('Plots saved\n');

%% Generate report
fid = fopen(fullfile(OUTPUT_DIR, 'comsol_report.txt'), 'w');
fprintf(fid, 'COMSOL DATA EXTRACTION REPORT\n');
fprintf(fid, '=============================\n\n');
fprintf(fid, 'Date: %s\n', datestr(now));
fprintf(fid, 'Model: %s\n\n', COMSOL_MODEL_PATH);

fprintf(fid, '--- TEST CASE RESULTS ---\n');
for i = 1:length(results)
    fprintf(fid, '%d. %s: T = %.1f°C (I=%.0fmA, t=%dum, q=%dW/m²)\n', ...
        i, results(i).name, results(i).T_C, ...
        results(i).params.I_current_A * 1000, ...
        results(i).params.thickness_um, ...
        results(i).params.q_flux_W_m2);
end

fprintf(fid, '\n--- OPTIMAL PARAMETERS ---\n');
fprintf(fid, 'Best current: %.0f mA (T = %.1f°C)\n', ...
    sweeps.current.I(find(sweeps.current.T == min(sweeps.current.T), 1)) * 1000, ...
    min(sweeps.current.T));
fprintf(fid, 'Best thickness: %d µm (T = %.1f°C)\n', ...
    sweeps.thickness.t(find(sweeps.thickness.T == min(sweeps.thickness.T), 1)), ...
    min(sweeps.thickness.T));
fprintf(fid, 'Max feasible flux: %.0f W/m² (T < 100°C)\n', ...
    sweeps.flux.q(find(sweeps.flux.T < 100, 1, 'last')));

fclose(fid);

fprintf('\nReport saved to: %s\n', fullfile(OUTPUT_DIR, 'comsol_report.txt'));

%% Done
fprintf('\n=== DATA EXTRACTION COMPLETE ===\n');
fprintf('All results saved to: %s\n', OUTPUT_DIR);
fprintf('\nTo view the plot, open: comsol_analysis.fig\n');
