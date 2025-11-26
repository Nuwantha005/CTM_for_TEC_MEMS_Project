% prepare_comsol_candidates.m
% Generate and evaluate candidate designs for COMSOL validation
% Run this BEFORE connecting to COMSOL to prepare the test cases
%
% This script creates candidate designs that are compatible with
% the 3-stage, TSV-in-stage-1-only COMSOL template.

clear; clc;
addpath(genpath('src'));

fprintf('=== PREPARING COMSOL CANDIDATES ===\n');
fprintf('Generating designs for 3-stage TEC with TSV in stage 1 only\n\n');

%% Output directory
timestamp = datestr(now, 'yyyy-mm-dd_HH-MM-SS');
OUTPUT_DIR = fullfile('output', 'comsol_candidates', timestamp);
if ~exist(OUTPUT_DIR, 'dir')
    mkdir(OUTPUT_DIR);
end

%% Configuration constraints (matching COMSOL template)
fprintf('Template constraints:\n');
fprintf('  - 3 stages\n');
fprintf('  - TSV only in stage 1\n');
fprintf('  - 30° wedge angle (12 wedges)\n');
fprintf('  - 10mm × 10mm chip\n\n');

% Fixed parameters
FIXED = struct();
FIXED.N_stages = 3;
FIXED.N_tsv_limit = 1;
FIXED.wedge_angle_deg = 30;
FIXED.w_chip_um = 10000;
FIXED.R_cyl_um = 1000;
FIXED.t_chip_um = 50;
FIXED.T_water_K = 300;
FIXED.h_conv_W_m2K = 1e6;

%% Design space for parameter sweep
% Based on our optimization analysis:

% Current: Found optimal around 20-100 mA
I_range = [0.02, 0.03, 0.05, 0.075, 0.1, 0.15];  % A

% TEC thickness: Thicker helps (200-300 um optimal)
t_range = [50, 100, 150, 200, 250, 300];  % um

% Radial expansion factor: 0.8-1.2 range
kr_range = [0.9, 1.0, 1.15];

% Fill factor: Higher is better
ff_range = [0.9, 0.95];

% Heat flux: Test feasible range
q_range = [500, 1000, 1500, 2000];  % W/m²

%% Generate all combinations
fprintf('Generating candidate combinations...\n');

candidates = [];
id = 0;

for I = I_range
    for t = t_range
        for kr = kr_range
            for ff = ff_range
                for q = q_range
                    id = id + 1;
                    
                    c = struct();
                    c.id = id;
                    
                    % Fixed params
                    c.geometry.N_stages = FIXED.N_stages;
                    c.geometry.N_tsv_limit = FIXED.N_tsv_limit;
                    c.geometry.wedge_angle_deg = FIXED.wedge_angle_deg;
                    c.geometry.w_chip_um = FIXED.w_chip_um;
                    c.geometry.R_cyl_um = FIXED.R_cyl_um;
                    c.geometry.t_chip_um = FIXED.t_chip_um;
                    
                    % Variable params
                    c.geometry.thickness_um = t;
                    c.geometry.radial_expansion_factor = kr;
                    c.geometry.fill_factor = ff;
                    c.geometry.interconnect_ratio = 0.15;
                    c.geometry.outerconnect_ratio = 0.15;
                    c.geometry.insulation_width_ratio = 0.04;
                    c.geometry.interconnect_angle_ratio = 0.16;
                    c.geometry.outerconnect_angle_ratio = 0.16;
                    c.geometry.interconnect_thickness_ratio = 1.0;
                    c.geometry.outerconnect_thickness_ratio = 1.0;
                    
                    % TSV parameters
                    c.geometry.tsv.R_TSV_um = 10;
                    c.geometry.tsv.P_TSV_um = 20;
                    c.geometry.tsv.g_rad_um = 10;
                    c.geometry.tsv.t_SOI_um = 100;
                    
                    % Operating conditions
                    c.operating_conditions.I_current_A = I;
                    
                    % Boundary conditions
                    c.boundary_conditions.q_flux_W_m2 = q;
                    c.boundary_conditions.T_water_K = FIXED.T_water_K;
                    c.boundary_conditions.h_conv_W_m2K = FIXED.h_conv_W_m2K;
                    
                    % Simulation parameters
                    c.simulation.max_iterations = 100;
                    c.simulation.tolerance = 1e-6;
                    c.simulation.T_ambient = 300;
                    c.simulation.T_initial_guess = 300;
                    
                    % Materials (simplified)
                    c.materials.Bi2Te3 = struct('k', 1.2, 'rho', 1e-5, 'S', 0.0002);
                    c.materials.Cu = struct('k', 400, 'rho', 1.7e-8);
                    c.materials.Si = struct('k', 150, 'rho', 0.01);
                    c.materials.AlN = struct('k', 170, 'rho', 1e10);
                    c.materials.SiO2 = struct('k', 1.4, 'rho', 1e14);
                    c.materials.Al2O3 = struct('k', 30, 'rho', 1e12);
                    
                    % Label
                    c.label = sprintf('I%03.0f_t%03d_kr%.2f_ff%.2f_q%04d', ...
                        I*1000, t, kr, ff, q);
                    
                    candidates = [candidates; c];
                end
            end
        end
    end
end

fprintf('Generated %d candidate designs\n\n', length(candidates));

%% Run MATLAB simulations for all candidates
fprintf('=== MATLAB PRE-SCREENING ===\n\n');

materials = MaterialProperties(candidates(1));
results = [];

progress_interval = ceil(length(candidates) / 20);

for i = 1:length(candidates)
    c = candidates(i);
    
    if mod(i, progress_interval) == 0 || i == 1
        fprintf('Processing candidate %d/%d (%.0f%%)...\n', i, length(candidates), 100*i/length(candidates));
    end
    
    try
        % Create geometry and network
        geometry = TECGeometry(c);
        network = ThermalNetwork(geometry, materials, c);
        
        % Solve iteratively
        N = geometry.N_stages;
        dim = 2*N + 1;
        T = ones(dim, 1) * 300;
        
        converged = false;
        for iter = 1:100
            T_old = T;
            [T, Q_out, Q_in] = network.solve(T);
            
            if max(abs(T - T_old)) < 1e-5
                converged = true;
                break;
            end
        end
        
        % Store results
        r = struct();
        r.id = c.id;
        r.label = c.label;
        r.I_mA = c.operating_conditions.I_current_A * 1000;
        r.t_um = c.geometry.thickness_um;
        r.k_r = c.geometry.radial_expansion_factor;
        r.ff = c.geometry.fill_factor;
        r.q_Wm2 = c.boundary_conditions.q_flux_W_m2;
        
        T_max = max(T);
        T_center = T(1);
        
        r.T_max_K = T_max;
        r.T_max_C = T_max - 273.15;
        r.T_center_K = T_center;
        r.T_center_C = T_center - 273.15;
        r.converged = converged;
        r.feasible = (r.T_max_C > 0) && (r.T_max_C < 300);  % Physical
        r.meets_target = r.T_max_C < 100;
        
        r.T_profile = T;
        r.Q_in = Q_in;
        r.Q_out = Q_out;
        
        % Store reference to full config
        r.config = c;
        
        results = [results; r];
        
    catch ME
        % Skip failed simulations
        r = struct();
        r.id = c.id;
        r.label = c.label;
        r.T_max_C = NaN;
        r.feasible = false;
        r.meets_target = false;
        r.error = ME.message;
        results = [results; r];
    end
end

fprintf('\nSimulation complete!\n\n');

%% Analyze results
fprintf('=== RESULTS ANALYSIS ===\n\n');

feasible = results([results.feasible]);
meeting_target = results([results.meets_target]);

fprintf('Total candidates: %d\n', length(results));
fprintf('Feasible (physical): %d (%.1f%%)\n', length(feasible), 100*length(feasible)/length(results));
fprintf('Meeting target (T<100°C): %d (%.1f%%)\n', length(meeting_target), 100*length(meeting_target)/length(results));

if ~isempty(meeting_target)
    [~, best_idx] = min([meeting_target.T_max_C]);
    best = meeting_target(best_idx);
    
    fprintf('\n--- BEST DESIGN ---\n');
    fprintf('Label: %s\n', best.label);
    fprintf('T_max: %.1f°C\n', best.T_max_C);
    fprintf('T_center: %.1f°C\n', best.T_center_C);
    fprintf('Current: %.0f mA\n', best.I_mA);
    fprintf('Thickness: %d µm\n', best.t_um);
    fprintf('k_r: %.2f\n', best.k_r);
    fprintf('Fill factor: %.2f\n', best.ff);
    fprintf('Heat flux: %d W/m²\n', best.q_Wm2);
end

%% Select top candidates for COMSOL
fprintf('\n=== TOP CANDIDATES FOR COMSOL ===\n\n');

% Sort feasible by temperature
if ~isempty(feasible)
    [~, sort_idx] = sort([feasible.T_max_C]);
    top_N = 15;  % Select top 15
    
    if length(feasible) < top_N
        top_candidates = feasible(sort_idx);
    else
        top_candidates = feasible(sort_idx(1:top_N));
    end
    
    fprintf('%-4s | %-25s | %8s | %6s | %5s | %5s | %6s\n', ...
        'Rank', 'Label', 'T_max(°C)', 'I(mA)', 't(um)', 'k_r', 'q(W/m²)');
    fprintf('%s\n', repmat('-', 1, 75));
    
    for i = 1:length(top_candidates)
        tc = top_candidates(i);
        fprintf('%4d | %-25s | %8.1f | %6.0f | %5d | %5.2f | %6d\n', ...
            i, tc.label, tc.T_max_C, tc.I_mA, tc.t_um, tc.k_r, tc.q_Wm2);
    end
    
    % Save top candidates
    save(fullfile(OUTPUT_DIR, 'top_candidates.mat'), 'top_candidates');
else
    fprintf('No feasible candidates found!\n');
    top_candidates = [];
end

%% Generate COMSOL parameter file for each top candidate
fprintf('\n=== GENERATING COMSOL PARAMETER FILES ===\n\n');

comsol_params_dir = fullfile(OUTPUT_DIR, 'comsol_params');
if ~exist(comsol_params_dir, 'dir')
    mkdir(comsol_params_dir);
end

for i = 1:length(top_candidates)
    tc = top_candidates(i);
    c = tc.config;
    
    filename = fullfile(comsol_params_dir, sprintf('params_%02d_%s.txt', i, tc.label));
    fid = fopen(filename, 'w');
    
    fprintf(fid, '%% COMSOL Parameters for Candidate %d\n', i);
    fprintf(fid, '%% Label: %s\n', tc.label);
    fprintf(fid, '%% MATLAB T_max: %.1f°C\n\n', tc.T_max_C);
    
    % Format for COMSOL LiveLink
    fprintf(fid, 'I0 %g[A]\n', c.operating_conditions.I_current_A);
    fprintf(fid, 'LL_t_TEC %g[um]\n', c.geometry.thickness_um);
    fprintf(fid, 'LL_k_r %g\n', c.geometry.radial_expansion_factor);
    fprintf(fid, 'LL_theta %g[deg]\n', c.geometry.wedge_angle_deg);
    fprintf(fid, 'q %g[W/m^2]\n', c.boundary_conditions.q_flux_W_m2);
    fprintf(fid, 'LL_R_cyl %g[um]\n', c.geometry.R_cyl_um);
    fprintf(fid, 'LL_t_chip %g[um]\n', c.geometry.t_chip_um);
    
    fclose(fid);
end

fprintf('Created %d parameter files in: %s\n', length(top_candidates), comsol_params_dir);

%% Save all results
save(fullfile(OUTPUT_DIR, 'all_candidates.mat'), 'candidates');
save(fullfile(OUTPUT_DIR, 'all_results.mat'), 'results');
save(fullfile(OUTPUT_DIR, 'feasible_results.mat'), 'feasible');

% Export to CSV
if ~isempty(feasible)
    T = table([feasible.id]', {feasible.label}', [feasible.T_max_C]', ...
        [feasible.I_mA]', [feasible.t_um]', [feasible.k_r]', ...
        [feasible.ff]', [feasible.q_Wm2]', [feasible.meets_target]', ...
        'VariableNames', {'ID', 'Label', 'T_max_C', 'I_mA', 't_um', 'k_r', 'ff', 'q_Wm2', 'MeetsTarget'});
    writetable(T, fullfile(OUTPUT_DIR, 'feasible_designs.csv'));
end

%% Generate plots
fprintf('\n=== GENERATING PLOTS ===\n\n');

figure('Position', [50, 50, 1400, 900], 'Name', 'Candidate Analysis');

% Temperature distribution
subplot(2,3,1);
if ~isempty(feasible)
    histogram([feasible.T_max_C], 20);
    xlabel('T_{max} (°C)');
    ylabel('Count');
    title('Temperature Distribution (Feasible Designs)');
    xline(100, 'r--', 'Target', 'LineWidth', 2);
end

% Temperature vs Current
subplot(2,3,2);
if ~isempty(feasible)
    scatter([feasible.I_mA], [feasible.T_max_C], 40, [feasible.q_Wm2], 'filled');
    xlabel('Current (mA)');
    ylabel('T_{max} (°C)');
    title('Temperature vs Current');
    colorbar; ylabel(colorbar, 'Heat Flux (W/m²)');
    yline(100, 'r--');
end

% Temperature vs Thickness
subplot(2,3,3);
if ~isempty(feasible)
    scatter([feasible.t_um], [feasible.T_max_C], 40, [feasible.I_mA], 'filled');
    xlabel('TEC Thickness (µm)');
    ylabel('T_{max} (°C)');
    title('Temperature vs Thickness');
    colorbar; ylabel(colorbar, 'Current (mA)');
    yline(100, 'r--');
end

% Temperature vs Heat Flux
subplot(2,3,4);
if ~isempty(feasible)
    scatter([feasible.q_Wm2], [feasible.T_max_C], 40, [feasible.t_um], 'filled');
    xlabel('Heat Flux (W/m²)');
    ylabel('T_{max} (°C)');
    title('Temperature vs Heat Flux');
    colorbar; ylabel(colorbar, 'Thickness (µm)');
    yline(100, 'r--');
end

% Optimal current by flux level
subplot(2,3,5);
for q = unique([feasible.q_Wm2])
    subset = feasible([feasible.q_Wm2] == q);
    if ~isempty(subset)
        % Find optimal current for this flux
        [T_min, idx] = min([subset.T_max_C]);
        I_opt = subset(idx).I_mA;
        scatter(q, I_opt, 100, 'filled');
        hold on;
        text(q + 50, I_opt, sprintf('%.0f°C', T_min), 'FontSize', 8);
    end
end
xlabel('Heat Flux (W/m²)');
ylabel('Optimal Current (mA)');
title('Optimal Current by Heat Flux');
grid on;

% Parameter space coverage
subplot(2,3,6);
if ~isempty(top_candidates)
    scatter3([top_candidates.I_mA], [top_candidates.t_um], [top_candidates.q_Wm2], ...
        80, [top_candidates.T_max_C], 'filled');
    xlabel('Current (mA)');
    ylabel('Thickness (µm)');
    zlabel('Heat Flux (W/m²)');
    title('Top Candidates in Parameter Space');
    colorbar; ylabel(colorbar, 'T_{max} (°C)');
end

saveas(gcf, fullfile(OUTPUT_DIR, 'candidate_analysis.png'));
saveas(gcf, fullfile(OUTPUT_DIR, 'candidate_analysis.fig'));

fprintf('Plots saved\n');

%% Summary
fprintf('\n=== SUMMARY ===\n');
fprintf('Output directory: %s\n', OUTPUT_DIR);
fprintf('Files created:\n');
fprintf('  - all_candidates.mat (all %d designs)\n', length(candidates));
fprintf('  - all_results.mat (MATLAB simulation results)\n');
fprintf('  - feasible_results.mat (%d feasible designs)\n', length(feasible));
fprintf('  - top_candidates.mat (%d for COMSOL)\n', length(top_candidates));
fprintf('  - feasible_designs.csv\n');
fprintf('  - comsol_params/*.txt (COMSOL parameter files)\n');
fprintf('  - candidate_analysis.png/fig\n');

fprintf('\n=== NEXT STEPS ===\n');
fprintf('1. Start COMSOL server:\n');
fprintf('   F:\\EngineeringSoftware\\COMSOL\\COMSOL63\\Multiphysics\\bin\\win64\\comsolmphserver -port 2036\n\n');
fprintf('2. Run extract_comsol_data.m to validate top candidates in COMSOL\n');
fprintf('   Or run run_comsol_validation.m for full automation\n');
