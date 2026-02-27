%% Verify Thermal Network Topology
% Check if stages are connected in series or parallel
% The 7.5x discrepancy might be due to incorrect topology

clear all; close all; clc;

script_dir = fileparts(mfilename('fullpath'));
project_root = fullfile(script_dir, '..', '..');
addpath(genpath(fullfile(project_root, 'src')));

fprintf('╔═══════════════════════════════════════════════════════════════╗\n');
fprintf('║     THERMAL NETWORK TOPOLOGY ANALYSIS                        ║\n');
fprintf('╚═══════════════════════════════════════════════════════════════╝\n\n');

%% The physical picture:
%
%  CHIP (heat input distributed across entire bottom surface)
%    │
%    ├──── Stage 1 TEC ──── radially ──→ cold side 1
%    │         (at r_in_1 to r_out_1)
%    │
%    ├──── Stage 2 TEC ──── radially ──→ cold side 2  
%    │         (at r_in_2 to r_out_2)
%    │
%    └──── Stage 3 TEC ──── radially ──→ cold side 3 ──→ water
%              (at r_in_3 to r_out_3)
%
% All stages receive heat from the chip DIRECTLY (in parallel)
% Then each stage conducts heat radially to its cold side
% The cold sides are connected to water

fprintf('=== PHYSICAL HEAT FLOW ===\n');
fprintf('Heat enters chip → spreads laterally → enters each TEC stage vertically\n');
fprintf('Each stage conducts heat RADIALLY (r direction) to cold side\n');
fprintf('All stages dump heat to water (in parallel)\n\n');

%% The MATLAB implementation seems to assume:
%
%  CHIP (heat input)
%    │
%    └── Stage 1 hot side ── K_stages(1) ──→ Stage 1 cold side
%                                             │
%                                             │ K_stages(1) connects to
%                                             ↓
%                              Stage 2 cold side ── K_stages(2) ──→ ...
%
% This is SERIES connection of stages!

fprintf('=== MATLAB IMPLEMENTATION (from ThermalNetwork.m) ===\n');
fprintf('Looking at matrix assembly:\n');
fprintf('  M(idx_c(i), idx_c(i+1)) = K_stages(i)  → stages connected in series!\n');
fprintf('  B(idx_c(N)) -= K_stages(N) * T_water   → only last stage to water!\n\n');

%% Let's compute what PARALLEL would give:
[~, ~, ~, ~, ~, CONFIG] = optimization_variables();

config = struct();
config.geometry.N_stages = 3;
config.geometry.wedge_angle_deg = 30;
config.geometry.w_chip_um = 10000;
config.geometry.t_chip_um = 50;
config.geometry.R_cyl_um = 1000;
config.geometry.thickness_um = 1000;
config.geometry.t_ins_um = 10;
config.geometry.radial_expansion_factor = 1.15;
config.geometry.azimuthal_gap_um = 20;
config.geometry.insulation_width_ratio = 0.05;
config.geometry.fill_factor = 0.9;
config.geometry.interconnect_ratio = 0.1;
config.geometry.interconnect_thickness_ratio = 0.6;
config.geometry.interconnect_angle_ratio = 0.5;
config.geometry.outerconnect_ratio = 0.1;
config.geometry.outerconnect_thickness_ratio = 0.6;
config.geometry.outerconnect_angle_ratio = 0.5;
config.materials = CONFIG.materials;

mat = MaterialProperties(config);
geom = TECGeometry(config);

k_TE = 1.6;  % W/mK
theta = geom.WedgeAngle;
t_tec = geom.Thickness;

fprintf('=== THERMAL RESISTANCE COMPARISON ===\n\n');

% Calculate resistance for each stage (radial conduction)
R_stages = zeros(3, 1);
for i = 1:3
    [r_in, L] = geom.get_stage_geometry(i);
    r_out = r_in + L;
    % Radial conduction: R = ln(r_out/r_in) / (k * θ * t)
    R_stages(i) = log(r_out/r_in) / (k_TE * theta * t_tec);
    fprintf('Stage %d: r_in=%.2f mm, r_out=%.2f mm, R=%.1f K/W\n', ...
            i, r_in*1e3, r_out*1e3, R_stages(i));
end

% SERIES: R_total = R_1 + R_2 + R_3
R_series = sum(R_stages);
fprintf('\nSERIES (current MATLAB): R_total = %.1f K/W\n', R_series);

% PARALLEL: 1/R_total = 1/R_1 + 1/R_2 + 1/R_3
R_parallel = 1 / (1/R_stages(1) + 1/R_stages(2) + 1/R_stages(3));
fprintf('PARALLEL (correct?): R_total = %.1f K/W\n', R_parallel);

fprintf('\nCOMSOL reference: R = 85.8 K/W\n');
fprintf('\nRatio (Series/COMSOL) = %.1f\n', R_series / 85.8);
fprintf('Ratio (Parallel/COMSOL) = %.1f\n', R_parallel / 85.8);

%% But wait - the stages don't have same heat input!
% Heat is distributed based on area:
% Q_1 = q * A_1, Q_2 = q * A_2, Q_3 = q * A_3

fprintf('\n=== HEAT DISTRIBUTION ===\n');
q_flux = 200e3;  % W/m²
A_stages = zeros(3, 1);
Q_stages = zeros(3, 1);

for i = 1:3
    [r_in, L] = geom.get_stage_geometry(i);
    r_out = r_in + L;
    A_stages(i) = (theta/2) * (r_out^2 - r_in^2);  % Annular sector area
    Q_stages(i) = q_flux * A_stages(i);
    fprintf('Stage %d: A=%.3e m², Q_in=%.3f W\n', i, A_stages(i), Q_stages(i));
end

Q_total = sum(Q_stages);
fprintf('Total Q_in (stages): %.3f W\n', Q_total);

% Add center cylinder contribution
A_cyl = (theta/2) * geom.R_cyl^2;
Q_cyl = q_flux * A_cyl;
Q_total_all = Q_total + Q_cyl;
fprintf('Center cylinder Q_in: %.3f W\n', Q_cyl);
fprintf('Total Q_in (all): %.3f W\n', Q_total_all);

%% Calculate effective R for each path
fprintf('\n=== EFFECTIVE THERMAL RESISTANCE BY PATH ===\n');
T_water = 300;
dT_stages = zeros(3, 1);
for i = 1:3
    dT_stages(i) = Q_stages(i) * R_stages(i);
    T_chip_i = T_water + dT_stages(i);
    fprintf('Stage %d path: Q=%.3f W, R=%.1f K/W, dT=%.1f K, T_chip=%.1f K\n', ...
            i, Q_stages(i), R_stages(i), dT_stages(i), T_chip_i);
end

% The chip temperature is determined by parallel combination
% T_chip - T_water = Q_total * R_equivalent
% where R_equivalent considers both heat distribution and resistance

% For parallel thermal network with distributed heat:
% Each stage: Q_i, R_i
% T_chip = T_water + Q_i * R_i (same T_chip for all if chip has infinite conductivity)
% But chip has finite conductivity, so T_chip varies!

fprintf('\n=== KEY INSIGHT ===\n');
fprintf('The stages are NOT simply in parallel or series!\n');
fprintf('They form a DISTRIBUTED thermal network with:\n');
fprintf('  - Chip layer providing lateral heat spreading\n');
fprintf('  - Each TEC stage receiving different Q_in\n');
fprintf('  - Each stage having different thermal resistance\n');
fprintf('\nThe MATLAB model DOES include chip spreading (R_lat_Si),\n');
fprintf('but the TEC stage connections might still be wrong.\n');
