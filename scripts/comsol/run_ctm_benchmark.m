%% CTM BENCHMARK RUNNER
% This script runs the CTM against the pre-calculated COMSOL benchmark values
% and evaluates the new accuracy.

clear; clc; close all;
% find root
root_dir = fullfile(fileparts(mfilename('fullpath')), '..', '..');
cd(root_dir);
addpath(genpath('src'));

fprintf('\n');
fprintf('================================================================================\n');
fprintf('           CTM BENCHMARK COMPARISON\n');
fprintf('================================================================================\n\n');

% Define 5 test cases from mathematical_modelling.tex
% Case 1
tc(1).name = 'Thin TEC Low Current';
tc(1).M = 12; tc(1).t_TEC_um = 100; tc(1).t_chip_um = 50; tc(1).t_SOI_um = 20;
tc(1).R_cyl_um = 1000; tc(1).f_L = 1.0; tc(1).fill_factor = 0.9; tc(1).w_is_um = 50;
tc(1).I_A = 0.05; tc(1).q_Wm2 = 500;
tc(1).ic_w_r = 0.1; tc(1).ic_t_r = 0.6; tc(1).ic_angle_r = 0.5;
tc(1).oc_w_r = 0.1; tc(1).oc_t_r = 0.6; tc(1).oc_angle_r = 0.5;
tc(1).Benchmark_T_max = 20.48;

% Test Case 2: Moderate Current
tc(2).name = 'Moderate Current';
tc(2).M = 12; tc(2).t_TEC_um = 150; tc(2).t_chip_um = 50; tc(2).t_SOI_um = 20;
tc(2).R_cyl_um = 1000; tc(2).f_L = 1.15; tc(2).fill_factor = 0.9; tc(2).w_is_um = 50;
tc(2).I_A = 0.10; tc(2).q_Wm2 = 500;
tc(2).ic_w_r = 0.1; tc(2).ic_t_r = 0.6; tc(2).ic_angle_r = 0.5;
tc(2).oc_w_r = 0.1; tc(2).oc_t_r = 0.6; tc(2).oc_angle_r = 0.5;
tc(2).Benchmark_T_max = 20.25;

% Test Case 3: Length Ratio Test
tc(3).name = 'Length Ratio Test';
tc(3).M = 12; tc(3).t_TEC_um = 150; tc(3).t_chip_um = 50; tc(3).t_SOI_um = 20;
tc(3).R_cyl_um = 1000; tc(3).f_L = 1.2; tc(3).fill_factor = 0.9; tc(3).w_is_um = 50;
tc(3).I_A = 0.08; tc(3).q_Wm2 = 500;
tc(3).ic_w_r = 0.1; tc(3).ic_t_r = 0.6; tc(3).ic_angle_r = 0.5;
tc(3).oc_w_r = 0.1; tc(3).oc_t_r = 0.6; tc(3).oc_angle_r = 0.5;
tc(3).Benchmark_T_max = 20.32;

% Test Case 4: High Heat Flux
tc(4).name = 'High Heat Flux';
tc(4).M = 12; tc(4).t_TEC_um = 200; tc(4).t_chip_um = 50; tc(4).t_SOI_um = 20;
tc(4).R_cyl_um = 1000; tc(4).f_L = 1.15; tc(4).fill_factor = 0.9; tc(4).w_is_um = 50;
tc(4).I_A = 0.12; tc(4).q_Wm2 = 2000;
tc(4).ic_w_r = 0.1; tc(4).ic_t_r = 0.6; tc(4).ic_angle_r = 0.5;
tc(4).oc_w_r = 0.1; tc(4).oc_t_r = 0.6; tc(4).oc_angle_r = 0.5;
tc(4).Benchmark_T_max = 22.42;

% Test Case 5: Wide Wedge Angle
tc(5).name = 'Wide Wedge Angle';
tc(5).M = 8; tc(5).t_TEC_um = 150; tc(5).t_chip_um = 50; tc(5).t_SOI_um = 20;
tc(5).R_cyl_um = 1000; tc(5).f_L = 1.15; tc(5).fill_factor = 0.9; tc(5).w_is_um = 50;
tc(5).I_A = 0.08; tc(5).q_Wm2 = 1000;
tc(5).ic_w_r = 0.1; tc(5).ic_t_r = 0.6; tc(5).ic_angle_r = 0.5;
tc(5).oc_w_r = 0.1; tc(5).oc_t_r = 0.6; tc(5).oc_angle_r = 0.5;
tc(5).Benchmark_T_max = 21.33;

fprintf('┌─────┬──────────────────────┬────────────┬────────────┬────────────┬──────────┐\n');
fprintf('│  #  │ Test Case            │ MATLAB [C] │ COMSOL [C] │ Error [C]  │ Error [%%] │\n');
fprintf('├─────┼──────────────────────┼────────────┼────────────┼────────────┼──────────┤\n');

for i = 1:length(tc)
    test = tc(i);
    theta_deg = 360 / test.M;
    
    config = struct();
    config.geometry.N_stages = 3;
    config.geometry.wedge_angle_deg = theta_deg;
    config.geometry.thickness_um = test.t_TEC_um;
    config.geometry.radial_expansion_factor = test.f_L;
    config.geometry.w_chip_um = 10000;
    config.geometry.R_cyl_um = test.R_cyl_um;
    config.geometry.t_chip_um = test.t_chip_um;
    config.geometry.t_ins_um = test.t_SOI_um;
    config.geometry.interconnect_ratio = test.ic_w_r;
    config.geometry.outerconnect_ratio = test.oc_w_r;
    config.geometry.insulation_width_um = test.w_is_um;
    config.geometry.interconnect_angle_ratio = test.ic_angle_r;
    config.geometry.outerconnect_angle_ratio = test.oc_angle_r;
    config.geometry.interconnect_thickness_ratio = test.ic_t_r;
    config.geometry.outerconnect_thickness_ratio = test.oc_t_r;
    config.geometry.fill_factor = test.fill_factor;
    
    config.operating_conditions.I_current_A = test.I_A;
    config.boundary_conditions.q_flux_W_m2 = test.q_Wm2;
    config.boundary_conditions.T_water_K = 293.15;
    config.boundary_conditions.h_conv_W_m2K = 1e6;
    
    config.materials.Bi2Te3 = struct('k', 1.6, 'rho', 1.15e-5, 'S', 210e-6);
    config.materials.Cu = struct('k', 400, 'rho', 1.667e-8);
    config.materials.Si = struct('k', 130, 'rho', 0.01);
    config.materials.AlN = struct('k', 170, 'rho', 1e10);
    config.materials.SiO2 = struct('k', 1.4, 'rho', 1e14);
    config.materials.Al2O3 = struct('k', 35, 'rho', 1e12);
    
    materials = MaterialProperties(config);
    geometry = TECGeometry(config);
    network = ThermalNetwork(geometry, materials, config);
    
    dim = 2 * config.geometry.N_stages + 2;
    T = ones(dim, 1) * 300;
    
    for iter = 1:100
        T_old = T;
        [T, Q_out, Q_in] = network.solve(T);
        if max(abs(T - T_old)) < 1e-6
            break;
        end
    end
    
    T_max = max(T) - 273.15;
    T_comsol = test.Benchmark_T_max;
    err_abs = T_max - T_comsol;
    err_pct = 100 * err_abs / T_comsol;
    
    fprintf('│  %d  │ %-20s │ %10.2f │ %10.2f │ %+10.2f │ %+8.1f │\n', ...
        i, test.name, T_max, T_comsol, err_abs, err_pct);
end
fprintf('└─────┴──────────────────────┴────────────┴────────────┴────────────┴──────────┘\n\n');
