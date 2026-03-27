clear; clc; close all;
% find root
root_dir = fullfile(fileparts(mfilename('fullpath')), '..', '..');
cd(root_dir);
addpath(genpath('src'));

config = struct();
config.geometry.N_stages = 3;
config.geometry.wedge_angle_deg = 30;
config.geometry.thickness_um = 100;
config.geometry.radial_expansion_factor = 1.0;
config.geometry.w_chip_um = 10000;
config.geometry.R_cyl_um = 1000;
config.geometry.t_chip_um = 50;
config.geometry.t_ins_um = 20;
config.geometry.interconnect_ratio = 0.1;
config.geometry.outerconnect_ratio = 0.1;
config.geometry.insulation_width_um = 50;
config.geometry.interconnect_angle_ratio = 0.5;
config.geometry.outerconnect_angle_ratio = 0.5;
config.geometry.interconnect_thickness_ratio = 0.6;
config.geometry.outerconnect_thickness_ratio = 0.6;
config.geometry.fill_factor = 0.9;

config.operating_conditions.I_current_A = 0.05;
config.boundary_conditions.q_flux_W_m2 = 500;
config.boundary_conditions.T_water_K = 300;
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

T_current = ones(8, 1) * 300;
[M, B] = network.assemble_system(T_current);
disp('Matrix M:');
disp(M);
disp('Vector B:');
disp(B);
T_new = M \ B
disp('Delta T:');
disp(T_new - 300)
