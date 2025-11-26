clc
clear all
addpath(genpath('src')); 
config = jsondecode(fileread('src/config/default_params.json'));

config.boundary_conditions.q_flux_W_m2 = 10000; 
config.geometry.N_tsv_limit = 5; 
config.operating_conditions.I_current_A = 0.02; 
fid = fopen('src/config/test_params.json', 'w'); 
fprintf(fid, '%s', jsonencode(config)); fclose(fid); solver = RadialTECSolver('src/config/test_params.json'); 
result = solver.run();