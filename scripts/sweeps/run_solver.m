% Run Solver Script

format longG;
clc;
clear all;

addpath(genpath('src'));
solver = RadialTECSolver('src/config/default_params.json');
solver.run();
