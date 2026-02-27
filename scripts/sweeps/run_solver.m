% Run Solver Script
% Uses default_params.json configuration

format longG;
clc;
clear all;

% Get the directory where this script is located
script_dir = fileparts(mfilename('fullpath'));
project_root = fullfile(script_dir, '..', '..');

addpath(genpath(fullfile(project_root, 'src')));

config_path = fullfile(project_root, 'src', 'config', 'default_params.json');
solver = RadialTECSolver(config_path);
solver.run();
