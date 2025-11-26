% Run Parametric Sweep Script

format longG;
clc;
clear all;

addpath(genpath('src'));

% Ensure output directory exists
if ~exist('output', 'dir')
    mkdir('output');
end

sweeper = ParametricSweeper('src/config/default_params.json');

% Example: Sweep Current from 0.1A to 5A in 20 steps
sweeper.run_sweep('current', 0.0, 0.2, 5);
%sweeper.run_sweep('interconnect_ratio', 0, 0.5, 5);
%sweeper.run_sweep('thickness_um',50,500,10)


% Example: Sweep Radial Expansion Factor
% sweeper.run_sweep('k_r', 1.0, 2.0, 10);
