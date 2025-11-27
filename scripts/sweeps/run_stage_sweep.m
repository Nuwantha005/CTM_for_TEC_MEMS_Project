% Run Stage Optimization Sweep

format longG;
clc;
clear all;

addpath(genpath('src'));

% Ensure output directory exists
if ~exist('output', 'dir')
    mkdir('output');
end

sweeper = StageSweeper('src/config/default_params.json');

% Sweep from 2 to 12 stages
sweeper.run_stage_sweep(2, 12);
