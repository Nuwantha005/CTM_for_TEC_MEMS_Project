% Run Optimization Script

format longG;
clc;
clear all;

addpath(genpath('src'));

% Ensure output directory exists
if ~exist('output', 'dir')
    mkdir('output');
end

optimizer = TECOptimizer('src/config/default_params.json');
optimizer.run_optimization();
