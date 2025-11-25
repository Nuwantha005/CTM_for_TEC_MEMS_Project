% Run Solver Script
addpath(genpath('src'));
solver = RadialTECSolver('src/config/default_params.json');
solver.run();
