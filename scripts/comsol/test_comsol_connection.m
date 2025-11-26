%% TEST COMSOL LIVELINK CONNECTION
% Run this script to verify COMSOL LiveLink is working
%
% Prerequisites:
%   1. COMSOL 6.3 installed with LiveLink for MATLAB
%   2. COMSOL server running OR will start automatically
%
% Steps to start COMSOL server manually:
%   Option 1: Command line
%     > comsolmphserver -port 2036
%   Option 2: From COMSOL Desktop
%     File > Client/Server > Start Server
%   Option 3: COMSOL with MATLAB (auto-connect)
%     > comsol mphserver matlab

clear; clc;

%% Add paths
addpath(genpath('src'));

%% Check installation
fprintf('=== CHECKING COMSOL LIVELINK INSTALLATION ===\n\n');
COMSOLInterface.checkInstallation();

%% Try to connect
fprintf('\n=== TESTING CONNECTION ===\n\n');

% Create interface
comsol = COMSOLInterface();

% Try to connect
if comsol.connect()
    fprintf('\n✓ CONNECTION SUCCESSFUL!\n\n');
    
    % If you have a model, try loading it
    % Uncomment and modify the path:
    % modelPath = 'path/to/your/model.mph';
    % if exist(modelPath, 'file')
    %     comsol.loadModel(modelPath);
    % end
    
    % Disconnect
    comsol.disconnect();
else
    fprintf('\n✗ CONNECTION FAILED\n\n');
    fprintf('Please ensure:\n');
    fprintf('1. COMSOL is installed\n');
    fprintf('2. LiveLink for MATLAB is licensed\n');
    fprintf('3. COMSOL server is running (port 2036)\n');
    fprintf('\nTo add LiveLink to path:\n');
    fprintf('  addpath(''C:\\Program Files\\COMSOL\\COMSOL63\\Multiphysics\\mli'')\n');
end
