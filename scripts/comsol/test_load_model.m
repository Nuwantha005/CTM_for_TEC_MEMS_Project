%% TEST_LOAD_MODEL - Simple test to diagnose COMSOL model loading issues
% Run this interactively in MATLAB (not batch mode) to see full error messages
%
% Steps:
%   1. Make sure SolidWorks is running with the correct assembly open
%   2. Start COMSOL server: comsolmphserver -port 2036
%   3. Run this script in MATLAB

clear; clc;

%% Configuration
addpath('F:\EngineeringSoftware\COMSOL\COMSOL63\Multiphysics\mli');

COMSOL_PORT = 2036;
COMSOL_MODEL_PATH = 'E:\Semester 7\ME4311 - MicroNano Electro Mechanical Systems and Nanotechnology\Project\COMSOL\Test_Wedge\asym.mph';

%% Check prerequisites
fprintf('=== COMSOL Model Load Test ===\n\n');

fprintf('1. Checking model file...\n');
if exist(COMSOL_MODEL_PATH, 'file')
    d = dir(COMSOL_MODEL_PATH);
    fprintf('   File exists: %.2f MB, modified %s\n', d.bytes/1e6, d.date);
else
    fprintf('   ERROR: File not found!\n');
    return;
end

%% Connect to COMSOL
fprintf('\n2. Connecting to COMSOL server...\n');

% Check if already connected
already_connected = false;
try
    import com.comsol.model.*
    import com.comsol.model.util.*
    ModelUtil.tags();  % This will throw if not connected
    already_connected = true;
    fprintf('   Already connected to COMSOL server.\n');
catch
    % Not connected, try to connect
    fprintf('   Connecting to port %d...\n', COMSOL_PORT);
    try
        mphstart(COMSOL_PORT);
        import com.comsol.model.*
        import com.comsol.model.util.*
        fprintf('   Connected successfully!\n');
    catch ME
        fprintf('   ERROR: %s\n', ME.message);
        fprintf('\n   Make sure COMSOL server is running:\n');
        fprintf('   "F:\\EngineeringSoftware\\COMSOL\\COMSOL63\\Multiphysics\\bin\\win64\\comsolmphserver.exe" -port 2036\n');
        return;
    end
end

%% Load model
fprintf('\n3. Loading COMSOL model (this may take 1-2 minutes with LiveLink)...\n');
fprintf('   Make sure SolidWorks has the linked assembly open!\n');

tic;
try
    model = mphload(COMSOL_MODEL_PATH);
    load_time = toc;
    fprintf('   Model loaded successfully in %.1f seconds!\n', load_time);
catch ME
    load_time = toc;
    fprintf('\n   ERROR after %.1f seconds:\n', load_time);
    fprintf('   ID: %s\n', ME.identifier);
    fprintf('   Message: %s\n', ME.message);
    
    % Print full stack trace
    fprintf('\n   Stack trace:\n');
    for i = 1:length(ME.stack)
        fprintf('   %d: %s (line %d)\n', i, ME.stack(i).name, ME.stack(i).line);
    end
    
    % Print causes if any
    if ~isempty(ME.cause)
        fprintf('\n   Causes:\n');
        for i = 1:length(ME.cause)
            fprintf('   %d: %s\n', i, ME.cause{i}.message);
        end
    end
    return;
end

%% Test parameter access
fprintf('\n4. Testing parameter access...\n');
try
    % List all parameters
    params = model.param.varnames();
    fprintf('   Found %d parameters\n', length(params));
    
    % Check a few key ones
    key_params = {'LL_t_TEC', 'LL_theta', 'LL_L_1', 'LL_k_r'};
    for i = 1:length(key_params)
        pname = key_params{i};
        try
            val = char(model.param.get(pname));
            fprintf('   %s = %s\n', pname, val);
        catch
            fprintf('   %s = (not found)\n', pname);
        end
    end
catch ME
    fprintf('   ERROR: %s\n', ME.message);
end

%% Test setting a parameter
fprintf('\n5. Testing parameter modification...\n');
try
    % Try setting t_TEC to 100 um
    model.param.set('LL_t_TEC', '100');
    fprintf('   Set LL_t_TEC = 100\n');
    
    % Read it back
    val = char(model.param.get('LL_t_TEC'));
    fprintf('   Readback: LL_t_TEC = %s\n', val);
catch ME
    fprintf('   ERROR: %s\n', ME.message);
end

%% Test geometry rebuild
fprintf('\n6. Testing geometry rebuild (triggers SolidWorks update)...\n');
fprintf('   Watch SolidWorks - it should update the model dimensions!\n');

tic;
try
    model.component('comp1').geom('geom1').run();
    geom_time = toc;
    fprintf('   Geometry rebuilt in %.1f seconds!\n', geom_time);
    fprintf('   Check SolidWorks - dimensions should have updated.\n');
catch ME
    geom_time = toc;
    fprintf('   ERROR after %.1f seconds: %s\n', geom_time, ME.message);
end

%% Summary
fprintf('\n=== Test Complete ===\n');
fprintf('If all steps passed, the model is ready for verification runs.\n');
fprintf('You can now run: run_verification_case.m\n');
