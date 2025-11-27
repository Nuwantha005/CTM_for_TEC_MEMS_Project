function config_path = get_config_path(filename)
% GET_CONFIG_PATH Return absolute path to config file under src/config
%   config_path = GET_CONFIG_PATH() returns path to default_params.json
%   config_path = GET_CONFIG_PATH(filename) returns path to the given file

if nargin < 1 || isempty(filename)
    filename = 'default_params.json';
end

% Determine location of this file so we can work out project root reliably
script_fullpath = which('get_config_path');
if isempty(script_fullpath)
    script_dir = pwd; % fallback to pwd
else
    script_dir = fileparts(script_fullpath);
end

% script_dir is .../src/utils -> project_root = fileparts(fileparts(script_dir))
project_root = fileparts(fileparts(script_dir));
src_dir = fullfile(project_root, 'src');

config_path = fullfile(src_dir, 'config', filename);
end
