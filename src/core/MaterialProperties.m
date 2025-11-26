classdef MaterialProperties
    % MATERIALPROPERTIES Handles material property lookups and temperature dependencies.
    %   Currently implements constant properties but designed for expansion.
    
    properties
        Library
        Tables
    end
    
    methods
        function obj = MaterialProperties(config)
            % Constructor: Loads material library from config struct
            if isfield(config, 'materials')
                obj.Library = config.materials;
            else
                error('Config must contain a "materials" field.');
            end
            obj.load_tables();
        end

        function load_tables(obj)
            obj.Tables = struct();
            % Define paths - assuming data is in project root/data
            % Get the directory of this file (src/core)
            [current_path, ~, ~] = fileparts(mfilename('fullpath'));
            % Go up two levels to project root (src/core -> src -> root)
            project_root = fileparts(fileparts(current_path));
            base_path = fullfile(project_root, 'data', 'material_props');
            
            % Map folder names to property keys
            prop_map = containers.Map(...
                {'thermal_conductivity', 'seebeck_coefficient', 'electrical_resistivity'}, ...
                {'k', 'S', 'rho'});
            
            folders = prop_map.keys;
            for i = 1:length(folders)
                folder = folders{i};
                prop_key = prop_map(folder);
                folder_path = fullfile(base_path, folder);
                
                if isfolder(folder_path)
                    files = dir(fullfile(folder_path, '*.txt'));
                    for j = 1:length(files)
                        fname = files(j).name;
                        [~, mat_name, ~] = fileparts(fname);
                        
                        % Load data (assuming space separated: Temp Value)
                        try
                            data = load(fullfile(folder_path, fname));
                            % Store as struct: Tables.Bi2Te3.k = [T, val]
                            if ~isfield(obj.Tables, mat_name)
                                obj.Tables.(mat_name) = struct();
                            end
                            obj.Tables.(mat_name).(prop_key) = data;
                            % Suppress verbose output for faster optimization
                            % fprintf('Loaded %s for %s\n', prop_key, mat_name);
                        catch ME
                            fprintf('Failed to load %s: %s\n', fname, ME.message);
                        end
                    end
                end
            end
        end
        
        function val = get_k(obj, material_name, T)
            % GET_K Returns thermal conductivity [W/mK]
            % T is temperature in Kelvin (unused for constant props)
            val = obj.get_property(material_name, 'k', T);
        end
        
        function val = get_rho(obj, material_name, T)
            % GET_RHO Returns electrical resistivity [Ohm*m]
            val = obj.get_property(material_name, 'rho', T);
        end

        function val = get_S(obj, material_name, T)
            % GET_S Returns Seebeck coefficient [V/K]
            val = obj.get_property(material_name, 'S', T);
        end
        
        function val = get_property(obj, material_name, prop_name, T)
            % Internal helper to fetch property
            
            % 1. Check for interpolation table
            if isfield(obj.Tables, material_name) && isfield(obj.Tables.(material_name), prop_name)
                table = obj.Tables.(material_name).(prop_name);
                % Linear interpolation, linear extrapolation
                % table(:,1) is Temperature, table(:,2) is Value
                val = interp1(table(:,1), table(:,2), T, 'linear', 'extrap');
                return;
            end

            % 2. Fallback to constant from Library
            if isfield(obj.Library, material_name)
                mat = obj.Library.(material_name);
                if isfield(mat, prop_name)
                    val = mat.(prop_name);
                    % TODO: Add temperature dependence logic here
                    % e.g., if val is a struct with coefficients, calculate polyval
                else
                    % Return 0 for missing properties (e.g. S for Copper)
                    val = 0; 
                end
            else
                error('Material "%s" not found in library.', material_name);
            end
        end
    end
end
