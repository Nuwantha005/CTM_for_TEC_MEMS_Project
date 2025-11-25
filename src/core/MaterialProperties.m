classdef MaterialProperties
    % MATERIALPROPERTIES Handles material property lookups and temperature dependencies.
    %   Currently implements constant properties but designed for expansion.
    
    properties
        Library
    end
    
    methods
        function obj = MaterialProperties(config)
            % Constructor: Loads material library from config struct
            if isfield(config, 'materials')
                obj.Library = config.materials;
            else
                error('Config must contain a "materials" field.');
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
