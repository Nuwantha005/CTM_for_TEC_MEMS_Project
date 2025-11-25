classdef RadialTECSolver
    % RADIALTECSOLVER Main orchestrator for the simulation.
    %   Manages the Picard iteration loop and coordinates modules.
    
    properties
        Config
        Geometry
        Materials
        Network
        Results
    end
    
    methods
        function obj = RadialTECSolver(config_path)
            % Constructor: Initializes modules
            text = fileread(config_path);
            obj.Config = jsondecode(text);
            
            obj.Materials = MaterialProperties(obj.Config);
            obj.Geometry = TECGeometry(obj.Config);
            obj.Network = ThermalNetwork(obj.Geometry, obj.Materials, obj.Config);
            obj.Results = ResultsManager(obj.Config);
        end
        
        function run(obj)
            % RUN Executes the simulation loop.
            
            fprintf('Starting Radial TEC Simulation...\n');
            
            % 1. Initialization
            N = obj.Geometry.N_stages;
            dim = 2*N + 1;
            T_guess = ones(dim, 1) * obj.Config.simulation.T_initial_guess;
            
            max_iter = obj.Config.simulation.max_iterations;
            tol = obj.Config.simulation.tolerance;
            
            % 2. Single Linear Solve (Debugging)
            fprintf('Running Linear Solve (No Iteration)...\n');
            
            % Solve linear system at initial guess (properties at T_initial)
            [T_new, Q_out] = obj.Network.solve(T_guess);
            
            % Print Temperatures
            fprintf('--- Temperature Results (K) ---\n');
            fprintf('Node 0 (Center): %.2f\n', T_new(1));
            for i = 1:N
                idx_Si = i + 1;
                idx_c = N + i + 1;
                fprintf('Stage %d: T_Si = %.2f, T_c = %.2f\n', i, T_new(idx_Si), T_new(idx_c));
            end
            fprintf('-------------------------------\n');
            fprintf('Q_out calculated: %.2e W\n', Q_out);
            
            % 3. Post-Processing
            obj.Results.save_results_text(T_new, Q_out, obj.Config); % New text logger
            obj.Results.save_results(T_new, Q_out, obj.Config); % Keep binary for now
            obj.Results.plot_temperature_profile(T_new, obj.Geometry);
            
            fprintf('Simulation Complete.\n');
        end
    end
end
