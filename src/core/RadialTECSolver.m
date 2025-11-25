classdef RadialTECSolver
    properties
        Config
        Geometry
        Materials
        Network
        Results
    end

    methods
        function obj = RadialTECSolver(config_path)
            text = fileread(config_path);
            obj.Config = jsondecode(text);

            obj.Materials = MaterialProperties(obj.Config);
            obj.Geometry = TECGeometry(obj.Config);
            obj.Network = ThermalNetwork(obj.Geometry, obj.Materials, obj.Config);
            obj.Results = ResultsManager(obj.Config);
        end

        function run(obj)
            fprintf('Starting Radial TEC Simulation...\n');
            fprintf('Running Linear Solve (No Iteration)...\n');

            N = obj.Geometry.N_stages;
            dim = 2*N + 1;
            T_guess = ones(dim, 1) * obj.Config.simulation.T_initial_guess;

            [T_new, Q_out, Q_in] = obj.Network.solve(T_guess);

            obj.Results.log_results(T_new, Q_out, Q_in, obj.Config);
            obj.Results.save_results_text(T_new, Q_out, Q_in, obj.Config);
            obj.Results.save_results(T_new, Q_out, obj.Config);
            obj.Results.plot_temperature_profile(T_new, obj.Geometry);

            fprintf('Simulation Complete.\n');
        end
    end
end
