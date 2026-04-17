classdef ThermalNetwork
    properties
        Geometry
        Materials
        Params
    end

    methods
        function obj = ThermalNetwork(geometry, materials, config)
            obj.Geometry = geometry;
            obj.Materials = materials;
            obj.Params = config;
        end

        function [T_new, Q_out, Q_in] = solve(obj, T_current)
            N = obj.Geometry.N_stages;
            dim = 2*N + 1;
            if length(T_current) < dim
                T_current_new = ones(dim, 1) * T_current(1);
                T_current_new(1:length(T_current)) = T_current;
                T_current = T_current_new;
            end
            
            [M, B] = obj.assemble_system(T_current);
            T_new = M \ B;
            
            idx_c_N = 2*N + 1;
            T_water = obj.Params.boundary_conditions.T_water_K;
            T_cold_N = T_new(idx_c_N);
            stage_N = obj.compute_stage_terms(N, T_cold_N, T_water);

            % Boundary heat rejection consistent with last TEC row coupling.
            Q_out = stage_N.K_stage * (T_cold_N - T_water);


            % n_wedges = 2*pi / obj.Geometry.WedgeAngle;
            % Q_out is already total for the wedge

            % Calculate total input heat (Simulated)
            % Note: The Latex model assigns heat generation based on TEC element areas.
            % This ignores flux falling on radial insulators.
            % To maintain energy balance consistency in results, we report the actual heat
            % injected into the nodes, not the theoretical full-area flux.

            q_flux = obj.Params.boundary_conditions.q_flux_W_m2;
            theta = obj.Geometry.WedgeAngle;
            R_cyl = obj.Geometry.R_cyl;

            r_in_arr = zeros(N, 1);
            r_out_arr = zeros(N, 1);
            for i = 1:N
                [r_in_st, L_st, ~, ~, ~, ~, ~, ~, ~, ~] = obj.Geometry.get_stage_geometry(i);
                r_in_arr(i) = r_in_st;
                r_out_arr(i) = r_in_st + L_st;
            end

            A_cv = zeros(N, 1);
            for i = 1:N
                if i == 1
                    r_inner_cv = R_cyl;
                else
                    r_inner_cv = (r_in_arr(i-1) + r_in_arr(i)) / 2;
                end

                if i == N
                    r_outer_cv = r_out_arr(N);
                else
                    r_outer_cv = (r_in_arr(i) + r_in_arr(i+1)) / 2;
                end
                A_cv(i) = 0.5 * theta * (r_outer_cv^2 - r_inner_cv^2);
            end

            A_cyl = 0.5 * theta * R_cyl^2;
            Q_in = q_flux * (A_cyl + sum(A_cv));
        end

        function [M, B] = assemble_system(obj, T_current)
            N = obj.Geometry.N_stages;
            dim = 2*N + 1;
            M = zeros(dim, dim); % spalloc(dim, dim, 5*dim);
            B = zeros(dim, 1);

            idx_0 = 1;
            idx_Si_start = 2;
            idx_c_start = N + 2;

            I = obj.Params.operating_conditions.I_current_A;
            T_water = obj.Params.boundary_conditions.T_water_K;
            q_flux = obj.Params.boundary_conditions.q_flux_W_m2;

            theta = obj.Geometry.WedgeAngle;
            t_chip = obj.Params.geometry.t_chip_um * 1e-6;
            t_tec = obj.Geometry.Thickness;

            % Get insulation width from geometry - support both ratio and absolute
            [~, ~, ~, ~, ~, ~, ~, ~, ~, w_is] = obj.Geometry.get_stage_geometry(1);

            R_cyl = obj.Geometry.R_cyl;

            if isfield(obj.Params.geometry, 't_ins_um')
                t_ins = obj.Params.geometry.t_ins_um * 1e-6;
            else
                t_ins = 10e-6; 
            end

            % Pre-compute Control Volume areas for nodes
            r_in_arr = zeros(N, 1);
            r_out_arr = zeros(N, 1);
            for i = 1:N
                [r_in_st, L_st, ~, ~, ~, ~, ~, ~, ~, ~] = obj.Geometry.get_stage_geometry(i);
                r_in_arr(i) = r_in_st;
                r_out_arr(i) = r_in_st + L_st;
            end
            
            A_cv = zeros(N, 1);
            for i = 1:N
                if i == 1
                    r_inner_cv = R_cyl;
                else
                    r_inner_cv = (r_in_arr(i-1) + r_in_arr(i)) / 2;
                end
                
                if i == N
                    r_outer_cv = r_out_arr(N);
                else
                    r_outer_cv = (r_in_arr(i) + r_in_arr(i+1)) / 2;
                end
                A_cv(i) = 0.5 * theta * (r_outer_cv^2 - r_inner_cv^2);
            end

            K_stages = zeros(N, 1);
            S_stages = zeros(N, 1);
            Re_stages_leg = zeros(N, 1);
            R_ic_stages = zeros(N, 1);
            R_oc_stages = zeros(N, 1);
            R_lat_Si = zeros(N, 1);
            R_vert = zeros(N, 1);
            Q_gen_nodes = zeros(N, 1);

            for i = 1:N
                [r_in, L, ~, ~, ~, ~, ~, ~, ~, ~] = obj.Geometry.get_stage_geometry(i);

                T_cold = T_current(idx_c_start + i - 1);
                if i < N
                    T_hot = T_current(idx_c_start + i);
                else
                    T_hot = T_water;
                end
                stage_i = obj.compute_stage_terms(i, T_cold, T_hot);

                K_stages(i) = stage_i.K_stage;
                S_stages(i) = stage_i.S_stage;
                Re_stages_leg(i) = stage_i.Re_leg;
                R_ic_stages(i) = stage_i.R_ic;
                R_oc_stages(i) = stage_i.R_oc;

                k_Si = obj.Materials.get_k('Si', T_current(idx_Si_start + i - 1));

                if i < N
                    [r_in_next, L_next] = obj.Geometry.get_stage_geometry(i+1);
                    r_mid_i = r_in + L/2;
                    r_mid_next =  r_in_next + L_next/2;
                    R_lat_Si(i) = log(r_mid_next/r_mid_i) / (k_Si * theta * t_chip);
                else
                    R_lat_Si(i) = inf;
                end

                A_top = A_cv(i);
                R_vert(i) = (2 * t_ins) / (stage_i.k_is * A_top);
                Q_gen_nodes(i) = q_flux * A_top;
            end

            A_cyl = 0.5 * theta * R_cyl^2;
            Q_gen_0 = q_flux * A_cyl;

            k_Si_0 = obj.Materials.get_k('Si', T_current(idx_0));
            R_Si_01 = 1 / (2 * theta * k_Si_0 * t_chip);

            k_is_0 = obj.Materials.get_k('AlN', T_current(idx_0));
            term1 = 1 / (2 * theta * k_Si_0 * t_tec);
            term2 = (1 / (k_is_0 * t_tec * theta)) * log((R_cyl + w_is) / R_cyl);
            R_TEC_01 = term1 + term2;

            M(idx_0, idx_0) = -(1/R_Si_01 + 1/R_TEC_01);
            M(idx_0, idx_Si_start) = 1/R_Si_01;
            M(idx_0, idx_c_start) = 1/R_TEC_01;
            B(idx_0) = -Q_gen_0;

            for i = 1:N
                idx_Si = idx_Si_start + i - 1;
                idx_c = idx_c_start + i - 1;

                if i == 1
                    M(idx_Si, idx_0) = 1/R_Si_01;
                    G_left = 1/R_Si_01;
                else
                    idx_Si_prev = idx_Si - 1;
                    M(idx_Si, idx_Si_prev) = +1/R_lat_Si(i-1); %% should be +
                    G_left = 1/R_lat_Si(i-1);
                end

                if i == N
                    G_right = 0;
                else
                    idx_Si_next = idx_Si + 1;
                    M(idx_Si, idx_Si_next) = +1/R_lat_Si(i); % should be +
                    G_right = 1/R_lat_Si(i);
                end

                G_vert = 1/R_vert(i);
                % FIXED: Was -G_vert, should be +G_vert
                M(idx_Si, idx_c) = G_vert;
                
                % FIXED: Added negative sign
                M(idx_Si, idx_Si) = -(G_left + G_right + G_vert);
                B(idx_Si) = B(idx_Si) - Q_gen_nodes(i);

                M(idx_c, idx_Si) = G_vert;

                if i == 1
                    M(idx_c, idx_0) = 1/R_TEC_01;
                end

                if i > 1
                    idx_c_prev = idx_c - 1;
                    M(idx_c, idx_c_prev) = K_stages(i-1);
                end

                if i < N
                    idx_c_next = idx_c + 1;
                    M(idx_c, idx_c_next) = K_stages(i); %+K
                end

                B(idx_c) = B(idx_c) - I^2 * (Re_stages_leg(i)/2 + R_ic_stages(i));
                if i > 1
                    B(idx_c) = B(idx_c) - I^2 * (Re_stages_leg(i-1)/2 + R_oc_stages(i-1));
                end

                sum_diag = G_vert + S_stages(i)*I + K_stages(i);
                
                if i > 1
                    sum_diag = sum_diag + K_stages(i-1) - S_stages(i-1)*I;
                end
                if i == 1
                    sum_diag = sum_diag + 1/R_TEC_01;
                end

                if i == N
                    B(idx_c) = B(idx_c) - K_stages(i) * T_water;
                end

                M(idx_c, idx_c) = -sum_diag;
            end

            % DEBUG: Display full matrix for N=3
            % if N == 3
            %     fprintf('\n=== DEBUG: Full Matrix M (7x7) ===\n');
            %     disp(full(M));
            %     fprintf('\n=== DEBUG: RHS Vector B ===\n');
            %     disp(B);
            % end
        end

        function stage = compute_stage_terms(obj, stage_idx, T_cold, T_hot)
            [r_in, L, w_ic, t_ic, beta_ic, w_oc, t_oc, beta_oc, w_az, w_is_stage] = obj.Geometry.get_stage_geometry(stage_idx);
            T_avg = (T_cold + T_hot) / 2;

            k_p = obj.Materials.get_k('Bi2Te3', T_avg);
            k_n = obj.Materials.get_k('Bi2Te3', T_avg);
            rho_p = obj.Materials.get_rho('Bi2Te3', T_avg);
            rho_n = obj.Materials.get_rho('Bi2Te3', T_avg);
            S_p = obj.Materials.get_S('Bi2Te3', T_avg);
            S_n = obj.Materials.get_S('Bi2Te3', T_avg);

            rho_c = obj.Materials.get_rho('Cu', T_avg);
            k_is = obj.Materials.get_k('AlN', T_avg);
            k_az = obj.Materials.get_k('SiO2', T_avg);
            k_Cu = obj.Materials.get_k('Cu', T_avg);

            R_TE_I_p = obj.Geometry.calculate_R_TE_I(r_in, w_ic, t_ic, beta_ic, w_az, k_p);
            R_TE_I_n = obj.Geometry.calculate_R_TE_I(r_in, w_ic, t_ic, beta_ic, w_az, k_n);
            R_TE_II_p = obj.Geometry.calculate_R_TE_II(r_in, L, w_ic, w_oc, w_az, k_p);
            R_TE_II_n = obj.Geometry.calculate_R_TE_II(r_in, L, w_ic, w_oc, w_az, k_n);
            R_TE_III_p = obj.Geometry.calculate_R_TE_III(r_in, L, w_oc, t_oc, beta_oc, w_az, k_p);
            R_TE_III_n = obj.Geometry.calculate_R_TE_III(r_in, L, w_oc, t_oc, beta_oc, w_az, k_n);
            [R_t_ic, R_t_oc] = obj.Geometry.calculate_R_thermal_interconnects(r_in, L, w_ic, t_ic, beta_ic, w_oc, t_oc, beta_oc, k_Cu, k_Cu);

            R_I_p_combined = 1 / (1/R_TE_I_p + 2/R_t_ic);
            R_I_n_combined = 1 / (1/R_TE_I_n + 2/R_t_ic);
            R_III_p_combined = 1 / (1/R_TE_III_p + 2/R_t_oc);
            R_III_n_combined = 1 / (1/R_TE_III_n + 2/R_t_oc);

            R_leg_p = R_I_p_combined + R_TE_II_p + R_III_p_combined;
            R_leg_n = R_I_n_combined + R_TE_II_n + R_III_n_combined;
            K_legs = 1/R_leg_p + 1/R_leg_n;

            Re_leg_p = R_TE_II_p * (k_p * rho_p);
            Re_leg_n = R_TE_II_n * (k_n * rho_n);
            Re_leg = Re_leg_p + Re_leg_n;

            [R_e_ic, R_e_oc] = obj.Geometry.calculate_R_electrical_interconnects(r_in, L, w_ic, t_ic, beta_ic, w_oc, t_oc, beta_oc, rho_c);

            R_is = obj.Geometry.calculate_R_thermal_insulator(r_in + L, w_is_stage, obj.Geometry.Thickness, obj.Geometry.WedgeAngle, k_is);
            K_az_val = obj.Geometry.calculate_K_azimuthal(r_in, L, w_az, obj.Geometry.Thickness, k_az, obj.Geometry.WedgeAngle);
            K_TE_combined = K_legs + K_az_val;
            K_stage = 1 / (1/K_TE_combined + R_is);

            stage = struct(...
                'K_stage', K_stage, ...
                'S_stage', S_p - (-abs(S_n)), ...
                'Re_leg', Re_leg, ...
                'R_ic', R_e_ic, ...
                'R_oc', R_e_oc, ...
                'k_is', k_is ...
            );
        end
    end
end
