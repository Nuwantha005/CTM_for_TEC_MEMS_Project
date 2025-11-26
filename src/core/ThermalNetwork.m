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
            [M, B] = obj.assemble_system(T_current);
            T_new = M \ B;

            N = obj.Geometry.N_stages;
            idx_c_N = 2*N + 1;
            T_water = obj.Params.boundary_conditions.T_water_K;

            T_cold = T_new(idx_c_N);
            T_avg = (T_cold + T_water) / 2;

            k_p = obj.Materials.get_k('Bi2Te3', T_avg);
            k_n = obj.Materials.get_k('Bi2Te3', T_avg);
            rho_p = obj.Materials.get_rho('Bi2Te3', T_avg);
            rho_n = obj.Materials.get_rho('Bi2Te3', T_avg);
            S_p = obj.Materials.get_S('Bi2Te3', T_avg);
            S_n = obj.Materials.get_S('Bi2Te3', T_avg);

            [r_in, L, w_ic, t_ic, beta_ic, w_oc, t_oc, beta_oc, w_az, w_is] = obj.Geometry.get_stage_geometry(N);
            G = obj.Geometry.calculate_G(r_in, L, w_ic, t_ic, beta_ic, w_oc, t_oc, beta_oc, w_az, w_is);
            
            K_legs = (k_p + k_n) / G;

            k_is = obj.Materials.get_k('AlN', T_avg);
            k_az = obj.Materials.get_k('SiO2', T_avg);

            r_end_leg = r_in + L - w_is;
            R_is = obj.Geometry.calculate_R_thermal_insulator(r_end_leg, w_is, obj.Geometry.Thickness, obj.Geometry.WedgeAngle, k_is);

            R_eff_series = R_is + 1/K_legs;
            K_eff_series = 1 / R_eff_series;

            K_az_val = obj.Geometry.calculate_K_azimuthal(r_in, L, w_az, obj.Geometry.Thickness, k_az, obj.Geometry.WedgeAngle);

            K_N = K_eff_series + K_az_val;

            Re_N = (rho_p + rho_n) * G;
            rho_c = obj.Materials.get_rho('Cu', T_avg);
            [R_ic, R_oc] = obj.Geometry.calculate_R_electrical_interconnects(r_in, L, w_ic, t_ic, beta_ic, w_oc, t_oc, beta_oc, rho_c);

            S_N = S_p - (-abs(S_n));
            I = obj.Params.operating_conditions.I_current_A;

            % Q_c equation: Heat INTO cold side of stage N
            % Q_c = S*I*T_c + K*(T_h - T_c) - [Joule at cold]
            % Note: Cold side Joule heating uses R_ic (interconnect), not R_oc (outerconnect)
            Q_c_N = S_N * I * T_cold + K_N * (T_water - T_cold) - (0.5 * I^2 * Re_N + I^2 * R_ic);

            %n_wedges = 2*pi / obj.Geometry.WedgeAngle;
            Q_out = Q_c_N ;%* n_wedges;

            % Calculate total input heat
            q_flux = obj.Params.boundary_conditions.q_flux_W_m2;
            theta = obj.Geometry.WedgeAngle;
            R_base = obj.Geometry.R_base;
            A_wedge = 0.5 * theta * R_base^2;
            Q_in = q_flux * A_wedge; %* n_wedges
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
            w_is = obj.Params.geometry.insulation_width_um * 1e-6;
            R_cyl = obj.Geometry.R_cyl;

            K_stages = zeros(N, 1);
            S_stages = zeros(N, 1);
            Re_stages_leg = zeros(N, 1);
            R_ic_stages = zeros(N, 1);
            R_oc_stages = zeros(N, 1);
            R_lat_Si = zeros(N, 1);
            R_vert = zeros(N, 1);
            Q_gen_nodes = zeros(N, 1);

            for i = 1:N
                [r_in, L, w_ic, t_ic, beta_ic, w_oc, t_oc, beta_oc, w_az, w_is_stage] = obj.Geometry.get_stage_geometry(i);
                r_out = r_in + L;

                T_cold = T_current(idx_c_start + i - 1);
                if i < N
                    T_hot = T_current(idx_c_start + i);
                else
                    T_hot = T_water;
                end
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
                k_TSV = obj.Materials.get_k('Cu', T_avg);  % Thermal conductivity for TSV

                G = obj.Geometry.calculate_G(r_in, L, w_ic, t_ic, beta_ic, w_oc, t_oc, beta_oc, w_az, w_is, w_is_stage);

                K_legs = (k_p + k_n) / G; % Thermal
                Re_stages_leg(i) = (rho_p + rho_n) * G; %% Electrical

                [R_ic, R_oc] = obj.Geometry.calculate_R_electrical_interconnects(r_in, L, w_ic, t_ic, beta_ic, w_oc, t_oc, beta_oc, rho_c);
                R_ic_stages(i) = R_ic;
                R_oc_stages(i) = R_oc;

                r_end_leg = r_in + L - w_is_stage;
                R_is = obj.Geometry.calculate_R_thermal_insulator(r_end_leg, w_is_stage, t_tec, theta, k_is);

                R_eff_series = R_is + 1/K_legs; % Thermal
                K_eff_series = 1 / R_eff_series; % Thermal

                K_az_val = obj.Geometry.calculate_K_azimuthal(r_in, L, w_az, t_tec, k_az, theta);

                K_stages(i) = K_eff_series + K_az_val; % Thermal
                S_stages(i) = S_p - (-abs(S_n));

                k_Si = obj.Materials.get_k('Si', T_current(idx_Si_start + i - 1));

                if i < N
                    [r_in_next, L_next] = obj.Geometry.get_stage_geometry(i+1);
                    r_mid_i = r_in + L/2;
                    r_mid_next =  r_in_next + L_next/2;
                    R_lat_Si(i) = log(r_mid_next/r_mid_i) / (k_Si * theta * t_chip);
                else
                    R_lat_Si(i) = inf;
                end

                [~, R_TSV_tot] = obj.Geometry.calculate_TSV_vertical_resistance(r_in, w_ic, beta_ic, k_TSV, i);
                R_vert(i) = R_TSV_tot;

                A_top = 0.5 * theta * (r_out^2 - r_in^2);
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

            M(idx_0, idx_0) = -(1/R_Si_01 + 1/R_TEC_01);  % Negative: sum of conductances leaving
            M(idx_0, idx_Si_start) = +1/R_Si_01;  % Positive: conductance to neighbor
            M(idx_0, idx_c_start) = +1/R_TEC_01;  % Positive: conductance to neighbor
            B(idx_0) = -Q_gen_0;

            for i = 1:N
                idx_Si = idx_Si_start + i - 1;
                idx_c = idx_c_start + i - 1;

                if i == 1
                    M(idx_Si, idx_0) = +1/R_Si_01;  % Positive: heat flows from T_0 to T_Si,1
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
                    M(idx_c, idx_0) = +1/R_TEC_01;  % Positive: heat flows from T_0 to T_c,1
                end

                if i > 1
                    idx_c_prev = idx_c - 1;
                    M(idx_c, idx_c_prev) = (S_stages(i-1)*I + K_stages(i-1));
                end

                if i < N
                    idx_c_next = idx_c + 1;
                    M(idx_c, idx_c_next) = -K_stages(i);  % Back conduction: -K_i from derivation
                else
                    B(idx_c) = B(idx_c) - K_stages(i) * T_water;  % Negative: K*T_water moved to RHS
                end
                
                B(idx_c) = B(idx_c) - I^2 * (Re_stages_leg(i)/2 + R_ic_stages(i));
                if i > 1
                    B(idx_c) = B(idx_c) - I^2 * (Re_stages_leg(i-1)/2 + R_oc_stages(i-1));
                end

                % Diagonal: -(G_vert + K_{i-1} + S_i*I + K_i) for boundary
                %           -(G_vert + K_{i-1} + S_i*I - K_i) for interior (K_i cancels with off-diag)
                % For i < N: -K_i in diagonal is cancelled by +K_i coefficient of T_{i+1}
                % For i = N: K_i goes to boundary term, so +K_i in diagonal (no cancellation)
                sum_diag = G_vert + S_stages(i)*I;
                if i < N
                    sum_diag = sum_diag - K_stages(i);  % Interior: subtract K_i (will add T_{i+1} term)
                else
                    sum_diag = sum_diag + K_stages(i);  % Boundary: add K_i (T_water is known)
                end
                if i > 1
                    sum_diag = sum_diag + K_stages(i-1);
                end
                if i == 1
                    sum_diag = sum_diag + 1/R_TEC_01;
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
    end
end
