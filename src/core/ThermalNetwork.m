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

            % Calculate thermal resistance using new approach (same as in assemble_system)
            k_is = obj.Materials.get_k('AlN', T_avg);
            k_az = obj.Materials.get_k('SiO2', T_avg);
            k_Cu = obj.Materials.get_k('Cu', T_avg);
            
            % Calculate TE region resistances
            R_TE_I_p = obj.Geometry.calculate_R_TE_I(r_in, w_ic, t_ic, beta_ic, w_az, k_p);
            R_TE_I_n = obj.Geometry.calculate_R_TE_I(r_in, w_ic, t_ic, beta_ic, w_az, k_n);
            R_TE_II_p = obj.Geometry.calculate_R_TE_II(r_in, L, w_ic, w_oc, w_az, k_p);
            R_TE_II_n = obj.Geometry.calculate_R_TE_II(r_in, L, w_ic, w_oc, w_az, k_n);
            R_TE_III_p = obj.Geometry.calculate_R_TE_III(r_in, L, w_oc, t_oc, beta_oc, w_az, k_p);
            R_TE_III_n = obj.Geometry.calculate_R_TE_III(r_in, L, w_oc, t_oc, beta_oc, w_az, k_n);
            
            % IC/OC thermal resistances
            [R_t_ic, R_t_oc] = obj.Geometry.calculate_R_thermal_interconnects(r_in, L, w_ic, t_ic, beta_ic, w_oc, t_oc, beta_oc, k_Cu, k_Cu);
            
            % Combine regions
            R_I_p_combined = 1 / (1/R_TE_I_p + 2/R_t_ic);
            R_I_n_combined = 1 / (1/R_TE_I_n + 2/R_t_ic);
            R_III_p_combined = 1 / (1/R_TE_III_p + 2/R_t_oc);
            R_III_n_combined = 1 / (1/R_TE_III_n + 2/R_t_oc);
            
            R_leg_p = R_I_p_combined + R_TE_II_p + R_III_p_combined;
            R_leg_n = R_I_n_combined + R_TE_II_n + R_III_n_combined;
            K_legs = 1/R_leg_p + 1/R_leg_n;

            r_out = r_in + L;
            r_end_leg = r_out;
            R_is = obj.Geometry.calculate_R_thermal_insulator(r_end_leg, w_is, obj.Geometry.Thickness, obj.Geometry.WedgeAngle, k_is);

            R_eff_series = R_is + 1/K_legs;
            K_eff_series = 1 / R_eff_series;

            K_az_val = obj.Geometry.calculate_K_azimuthal(r_in, L, w_az, obj.Geometry.Thickness, k_az, obj.Geometry.WedgeAngle);

            K_N = K_eff_series + 5 * K_az_val;

            % Electrical resistance (Region II only)
            Re_leg_p = R_TE_II_p * (k_p * rho_p);
            Re_leg_n = R_TE_II_n * (k_n * rho_n);
            Re_N = Re_leg_p + Re_leg_n;
            
            rho_c = obj.Materials.get_rho('Cu', T_avg);
            [R_ic, R_oc] = obj.Geometry.calculate_R_electrical_interconnects(r_in, L, w_ic, t_ic, beta_ic, w_oc, t_oc, beta_oc, rho_c);

            S_N = S_p - (-abs(S_n));
            I = obj.Params.operating_conditions.I_current_A;

            % Q_c equation: Heat INTO cold side of stage N
            % Q_c = S*I*T_c + K*(T_h - T_c) - [Joule at cold]

            % We want Q_out (Heat rejected to water at hot side)
            % Q_h = S*I*T_h - k*(T_h - T_c) + 0.5*I^2*R + I^2*R_oc
            % Note: Check directions. K*(T_c - T_h) flows TO water.

            % Correct Q_out calculation:
            % S term: Carries heat to hot side
            % K term: Conducts heat from hot to cold (backflow). We want net flow OUT.
            % Net flow OUT = Peltier_at_hot - Back_Conduction + Joule_Heating

            Q_Peltier_hot = S_N * I * T_water;
            Q_Back_Cond = K_N * (T_water - T_cold); % Heat flowing flowing from Water to Cold
            % If T_cold > T_water, this term is negative, meaning heat flows TO water.
            % But usually formulated as: Q_cond = K(dT).
            % Let's use: Q_out = S*I*T_h + K*(T_c - T_h) + Joule

            Q_Conduction_to_Water = K_N * (T_cold - T_water);
            Q_Joule_Hot = 0.5 * I^2 * Re_N + I^2 * R_oc;

            Q_out = Q_Peltier_hot + Q_Conduction_to_Water + Q_Joule_Hot;


            % n_wedges = 2*pi / obj.Geometry.WedgeAngle;
            % Q_out is already total for the wedge

            % Calculate total input heat (Simulated)
            % Note: The Latex model assigns heat generation based on TEC element areas.
            % This ignores flux falling on radial insulators.
            % To maintain energy balance consistency in results, we report the actual heat
            % injected into the nodes, not the theoretical full-area flux.

            % We need to access Q_gen terms calculated in assemble_system.
            % Since they aren't stored, we re-calculate or sum them here.
            % Actually, let's just recalculate the effective area sum.

            q_flux = obj.Params.boundary_conditions.q_flux_W_m2;
            theta = obj.Geometry.WedgeAngle;
            A_active_sum = 0;
            % Cylinder
            A_active_sum = A_active_sum + 0.5 * theta * obj.Geometry.R_cyl^2;

            % Stages
            for i = 1:N
                [r_in, L] = obj.Geometry.get_stage_geometry(i);
                r_out = r_in + L;
                A_active_sum = A_active_sum + 0.5 * theta * (r_out^2 - r_in^2);
            end

            Q_in = q_flux * A_active_sum;
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
            [~, L1, ~, ~, ~, ~, ~, ~, ~, w_is] = obj.Geometry.get_stage_geometry(1);

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
                k_Cu = obj.Materials.get_k('Cu', T_avg);  % IC/OC thermal conductivity

                % Calculate thermal resistances for the three TE regions
                % Region I: IC region (IC and TE in parallel)
                R_TE_I_p = obj.Geometry.calculate_R_TE_I(r_in, w_ic, t_ic, beta_ic, w_az, k_p);
                R_TE_I_n = obj.Geometry.calculate_R_TE_I(r_in, w_ic, t_ic, beta_ic, w_az, k_n);
                
                % Region II: Middle region (no IC/OC)
                R_TE_II_p = obj.Geometry.calculate_R_TE_II(r_in, L, w_ic, w_oc, w_az, k_p);
                R_TE_II_n = obj.Geometry.calculate_R_TE_II(r_in, L, w_ic, w_oc, w_az, k_n);
                
                % Region III: OC region (OC and TE in parallel)
                R_TE_III_p = obj.Geometry.calculate_R_TE_III(r_in, L, w_oc, t_oc, beta_oc, w_az, k_p);
                R_TE_III_n = obj.Geometry.calculate_R_TE_III(r_in, L, w_oc, t_oc, beta_oc, w_az, k_n);
                
                % Calculate IC/OC thermal resistances (radial conduction through Cu)
                [R_t_ic, R_t_oc] = obj.Geometry.calculate_R_thermal_interconnects(r_in, L, w_ic, t_ic, beta_ic, w_oc, t_oc, beta_oc, k_Cu, k_Cu);
                
                % Combine Region I: TE and IC in parallel (IC serves both P and N legs)
                % 1/R_I_combined = 1/R_TE_I + 1/(R_t_ic/2) where R_t_ic/2 for each leg
                R_I_p_combined = 1 / (1/R_TE_I_p + 2/R_t_ic);
                R_I_n_combined = 1 / (1/R_TE_I_n + 2/R_t_ic);
                
                % Combine Region III: TE and OC in parallel
                R_III_p_combined = 1 / (1/R_TE_III_p + 2/R_t_oc);
                R_III_n_combined = 1 / (1/R_TE_III_n + 2/R_t_oc);
                
                % Total resistance for each leg: series combination of three regions
                R_leg_p = R_I_p_combined + R_TE_II_p + R_III_p_combined;
                R_leg_n = R_I_n_combined + R_TE_II_n + R_III_n_combined;
                
                % Two legs in parallel
                K_legs = 1/R_leg_p + 1/R_leg_n;
                
                % Electrical resistance (only through Region II, no lateral conduction)
                % Re = rho * G_II where G_II is the geometric factor for Region II
                % For region II: Use same formula as R_TE_II but with 1/rho instead of k
                Re_leg_p = R_TE_II_p * (k_p * rho_p);  % Re = (rho/k) * R_thermal
                Re_leg_n = R_TE_II_n * (k_n * rho_n);
                Re_stages_leg(i) = Re_leg_p + Re_leg_n;  % P and N in series electrically

                % Electrical resistances of IC/OC
                [R_e_ic, R_e_oc] = obj.Geometry.calculate_R_electrical_interconnects(r_in, L, w_ic, t_ic, beta_ic, w_oc, t_oc, beta_oc, rho_c);
                R_ic_stages(i) = R_e_ic;
                R_oc_stages(i) = R_e_oc;

                % Radial insulator thermal resistance
                r_end_leg = r_out;
                R_is = obj.Geometry.calculate_R_thermal_insulator(r_end_leg, w_is_stage, t_tec, theta, k_is);

                % Combine with radial insulator (series)
                R_eff_series = R_is + 1/K_legs;
                K_eff_series = 1 / R_eff_series;

                % Azimuthal conductance (5 parallel paths)
                K_az_val = obj.Geometry.calculate_K_azimuthal(r_in, L, w_az, t_tec, k_az, theta);

                % Total stage conductance: TE legs + 5 azimuthal paths in parallel
                K_stages(i) = K_eff_series + 5 * K_az_val;
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

                R_vert(i) = obj.Geometry.calculate_vertical_resistance(r_in, r_out, k_is);

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
                    M(idx_c, idx_c_prev) = (S_stages(i-1)*I + K_stages(i-1));
                end

                if i < N
                    idx_c_next = idx_c + 1;
                    M(idx_c, idx_c_next) = K_stages(i); %+K
                else
                    B(idx_c) = B(idx_c) - K_stages(i) * T_water; %-K
                end

                B(idx_c) = B(idx_c) - I^2 * (Re_stages_leg(i)/2 + R_ic_stages(i));
                if i > 1
                    B(idx_c) = B(idx_c) - I^2 * (Re_stages_leg(i-1)/2 + R_oc_stages(i-1));
                end

                sum_diag = G_vert + S_stages(i)*I + K_stages(i); % -K
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
