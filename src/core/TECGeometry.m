classdef TECGeometry < handle
    % TECGEOMETRY Handles geometric calculations for the Radial TEC.

    properties
        Params
        N_stages
        WedgeAngle
        Thickness
        L_1
        R_cyl
        R_base
    end

    methods
        function obj = TECGeometry(config)
            if isfield(config, 'geometry')
                obj.Params = config.geometry;
                obj.N_stages = obj.Params.N_stages;
                obj.WedgeAngle = deg2rad(obj.Params.wedge_angle_deg);
                obj.Thickness = obj.Params.thickness_um * 1e-6;
                obj.R_cyl = obj.Params.R_cyl_um * 1e-6;
                w_chip = obj.Params.w_chip_um * 1e-6;
                obj.R_base = w_chip / sqrt(2);
                obj.calculate_L1();
            else
                error('Config must contain a "geometry" field.');
            end
        end

        function calculate_L1(obj)
            w_is = obj.Params.insulation_width_um * 1e-6;
            k_r = obj.Params.radial_expansion_factor;
            N = obj.N_stages;

            L_total_active = (obj.R_base - obj.R_cyl) - (N + 1) * w_is;

            if L_total_active <= 0
                error('Geometry Error: No space for TE material. Check dimensions.');
            end

            if k_r == 1
                obj.L_1 = L_total_active / N;
            else
                obj.L_1 = L_total_active * (1 - k_r) / (1 - k_r^N);
            end

            fprintf('L_1 calculated = %e m\n', obj.L_1);
        end

        function [r_in, L, w_ic, t_ic, beta_ic, w_oc, t_oc, beta_oc, w_az, w_is] = get_stage_geometry(obj, i)
            k_r = obj.Params.radial_expansion_factor;
            L = obj.L_1 * k_r^(i-1);

            w_is = obj.Params.insulation_width_um * 1e-6;

            if k_r == 1
                sum_L_prev = (i-1) * obj.L_1;
            else
                sum_L_prev = obj.L_1 * (1 - k_r^(i-1)) / (1 - k_r);
            end

            r_in = obj.R_cyl + sum_L_prev; % + i * w_is becasue in the original formulation, w_is is inside the L

            w_ic = L * obj.Params.interconnect_ratio;
            w_oc = L * obj.Params.outerconnect_ratio;

            t_ic = obj.Thickness;
            t_oc = obj.Thickness;

            if isfield(obj.Params, 'interconnect_angle_deg')
                beta_ic = deg2rad(obj.Params.interconnect_angle_deg);
            else
                beta_ic = deg2rad(5);
            end

            if isfield(obj.Params, 'outerconnect_angle_deg')
                beta_oc = deg2rad(obj.Params.outerconnect_angle_deg);
            else
                beta_oc = deg2rad(5);
            end

            w_az = obj.Params.azimuthal_gap_um * 1e-6;
        end

        function G = calculate_G(obj, r1, L, w_ic, t_ic, beta_ic, w_oc, t_oc, beta_oc, w_az, w_is, ~)
            t = obj.Thickness;
            theta = obj.WedgeAngle;

            C1 = t*theta/2 - t_ic*beta_ic/2;
            D = t*w_az;

            r_start = r1 + w_is/2;
            r_limit_1 = r1 + w_ic;

            term1 = (1/C1) * log( abs( ( C1*r_limit_1 - D ) / ( C1*r_start - D ) ) );

            C2 = t*theta/2;
            r_limit_2 = r1 + L - w_oc;

            term2 = (1/C2) * log( abs( ( C2*r_limit_2 - D ) / ( C2*r_limit_1 - D ) ) );

            C3 = t*theta/2 - t_oc*beta_oc/2;
            r_end = r1 + L - w_is/2;

            term3 = (1/C3) * log( abs( ( C3*r_end - D ) / ( C3*r_limit_2 - D ) ) );

            G = term1 + term2 + term3;
        end

        function [R_ic, R_oc] = calculate_R_electrical_interconnects(obj, r1, L, w_ic, t_ic, beta_ic, w_oc, t_oc, beta_oc, rho_c)
            term_ic = log((r1 + w_ic) / r1);
            if term_ic == 0
                R_ic = 0;
            else
                R_ic = (rho_c * beta_ic) / (t_ic * term_ic);
            end

            term_oc = log((r1 + L) / (r1 + L - w_oc));
            if term_oc == 0
                R_oc = 0;
            else
                R_oc = (rho_c * beta_oc) / (2 * t_oc * term_oc);
            end
        end

        function R_is = calculate_R_thermal_insulator(obj, r_end_leg, w_is, t, theta, k_is)
            r_outer = r_end_leg + w_is;
            r_inner = r_end_leg;

            R_is = (1 / (k_is * t * theta)) * log(r_outer / r_inner);
        end

        function K_az = calculate_K_azimuthal(obj, r1, L, w_az, t, k_az, theta)
            K_az = k_az * (w_az * t) / L;
        end

        function [N_TSV, R_TSV_tot] = calculate_TSV_vertical_resistance(obj, r, w_ic, beta_ic, rho_TSV, stage_idx)
            if ~isfield(obj.Params, 'tsv')
                N_TSV = 1;
                R_TSV_tot = 1e-4;
                return;
            end

            N_tsv_limit = 1;
            if isfield(obj.Params, 'N_tsv_limit')
                N_tsv_limit = obj.Params.N_tsv_limit;
            end

            if stage_idx > N_tsv_limit
                N_TSV = 0;
                R_TSV_tot = 1e9;
                return;
            end

            tsv = obj.Params.tsv;
            R_TSV_rad = tsv.R_TSV_um * 1e-6;
            P_TSV = tsv.P_TSV_um * 1e-6;
            g_rad = tsv.g_rad_um * 1e-6;
            t_SOI = tsv.t_SOI_um * 1e-6;

            N_row = floor(w_ic / (2 * R_TSV_rad + g_rad));
            N_per_row = floor((r * beta_ic) / P_TSV);
            N_TSV = N_row * N_per_row;

            if N_TSV < 1
                N_TSV = 0;
                R_TSV_tot = 1e9;
            else
                R_single = rho_TSV * t_SOI / (pi * R_TSV_rad^2);
                R_TSV_tot = R_single / N_TSV;
            end
        end

        function export_comsol_params(obj, filename)
            fid = fopen(filename, 'w');
            if fid == -1
                error('Could not open file %s for writing.', filename);
            end

            fprintf(fid, '%% Radial TEC Parameters\n');
            fprintf(fid, 'N_stages %d\n', obj.N_stages);
            fprintf(fid, 'theta %f[deg]\n', rad2deg(obj.WedgeAngle));
            fprintf(fid, 't_tec %e[m]\n', obj.Thickness);
            fprintf(fid, 'R_cyl %e[m]\n', obj.R_cyl);
            fprintf(fid, 'L_1 %e[m]\n', obj.L_1);

            fclose(fid);
            fprintf('COMSOL parameters exported to %s\n', filename);
        end
    end
end
