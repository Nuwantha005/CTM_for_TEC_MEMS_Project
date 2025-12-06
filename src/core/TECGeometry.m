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
            % TECGeometry constructor
            % Supports both old notation and new paper notation (Thermal_Network_For_Radial_TEC.tex)
            if isfield(config, 'geometry')
                obj.Params = config.geometry;
                obj.N_stages = obj.Params.N_stages;

                % Wedge angle: theta_deg (new) or wedge_angle_deg (old)
                if isfield(obj.Params, 'theta_deg')
                    obj.WedgeAngle = deg2rad(obj.Params.theta_deg);
                elseif isfield(obj.Params, 'wedge_angle_deg')
                    obj.WedgeAngle = deg2rad(obj.Params.wedge_angle_deg);
                else
                    obj.WedgeAngle = deg2rad(30);
                end

                % TEC thickness: t_TEC_um (new) or thickness_um (old)
                if isfield(obj.Params, 't_TEC_um')
                    obj.Thickness = obj.Params.t_TEC_um * 1e-6;
                elseif isfield(obj.Params, 'thickness_um')
                    obj.Thickness = obj.Params.thickness_um * 1e-6;
                else
                    obj.Thickness = 200e-6;
                end

                % Cylinder radius: r_cyl_um (new) or R_cyl_um (old)
                if isfield(obj.Params, 'r_cyl_um')
                    obj.R_cyl = obj.Params.r_cyl_um * 1e-6;
                elseif isfield(obj.Params, 'R_cyl_um')
                    obj.R_cyl = obj.Params.R_cyl_um * 1e-6;
                else
                    obj.R_cyl = 1000e-6;
                end

                % Chip width: W_chip_um (new) or w_chip_um (old)
                if isfield(obj.Params, 'W_chip_um')
                    w_chip = obj.Params.W_chip_um * 1e-6;
                elseif isfield(obj.Params, 'w_chip_um')
                    w_chip = obj.Params.w_chip_um * 1e-6;
                else
                    w_chip = 10000e-6;
                end
                obj.R_base = w_chip / sqrt(2);
                obj.calculate_L1();
            else
                error('Config must contain a "geometry" field.');
            end
        end

        function calculate_L1(obj)
            % Calculate L_1 using Eq. 238 from paper
            % Supports both old and new paper notation

            % Get insulation width ratio: W_is_ratio (new) or insulation_width_ratio (old)
            if isfield(obj.Params, 'W_is_ratio')
                L_avg_approx = (obj.R_base - obj.R_cyl) / obj.N_stages;
                w_is = obj.Params.W_is_ratio * L_avg_approx;
            elseif isfield(obj.Params, 'insulation_width_ratio')
                L_avg_approx = (obj.R_base - obj.R_cyl) / obj.N_stages;
                w_is = obj.Params.insulation_width_ratio * L_avg_approx;
            elseif isfield(obj.Params, 'insulation_width_um')
                w_is = obj.Params.insulation_width_um * 1e-6;
            else
                w_is = 40e-6;  % Default 40 um
            end

            % Length ratio: f_L (new) or radial_expansion_factor (old)
            if isfield(obj.Params, 'f_L')
                k_r = obj.Params.f_L;
            elseif isfield(obj.Params, 'radial_expansion_factor')
                k_r = obj.Params.radial_expansion_factor;
            else
                k_r = 1.15;
            end

            N = obj.N_stages;

            % Total length available for TEC material (Eq. 238)
            % Subtract (N+1) insulators
            L_total_active = (obj.R_base - obj.R_cyl) - (N + 1) * w_is;

            if L_total_active <= 0
                error('Geometry Error: No space for TE material. Check dimensions and insulator widths.');
            end

            if k_r == 1
                obj.L_1 = L_total_active / N;
            else
                obj.L_1 = L_total_active * (1 - k_r) / (1 - k_r^N);
            end

        end

        function [r_in, L, w_ic, t_ic, beta_ic, w_oc, t_oc, beta_oc, w_az, w_is] = get_stage_geometry(obj, i)
            % Get stage geometry with paper notation support
            % Supports both old and new (Thermal_Network_For_Radial_TEC.tex) notation

            % f_L (new) or radial_expansion_factor (old)
            if isfield(obj.Params, 'f_L')
                k_r = obj.Params.f_L;
            elseif isfield(obj.Params, 'radial_expansion_factor')
                k_r = obj.Params.radial_expansion_factor;
            else
                k_r = 1.15;
            end

            L = obj.L_1 * k_r^(i-1);

            % Insulation width: W_is_ratio (new) or insulation_width_ratio (old)
            if isfield(obj.Params, 'W_is_ratio')
                w_is = obj.Params.W_is_ratio * L;
            elseif isfield(obj.Params, 'insulation_width_ratio')
                w_is = obj.Params.insulation_width_ratio * L;
            elseif isfield(obj.Params, 'insulation_width_um')
                w_is = obj.Params.insulation_width_um * 1e-6;
            else
                w_is = 40e-6;  % Default 40 um
            end

            if k_r == 1
                sum_L_prev = (i-1) * obj.L_1;
            else
                sum_L_prev = obj.L_1 * (1 - k_r^(i-1)) / (1 - k_r);
            end

            % Radius start calculation (Eq. from paper)
            r_in = obj.R_cyl + sum_L_prev + i * w_is;

            % Interconnector width: f_ic_W (new) or interconnect_ratio (old)
            if isfield(obj.Params, 'f_ic_W')
                w_ic = L * obj.Params.f_ic_W;
            elseif isfield(obj.Params, 'interconnect_ratio')
                w_ic = L * obj.Params.interconnect_ratio;
            else
                w_ic = L * 0.15;
            end

            % Outerconnector width: f_oc_W (new) or outerconnect_ratio (old)
            if isfield(obj.Params, 'f_oc_W')
                w_oc = L * obj.Params.f_oc_W;
            elseif isfield(obj.Params, 'outerconnect_ratio')
                w_oc = L * obj.Params.outerconnect_ratio;
            else
                w_oc = L * 0.15;
            end

            % Interconnector thickness: f_ic_t (new) or interconnect_thickness_ratio (old)
            if isfield(obj.Params, 'f_ic_t')
                t_ic = obj.Thickness * obj.Params.f_ic_t;
            elseif isfield(obj.Params, 'interconnect_thickness_ratio')
                t_ic = obj.Thickness * obj.Params.interconnect_thickness_ratio;
            else
                t_ic = obj.Thickness;
            end

            % Outerconnector thickness: f_oc_t (new) or outerconnect_thickness_ratio (old)
            if isfield(obj.Params, 'f_oc_t')
                t_oc = obj.Thickness * obj.Params.f_oc_t;
            elseif isfield(obj.Params, 'outerconnect_thickness_ratio')
                t_oc = obj.Thickness * obj.Params.outerconnect_thickness_ratio;
            else
                t_oc = obj.Thickness;
            end

            % Interconnector angle: f_ic_beta (new) or interconnect_angle_ratio (old)
            if isfield(obj.Params, 'f_ic_beta')
                beta_ic = obj.WedgeAngle * obj.Params.f_ic_beta;
            elseif isfield(obj.Params, 'interconnect_angle_ratio')
                beta_ic = obj.WedgeAngle * obj.Params.interconnect_angle_ratio;
            elseif isfield(obj.Params, 'interconnect_angle_deg')
                beta_ic = deg2rad(obj.Params.interconnect_angle_deg);
            else
                beta_ic = deg2rad(5);
            end

            % Outerconnector angle: f_oc_beta (new) or outerconnect_angle_ratio (old)
            if isfield(obj.Params, 'f_oc_beta')
                beta_oc = obj.WedgeAngle * obj.Params.f_oc_beta;
            elseif isfield(obj.Params, 'outerconnect_angle_ratio')
                beta_oc = obj.WedgeAngle * obj.Params.outerconnect_angle_ratio;
            elseif isfield(obj.Params, 'outerconnect_angle_deg')
                beta_oc = deg2rad(obj.Params.outerconnect_angle_deg);
            else
                beta_oc = deg2rad(5);
            end

            % Azimuthal Gap: W_az_um (new) or fill_factor (old)
            if isfield(obj.Params, 'W_az_um')
                w_az = obj.Params.W_az_um * 1e-6;
            elseif isfield(obj.Params, 'azimuthal_gap_um')
                w_az = obj.Params.azimuthal_gap_um * 1e-6;
            elseif isfield(obj.Params, 'fill_factor')
                r_mid = r_in + L/2;
                arc_length = r_mid * obj.WedgeAngle;
                w_az = (1 - obj.Params.fill_factor) * arc_length;
            else
                w_az = 20e-6;
            end
        end

        function G = calculate_G(obj, r1, L, w_ic, t_ic, beta_ic, w_oc, t_oc, beta_oc, w_az, w_is, ~)
            t = obj.Thickness;
            theta = obj.WedgeAngle;

            C1 = t*theta/2 - t_ic*beta_ic/2;
            D = t*w_az;

            r_start = r1;
            r_limit_1 = r1 + w_ic;

            term1 = (1/C1) * log( abs( ( C1*r_limit_1 - D ) / ( C1*r_start - D ) ) );

            C2 = t*theta/2;
            r_limit_2 = r1 + L - w_oc;

            term2 = (1/C2) * log( abs( ( C2*r_limit_2 - D ) / ( C2*r_limit_1 - D ) ) );

            C3 = t*theta/2 - t_oc*beta_oc/2;
            r_end = r1 + L;

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

        function R_ve = calculate_vertical_resistance(obj, r_in, r_out, k_ins)
            % Calculates Vertical thermal resistance through insulator layer
            % Eq 785: R_ve = 2*t_ins / (k_ins * theta * (r_out^2 - r_in^2))

            if isfield(obj.Params, 't_ins_um')
                t_ins = obj.Params.t_ins_um * 1e-6;
            else
                t_ins = 10e-6; % Default
            end

            theta = obj.WedgeAngle;

            % Area of the ANNULAR SEGMENT under the TEC element
            % A = (theta/2) * (r_out^2 - r_in^2)
            % Note: Formula in Latex uses A_TEC, which is correct.

            R_ve = (2 * t_ins) / (k_ins * theta * (r_out^2 - r_in^2));
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
