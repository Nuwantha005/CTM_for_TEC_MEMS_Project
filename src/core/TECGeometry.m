classdef TECGeometry
    % TECGEOMETRY Handles geometric calculations for the Radial TEC.
    %   Calculates the Geometric Factor (G) and resistance components.
    
    properties
        Params
        N_stages
        WedgeAngle % theta
        Thickness  % t
    end
    
    methods
        function obj = TECGeometry(config)
            % Constructor
            if isfield(config, 'geometry')
                obj.Params = config.geometry;
                obj.N_stages = obj.Params.N_stages;
                obj.WedgeAngle = deg2rad(obj.Params.wedge_angle_deg);
                obj.Thickness = obj.Params.thickness_um * 1e-6; % Convert to meters
            else
                error('Config must contain a "geometry" field.');
            end
        end
        
        function [r_in, L, w_ic, t_ic, beta_ic, w_oc, t_oc, beta_oc, w_az, w_is] = get_stage_geometry(obj, i)
            % GET_STAGE_GEOMETRY Returns geometric params for stage i.
            %   Assumes linear radial progression for now.
            
            % Base params
            L_base = obj.Params.L_element_um * 1e-6;
            r_start = obj.Params.r_start_um * 1e-6;
            k_r = obj.Params.radial_expansion_factor;
            
            % Calculate L_i based on radial expansion
            % L_i = k_r^(i-1) * L_1
            L = L_base * k_r^(i-1);
            
            % Calculate r_in for stage i
            % r_in = r_start + sum(L_j) for j=1 to i-1
            % For constant k_r, this is a geometric series sum
            if k_r == 1
                r_in = r_start + (i-1) * L_base;
            else
                r_in = r_start + L_base * (1 - k_r^(i-1)) / (1 - k_r);
            end
            
            % Interconnects (using ratios from config)
            w_ic = L * obj.Params.interconnect_ratio;
            w_oc = L * obj.Params.outerconnect_ratio;
            
            % Thicknesses (assume same as element for now, or from params)
            t_ic = obj.Thickness; 
            t_oc = obj.Thickness;
            
            % Angles
            if isfield(obj.Params, 'interconnect_angle_deg')
                beta_ic = deg2rad(obj.Params.interconnect_angle_deg);
            else
                beta_ic = deg2rad(5); % Default from notes
            end
            
            if isfield(obj.Params, 'outerconnect_angle_deg')
                beta_oc = deg2rad(obj.Params.outerconnect_angle_deg);
            else
                beta_oc = deg2rad(5); % Default from notes
            end
            
            % Gaps
            w_az = obj.Params.azimuthal_gap_um * 1e-6;
            w_is = obj.Params.insulation_width_um * 1e-6;
        end
        
        function G = calculate_G(obj, r1, L, w_ic, t_ic, beta_ic, w_oc, t_oc, beta_oc, w_az, w_is)
            % CALCULATE_G Calculates the Geometric Factor G for a single leg.
            %   All inputs should be in SI units (meters, radians).
            %   Formula derived in "Variable Cross section area TEC.md".
            
            t = obj.Thickness;
            theta = obj.WedgeAngle;
            
            % Derived boundaries
            r_start = r1 + w_is/2;
            r_end = r1 + L - w_is/2;
            
            % Term 1: Inner Connect Region
            % Range: [r_start, r1 + w_ic]
            % A1(r) = r(t*theta/2 - t_ic*beta_ic/2) - t*w_az
            C1 = t*theta/2 - t_ic*beta_ic/2;
            D = t*w_az;
            
            r_limit_1 = r1 + w_ic;
            
            % Check for validity (log argument must be positive)
            % (C1*r - D) > 0
            
            term1 = (1/C1) * log( abs( ( C1*r_limit_1 - D ) / ( C1*r_start - D ) ) );
            
            % Term 2: Full TE Region
            % Range: [r1 + w_ic, r1 + L - w_oc]
            % A2(r) = r(t*theta/2) - t*w_az
            C2 = t*theta/2;
            % D is same
            
            r_limit_2 = r1 + L - w_oc;
            
            term2 = (1/C2) * log( abs( ( C2*r_limit_2 - D ) / ( C2*r_limit_1 - D ) ) );
            
            % Term 3: Outer Connect Region
            % Range: [r1 + L - w_oc, r_end]
            % A3(r) = r(t*theta/2 - t_oc*beta_oc/2) - t*w_az
            C3 = t*theta/2 - t_oc*beta_oc/2;
            % D is same
            
            term3 = (1/C3) * log( abs( ( C3*r_end - D ) / ( C3*r_limit_2 - D ) ) );
            
            G = term1 + term2 + term3;
        end
        
        function [R_ic, R_oc] = calculate_R_electrical_interconnects(obj, r1, L, w_ic, t_ic, beta_ic, w_oc, t_oc, beta_oc, rho_c)
            % CALCULATE_R_ELECTRICAL_INTERCONNECTS Calculates electrical resistance of interconnects.
            %   Eq. from "Variable Cross section area TEC.md" Section 2.
            
            % Inner Interconnect
            % R_ic = (rho_c * beta_ic) / (t_ic * ln((r1+w_ic)/r1))
            term_ic = log((r1 + w_ic) / r1);
            if term_ic == 0
                R_ic = 0;
            else
                R_ic = (rho_c * beta_ic) / (t_ic * term_ic);
            end
            
            % Outer Interconnect
            % R_oc = (rho_c * beta_oc) / (2 * t_oc * ln((r1+L)/(r1+L-w_oc)))
            % Note: The factor 2 is from the derivation (parallel strips).
            term_oc = log((r1 + L) / (r1 + L - w_oc));
            if term_oc == 0
                R_oc = 0;
            else
                R_oc = (rho_c * beta_oc) / (2 * t_oc * term_oc);
            end
        end
        
        function R_is = calculate_R_thermal_insulator(obj, r_end_leg, w_is, t, theta, k_is)
            % CALCULATE_R_THERMAL_INSULATOR Calculates radial thermal resistance of insulator.
            %   Eq. from "TEC Thermal Network.md" Section 2.1.
            %   R_is = (1 / (k_is * t * theta)) * ln((r_end_leg + w_is) / r_end_leg)
            
            r_outer = r_end_leg + w_is;
            r_inner = r_end_leg;
            
            R_is = (1 / (k_is * t * theta)) * log(r_outer / r_inner);
        end
        
        function K_az = calculate_K_azimuthal(obj, r1, L, w_az, t, k_az, theta)
            % CALCULATE_K_AZIMUTHAL Calculates azimuthal leakage conductance.
            %   Eq. from "TEC Thermal Network.md" Section 2.2.
            %   K_az = k_az * (w_az * t) / L
            %   Wait, the note says: "Since this is characterized by the arc length... A(r)=w_az*t... K_az = k_az * w_az * t / L"
            %   This assumes w_az is constant width.
            %   Yes, "w_az is the Azimuthal insulation width".
            
            K_az = k_az * (w_az * t) / L;
        end
        
        function [N_TSV, R_TSV_tot] = calculate_TSV_vertical_resistance(obj, r, w_ic, beta_ic, rho_TSV)
            % CALCULATE_TSV_VERTICAL_RESISTANCE Calculates vertical resistance.
            %   Eq. 306-318 from notes.
            
            if ~isfield(obj.Params, 'tsv')
                % Fallback if config not updated
                N_TSV = 1;
                R_TSV_tot = 1e-4; 
                return;
            end
            
            tsv = obj.Params.tsv;
            R_TSV_rad = tsv.R_TSV_um * 1e-6;
            P_TSV = tsv.P_TSV_um * 1e-6;
            g_rad = tsv.g_rad_um * 1e-6;
            t_SOI = tsv.t_SOI_um * 1e-6;
            
            % Eq 316
            N_row = floor(w_ic / (2 * R_TSV_rad + g_rad));
            
            % Eq 317
            % Arc length available = r * beta_ic
            N_per_row = floor((r * beta_ic) / P_TSV);
            
            % Eq 318
            N_TSV = N_row * N_per_row;
            
            if N_TSV < 1
                N_TSV = 0;
                R_TSV_tot = 1e9; % High resistance if no TSVs fit
            else
                % Eq 306
                R_single = rho_TSV * t_SOI / (pi * R_TSV_rad^2);
                
                % Eq 311
                R_TSV_tot = R_single / N_TSV;
            end
        end
        
        function export_comsol_params(obj, filename)
            % EXPORT_COMSOL_PARAMS Exports parameters to a text file for COMSOL.
            
            fid = fopen(filename, 'w');
            if fid == -1
                error('Could not open file %s for writing.', filename);
            end
            
            fprintf(fid, '%% Radial TEC Parameters exported from MATLAB Solver\n');
            fprintf(fid, 'N_stages %d "Number of stages"\n', obj.N_stages);
            fprintf(fid, 'theta %f[deg] "Wedge angle"\n', rad2deg(obj.WedgeAngle));
            fprintf(fid, 't_tec %e[m] "Thickness"\n', obj.Thickness);
            fprintf(fid, 'r_start %e[m] "Start Radius"\n', obj.Params.r_start_um * 1e-6);
            fprintf(fid, 'L_element %e[m] "Element Length"\n', obj.Params.L_element_um * 1e-6);
            fprintf(fid, 'ratio_ic %f "Interconnect Ratio"\n', obj.Params.interconnect_ratio);
            fprintf(fid, 'ratio_oc %f "Outerconnect Ratio"\n', obj.Params.outerconnect_ratio);
            fprintf(fid, 'w_az %e[m] "Azimuthal Gap"\n', obj.Params.azimuthal_gap_um * 1e-6);
            fprintf(fid, 'w_is %e[m] "Insulation Width"\n', obj.Params.insulation_width_um * 1e-6);
            
            fclose(fid);
            fprintf('COMSOL parameters exported to %s\n', filename);
        end
    end
end
