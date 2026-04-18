function results = run_2d_validity_sweep(sweep_name, var_name1, values1, var_name2, values2, base_cases, output_dir, comsol_model_path)
% run_2d_validity_sweep Executes a 1D vs 3D sweep for two variables on a grid
%   sweep_name - String label for the sweep (e.g. 'HeatFlux_vs_Thickness')
%   var_name1 - Property 1 ('q_Wm2', 'M', or 't_TEC_um')
%   values1 - Array of values for Property 1
%   var_name2 - Property 2 ('q_Wm2', 'M', or 't_TEC_um')
%   values2 - Array of values for Property 2
%   base_cases - Base configuration structure
%   output_dir - Directory to save CSV and Plots
%   comsol_model_path - Path to the mph file

    results = table();
    run_idx = 1;
    total_runs = length(values1) * length(values2);
    
    %% Helper to calculate L_1
    calc_L1 = @(R_cyl, w_is, f_L) ((10000e-6 / sqrt(2)) - (R_cyl*1e-6) - (3+1)*(w_is*1e-6)) * (1-f_L)/(1-f_L^3) * 1e6;

    % Prepare matrices for 2D plotting. Rows correspond to values2 (Y-axis), cols correspond to values1 (X-axis)
    [V1, V2] = meshgrid(values1, values2);
    matlab_T_matrix = zeros(length(values2), length(values1));
    comsol_T_matrix = zeros(length(values2), length(values1));
    error_matrix = zeros(length(values2), length(values1));

    for v1_idx = 1:length(values1)
        for v2_idx = 1:length(values2)
            val1 = values1(v1_idx);
            val2 = values2(v2_idx);
            
            fprintf('\n[%d/%d] 2D %s: %s = %.2f, %s = %.2f\n', ...
                run_idx, total_runs, sweep_name, var_name1, val1, var_name2, val2);

            % Start with baseline, override current sweep variables
            q_now = base_cases.base_q;
            M_now = base_cases.base_M;
            t_now = base_cases.base_t_TEC;
            f_now = base_cases.f_L;
            w_is_now = base_cases.w_is_um;
            fill_factor_now = base_cases.fill_factor;

            if strcmp(var_name1, 'q_Wm2'), q_now = val1; end
            if strcmp(var_name1, 'M'), M_now = val1; end
            if strcmp(var_name1, 't_TEC_um'), t_now = val1; end
            if strcmp(var_name1, 'f_L'), f_now = val1; end
            if strcmp(var_name1, 'w_is_um'), w_is_now = val1; end
            if strcmp(var_name1, 'fill_factor'), fill_factor_now = val1; end

            if strcmp(var_name2, 'q_Wm2'), q_now = val2; end
            if strcmp(var_name2, 'M'), M_now = val2; end
            if strcmp(var_name2, 't_TEC_um'), t_now = val2; end
            if strcmp(var_name2, 'f_L'), f_now = val2; end
            if strcmp(var_name2, 'w_is_um'), w_is_now = val2; end
            if strcmp(var_name2, 'fill_factor'), fill_factor_now = val2; end

            theta_deg = 360 / M_now;
            L1_now = calc_L1(base_cases.R_cyl_um, w_is_now, f_now);

            %% Run MATLAB CTM
            try
                config = struct();
                config.geometry.N_stages = 3;
                config.geometry.wedge_angle_deg = theta_deg;
                config.geometry.thickness_um = t_now;
                config.geometry.radial_expansion_factor = f_now;
                config.geometry.w_chip_um = 10000;
                config.geometry.R_cyl_um = base_cases.R_cyl_um;
                config.geometry.t_chip_um = base_cases.t_chip_um;
                config.geometry.t_ins_um = base_cases.t_SOI_um;
                config.geometry.interconnect_ratio = base_cases.ic_w_r;
                config.geometry.outerconnect_ratio = base_cases.oc_w_r;
                config.geometry.insulation_width_um = w_is_now;
                config.geometry.interconnect_angle_ratio = base_cases.ic_angle_r;
                config.geometry.outerconnect_angle_ratio = base_cases.oc_angle_r;
                config.geometry.interconnect_thickness_ratio = base_cases.ic_t_r;
                config.geometry.outerconnect_thickness_ratio = base_cases.oc_t_r;
                config.geometry.fill_factor = fill_factor_now;

                config.operating_conditions.I_current_A = base_cases.I_A;
                config.boundary_conditions.q_flux_W_m2 = q_now;
                config.boundary_conditions.T_water_K = base_cases.T_water;
                config.boundary_conditions.h_conv_W_m2K = 1e6;

                config.materials.Bi2Te3 = struct('k', 1.6, 'rho', 1.15e-5, 'S', 210e-6);
                config.materials.Cu    = struct('k', 400, 'rho', 1.667e-8);
                config.materials.Si    = struct('k', 130, 'rho', 0.01);
                config.materials.AlN   = struct('k', 170, 'rho', 1e10);
                config.materials.SiO2  = struct('k', 1.4, 'rho', 1e14);
                config.materials.Al2O3 = struct('k', 35, 'rho', 1e12);

                materials = MaterialProperties(config);
                geometry = TECGeometry(config);
                network = ThermalNetwork(geometry, materials, config);

                dim = 2 * config.geometry.N_stages + 1;
                T = ones(dim, 1) * 300;
                for iter = 1:100
                    T_old = T(1:length(network.solve(T)));
                    [T, Q_out, Q_in] = network.solve(T);
                    if max(abs(T - T_old)) < 1e-6
                        break;
                    end
                end
                matlab_T = max(T) - 273.15;
                fprintf('  MATLAB T_max = %.2f C\n', matlab_T);
            catch ME
                fprintf('  MATLAB Error: %s\n', ME.message);
                matlab_T = NaN;
            end

            %% Run COMSOL
            try
                import com.comsol.model.*
                import com.comsol.model.util.*
                try
                    model = mphload(comsol_model_path);
                catch
                    pause(5);
                    model = mphload(comsol_model_path);
                end

                model.param.set('LL_k_r', sprintf('%g', f_now));
                model.param.set('LL_L_1', sprintf('%g', L1_now));
                model.param.set('LL_R_cyl', sprintf('%g', base_cases.R_cyl_um));
                model.param.set('LL_t_chip', sprintf('%g', base_cases.t_chip_um));
                model.param.set('LL_t_SOI', sprintf('%g', base_cases.t_SOI_um));
                model.param.set('LL_t_TEC', sprintf('%g', t_now));
                model.param.set('LL_theta', sprintf('%g', theta_deg));
                model.param.set('LL_w_is', sprintf('%g', w_is_now));
                model.param.set('q_i', sprintf('%g[W/m^2]', q_now));
                model.param.set('I_0', sprintf('%g[A]', base_cases.I_A));
                model.param.set('LL_fill_factor', sprintf('%g', fill_factor_now));
                model.param.set('LL_ic_angle_r', sprintf('%g', base_cases.ic_angle_r));
                model.param.set('LL_ic_t_r', sprintf('%g', base_cases.ic_t_r));
                model.param.set('LL_ic_w_r', sprintf('%g', base_cases.ic_w_r));
                model.param.set('LL_oc_angle_r', sprintf('%g', base_cases.oc_angle_r));
                model.param.set('LL_oc_t_r', sprintf('%g', base_cases.oc_t_r));
                model.param.set('LL_oc_w_r', sprintf('%g', base_cases.oc_w_r));

                model.study('std1').run();
                try
                    T_all = mpheval(model, 'T', 'dataset', 'dset1');
                    comsol_T = max(T_all.d1) - 273.15;
                catch
                    comsol_T = NaN;
                end
                fprintf('  COMSOL T_max = %.2f C\n', comsol_T);

                ModelUtil.remove('Model'); 

            catch ME
                fprintf('  COMSOL Error: %s\n', ME.message);
                comsol_T = NaN;
                try, ModelUtil.remove('Model'); catch, end 
                
                fprintf('  Restarting COMSOL connection due to error...\n');
                system('taskkill /F /IM comsolmphserver.exe 2>nul', '-echo');
                pause(2);
                system(sprintf('start "" "%s" -port %d', 'F:\EngineeringSoftware\COMSOL\COMSOL63\Multiphysics\bin\win64\comsolmphserver.exe', 2036));
                pause(15);
                try, mphstart(2036); catch, end
            end

            %% Save to results & matrix
            err_abs = comsol_T - matlab_T;
            err_pct = 100 * abs(err_abs) / abs(comsol_T);

            matlab_T_matrix(v2_idx, v1_idx) = matlab_T;
            comsol_T_matrix(v2_idx, v1_idx) = comsol_T;
            error_matrix(v2_idx, v1_idx) = err_pct;

            new_row = table(run_idx, categorical({sweep_name}), val1, val2, ...
                q_now, M_now, theta_deg, t_now, w_is_now, fill_factor_now, ...
                matlab_T, comsol_T, err_abs, err_pct, ...
                'VariableNames', {'RunID', 'SweepType', ['Var1_', var_name1], ['Var2_', var_name2], ...
                'q_Wm2', 'M', 'theta_deg', 't_TEC_um', 'w_is_um', 'fill_factor', ...
                'MATLAB_Tmax', 'COMSOL_Tmax', 'Error_Abs', 'Error_Pct'});

            results = [results; new_row];
            
            % Save incremental CSV
            file_path = fullfile(output_dir, sprintf('%s_2D_results.csv', sweep_name));
            writetable(results, file_path);
            
            run_idx = run_idx + 1;
        end
    end
    
    %% Generate 2D Plots
    label1 = get_latex_param_label(var_name1);
    label2 = get_latex_param_label(var_name2);

    try
        % 1. Error Plot
        fig_err = figure('Name', sprintf('%s 2D Error', sweep_name), 'Visible', 'off');
        contourf(V1, V2, error_matrix, 20, 'LineStyle', 'none'); 
        colorbar; colormap jet;
        xlabel(label1, 'Interpreter', 'latex');
        ylabel(label2, 'Interpreter', 'latex');
        title(sprintf('Relative Error (\\%%): %s', strrep(sweep_name, '_', ' ')), 'Interpreter', 'none');
        saveas(fig_err, fullfile(output_dir, sprintf('%s_Error_plot.png', sweep_name)));
        savefig(fig_err, fullfile(output_dir, sprintf('%s_Error_plot.fig', sweep_name)));
        close(fig_err);

        % 2. Side-by-side temperature comparison (CTM vs COMSOL)
        fig_temp = figure('Name', sprintf('%s 2D Temperature Comparison', sweep_name), 'Visible', 'off');

        tiledlayout(1, 2, 'Padding', 'compact', 'TileSpacing', 'compact');

        ax1 = nexttile;
        contourf(ax1, V1, V2, matlab_T_matrix, 20, 'LineStyle', 'none');
        colormap(ax1, hot);
        cb1 = colorbar(ax1);
        cb1.Label.String = 'CTM $T_{\\max}$ ($^\\circ$C)';
        cb1.Label.Interpreter = 'latex';
        xlabel(ax1, label1, 'Interpreter', 'latex');
        ylabel(ax1, label2, 'Interpreter', 'latex');
        title(ax1, 'CTM $T_{\\max}$', 'Interpreter', 'latex');
        grid(ax1, 'on');

        ax2 = nexttile;
        contourf(ax2, V1, V2, comsol_T_matrix, 20, 'LineStyle', 'none');
        colormap(ax2, hot);
        cb2 = colorbar(ax2);
        cb2.Label.String = 'COMSOL $T_{\\max}$ ($^\\circ$C)';
        cb2.Label.Interpreter = 'latex';
        xlabel(ax2, label1, 'Interpreter', 'latex');
        ylabel(ax2, label2, 'Interpreter', 'latex');
        title(ax2, 'COMSOL $T_{\\max}$', 'Interpreter', 'latex');
        grid(ax2, 'on');

        sgtitle(sprintf('2D Temperature Comparison: %s', strrep(sweep_name, '_', ' ')), 'Interpreter', 'none');

        saveas(fig_temp, fullfile(output_dir, sprintf('%s_TemperatureComparison_plot.png', sweep_name)));
        savefig(fig_temp, fullfile(output_dir, sprintf('%s_TemperatureComparison_plot.fig', sweep_name)));
        close(fig_temp);

        fprintf('\n  Saved 2D plots for %s\n', sweep_name);
    catch ME
        fprintf('\n  Warning: Could not save 2D plots: %s\n', ME.message);
    end
end

function label = get_latex_param_label(var_name)
% Returns publication-ready LaTeX axis labels for sweep variables.
    switch var_name
        case 'q_Wm2'
            label = '$q_{\mathrm{i}}\;[\mathrm{W/m^2}]$';
        case 'M'
            label = '$N$';
        case 't_TEC_um'
            label = '$t_{\mathrm{TEC}}\;[\mu\mathrm{m}]$';
        case 'f_L'
            label = '$f_{L}$';
        case 'fill_factor'
            label = '$f_{\mathrm{f}}$';
        case 'w_is_um'
            label = '$W_{\mathrm{is}}\;[\mu\mathrm{m}]$';
        case 't_chip_um'
            label = '$t_{\mathrm{gen}}\;[\mu\mathrm{m}]$';
        otherwise
            label = strrep(var_name, '_', '\_');
    end
end