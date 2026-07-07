function [results, next_run_idx] = run_single_validity_sweep(sweep_name, var_name, values, base_cases, output_dir, run_idx_offset, comsol_model_path)
% run_single_validity_sweep Executes a 1D vs 3D sweep for a single variable
%   sweep_name - String label for the sweep (e.g. 'HeatFlux')
%   var_name - Property to vary ('q_Wm2', 'M', or 't_TEC_um')
%   values - Array of values to sweep over
%   base_cases - Base configuration structure
%   output_dir - Directory to save CSV incrementally
%   run_idx_offset - Integer for global run numbering
%   comsol_model_path - Path to the mph file

    results = table();
    run_idx = run_idx_offset;
    total_runs = length(values);
    
    %% Helper to calculate L_1
    calc_L1 = @(R_cyl, w_is, f_L) ((10000e-6 / sqrt(2)) - (R_cyl*1e-6) - (3+1)*(w_is*1e-6)) * (1-f_L)/(1-f_L^3) * 1e6;

    for vi = 1:length(values)
        val = values(vi);
        fprintf('\n[%d/%d] %s Sweep: %s = %.2f\n', vi, total_runs, sweep_name, var_name, val);

        % Start with baseline, override current sweep variable
        q_now = base_cases.base_q;
        M_now = base_cases.base_M;
        t_now = base_cases.base_t_TEC;
        f_now = base_cases.f_L;
        w_is_now = base_cases.w_is_um;
        fill_factor_now = base_cases.fill_factor;
        t_chip_now = base_cases.t_chip_um;
        I_A_now = base_cases.I_A;

        if strcmp(var_name, 'q_Wm2'), q_now = val; end
        if strcmp(var_name, 'M'), M_now = val; end
        if strcmp(var_name, 't_TEC_um'), t_now = val; end
        if strcmp(var_name, 'f_L'), f_now = val; end
        if strcmp(var_name, 'w_is_um'), w_is_now = val; end
        if strcmp(var_name, 'fill_factor'), fill_factor_now = val; end
        if strcmp(var_name, 't_chip_um'), t_chip_now = val; end
        if strcmp(var_name, 'I_A'), I_A_now = val; end

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
            config.geometry.t_chip_um = t_chip_now;
            config.geometry.t_ins_um = base_cases.t_SOI_um;
            config.geometry.interconnect_ratio = base_cases.ic_w_r;
            config.geometry.outerconnect_ratio = base_cases.oc_w_r;
            config.geometry.insulation_width_um = w_is_now;
            config.geometry.interconnect_angle_ratio = base_cases.ic_angle_r;
            config.geometry.outerconnect_angle_ratio = base_cases.oc_angle_r;
            config.geometry.interconnect_thickness_ratio = base_cases.ic_t_r;
            config.geometry.outerconnect_thickness_ratio = base_cases.oc_t_r;
            config.geometry.fill_factor = fill_factor_now;

            config.operating_conditions.I_current_A = I_A_now;
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
            %% For older models this might have been +2, using +1 safely based on solver
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
            for k=1:min(3, length(ME.stack))
                fprintf('    In %s (line %d)\n', ME.stack(k).name, ME.stack(k).line);
            end
            matlab_T = NaN;
        end

        %% Run COMSOL
        try
            import com.comsol.model.*
            import com.comsol.model.util.*
            try
                model = mphload(comsol_model_path);
            catch
                fprintf('  Error loading model. Ensuring mphserver is ready...\n');
                % Give the server a moment if it just restarted
                pause(5);
                model = mphload(comsol_model_path);
            end

            model.param.set('LL_k_r', sprintf('%g', f_now));
            model.param.set('LL_L_1', sprintf('%g', L1_now));
            model.param.set('LL_R_cyl', sprintf('%g', base_cases.R_cyl_um));
            model.param.set('LL_t_chip', sprintf('%g', t_chip_now));
            model.param.set('LL_t_SOI', sprintf('%g', base_cases.t_SOI_um));
            model.param.set('LL_t_TEC', sprintf('%g', t_now));
            model.param.set('LL_theta', sprintf('%g', theta_deg));
            model.param.set('LL_w_is', sprintf('%g', w_is_now));
            model.param.set('q_i', sprintf('%g[W/m^2]', q_now));
            model.param.set('I_0', sprintf('%g[A]', I_A_now));
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

            ModelUtil.remove('Model'); % Clean up RAM

        catch ME
            fprintf('  COMSOL Error: %s\n', ME.message);
            comsol_T = NaN;
            
            % Attempt to recover Model if loaded but failed
            try, ModelUtil.remove('Model'); catch, end 
            
            fprintf('  Restarting COMSOL connection due to error...\n');
            system('taskkill /F /IM comsolmphserver.exe 2>nul', '-echo');
            pause(2);
            system(sprintf('start "" "%s" -port %d', 'F:\EngineeringSoftware\COMSOL\COMSOL63\Multiphysics\bin\win64\comsolmphserver.exe', 2036));
            pause(15);
            try
                mphstart(2036);
            catch
            end
        end

        %% Save to results
        err_abs = comsol_T - matlab_T;
        err_pct = 100 * abs(err_abs) / abs(comsol_T);

        new_row = table(run_idx, categorical({sweep_name}), categorical({var_name}), val, ...
            q_now, M_now, theta_deg, t_now, w_is_now, fill_factor_now, t_chip_now, ...
            matlab_T, comsol_T, err_abs, err_pct, ...
            'VariableNames', {'RunID', 'SweepType', 'Variable', 'Value', ...
            'q_Wm2', 'M', 'theta_deg', 't_TEC_um', 'w_is_um', 'fill_factor', 't_chip_um', ...
            'MATLAB_Tmax', 'COMSOL_Tmax', 'Error_Abs', 'Error_Pct'});

        results = [results; new_row];
        
        file_path = fullfile(output_dir, sprintf('%s_sweep_results.csv', sweep_name));
        writetable(results, file_path);
        
        run_idx = run_idx + 1;
    end
    
    %% Generate and Save Plots for the Sweep
    try
        x_label_latex = get_latex_x_label(var_name);

        % -------- Error-only plot --------
        fig_err = figure('Name', sprintf('%s Sweep Error', sweep_name), 'Visible', 'off');
        plot(results.Value, results.Error_Pct, '-o', 'LineWidth', 2, 'MarkerSize', 6, ...
            'Color', [0.85 0.33 0.10], 'MarkerFaceColor', [0.85 0.33 0.10]);
        ylabel('Relative Error (\%)', 'Interpreter', 'latex');
        xlabel(x_label_latex, 'Interpreter', 'latex');
        title(sprintf('1D vs 3D Validation Error: %s', strrep(sweep_name, '_', ' ')), ...
            'Interpreter', 'none');
        grid on;

        err_plot_path_png = fullfile(output_dir, sprintf('%s_sweep_error_plot.png', sweep_name));
        err_plot_path_fig = fullfile(output_dir, sprintf('%s_sweep_error_plot.fig', sweep_name));
        saveas(fig_err, err_plot_path_png);
        savefig(fig_err, err_plot_path_fig);
        close(fig_err);

        % -------- Temperature-only plot (left: MATLAB linear, right: COMSOL log) --------
        fig_temp = figure('Name', sprintf('%s Sweep Temperatures', sweep_name), 'Visible', 'off');
        ax = gca;

        yyaxis left
        plot(results.Value, results.MATLAB_Tmax, '-s', 'LineWidth', 1.8, 'MarkerSize', 6, ...
            'Color', [0.00 0.45 0.74], 'MarkerFaceColor', [0.00 0.45 0.74]);
        ylabel('MATLAB $T_{\max}$ ($^\circ$C)', 'Interpreter', 'latex', 'Color', [0.00 0.45 0.74]);
        ax.YAxis(1).Color = [0.00 0.45 0.74];
        ax.YAxis(1).Scale = 'linear';

        yyaxis right
        plot(results.Value, results.COMSOL_Tmax, '-d', 'LineWidth', 1.8, 'MarkerSize', 6, ...
            'Color', [0.85 0.33 0.10], 'MarkerFaceColor', [0.85 0.33 0.10]);
        ylabel('COMSOL $T_{\max}$ ($^\circ$C, log scale)', 'Interpreter', 'latex', 'Color', [0.85 0.33 0.10]);
        ax.YAxis(2).Color = [0.85 0.33 0.10];
        ax.YAxis(2).Scale = 'log';

        xlabel(x_label_latex, 'Interpreter', 'latex');
        title(sprintf('1D vs 3D Validation Temperatures: %s', strrep(sweep_name, '_', ' ')), ...
            'Interpreter', 'none');
        legend({'MATLAB $T_{\max}$', 'COMSOL $T_{\max}$'}, 'Interpreter', 'latex', 'Location', 'best');
        grid on;

        temp_plot_path_png = fullfile(output_dir, sprintf('%s_sweep_temperature_plot.png', sweep_name));
        temp_plot_path_fig = fullfile(output_dir, sprintf('%s_sweep_temperature_plot.fig', sweep_name));
        saveas(fig_temp, temp_plot_path_png);
        savefig(fig_temp, temp_plot_path_fig);
        close(fig_temp);

        fprintf('  Saved error plot to %s\n', err_plot_path_png);
        fprintf('  Saved temperature plot to %s\n', temp_plot_path_png);
    catch ME
        fprintf('  Warning: Could not save sweep plots: %s\n', ME.message);
    end
    
    next_run_idx = run_idx;
end

function x_label_latex = get_latex_x_label(var_name)
% Maps sweep variable names to publication-quality LaTeX axis labels.
    switch var_name
        case 'q_Wm2'
            x_label_latex = '$q_{\mathrm{i}}\;[\mathrm{W/m^2}]$';
        case 'M'
            x_label_latex = '$N$';
        case 't_TEC_um'
            x_label_latex = '$t_{\mathrm{TEC}}\;[\mu\mathrm{m}]$';
        case 'f_L'
            x_label_latex = '$f_L$';
        case 'w_is_um'
            x_label_latex = '$W_{\mathrm{is}}\;[\mu\mathrm{m}]$';
        case 'fill_factor'
            x_label_latex = '$f_{\mathrm{f}}$';
        case 't_chip_um'
            x_label_latex = '$t_{\mathrm{gen}}\;[\mu\mathrm{m}]$';
        otherwise
            x_label_latex = strrep(var_name, '_', '\_');
    end
end