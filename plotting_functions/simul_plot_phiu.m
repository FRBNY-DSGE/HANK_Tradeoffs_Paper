% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % simul_plot_phiu.m
% %
% % Description: Given simulation results, file recreates figure A4 in the
% % paper
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

close all;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% (1.2) Set up output file paths
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

plt_subfolder = strcat('ELB=',simul_ELB);

pltroot.frontierplot_phiu = plt_folder + "/Simulations/" + plt_subfolder + "/Phi_u/FrontierPlots/";


manual_std = @(x) sqrt(sum(x .^2)/length(x));



%Load in Taylor
for vv=1:length(phipi1)
    pp=phipi1(vv);
    filepath1 = results_folder + "/precovid_ugap_v4_pol=" + pol1 + "_shocks=" + simul_str +"_phiu=0_2_ELB=" + simul_ELB+ "_all_phipi=" + pp + "_kappa=" + kappa_val1 + ".mat";
    filepath2 = results_folder + "/precovid_ugap_v4_pol=" + pol1 + "_shocks=" + simul_str +"_ELB=" + simul_ELB+ "_all_phipi=" + pp + "_kappa=" + kappa_val1 + ".mat";
    filepath3 = results_folder + "/precovid_ugap_v4_pol=" + pol1 + "_shocks=" + simul_str +"_phiu=0_5_ELB=" + simul_ELB+ "_all_phipi=" + pp + "_kappa=" + kappa_val1 + ".mat";
    var_series = 'X_series';

    Pol1_results0_1(vv) = load(filepath1).results;
    Pol1_results0_3(vv) = load(filepath2).results;
    Pol1_results0_5(vv) = load(filepath3).results;


    pol1_phi_u_calculations_helper;
end

%%

frontier_plots_phi_u(pol1_01_stds, pol1_03_stds, pol1_05_stds, pol1 ,phipi1_vals, baseline_phi_ind, string(fieldnames(pol1_03_stds)), "SS_Dev", pltroot);


%%%%%% Defining Functions %%%%%%%%%%%%%

function frontier_plots_phi_u(pol1_01_std, pol1_03_std, pol1_05_std, pol1,phi_pi_taylor, baseline_phi_ind, fields, type, pltroot)
%asserting that fields is of type string(array)

%Inflation relative to Taylor Baseline:
%Do this outside of the loop, it doesn't change:
if strcmp(type, "Mean")
    pol1_01_infl = getfield(pol1_01_std, "pi_mean");
    pol1_03_infl = getfield(pol1_03_std, "pi_mean");
    pol1_05_infl = getfield(pol1_05_std, "pi_mean");
else
    pol1_01_infl = getfield(pol1_01_std, "pi_std");
    pol1_03_infl = getfield(pol1_03_std, "pi_std");
    pol1_05_infl = getfield(pol1_05_std, "pi_std");
end


pol1_01_infl_norm = 100 * (pol1_01_infl - pol1_03_infl(baseline_phi_ind(1)));
pol1_03_infl_norm = 100 * (pol1_03_infl - pol1_03_infl(baseline_phi_ind(1)));
pol1_05_infl_norm = 100 * (pol1_05_infl - pol1_03_infl(baseline_phi_ind(1)));



% Point label offset for plotting:
dx = 0.0001;
dy = 0.0001;


for f = 1:length(fields)


    char_f = char(fields(f));
    para_name = string(char_f(1:end-4)); %removes '_std' from string name
    x_lab = "pi";


    %Variable relative to Taylor Baseline:
    pol1_01_vals = getfield(pol1_01_std, fields(f));
    pol1_03_vals = getfield(pol1_03_std, fields(f));
    pol1_05_vals = getfield(pol1_05_std, fields(f));


    pol1_01_vals_norm = 100 * (pol1_01_vals - pol1_03_vals(baseline_phi_ind(1)));
    pol1_03_vals_norm = 100 * (pol1_03_vals - pol1_03_vals(baseline_phi_ind(1)));
    pol1_05_vals_norm = 100 * (pol1_05_vals - pol1_03_vals(baseline_phi_ind(1)));
    %%%%%% Actual Plotting %%%%%%%%%
    figure('name', strcat(para_name,"_frontier"));

    hold on
    plot(pol1_01_infl_norm, pol1_01_vals_norm, "--o", 'Color', "#AED6F1", "LineWidth", 1.5);
    plot(pol1_03_infl_norm, pol1_03_vals_norm, "--o", 'Color', "#3498DB", "LineWidth", 1.5);
    plot(pol1_05_infl_norm, pol1_05_vals_norm, "--o", 'Color', "#2471A3", "LineWidth", 1.5);

    yline(0);
    xline(0);

    text(pol1_01_infl_norm+dx, pol1_01_vals_norm+dy, string(phi_pi_taylor));
    text(pol1_03_infl_norm+dx, pol1_03_vals_norm+dy, string(phi_pi_taylor));
    text(pol1_05_infl_norm+dx, pol1_05_vals_norm+dy, string(phi_pi_taylor));

    hold off

    legend(["\phi_{u} = 0.2",  "\phi_{u} = 0.36", "\phi_{u} = 0.5"])
    if strcmp(type, "Mean")
        title(strrep(para_name, "_", "-") + " v.s " + x_lab + ' (p.p )' );
        xlabel('Mean of ' + x_lab + ' (in p.p)');
        ylabel('Mean (in p.p)');
    else
        title(strrep(para_name, "_", "-") + " v.s " + x_lab + ' (p.p deviation from baseline volatility)' );
        xlabel('Deviation of ' + x_lab + ' from baseline (in p.p)');
        ylabel('Deviation from baseline (in p.p)');
    end


    fi = gcf;

    if strcmp(type, "Mean")
        base_fl = pltroot.frontierplots_mean;
    else
        base_fl = pltroot.frontierplot_phiu;
    end
    if ~exist(base_fl, 'dir')
        mkdir(base_fl)
    end
    exportgraphics(fi, base_fl + type + '_phi_u_' + para_name + '.png', 'Resolution', 300);



end

end





