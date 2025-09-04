% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % simul_plot_tfp_mk.m
% %
% % Description: Given simulation results, file recreates
% % (1) disaggregation of aggregate and distributional volatility vs phi_pi
% % parameter by markup and tfp shocks
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% (1.1) Set up output file paths
%%%%%%%%%%

% Output file folders: Root folder/sub-root folder/[changable final path]
plt_subfolder = strcat('ELB=',simul_ELB);
% For simulations run with baseline kappa parameter value
pltroot.disagg = plt_folder + "/Simulations/" + plt_subfolder + "/TFP_MK_Disagg/";


%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% (1.2) Volatility calculation
%%%%%%%%%%

% For each variable, volatility (labeled 'SS_Dev') is calculated either
% using the standard MATLAB std() function or a root mean-squared
% calculation, referred to henceforth as "manual." For variables with
% mean-zero deviation from their steady state value, both calculations
% should yield similar plots. Otherwise, volatility plots are made using
% the manual calculation.

manual_std = @(x) sqrt(sum(x .^2)/length(x));

%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% (2) Loading simulation results: Non-Disaggregated
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% (2.1) Loads results using Taylor policy rule based on settings
%%%%%%%%%%

for vv=1:length(phipi1)
    pp=phipi1(vv);
    filepath1 = results_folder + "/precovid_ugap_v4_pol=" + pol1 + "_shocks=" + simul_str +"_ELB=" + simul_ELB+ "_supply_phipi=" + pp + "_kappa=" + kappa_val1 + ".mat";
    var_series = 'X_series';


    Pol1_results(vv) = load(filepath1).results;


    pol1_ELB_freq(vv) = sum(exp(Pol1_results(vv).(var_series)(env.grid.numstates - ...
        env.grid.os + ...
        env.idx.R_cb_ind,:))*1.01 ...
        <=1+1e-8) / length(Pol1_results(vv).(var_series)(env.grid.numstates - ...
        env.grid.os + ...
        env.idx.R_cb_ind,:));

    % Calculate variable means:
    pol1_means.C_10_mean(vv) = mean(Pol1_results(vv).(var_series)(end - env.grid.oc + env.idx.C_10_ind,:));
    pol1_means.C_Q1_mean(vv) = mean(Pol1_results(vv).(var_series)(end - env.grid.oc + env.idx.C_Q1_ind,:));
    pol1_means.C_Q2_mean(vv) = mean(Pol1_results(vv).(var_series)(end - env.grid.oc + env.idx.C_Q2_ind,:));
    pol1_means.C_Q3_mean(vv) = mean(Pol1_results(vv).(var_series)(end - env.grid.oc + env.idx.C_Q3_ind,:));
    pol1_means.C_Q4_mean(vv) = mean(Pol1_results(vv).(var_series)(end - env.grid.oc + env.idx.C_Q4_ind,:));
    pol1_means.C_Q5_mean(vv) = mean(Pol1_results(vv).(var_series)(end - env.grid.oc + env.idx.C_Q5_ind,:));
    pol1_means.C_90_mean(vv) = mean(Pol1_results(vv).(var_series)(end - env.grid.oc + env.idx.C_90_ind,:));


    pol1_means.W_10_mean(vv) = mean(Pol1_results(vv).(var_series)(end - env.grid.oc + env.idx.w_10_ind,:));
    pol1_means.W_Q1_mean(vv) = mean(Pol1_results(vv).(var_series)(end - env.grid.oc + env.idx.w_Q1_ind,:));
    pol1_means.W_Q2_mean(vv) = mean(Pol1_results(vv).(var_series)(end - env.grid.oc + env.idx.w_Q2_ind,:));
    pol1_means.W_Q3_mean(vv) = mean(Pol1_results(vv).(var_series)(end - env.grid.oc + env.idx.w_Q3_ind,:));
    pol1_means.W_Q4_mean(vv) = mean(Pol1_results(vv).(var_series)(end - env.grid.oc + env.idx.w_Q4_ind,:));
    pol1_means.W_Q5_mean(vv) = mean(Pol1_results(vv).(var_series)(end - env.grid.oc + env.idx.w_Q5_ind,:));
    pol1_means.W_90_mean(vv) = mean(Pol1_results(vv).(var_series)(end - env.grid.oc + env.idx.w_90_ind,:));

    pol1_means.C_I_10_mean(vv) = mean(Pol1_results(vv).(var_series)(end - env.grid.oc + env.idx.C_I_10_ind,:));
    pol1_means.C_I_Q1_mean(vv) = mean(Pol1_results(vv).(var_series)(end - env.grid.oc + env.idx.C_I_Q1_ind,:));
    pol1_means.C_I_Q2_mean(vv) = mean(Pol1_results(vv).(var_series)(end - env.grid.oc + env.idx.C_I_Q2_ind,:));
    pol1_means.C_I_Q3_mean(vv) = mean(Pol1_results(vv).(var_series)(end - env.grid.oc + env.idx.C_I_Q3_ind,:));
    pol1_means.C_I_Q4_mean(vv) = mean(Pol1_results(vv).(var_series)(end - env.grid.oc + env.idx.C_I_Q4_ind,:));
    pol1_means.C_I_Q5_mean(vv) = mean(Pol1_results(vv).(var_series)(end - env.grid.oc + env.idx.C_I_Q5_ind,:));
    pol1_means.C_I_90_mean(vv) = mean(Pol1_results(vv).(var_series)(end - env.grid.oc + env.idx.C_I_90_ind,:));
    pol1_means.C_mean(vv)         = mean(Pol1_results(vv).(var_series)(end - env.grid.oc + env.idx.C_ind,:) );
    pol1_means.Y_mean(vv)         = mean(Pol1_results(vv).(var_series)(end - env.grid.oc + env.idx.Y_ind,:) );
    pol1_means.pi_mean(vv)        = mean(Pol1_results(vv).(var_series)(end - env.grid.oc + env.idx.pi_ind,:) );
    pol1_means.unemp_mean(vv)     = mean(Pol1_results(vv).(var_series)(end - env.grid.oc + env.idx.unemp_ind,:));
    pol1_means.RR_mean(vv)        = mean(Pol1_results(vv).(var_series)(end - env.grid.oc + env.idx.RR_ind,:));
    pol1_means.Profit_mean(vv)    = mean(Pol1_results(vv).(var_series)(end - env.grid.oc + env.idx.Profit_ind,:));
    pol1_means.Profit_FI_mean(vv) = mean(Pol1_results(vv).(var_series)(end - env.grid.oc + env.idx.unemp_ind,:));

    %scalars:
    pol1_means.Q_mean(vv) = mean(Pol1_results(vv).(var_series)(env.grid.numstates - env.grid.os + env.idx.Q_ind,:));
    pol1_means.R_cb_mean(vv) = mean(Pol1_results(vv).(var_series)(env.grid.numstates - env.grid.os + env.idx.R_cb_ind,:));
    pol1_means.W_mean(vv) = mean(Pol1_results(vv).(var_series)(env.grid.numstates - env.grid.os + env.idx.w_ind,:));


    % Calculate variable standard deviations
    pol1_stds.C_10_std(vv) = std(Pol1_results(vv).(var_series)(end - env.grid.oc + env.idx.C_10_ind,:));
    pol1_stds.C_Q1_std(vv) = std(Pol1_results(vv).(var_series)(end - env.grid.oc + env.idx.C_Q1_ind,:));
    pol1_stds.C_Q2_std(vv) = std(Pol1_results(vv).(var_series)(end - env.grid.oc + env.idx.C_Q2_ind,:));
    pol1_stds.C_Q3_std(vv) = std(Pol1_results(vv).(var_series)(end - env.grid.oc + env.idx.C_Q3_ind,:));
    pol1_stds.C_Q4_std(vv) = std(Pol1_results(vv).(var_series)(end - env.grid.oc + env.idx.C_Q4_ind,:));
    pol1_stds.C_Q5_std(vv) = std(Pol1_results(vv).(var_series)(end - env.grid.oc + env.idx.C_Q5_ind,:));
    pol1_stds.C_90_std(vv) = std(Pol1_results(vv).(var_series)(end - env.grid.oc + env.idx.C_90_ind,:));

    pol1_stds.W_10_std(vv) = std(Pol1_results(vv).(var_series)(end - env.grid.oc + env.idx.w_10_ind,:));
    pol1_stds.W_Q1_std(vv) = std(Pol1_results(vv).(var_series)(end - env.grid.oc + env.idx.w_Q1_ind,:));
    pol1_stds.W_Q2_std(vv) = std(Pol1_results(vv).(var_series)(end - env.grid.oc + env.idx.w_Q2_ind,:));
    pol1_stds.W_Q3_std(vv) = std(Pol1_results(vv).(var_series)(end - env.grid.oc + env.idx.w_Q3_ind,:));
    pol1_stds.W_Q4_std(vv) = std(Pol1_results(vv).(var_series)(end - env.grid.oc + env.idx.w_Q4_ind,:));
    pol1_stds.W_Q5_std(vv) = std(Pol1_results(vv).(var_series)(end - env.grid.oc + env.idx.w_Q5_ind,:));
    pol1_stds.W_90_std(vv) = std(Pol1_results(vv).(var_series)(end - env.grid.oc + env.idx.w_90_ind,:));



    pol1_stds.C_I_10_std(vv) = std(Pol1_results(vv).(var_series)(end - env.grid.oc + env.idx.C_I_10_ind,:));
    pol1_stds.C_I_Q1_std(vv) = std(Pol1_results(vv).(var_series)(end - env.grid.oc + env.idx.C_I_Q1_ind,:));
    pol1_stds.C_I_Q2_std(vv) = std(Pol1_results(vv).(var_series)(end - env.grid.oc + env.idx.C_I_Q2_ind,:));
    pol1_stds.C_I_Q3_std(vv) = std(Pol1_results(vv).(var_series)(end - env.grid.oc + env.idx.C_I_Q3_ind,:));
    pol1_stds.C_I_Q4_std(vv) = std(Pol1_results(vv).(var_series)(end - env.grid.oc + env.idx.C_I_Q4_ind,:));
    pol1_stds.C_I_Q5_std(vv) = std(Pol1_results(vv).(var_series)(end - env.grid.oc + env.idx.C_I_Q5_ind,:));
    pol1_stds.C_I_90_std(vv) = std(Pol1_results(vv).(var_series)(end - env.grid.oc + env.idx.C_I_90_ind,:));
    pol1_stds.C_std(vv)         = std(Pol1_results(vv).(var_series)(end - env.grid.oc + env.idx.C_ind,:) );
    pol1_stds.Y_std(vv)         = std(Pol1_results(vv).(var_series)(end - env.grid.oc + env.idx.Y_ind,:) );
    pol1_stds.pi_std(vv)        = std(Pol1_results(vv).(var_series)(end - env.grid.oc + env.idx.pi_ind,:) );
    pol1_stds.unemp_std(vv)     = std(Pol1_results(vv).(var_series)(end - env.grid.oc + env.idx.unemp_ind,:));
    pol1_stds.RR_std(vv)        = std(Pol1_results(vv).(var_series)(end - env.grid.oc + env.idx.RR_ind,:));
    pol1_stds.Profit_std(vv)    = std(Pol1_results(vv).(var_series)(end - env.grid.oc + env.idx.Profit_ind,:));
    pol1_stds.Profit_FI_std(vv) = std(Pol1_results(vv).(var_series)(end - env.grid.oc + env.idx.unemp_ind,:));

    %scalars:
    pol1_stds.Q_std(vv) = std(Pol1_results(vv).(var_series)(env.grid.numstates - env.grid.os + env.idx.Q_ind,:));
    pol1_stds.R_cb_std(vv) = std(Pol1_results(vv).(var_series)(env.grid.numstates - env.grid.os + env.idx.R_cb_ind,:));
    pol1_stds.W_std(vv) = std(Pol1_results(vv).(var_series)(env.grid.numstates - env.grid.os + env.idx.w_ind,:));


    %Manual Standard Deviations:
    pol1_stds2.C_10_std(vv) = manual_std(Pol1_results(vv).(var_series)(end - env.grid.oc + env.idx.C_10_ind,:));
    pol1_stds2.C_Q1_std(vv) = manual_std(Pol1_results(vv).(var_series)(end - env.grid.oc + env.idx.C_Q1_ind,:));
    pol1_stds2.C_Q2_std(vv) = manual_std(Pol1_results(vv).(var_series)(end - env.grid.oc + env.idx.C_Q2_ind,:));
    pol1_stds2.C_Q3_std(vv) = manual_std(Pol1_results(vv).(var_series)(end - env.grid.oc + env.idx.C_Q3_ind,:));
    pol1_stds2.C_Q4_std(vv) = manual_std(Pol1_results(vv).(var_series)(end - env.grid.oc + env.idx.C_Q4_ind,:));
    pol1_stds2.C_Q5_std(vv) = manual_std(Pol1_results(vv).(var_series)(end - env.grid.oc + env.idx.C_Q5_ind,:));
    pol1_stds2.C_90_std(vv) = manual_std(Pol1_results(vv).(var_series)(end - env.grid.oc + env.idx.C_90_ind,:));


    pol1_stds2.W_10_std(vv) = manual_std(Pol1_results(vv).(var_series)(end - env.grid.oc + env.idx.w_10_ind,:));
    pol1_stds2.W_Q1_std(vv) = manual_std(Pol1_results(vv).(var_series)(end - env.grid.oc + env.idx.w_Q1_ind,:));
    pol1_stds2.W_Q2_std(vv) = manual_std(Pol1_results(vv).(var_series)(end - env.grid.oc + env.idx.w_Q2_ind,:));
    pol1_stds2.W_Q3_std(vv) = manual_std(Pol1_results(vv).(var_series)(end - env.grid.oc + env.idx.w_Q3_ind,:));
    pol1_stds2.W_Q4_std(vv) = manual_std(Pol1_results(vv).(var_series)(end - env.grid.oc + env.idx.w_Q4_ind,:));
    pol1_stds2.W_Q5_std(vv) = manual_std(Pol1_results(vv).(var_series)(end - env.grid.oc + env.idx.w_Q5_ind,:));
    pol1_stds2.W_90_std(vv) = manual_std(Pol1_results(vv).(var_series)(end - env.grid.oc + env.idx.w_90_ind,:));

    pol1_stds2.C_I_10_std(vv) = manual_std(Pol1_results(vv).(var_series)(end - env.grid.oc + env.idx.C_I_10_ind,:));
    pol1_stds2.C_I_Q1_std(vv) = manual_std(Pol1_results(vv).(var_series)(end - env.grid.oc + env.idx.C_I_Q1_ind,:));
    pol1_stds2.C_I_Q2_std(vv) = manual_std(Pol1_results(vv).(var_series)(end - env.grid.oc + env.idx.C_I_Q2_ind,:));
    pol1_stds2.C_I_Q3_std(vv) = manual_std(Pol1_results(vv).(var_series)(end - env.grid.oc + env.idx.C_I_Q3_ind,:));
    pol1_stds2.C_I_Q4_std(vv) = manual_std(Pol1_results(vv).(var_series)(end - env.grid.oc + env.idx.C_I_Q4_ind,:));
    pol1_stds2.C_I_Q5_std(vv) = manual_std(Pol1_results(vv).(var_series)(end - env.grid.oc + env.idx.C_I_Q5_ind,:));
    pol1_stds2.C_I_90_std(vv) = manual_std(Pol1_results(vv).(var_series)(end - env.grid.oc + env.idx.C_I_90_ind,:));
    pol1_stds2.C_std(vv)         = manual_std(Pol1_results(vv).(var_series)(end - env.grid.oc + env.idx.C_ind,:) );
    pol1_stds2.Y_std(vv)         = manual_std(Pol1_results(vv).(var_series)(end - env.grid.oc + env.idx.Y_ind,:) );
    pol1_stds2.pi_std(vv)        = manual_std(Pol1_results(vv).(var_series)(end - env.grid.oc + env.idx.pi_ind,:) );
    pol1_stds2.unemp_std(vv)     = manual_std(Pol1_results(vv).(var_series)(end - env.grid.oc + env.idx.unemp_ind,:));
    pol1_stds2.RR_std(vv)        = manual_std(Pol1_results(vv).(var_series)(end - env.grid.oc + env.idx.RR_ind,:));
    pol1_stds2.Profit_std(vv)    = manual_std(Pol1_results(vv).(var_series)(end - env.grid.oc + env.idx.Profit_ind,:));
    pol1_stds2.Profit_FI_std(vv) = manual_std(Pol1_results(vv).(var_series)(end - env.grid.oc + env.idx.unemp_ind,:));

    %scalars:
    pol1_stds2.Q_std(vv) = manual_std(Pol1_results(vv).(var_series)(env.grid.numstates - env.grid.os + env.idx.Q_ind,:));
    pol1_stds2.R_cb_std(vv) = manual_std(Pol1_results(vv).(var_series)(env.grid.numstates - env.grid.os + env.idx.R_cb_ind,:));
    pol1_stds2.W_std(vv) = manual_std(Pol1_results(vv).(var_series)(env.grid.numstates - env.grid.os + env.idx.w_ind,:));

    pol1_means.R_a_mean(vv)     = mean(Pol1_results(vv).(var_series)(end - env.grid.oc + env.idx.r_a_ind,:));
    pol1_stds.R_a_std(vv)     = std(Pol1_results(vv).(var_series)(end - env.grid.oc + env.idx.r_a_ind,:));
    pol1_stds2.R_a_std(vv)     = manual_std(Pol1_results(vv).(var_series)(end - env.grid.oc + env.idx.r_a_ind,:));


    pol1_means.RR_new_mean(vv)        = mean(Pol1_results(vv).(var_series)(env.grid.numstates - env.grid.os + env.idx.R_cb_ind,:) - Pol1_results(vv).(var_series)(end - env.grid.oc + env.idx.pi_ind,:) );
    pol1_stds.RR_new_std(vv)        = std(Pol1_results(vv).(var_series)(env.grid.numstates - env.grid.os + env.idx.R_cb_ind,:) - Pol1_results(vv).(var_series)(end - env.grid.oc + env.idx.pi_ind,:) );
    pol1_stds2.RR_new_std(vv)        = manual_std(Pol1_results(vv).(var_series)(env.grid.numstates - env.grid.os + env.idx.R_cb_ind,:) - Pol1_results(vv).(var_series)(end - env.grid.oc + env.idx.pi_ind,:) );
end


%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% (3) Loading simulation results: Disaggregated
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    for vv=1:length(phipi1)
        pp=phipi1(vv);
        filepath_pol1_tfp = results_folder + "/precovid_ugap_v4_pol=" + pol1 + "_shocks=" + simul_str_2 +"_ELB=" + simul_ELB+ "_all_phipi=" + pp + "_kappa=0_054.mat";
        filepath_pol1_mk = results_folder + "/precovid_ugap_v4_pol=" + pol1 + "_shocks=" + simul_str_3 +"_ELB=" + simul_ELB+ "_all_phipi=" + pp + "_kappa=0_054.mat";

        vers_pol1_tfp(vv) = load(filepath_pol1_tfp).results;
        vers_pol1_mk(vv) = load(filepath_pol1_mk).results;

        %%% Policy 1:

        pol1_tfp_means.C_10_mean(vv) = mean(vers_pol1_tfp(vv).(var_series)(end - env.grid.oc + env.idx.C_10_ind,:));
        pol1_tfp_means.C_Q1_mean(vv) = mean(vers_pol1_tfp(vv).(var_series)(end - env.grid.oc + env.idx.C_Q1_ind,:));
        pol1_tfp_means.C_Q2_mean(vv) = mean(vers_pol1_tfp(vv).(var_series)(end - env.grid.oc + env.idx.C_Q2_ind,:));
        pol1_tfp_means.C_Q3_mean(vv) = mean(vers_pol1_tfp(vv).(var_series)(end - env.grid.oc + env.idx.C_Q3_ind,:));
        pol1_tfp_means.C_Q4_mean(vv) = mean(vers_pol1_tfp(vv).(var_series)(end - env.grid.oc + env.idx.C_Q4_ind,:));
        pol1_tfp_means.C_Q5_mean(vv) = mean(vers_pol1_tfp(vv).(var_series)(end - env.grid.oc + env.idx.C_Q5_ind,:));
        pol1_tfp_means.C_90_mean(vv) = mean(vers_pol1_tfp(vv).(var_series)(end - env.grid.oc + env.idx.C_90_ind,:));

        pol1_tfp_means.C_I_10_mean(vv) = mean(vers_pol1_tfp(vv).(var_series)(end - env.grid.oc + env.idx.C_I_10_ind,:));
        pol1_tfp_means.C_I_Q1_mean(vv) = mean(vers_pol1_tfp(vv).(var_series)(end - env.grid.oc + env.idx.C_I_Q1_ind,:));
        pol1_tfp_means.C_I_Q2_mean(vv) = mean(vers_pol1_tfp(vv).(var_series)(end - env.grid.oc + env.idx.C_I_Q2_ind,:));
        pol1_tfp_means.C_I_Q3_mean(vv) = mean(vers_pol1_tfp(vv).(var_series)(end - env.grid.oc + env.idx.C_I_Q3_ind,:));
        pol1_tfp_means.C_I_Q4_mean(vv) = mean(vers_pol1_tfp(vv).(var_series)(end - env.grid.oc + env.idx.C_I_Q4_ind,:));
        pol1_tfp_means.C_I_Q5_mean(vv) = mean(vers_pol1_tfp(vv).(var_series)(end - env.grid.oc + env.idx.C_I_Q5_ind,:));
        pol1_tfp_means.C_I_90_mean(vv) = mean(vers_pol1_tfp(vv).(var_series)(end - env.grid.oc + env.idx.C_I_90_ind,:));
        pol1_tfp_means.C_mean(vv)         = mean(vers_pol1_tfp(vv).(var_series)(end - env.grid.oc + env.idx.C_ind,:) );
        pol1_tfp_means.Y_mean(vv)         = mean(vers_pol1_tfp(vv).(var_series)(end - env.grid.oc + env.idx.Y_ind,:) );
        pol1_tfp_means.pi_mean(vv)        = mean(vers_pol1_tfp(vv).(var_series)(end - env.grid.oc + env.idx.pi_ind,:) );
        pol1_tfp_means.unemp_mean(vv)     = mean(vers_pol1_tfp(vv).(var_series)(end - env.grid.oc + env.idx.unemp_ind,:));
        pol1_tfp_means.RR_mean(vv)        = mean(vers_pol1_tfp(vv).(var_series)(end - env.grid.oc + env.idx.RR_ind,:));
        pol1_tfp_means.Profit_mean(vv)    = mean(vers_pol1_tfp(vv).(var_series)(end - env.grid.oc + env.idx.Profit_ind,:));
        pol1_tfp_means.Profit_FI_mean(vv) = mean(vers_pol1_tfp(vv).(var_series)(end - env.grid.oc + env.idx.unemp_ind,:));

        %scalars:
        pol1_tfp_means.Q_mean(vv) = mean(vers_pol1_tfp(vv).(var_series)(env.grid.numstates - env.grid.os + env.idx.Q_ind,:));
        pol1_tfp_means.R_cb_mean(vv) = mean(vers_pol1_tfp(vv).(var_series)(env.grid.numstates - env.grid.os + env.idx.R_cb_ind,:));
        pol1_tfp_means.W_mean(vv) = mean(vers_pol1_tfp(vv).(var_series)(env.grid.numstates - env.grid.os + env.idx.w_ind,:));

        %Standard Deviations
        pol1_tfp_stds.C_10_std(vv) = std(vers_pol1_tfp(vv).(var_series)(end - env.grid.oc + env.idx.C_10_ind,:));
        pol1_tfp_stds.C_Q1_std(vv) = std(vers_pol1_tfp(vv).(var_series)(end - env.grid.oc + env.idx.C_Q1_ind,:));
        pol1_tfp_stds.C_Q2_std(vv) = std(vers_pol1_tfp(vv).(var_series)(end - env.grid.oc + env.idx.C_Q2_ind,:));
        pol1_tfp_stds.C_Q3_std(vv) = std(vers_pol1_tfp(vv).(var_series)(end - env.grid.oc + env.idx.C_Q3_ind,:));
        pol1_tfp_stds.C_Q4_std(vv) = std(vers_pol1_tfp(vv).(var_series)(end - env.grid.oc + env.idx.C_Q4_ind,:));
        pol1_tfp_stds.C_Q5_std(vv) = std(vers_pol1_tfp(vv).(var_series)(end - env.grid.oc + env.idx.C_Q5_ind,:));
        pol1_tfp_stds.C_90_std(vv) = std(vers_pol1_tfp(vv).(var_series)(end - env.grid.oc + env.idx.C_90_ind,:));
        pol1_tfp_stds.C_I_10_std(vv) = std(vers_pol1_tfp(vv).(var_series)(end - env.grid.oc + env.idx.C_I_10_ind,:));
        pol1_tfp_stds.C_I_Q1_std(vv) = std(vers_pol1_tfp(vv).(var_series)(end - env.grid.oc + env.idx.C_I_Q1_ind,:));
        pol1_tfp_stds.C_I_Q2_std(vv) = std(vers_pol1_tfp(vv).(var_series)(end - env.grid.oc + env.idx.C_I_Q2_ind,:));
        pol1_tfp_stds.C_I_Q3_std(vv) = std(vers_pol1_tfp(vv).(var_series)(end - env.grid.oc + env.idx.C_I_Q3_ind,:));
        pol1_tfp_stds.C_I_Q4_std(vv) = std(vers_pol1_tfp(vv).(var_series)(end - env.grid.oc + env.idx.C_I_Q4_ind,:));
        pol1_tfp_stds.C_I_Q5_std(vv) = std(vers_pol1_tfp(vv).(var_series)(end - env.grid.oc + env.idx.C_I_Q5_ind,:));
        pol1_tfp_stds.C_I_90_std(vv) = std(vers_pol1_tfp(vv).(var_series)(end - env.grid.oc + env.idx.C_I_90_ind,:));
        pol1_tfp_stds.C_std(vv)         = std(vers_pol1_tfp(vv).(var_series)(end - env.grid.oc + env.idx.C_ind,:) );
        pol1_tfp_stds.Y_std(vv)         = std(vers_pol1_tfp(vv).(var_series)(end - env.grid.oc + env.idx.Y_ind,:) );
        pol1_tfp_stds.pi_std(vv)        = std(vers_pol1_tfp(vv).(var_series)(end - env.grid.oc + env.idx.pi_ind,:) );
        pol1_tfp_stds.unemp_std(vv)     = std(vers_pol1_tfp(vv).(var_series)(end - env.grid.oc + env.idx.unemp_ind,:));
        pol1_tfp_stds.RR_std(vv)        = std(vers_pol1_tfp(vv).(var_series)(end - env.grid.oc + env.idx.RR_ind,:));
        pol1_tfp_stds.Profit_std(vv)    = std(vers_pol1_tfp(vv).(var_series)(end - env.grid.oc + env.idx.Profit_ind,:));
        pol1_tfp_stds.Profit_FI_std(vv) = std(vers_pol1_tfp(vv).(var_series)(end - env.grid.oc + env.idx.unemp_ind,:));

        %scalars:
        pol1_tfp_stds.Q_std(vv) = std(vers_pol1_tfp(vv).(var_series)(env.grid.numstates - env.grid.os + env.idx.Q_ind,:));
        pol1_tfp_stds.R_cb_std(vv) = std(vers_pol1_tfp(vv).(var_series)(env.grid.numstates - env.grid.os + env.idx.R_cb_ind,:));
        pol1_tfp_stds.W_std(vv) = std(vers_pol1_tfp(vv).(var_series)(env.grid.numstates - env.grid.os + env.idx.w_ind,:));


        % Manual Standard Deviations:
        pol1_tfp_stds2.C_10_std(vv) = manual_std(vers_pol1_tfp(vv).(var_series)(end - env.grid.oc + env.idx.C_10_ind,:));
        pol1_tfp_stds2.C_Q1_std(vv) = manual_std(vers_pol1_tfp(vv).(var_series)(end - env.grid.oc + env.idx.C_Q1_ind,:));
        pol1_tfp_stds2.C_Q2_std(vv) = manual_std(vers_pol1_tfp(vv).(var_series)(end - env.grid.oc + env.idx.C_Q2_ind,:));
        pol1_tfp_stds2.C_Q3_std(vv) = manual_std(vers_pol1_tfp(vv).(var_series)(end - env.grid.oc + env.idx.C_Q3_ind,:));
        pol1_tfp_stds2.C_Q4_std(vv) = manual_std(vers_pol1_tfp(vv).(var_series)(end - env.grid.oc + env.idx.C_Q4_ind,:));
        pol1_tfp_stds2.C_Q5_std(vv) = manual_std(vers_pol1_tfp(vv).(var_series)(end - env.grid.oc + env.idx.C_Q5_ind,:));
        pol1_tfp_stds2.C_90_std(vv) = manual_std(vers_pol1_tfp(vv).(var_series)(end - env.grid.oc + env.idx.C_90_ind,:));
        pol1_tfp_stds2.C_I_10_std(vv) = manual_std(vers_pol1_tfp(vv).(var_series)(end - env.grid.oc + env.idx.C_I_10_ind,:));
        pol1_tfp_stds2.C_I_Q1_std(vv) = manual_std(vers_pol1_tfp(vv).(var_series)(end - env.grid.oc + env.idx.C_I_Q1_ind,:));
        pol1_tfp_stds2.C_I_Q2_std(vv) = manual_std(vers_pol1_tfp(vv).(var_series)(end - env.grid.oc + env.idx.C_I_Q2_ind,:));
        pol1_tfp_stds2.C_I_Q3_std(vv) = manual_std(vers_pol1_tfp(vv).(var_series)(end - env.grid.oc + env.idx.C_I_Q3_ind,:));
        pol1_tfp_stds2.C_I_Q4_std(vv) = manual_std(vers_pol1_tfp(vv).(var_series)(end - env.grid.oc + env.idx.C_I_Q4_ind,:));
        pol1_tfp_stds2.C_I_Q5_std(vv) = manual_std(vers_pol1_tfp(vv).(var_series)(end - env.grid.oc + env.idx.C_I_Q5_ind,:));
        pol1_tfp_stds2.C_I_90_std(vv) = manual_std(vers_pol1_tfp(vv).(var_series)(end - env.grid.oc + env.idx.C_I_90_ind,:));
        pol1_tfp_stds2.C_std(vv)         = manual_std(vers_pol1_tfp(vv).(var_series)(end - env.grid.oc + env.idx.C_ind,:) );
        pol1_tfp_stds2.Y_std(vv)         = manual_std(vers_pol1_tfp(vv).(var_series)(end - env.grid.oc + env.idx.Y_ind,:) );
        pol1_tfp_stds2.pi_std(vv)        = manual_std(vers_pol1_tfp(vv).(var_series)(end - env.grid.oc + env.idx.pi_ind,:) );
        pol1_tfp_stds2.unemp_std(vv)     = manual_std(vers_pol1_tfp(vv).(var_series)(end - env.grid.oc + env.idx.unemp_ind,:));
        pol1_tfp_stds2.RR_std(vv)        = manual_std(vers_pol1_tfp(vv).(var_series)(end - env.grid.oc + env.idx.RR_ind,:));
        pol1_tfp_stds2.Profit_std(vv)    = manual_std(vers_pol1_tfp(vv).(var_series)(end - env.grid.oc + env.idx.Profit_ind,:));
        pol1_tfp_stds2.Profit_FI_std(vv) = manual_std(vers_pol1_tfp(vv).(var_series)(end - env.grid.oc + env.idx.unemp_ind,:));

        %scalars:
        pol1_tfp_stds2.Q_std(vv) = manual_std(vers_pol1_tfp(vv).(var_series)(env.grid.numstates - env.grid.os + env.idx.Q_ind,:));
        pol1_tfp_stds2.R_cb_std(vv) = manual_std(vers_pol1_tfp(vv).(var_series)(env.grid.numstates - env.grid.os + env.idx.R_cb_ind,:));
        pol1_tfp_stds2.W_std(vv) = manual_std(vers_pol1_tfp(vv).(var_series)(env.grid.numstates - env.grid.os + env.idx.w_ind,:));

        pol1_tfp_means.W_10_std(vv) = mean(vers_pol1_tfp(vv).(var_series)(end - env.grid.oc + env.idx.w_10_ind,:));
        pol1_tfp_means.W_Q1_std(vv) = mean(vers_pol1_tfp(vv).(var_series)(end - env.grid.oc + env.idx.w_Q1_ind,:));
        pol1_tfp_means.W_Q2_std(vv) = mean(vers_pol1_tfp(vv).(var_series)(end - env.grid.oc + env.idx.w_Q2_ind,:));
        pol1_tfp_means.W_Q3_std(vv) = mean(vers_pol1_tfp(vv).(var_series)(end - env.grid.oc + env.idx.w_Q3_ind,:));
        pol1_tfp_means.W_Q4_std(vv) = mean(vers_pol1_tfp(vv).(var_series)(end - env.grid.oc + env.idx.w_Q4_ind,:));
        pol1_tfp_means.W_Q5_std(vv) = mean(vers_pol1_tfp(vv).(var_series)(end - env.grid.oc + env.idx.w_Q5_ind,:));
        pol1_tfp_means.W_90_std(vv) = mean(vers_pol1_tfp(vv).(var_series)(end - env.grid.oc + env.idx.w_90_ind,:));
        pol1_tfp_stds.W_10_std(vv) = std(vers_pol1_tfp(vv).(var_series)(end - env.grid.oc + env.idx.w_10_ind,:));
        pol1_tfp_stds.W_Q1_std(vv) = std(vers_pol1_tfp(vv).(var_series)(end - env.grid.oc + env.idx.w_Q1_ind,:));
        pol1_tfp_stds.W_Q2_std(vv) = std(vers_pol1_tfp(vv).(var_series)(end - env.grid.oc + env.idx.w_Q2_ind,:));
        pol1_tfp_stds.W_Q3_std(vv) = std(vers_pol1_tfp(vv).(var_series)(end - env.grid.oc + env.idx.w_Q3_ind,:));
        pol1_tfp_stds.W_Q4_std(vv) = std(vers_pol1_tfp(vv).(var_series)(end - env.grid.oc + env.idx.w_Q4_ind,:));
        pol1_tfp_stds.W_Q5_std(vv) = std(vers_pol1_tfp(vv).(var_series)(end - env.grid.oc + env.idx.w_Q5_ind,:));
        pol1_tfp_stds.W_90_std(vv) = std(vers_pol1_tfp(vv).(var_series)(end - env.grid.oc + env.idx.w_90_ind,:));
        pol1_tfp_stds2.W_10_std(vv) = manual_std(vers_pol1_tfp(vv).(var_series)(end - env.grid.oc + env.idx.w_10_ind,:));
        pol1_tfp_stds2.W_Q1_std(vv) = manual_std(vers_pol1_tfp(vv).(var_series)(end - env.grid.oc + env.idx.w_Q1_ind,:));
        pol1_tfp_stds2.W_Q2_std(vv) = manual_std(vers_pol1_tfp(vv).(var_series)(end - env.grid.oc + env.idx.w_Q2_ind,:));
        pol1_tfp_stds2.W_Q3_std(vv) = manual_std(vers_pol1_tfp(vv).(var_series)(end - env.grid.oc + env.idx.w_Q3_ind,:));
        pol1_tfp_stds2.W_Q4_std(vv) = manual_std(vers_pol1_tfp(vv).(var_series)(end - env.grid.oc + env.idx.w_Q4_ind,:));
        pol1_tfp_stds2.W_Q5_std(vv) = manual_std(vers_pol1_tfp(vv).(var_series)(end - env.grid.oc + env.idx.w_Q5_ind,:));
        pol1_tfp_stds2.W_90_std(vv) = manual_std(vers_pol1_tfp(vv).(var_series)(end - env.grid.oc + env.idx.w_90_ind,:));

        pol1_tfp_means.R_a_mean(vv)     = mean(vers_pol1_tfp(vv).(var_series)(end - env.grid.oc + env.idx.r_a_ind,:));
        pol1_tfp_stds.R_a_std(vv)     = std(vers_pol1_tfp(vv).(var_series)(end - env.grid.oc + env.idx.r_a_ind,:));
        pol1_tfp_stds2.R_a_std(vv)     = manual_std(vers_pol1_tfp(vv).(var_series)(end - env.grid.oc + env.idx.r_a_ind,:));
        pol1_tfp_means.RR_new_mean(vv)        = mean(vers_pol1_tfp(vv).(var_series)(env.grid.numstates - env.grid.os + env.idx.R_cb_ind,:) - vers_pol1_tfp(vv).(var_series)(end - env.grid.oc + env.idx.pi_ind,:) );
        pol1_tfp_stds.RR_new_std(vv)        = std(vers_pol1_tfp(vv).(var_series)(env.grid.numstates - env.grid.os + env.idx.R_cb_ind,:) - vers_pol1_tfp(vv).(var_series)(end - env.grid.oc + env.idx.pi_ind,:) );
        pol1_tfp_stds2.RR_new_std(vv)        = manual_std(vers_pol1_tfp(vv).(var_series)(env.grid.numstates - env.grid.os + env.idx.R_cb_ind,:) - vers_pol1_tfp(vv).(var_series)(end - env.grid.oc + env.idx.pi_ind,:) );

        pol1_mk_means.C_10_mean(vv) = mean(vers_pol1_mk(vv).(var_series)(end - env.grid.oc + env.idx.C_10_ind,:));
        pol1_mk_means.C_Q1_mean(vv) = mean(vers_pol1_mk(vv).(var_series)(end - env.grid.oc + env.idx.C_Q1_ind,:));
        pol1_mk_means.C_Q2_mean(vv) = mean(vers_pol1_mk(vv).(var_series)(end - env.grid.oc + env.idx.C_Q2_ind,:));
        pol1_mk_means.C_Q3_mean(vv) = mean(vers_pol1_mk(vv).(var_series)(end - env.grid.oc + env.idx.C_Q3_ind,:));
        pol1_mk_means.C_Q4_mean(vv) = mean(vers_pol1_mk(vv).(var_series)(end - env.grid.oc + env.idx.C_Q4_ind,:));
        pol1_mk_means.C_Q5_mean(vv) = mean(vers_pol1_mk(vv).(var_series)(end - env.grid.oc + env.idx.C_Q5_ind,:));
        pol1_mk_means.C_90_mean(vv) = mean(vers_pol1_mk(vv).(var_series)(end - env.grid.oc + env.idx.C_90_ind,:));
        pol1_mk_means.C_I_10_mean(vv) = mean(vers_pol1_mk(vv).(var_series)(end - env.grid.oc + env.idx.C_I_10_ind,:));
        pol1_mk_means.C_I_Q1_mean(vv) = mean(vers_pol1_mk(vv).(var_series)(end - env.grid.oc + env.idx.C_I_Q1_ind,:));
        pol1_mk_means.C_I_Q2_mean(vv) = mean(vers_pol1_mk(vv).(var_series)(end - env.grid.oc + env.idx.C_I_Q2_ind,:));
        pol1_mk_means.C_I_Q3_mean(vv) = mean(vers_pol1_mk(vv).(var_series)(end - env.grid.oc + env.idx.C_I_Q3_ind,:));
        pol1_mk_means.C_I_Q4_mean(vv) = mean(vers_pol1_mk(vv).(var_series)(end - env.grid.oc + env.idx.C_I_Q4_ind,:));
        pol1_mk_means.C_I_Q5_mean(vv) = mean(vers_pol1_mk(vv).(var_series)(end - env.grid.oc + env.idx.C_I_Q5_ind,:));
        pol1_mk_means.C_I_90_mean(vv) = mean(vers_pol1_mk(vv).(var_series)(end - env.grid.oc + env.idx.C_I_90_ind,:));
        pol1_mk_means.C_mean(vv)         = mean(vers_pol1_mk(vv).(var_series)(end - env.grid.oc + env.idx.C_ind,:) );
        pol1_mk_means.Y_mean(vv)         = mean(vers_pol1_mk(vv).(var_series)(end - env.grid.oc + env.idx.Y_ind,:) );
        pol1_mk_means.pi_mean(vv)        = mean(vers_pol1_mk(vv).(var_series)(end - env.grid.oc + env.idx.pi_ind,:) );
        pol1_mk_means.unemp_mean(vv)     = mean(vers_pol1_mk(vv).(var_series)(end - env.grid.oc + env.idx.unemp_ind,:));
        pol1_mk_means.RR_mean(vv)        = mean(vers_pol1_mk(vv).(var_series)(end - env.grid.oc + env.idx.RR_ind,:));
        pol1_mk_means.Profit_mean(vv)    = mean(vers_pol1_mk(vv).(var_series)(end - env.grid.oc + env.idx.Profit_ind,:));
        pol1_mk_means.Profit_FI_mean(vv) = mean(vers_pol1_mk(vv).(var_series)(end - env.grid.oc + env.idx.unemp_ind,:));

        %scalars:
        pol1_mk_means.Q_mean(vv) = mean(vers_pol1_mk(vv).(var_series)(env.grid.numstates - env.grid.os + env.idx.Q_ind,:));
        pol1_mk_means.R_cb_mean(vv) = mean(vers_pol1_mk(vv).(var_series)(env.grid.numstates - env.grid.os + env.idx.R_cb_ind,:));
        pol1_mk_means.W_mean(vv) = mean(vers_pol1_mk(vv).(var_series)(env.grid.numstates - env.grid.os + env.idx.w_ind,:));


        %Standard Deviations
        pol1_mk_stds.C_10_std(vv) = std(vers_pol1_mk(vv).(var_series)(end - env.grid.oc + env.idx.C_10_ind,:));
        pol1_mk_stds.C_Q1_std(vv) = std(vers_pol1_mk(vv).(var_series)(end - env.grid.oc + env.idx.C_Q1_ind,:));
        pol1_mk_stds.C_Q2_std(vv) = std(vers_pol1_mk(vv).(var_series)(end - env.grid.oc + env.idx.C_Q2_ind,:));
        pol1_mk_stds.C_Q3_std(vv) = std(vers_pol1_mk(vv).(var_series)(end - env.grid.oc + env.idx.C_Q3_ind,:));
        pol1_mk_stds.C_Q4_std(vv) = std(vers_pol1_mk(vv).(var_series)(end - env.grid.oc + env.idx.C_Q4_ind,:));
        pol1_mk_stds.C_Q5_std(vv) = std(vers_pol1_mk(vv).(var_series)(end - env.grid.oc + env.idx.C_Q5_ind,:));
        pol1_mk_stds.C_90_std(vv) = std(vers_pol1_mk(vv).(var_series)(end - env.grid.oc + env.idx.C_90_ind,:));
        pol1_mk_stds.C_I_10_std(vv) = std(vers_pol1_mk(vv).(var_series)(end - env.grid.oc + env.idx.C_I_10_ind,:));
        pol1_mk_stds.C_I_Q1_std(vv) = std(vers_pol1_mk(vv).(var_series)(end - env.grid.oc + env.idx.C_I_Q1_ind,:));
        pol1_mk_stds.C_I_Q2_std(vv) = std(vers_pol1_mk(vv).(var_series)(end - env.grid.oc + env.idx.C_I_Q2_ind,:));
        pol1_mk_stds.C_I_Q3_std(vv) = std(vers_pol1_mk(vv).(var_series)(end - env.grid.oc + env.idx.C_I_Q3_ind,:));
        pol1_mk_stds.C_I_Q4_std(vv) = std(vers_pol1_mk(vv).(var_series)(end - env.grid.oc + env.idx.C_I_Q4_ind,:));
        pol1_mk_stds.C_I_Q5_std(vv) = std(vers_pol1_mk(vv).(var_series)(end - env.grid.oc + env.idx.C_I_Q5_ind,:));
        pol1_mk_stds.C_I_90_std(vv) = std(vers_pol1_mk(vv).(var_series)(end - env.grid.oc + env.idx.C_I_90_ind,:));
        pol1_mk_stds.C_std(vv)         = std(vers_pol1_mk(vv).(var_series)(end - env.grid.oc + env.idx.C_ind,:) );
        pol1_mk_stds.Y_std(vv)         = std(vers_pol1_mk(vv).(var_series)(end - env.grid.oc + env.idx.Y_ind,:) );
        pol1_mk_stds.pi_std(vv)        = std(vers_pol1_mk(vv).(var_series)(end - env.grid.oc + env.idx.pi_ind,:) );
        pol1_mk_stds.unemp_std(vv)     = std(vers_pol1_mk(vv).(var_series)(end - env.grid.oc + env.idx.unemp_ind,:));
        pol1_mk_stds.RR_std(vv)        = std(vers_pol1_mk(vv).(var_series)(end - env.grid.oc + env.idx.RR_ind,:));
        pol1_mk_stds.Profit_std(vv)    = std(vers_pol1_mk(vv).(var_series)(end - env.grid.oc + env.idx.Profit_ind,:));
        pol1_mk_stds.Profit_FI_std(vv) = std(vers_pol1_mk(vv).(var_series)(end - env.grid.oc + env.idx.unemp_ind,:));

        %scalars:
        pol1_mk_stds.Q_std(vv) = std(vers_pol1_mk(vv).(var_series)(env.grid.numstates - env.grid.os + env.idx.Q_ind,:));
        pol1_mk_stds.R_cb_std(vv) = std(vers_pol1_mk(vv).(var_series)(env.grid.numstates - env.grid.os + env.idx.R_cb_ind,:));
        pol1_mk_stds.W_std(vv) = std(vers_pol1_mk(vv).(var_series)(env.grid.numstates - env.grid.os + env.idx.w_ind,:));


        % Manual Standard Deviations:
        pol1_mk_stds2.C_10_std(vv) = manual_std(vers_pol1_mk(vv).(var_series)(end - env.grid.oc + env.idx.C_10_ind,:));
        pol1_mk_stds2.C_Q1_std(vv) = manual_std(vers_pol1_mk(vv).(var_series)(end - env.grid.oc + env.idx.C_Q1_ind,:));
        pol1_mk_stds2.C_Q2_std(vv) = manual_std(vers_pol1_mk(vv).(var_series)(end - env.grid.oc + env.idx.C_Q2_ind,:));
        pol1_mk_stds2.C_Q3_std(vv) = manual_std(vers_pol1_mk(vv).(var_series)(end - env.grid.oc + env.idx.C_Q3_ind,:));
        pol1_mk_stds2.C_Q4_std(vv) = manual_std(vers_pol1_mk(vv).(var_series)(end - env.grid.oc + env.idx.C_Q4_ind,:));
        pol1_mk_stds2.C_Q5_std(vv) = manual_std(vers_pol1_mk(vv).(var_series)(end - env.grid.oc + env.idx.C_Q5_ind,:));
        pol1_mk_stds2.C_90_std(vv) = manual_std(vers_pol1_mk(vv).(var_series)(end - env.grid.oc + env.idx.C_90_ind,:));
        pol1_mk_stds2.C_I_10_std(vv) = manual_std(vers_pol1_mk(vv).(var_series)(end - env.grid.oc + env.idx.C_I_10_ind,:));
        pol1_mk_stds2.C_I_Q1_std(vv) = manual_std(vers_pol1_mk(vv).(var_series)(end - env.grid.oc + env.idx.C_I_Q1_ind,:));
        pol1_mk_stds2.C_I_Q2_std(vv) = manual_std(vers_pol1_mk(vv).(var_series)(end - env.grid.oc + env.idx.C_I_Q2_ind,:));
        pol1_mk_stds2.C_I_Q3_std(vv) = manual_std(vers_pol1_mk(vv).(var_series)(end - env.grid.oc + env.idx.C_I_Q3_ind,:));
        pol1_mk_stds2.C_I_Q4_std(vv) = manual_std(vers_pol1_mk(vv).(var_series)(end - env.grid.oc + env.idx.C_I_Q4_ind,:));
        pol1_mk_stds2.C_I_Q5_std(vv) = manual_std(vers_pol1_mk(vv).(var_series)(end - env.grid.oc + env.idx.C_I_Q5_ind,:));
        pol1_mk_stds2.C_I_90_std(vv) = manual_std(vers_pol1_mk(vv).(var_series)(end - env.grid.oc + env.idx.C_I_90_ind,:));
        pol1_mk_stds2.C_std(vv)         = manual_std(vers_pol1_mk(vv).(var_series)(end - env.grid.oc + env.idx.C_ind,:) );
        pol1_mk_stds2.Y_std(vv)         = manual_std(vers_pol1_mk(vv).(var_series)(end - env.grid.oc + env.idx.Y_ind,:) );
        pol1_mk_stds2.pi_std(vv)        = manual_std(vers_pol1_mk(vv).(var_series)(end - env.grid.oc + env.idx.pi_ind,:) );
        pol1_mk_stds2.unemp_std(vv)     = manual_std(vers_pol1_mk(vv).(var_series)(end - env.grid.oc + env.idx.unemp_ind,:));
        pol1_mk_stds2.RR_std(vv)        = manual_std(vers_pol1_mk(vv).(var_series)(end - env.grid.oc + env.idx.RR_ind,:));
        pol1_mk_stds2.Profit_std(vv)    = manual_std(vers_pol1_mk(vv).(var_series)(end - env.grid.oc + env.idx.Profit_ind,:));
        pol1_mk_stds2.Profit_FI_std(vv) = manual_std(vers_pol1_mk(vv).(var_series)(end - env.grid.oc + env.idx.unemp_ind,:));

        %scalars:
        pol1_mk_stds2.Q_std(vv) = manual_std(vers_pol1_mk(vv).(var_series)(env.grid.numstates - env.grid.os + env.idx.Q_ind,:));
        pol1_mk_stds2.R_cb_std(vv) = manual_std(vers_pol1_mk(vv).(var_series)(env.grid.numstates - env.grid.os + env.idx.R_cb_ind,:));
        pol1_mk_stds2.W_std(vv) = manual_std(vers_pol1_mk(vv).(var_series)(env.grid.numstates - env.grid.os + env.idx.w_ind,:));

        pol1_mk_means.W_10_std(vv) = mean(vers_pol1_mk(vv).(var_series)(end - env.grid.oc + env.idx.w_10_ind,:));
        pol1_mk_means.W_Q1_std(vv) = mean(vers_pol1_mk(vv).(var_series)(end - env.grid.oc + env.idx.w_Q1_ind,:));
        pol1_mk_means.W_Q2_std(vv) = mean(vers_pol1_mk(vv).(var_series)(end - env.grid.oc + env.idx.w_Q2_ind,:));
        pol1_mk_means.W_Q3_std(vv) = mean(vers_pol1_mk(vv).(var_series)(end - env.grid.oc + env.idx.w_Q3_ind,:));
        pol1_mk_means.W_Q4_std(vv) = mean(vers_pol1_mk(vv).(var_series)(end - env.grid.oc + env.idx.w_Q4_ind,:));
        pol1_mk_means.W_Q5_std(vv) = mean(vers_pol1_mk(vv).(var_series)(end - env.grid.oc + env.idx.w_Q5_ind,:));
        pol1_mk_means.W_90_std(vv) = mean(vers_pol1_mk(vv).(var_series)(end - env.grid.oc + env.idx.w_90_ind,:));
        pol1_mk_stds.W_10_std(vv) = std(vers_pol1_mk(vv).(var_series)(end - env.grid.oc + env.idx.w_10_ind,:));
        pol1_mk_stds.W_Q1_std(vv) = std(vers_pol1_mk(vv).(var_series)(end - env.grid.oc + env.idx.w_Q1_ind,:));
        pol1_mk_stds.W_Q2_std(vv) = std(vers_pol1_mk(vv).(var_series)(end - env.grid.oc + env.idx.w_Q2_ind,:));
        pol1_mk_stds.W_Q3_std(vv) = std(vers_pol1_mk(vv).(var_series)(end - env.grid.oc + env.idx.w_Q3_ind,:));
        pol1_mk_stds.W_Q4_std(vv) = std(vers_pol1_mk(vv).(var_series)(end - env.grid.oc + env.idx.w_Q4_ind,:));
        pol1_mk_stds.W_Q5_std(vv) = std(vers_pol1_mk(vv).(var_series)(end - env.grid.oc + env.idx.w_Q5_ind,:));
        pol1_mk_stds.W_90_std(vv) = std(vers_pol1_mk(vv).(var_series)(end - env.grid.oc + env.idx.w_90_ind,:));
        pol1_mk_stds2.W_10_std(vv) = manual_std(vers_pol1_mk(vv).(var_series)(end - env.grid.oc + env.idx.w_10_ind,:));
        pol1_mk_stds2.W_Q1_std(vv) = manual_std(vers_pol1_mk(vv).(var_series)(end - env.grid.oc + env.idx.w_Q1_ind,:));
        pol1_mk_stds2.W_Q2_std(vv) = manual_std(vers_pol1_mk(vv).(var_series)(end - env.grid.oc + env.idx.w_Q2_ind,:));
        pol1_mk_stds2.W_Q3_std(vv) = manual_std(vers_pol1_mk(vv).(var_series)(end - env.grid.oc + env.idx.w_Q3_ind,:));
        pol1_mk_stds2.W_Q4_std(vv) = manual_std(vers_pol1_mk(vv).(var_series)(end - env.grid.oc + env.idx.w_Q4_ind,:));
        pol1_mk_stds2.W_Q5_std(vv) = manual_std(vers_pol1_mk(vv).(var_series)(end - env.grid.oc + env.idx.w_Q5_ind,:));
        pol1_mk_stds2.W_90_std(vv) = manual_std(vers_pol1_mk(vv).(var_series)(end - env.grid.oc + env.idx.w_90_ind,:));

        pol1_mk_means.R_a_mean(vv)    = mean(vers_pol1_mk(vv).(var_series)(end - env.grid.oc + env.idx.r_a_ind,:));
        pol1_mk_stds.R_a_std(vv)      = std(vers_pol1_mk(vv).(var_series)(end - env.grid.oc + env.idx.r_a_ind,:));
        pol1_mk_stds2.R_a_std(vv)     = manual_std(vers_pol1_mk(vv).(var_series)(end - env.grid.oc + env.idx.r_a_ind,:));
        pol1_mk_means.RR_new_mean(vv) = mean(vers_pol1_mk(vv).(var_series)(env.grid.numstates - env.grid.os + env.idx.R_cb_ind,:) - vers_pol1_mk(vv).(var_series)(end - env.grid.oc + env.idx.pi_ind,:) );
        pol1_mk_stds.RR_new_std(vv)   = std(vers_pol1_mk(vv).(var_series)(env.grid.numstates - env.grid.os + env.idx.R_cb_ind,:) - vers_pol1_mk(vv).(var_series)(end - env.grid.oc + env.idx.pi_ind,:) );
        pol1_mk_stds2.RR_new_std(vv)  = manual_std(vers_pol1_mk(vv).(var_series)(env.grid.numstates - env.grid.os + env.idx.R_cb_ind,:) - vers_pol1_mk(vv).(var_series)(end - env.grid.oc + env.idx.pi_ind,:) );

    end

    %%

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % (4) Plotting Function Calls: Disaggregated Plots
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % (4.1) Create  plots for disaggregated policy-rule plots
    %%%%%%%%%%

    disagg_plots(pol1_stds, pol1_mk_stds, pol1_tfp_stds, phipi1_vals, string(fieldnames(pol1_stds)), "SS_Dev", pltroot)


%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% (5) Plotting functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% (5.1)
% function disagg_plots()
%
% Arguments:
% (1) all_std: Struct of volatilities for each variable by phi_pi value
% (under simulation where both markup and tfp shocks are active)
% (2) markup_std: Struct of volatilities for each variable by phi_pi value
% (under simulation with isolated markup shocks)
% (3) tfp_std: Struct of volatilities for each variable by phi_pi value
% (under simulation with isolated tfp shocks)
% (4) phi_pi_s: Iterables for phi_pi values under specified policy rule
% (5) type: String indentifier for plots (MATLAB std() = "SS_Dev, manual =
% "SS_Dev_manual)
% (6) pltroot: File path to save plots
%
% Description: Produces phi-pi disaggregation plots
%%%%%%%%%%

function disagg_plots(all_std, markup_std, tfp_std, phi_pis, fields, type, pltroot)
for f = 1:length(fields) %% Assuming that the two have same fields

    char_f = char(fields(f));
    para_name = string(char_f(1:end-4)); %removes '_std' from string name
    x_lab = "\phi_{\pi}";

    figure('name', strcat(para_name,"_phi_pi"));
    hold on
    plot(phi_pis, 100 * getfield(all_std, fields(f)), "--o", 'Color', "black", "LineWidth", 1.5)
    plot(phi_pis, 100 * getfield(markup_std, fields(f)), "--o", 'Color', "red", "LineWidth", 1.5)
    plot(phi_pis, 100 * getfield(tfp_std, fields(f)), "--o", 'Color', "blue", "LineWidth", 1.5)
    hold off
    title('Deviation of ' + strrep(para_name, "_", "-") + " v.s " + x_lab);
    xlabel(x_lab);
    ylabel('Deviation of ' + strrep(para_name, "_", "-") + ' (in %)');

    l = legend('Supply Shocks', 'Markup Shocks', 'Productivity Shocks');
    legend('boxoff')

    fi = gcf;
    base = pltroot.disagg;

    plt_type = type;
    adl_fld = "Taylor/";
    base_2 = base + adl_fld;

    if ~exist(base_2, 'dir')
        mkdir(base_2)
    end
    exportgraphics(fi, base_2 + plt_type+ '_disagg_' + para_name + '.png', 'Resolution', 300);

    close all;

end


end

