% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % simul_plot.m
% %
% % Description: Given simulation results, file recreates
% % (1) aggregate and distributional variable volatility vs phi_pi values
% % values ('standard deviation plots' in paper|'phi-pi' plots in our code)
% % (2) disaggregation of aggregate and distributional volatility vs phi_pi
% % parameter by supply and demand shocks
% % (3) aggregate and distributional variable volatility vs inflation
% % volatility ('frontier plots' in the paper)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% (1.2) Set up output file paths
%%%%%%%%%%

% Output file folders: Root folder/sub-root folder/[changable final path]
plt_subfolder = strcat('ELB=',simul_ELB);
% For simulations run with baseline kappa parameter value
if strcmp(simul_str, 'shock_series_rand_new')
    pltroot.onefrontierplot = plt_folder + "/Simulations/" + plt_subfolder + "/nocovidGR/TaylorFrontierPlots/";
else
    pltroot.kappafrontierplots = plt_folder + "/Simulations/" + plt_subfolder + "/FrontierPlots/HighKappa/";
    pltroot.SDfrontierplots = plt_folder + "/Simulations/" + plt_subfolder + "/FrontierPlots/SupplyDemand/";
    pltroot.supply_demand_plots = plt_folder + "/Simulations/" + plt_subfolder + "/Supply_Demand_Disagg/";
    pltroot.ZLBplots = plt_folder + "/Simulations/" + plt_subfolder + "/";
    pltroot.onefrontierplot = plt_folder + "/Simulations/" + plt_subfolder + "/TaylorFrontierPlots/";
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% (1.3) Volatility calculation
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
    filepath1 = results_folder + "/precovid_ugap_v4_pol=" + pol1 + "_shocks=" + simul_str +"_ELB=" + simul_ELB+ "_all_phipi=" + pp + "_kappa=" + kappa_val1 + ".mat";
    if ~strcmp(simul_str, 'shock_series_rand_new')
        filepath2 =  results_folder + "/precovid_ugap_v4_pol=" + pol1 + "_shocks=" + simul_str +"_ELB=" + simul_ELB+ "_all_phipi=" + pp + "_kappa=" + kappa_val2 + ".mat";
        Pol2_results(vv) = load(filepath2).results;

    end
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


    %%%%%%%%% Repeat for high kappa:
    if ~strcmp(simul_str, 'shock_series_rand_new')

        pol2_ELB_freq(vv) = sum(exp(Pol2_results(vv).(var_series)(env.grid.numstates - ...
            env.grid.os + ...
            env.idx.R_cb_ind,:))*1.01 ...
            <=1+1e-8) / length(Pol2_results(vv).(var_series)(env.grid.numstates - ...
            env.grid.os + ...
            env.idx.R_cb_ind,:));
        % Calculate variable means
        pol2_means.C_10_mean(vv) = mean(Pol2_results(vv).(var_series)(end - env.grid.oc + env.idx.C_10_ind,:));
        pol2_means.C_Q1_mean(vv) = mean(Pol2_results(vv).(var_series)(end - env.grid.oc + env.idx.C_Q1_ind,:));
        pol2_means.C_Q2_mean(vv) = mean(Pol2_results(vv).(var_series)(end - env.grid.oc + env.idx.C_Q2_ind,:));
        pol2_means.C_Q3_mean(vv) = mean(Pol2_results(vv).(var_series)(end - env.grid.oc + env.idx.C_Q3_ind,:));
        pol2_means.C_Q4_mean(vv) = mean(Pol2_results(vv).(var_series)(end - env.grid.oc + env.idx.C_Q4_ind,:));
        pol2_means.C_Q5_mean(vv) = mean(Pol2_results(vv).(var_series)(end - env.grid.oc + env.idx.C_Q5_ind,:));
        pol2_means.C_90_mean(vv) = mean(Pol2_results(vv).(var_series)(end - env.grid.oc + env.idx.C_90_ind,:));

        pol2_means.W_10_mean(vv) = mean(Pol2_results(vv).(var_series)(end - env.grid.oc + env.idx.w_10_ind,:));
        pol2_means.W_Q1_mean(vv) = mean(Pol2_results(vv).(var_series)(end - env.grid.oc + env.idx.w_Q1_ind,:));
        pol2_means.W_Q2_mean(vv) = mean(Pol2_results(vv).(var_series)(end - env.grid.oc + env.idx.w_Q2_ind,:));
        pol2_means.W_Q3_mean(vv) = mean(Pol2_results(vv).(var_series)(end - env.grid.oc + env.idx.w_Q3_ind,:));
        pol2_means.W_Q4_mean(vv) = mean(Pol2_results(vv).(var_series)(end - env.grid.oc + env.idx.w_Q4_ind,:));
        pol2_means.W_Q5_mean(vv) = mean(Pol2_results(vv).(var_series)(end - env.grid.oc + env.idx.w_Q5_ind,:));
        pol2_means.W_90_mean(vv) = mean(Pol2_results(vv).(var_series)(end - env.grid.oc + env.idx.w_90_ind,:));

        pol2_means.C_I_10_mean(vv) = mean(Pol2_results(vv).(var_series)(end - env.grid.oc + env.idx.C_I_10_ind,:));
        pol2_means.C_I_Q1_mean(vv) = mean(Pol2_results(vv).(var_series)(end - env.grid.oc + env.idx.C_I_Q1_ind,:));
        pol2_means.C_I_Q2_mean(vv) = mean(Pol2_results(vv).(var_series)(end - env.grid.oc + env.idx.C_I_Q2_ind,:));
        pol2_means.C_I_Q3_mean(vv) = mean(Pol2_results(vv).(var_series)(end - env.grid.oc + env.idx.C_I_Q3_ind,:));
        pol2_means.C_I_Q4_mean(vv) = mean(Pol2_results(vv).(var_series)(end - env.grid.oc + env.idx.C_I_Q4_ind,:));
        pol2_means.C_I_Q5_mean(vv) = mean(Pol2_results(vv).(var_series)(end - env.grid.oc + env.idx.C_I_Q5_ind,:));
        pol2_means.C_I_90_mean(vv) = mean(Pol2_results(vv).(var_series)(end - env.grid.oc + env.idx.C_I_90_ind,:));

        pol2_means.C_mean(vv)         = mean(Pol2_results(vv).(var_series)(end - env.grid.oc + env.idx.C_ind,:) );
        pol2_means.Y_mean(vv)         = mean(Pol2_results(vv).(var_series)(end - env.grid.oc + env.idx.Y_ind,:) );
        pol2_means.pi_mean(vv)        = mean(Pol2_results(vv).(var_series)(end - env.grid.oc + env.idx.pi_ind,:) );
        pol2_means.unemp_mean(vv)     = mean(Pol2_results(vv).(var_series)(end - env.grid.oc + env.idx.unemp_ind,:));
        pol2_means.RR_mean(vv)        = mean(Pol2_results(vv).(var_series)(end - env.grid.oc + env.idx.RR_ind,:));
        pol2_means.Profit_mean(vv)    = mean(Pol2_results(vv).(var_series)(end - env.grid.oc + env.idx.Profit_ind,:));
        pol2_means.Profit_FI_mean(vv) = mean(Pol2_results(vv).(var_series)(end - env.grid.oc + env.idx.unemp_ind,:));

        %scalars:
        pol2_means.Q_mean(vv) = mean(Pol2_results(vv).(var_series)(env.grid.numstates - env.grid.os + env.idx.Q_ind,:));
        pol2_means.R_cb_mean(vv) = mean(Pol2_results(vv).(var_series)(env.grid.numstates - env.grid.os + env.idx.R_cb_ind,:));
        pol2_means.W_mean(vv) = mean(Pol2_results(vv).(var_series)(env.grid.numstates - env.grid.os + env.idx.w_ind,:));


        %Standard Deviations
        pol2_stds.C_10_std(vv) = std(Pol2_results(vv).(var_series)(end - env.grid.oc + env.idx.C_10_ind,:));
        pol2_stds.C_Q1_std(vv) = std(Pol2_results(vv).(var_series)(end - env.grid.oc + env.idx.C_Q1_ind,:));
        pol2_stds.C_Q2_std(vv) = std(Pol2_results(vv).(var_series)(end - env.grid.oc + env.idx.C_Q2_ind,:));
        pol2_stds.C_Q3_std(vv) = std(Pol2_results(vv).(var_series)(end - env.grid.oc + env.idx.C_Q3_ind,:));
        pol2_stds.C_Q4_std(vv) = std(Pol2_results(vv).(var_series)(end - env.grid.oc + env.idx.C_Q4_ind,:));
        pol2_stds.C_Q5_std(vv) = std(Pol2_results(vv).(var_series)(end - env.grid.oc + env.idx.C_Q5_ind,:));
        pol2_stds.C_90_std(vv) = std(Pol2_results(vv).(var_series)(end - env.grid.oc + env.idx.C_90_ind,:));

        pol2_stds.W_10_std(vv) = std(Pol2_results(vv).(var_series)(end - env.grid.oc + env.idx.w_10_ind,:));
        pol2_stds.W_Q1_std(vv) = std(Pol2_results(vv).(var_series)(end - env.grid.oc + env.idx.w_Q1_ind,:));
        pol2_stds.W_Q2_std(vv) = std(Pol2_results(vv).(var_series)(end - env.grid.oc + env.idx.w_Q2_ind,:));
        pol2_stds.W_Q3_std(vv) = std(Pol2_results(vv).(var_series)(end - env.grid.oc + env.idx.w_Q3_ind,:));
        pol2_stds.W_Q4_std(vv) = std(Pol2_results(vv).(var_series)(end - env.grid.oc + env.idx.w_Q4_ind,:));
        pol2_stds.W_Q5_std(vv) = std(Pol2_results(vv).(var_series)(end - env.grid.oc + env.idx.w_Q5_ind,:));
        pol2_stds.W_90_std(vv) = std(Pol2_results(vv).(var_series)(end - env.grid.oc + env.idx.w_90_ind,:));

        pol2_stds.C_I_10_std(vv) = std(Pol2_results(vv).(var_series)(end - env.grid.oc + env.idx.C_I_10_ind,:));
        pol2_stds.C_I_Q1_std(vv) = std(Pol2_results(vv).(var_series)(end - env.grid.oc + env.idx.C_I_Q1_ind,:));
        pol2_stds.C_I_Q2_std(vv) = std(Pol2_results(vv).(var_series)(end - env.grid.oc + env.idx.C_I_Q2_ind,:));
        pol2_stds.C_I_Q3_std(vv) = std(Pol2_results(vv).(var_series)(end - env.grid.oc + env.idx.C_I_Q3_ind,:));
        pol2_stds.C_I_Q4_std(vv) = std(Pol2_results(vv).(var_series)(end - env.grid.oc + env.idx.C_I_Q4_ind,:));
        pol2_stds.C_I_Q5_std(vv) = std(Pol2_results(vv).(var_series)(end - env.grid.oc + env.idx.C_I_Q5_ind,:));
        pol2_stds.C_I_90_std(vv) = std(Pol2_results(vv).(var_series)(end - env.grid.oc + env.idx.C_I_90_ind,:));

        pol2_stds.C_std(vv)         = std(Pol2_results(vv).(var_series)(end - env.grid.oc + env.idx.C_ind,:) );
        pol2_stds.Y_std(vv)         = std(Pol2_results(vv).(var_series)(end - env.grid.oc + env.idx.Y_ind,:) );
        pol2_stds.pi_std(vv)        = std(Pol2_results(vv).(var_series)(end - env.grid.oc + env.idx.pi_ind,:) );
        pol2_stds.unemp_std(vv)     = std(Pol2_results(vv).(var_series)(end - env.grid.oc + env.idx.unemp_ind,:));
        pol2_stds.RR_std(vv)        = std(Pol2_results(vv).(var_series)(end - env.grid.oc + env.idx.RR_ind,:));
        pol2_stds.Profit_std(vv)    = std(Pol2_results(vv).(var_series)(end - env.grid.oc + env.idx.Profit_ind,:));
        pol2_stds.Profit_FI_std(vv) = std(Pol2_results(vv).(var_series)(end - env.grid.oc + env.idx.unemp_ind,:));

        %scalars:
        pol2_stds.Q_std(vv) = std(Pol2_results(vv).(var_series)(env.grid.numstates - env.grid.os + env.idx.Q_ind,:));
        pol2_stds.R_cb_std(vv) = std(Pol2_results(vv).(var_series)(env.grid.numstates - env.grid.os + env.idx.R_cb_ind,:));
        pol2_stds.W_std(vv) = std(Pol2_results(vv).(var_series)(env.grid.numstates - env.grid.os + env.idx.w_ind,:));


        % Manual Standard Deviations
        pol2_stds2.C_10_std(vv) = manual_std(Pol2_results(vv).(var_series)(end - env.grid.oc + env.idx.C_10_ind,:));
        pol2_stds2.C_Q1_std(vv) = manual_std(Pol2_results(vv).(var_series)(end - env.grid.oc + env.idx.C_Q1_ind,:));
        pol2_stds2.C_Q2_std(vv) = manual_std(Pol2_results(vv).(var_series)(end - env.grid.oc + env.idx.C_Q2_ind,:));
        pol2_stds2.C_Q3_std(vv) = manual_std(Pol2_results(vv).(var_series)(end - env.grid.oc + env.idx.C_Q3_ind,:));
        pol2_stds2.C_Q4_std(vv) = manual_std(Pol2_results(vv).(var_series)(end - env.grid.oc + env.idx.C_Q4_ind,:));
        pol2_stds2.C_Q5_std(vv) = manual_std(Pol2_results(vv).(var_series)(end - env.grid.oc + env.idx.C_Q5_ind,:));
        pol2_stds2.C_90_std(vv) = manual_std(Pol2_results(vv).(var_series)(end - env.grid.oc + env.idx.C_90_ind,:));

        pol2_stds2.W_10_std(vv) = manual_std(Pol2_results(vv).(var_series)(end - env.grid.oc + env.idx.w_10_ind,:));
        pol2_stds2.W_Q1_std(vv) = manual_std(Pol2_results(vv).(var_series)(end - env.grid.oc + env.idx.w_Q1_ind,:));
        pol2_stds2.W_Q2_std(vv) = manual_std(Pol2_results(vv).(var_series)(end - env.grid.oc + env.idx.w_Q2_ind,:));
        pol2_stds2.W_Q3_std(vv) = manual_std(Pol2_results(vv).(var_series)(end - env.grid.oc + env.idx.w_Q3_ind,:));
        pol2_stds2.W_Q4_std(vv) = manual_std(Pol2_results(vv).(var_series)(end - env.grid.oc + env.idx.w_Q4_ind,:));
        pol2_stds2.W_Q5_std(vv) = manual_std(Pol2_results(vv).(var_series)(end - env.grid.oc + env.idx.w_Q5_ind,:));
        pol2_stds2.W_90_std(vv) = manual_std(Pol2_results(vv).(var_series)(end - env.grid.oc + env.idx.w_90_ind,:));

        pol2_stds2.C_I_10_std(vv) = manual_std(Pol2_results(vv).(var_series)(end - env.grid.oc + env.idx.C_I_10_ind,:));
        pol2_stds2.C_I_Q1_std(vv) = manual_std(Pol2_results(vv).(var_series)(end - env.grid.oc + env.idx.C_I_Q1_ind,:));
        pol2_stds2.C_I_Q2_std(vv) = manual_std(Pol2_results(vv).(var_series)(end - env.grid.oc + env.idx.C_I_Q2_ind,:));
        pol2_stds2.C_I_Q3_std(vv) = manual_std(Pol2_results(vv).(var_series)(end - env.grid.oc + env.idx.C_I_Q3_ind,:));
        pol2_stds2.C_I_Q4_std(vv) = manual_std(Pol2_results(vv).(var_series)(end - env.grid.oc + env.idx.C_I_Q4_ind,:));
        pol2_stds2.C_I_Q5_std(vv) = manual_std(Pol2_results(vv).(var_series)(end - env.grid.oc + env.idx.C_I_Q5_ind,:));
        pol2_stds2.C_I_90_std(vv) = manual_std(Pol2_results(vv).(var_series)(end - env.grid.oc + env.idx.C_I_90_ind,:));

        pol2_stds2.C_std(vv)         = manual_std(Pol2_results(vv).(var_series)(end - env.grid.oc + env.idx.C_ind,:) );
        pol2_stds2.Y_std(vv)         = manual_std(Pol2_results(vv).(var_series)(end - env.grid.oc + env.idx.Y_ind,:) );
        pol2_stds2.pi_std(vv)        = manual_std(Pol2_results(vv).(var_series)(end - env.grid.oc + env.idx.pi_ind,:) );
        pol2_stds2.unemp_std(vv)     = manual_std(Pol2_results(vv).(var_series)(end - env.grid.oc + env.idx.unemp_ind,:));
        pol2_stds2.RR_std(vv)        = manual_std(Pol2_results(vv).(var_series)(end - env.grid.oc + env.idx.RR_ind,:));
        pol2_stds2.Profit_std(vv)    = manual_std(Pol2_results(vv).(var_series)(end - env.grid.oc + env.idx.Profit_ind,:));
        pol2_stds2.Profit_FI_std(vv) = manual_std(Pol2_results(vv).(var_series)(end - env.grid.oc + env.idx.unemp_ind,:));

        %scalars:
        pol2_stds2.Q_std(vv) = manual_std(Pol2_results(vv).(var_series)(env.grid.numstates - env.grid.os + env.idx.Q_ind,:));
        pol2_stds2.R_cb_std(vv) = manual_std(Pol2_results(vv).(var_series)(env.grid.numstates - env.grid.os + env.idx.R_cb_ind,:));
        pol2_stds2.W_std(vv) = manual_std(Pol2_results(vv).(var_series)(env.grid.numstates - env.grid.os + env.idx.w_ind,:));

        pol2_means.R_a_mean(vv) = mean(Pol2_results(vv).(var_series)(end - env.grid.oc + env.idx.r_a_ind,:));
        pol2_stds.R_a_std(vv) = std(Pol2_results(vv).(var_series)(end - env.grid.oc + env.idx.r_a_ind,:));
        pol2_stds2.R_a_std(vv) = manual_std(Pol2_results(vv).(var_series)(end - env.grid.oc + env.idx.r_a_ind,:));

        pol2_means.RR_new_mean(vv) = mean(Pol2_results(vv).(var_series)(env.grid.numstates - env.grid.os + env.idx.R_cb_ind,:) - Pol2_results(vv).(var_series)(end - env.grid.oc + env.idx.pi_ind,:) );
        pol2_stds.RR_new_std(vv) = std(Pol2_results(vv).(var_series)(env.grid.numstates - env.grid.os + env.idx.R_cb_ind,:) - Pol2_results(vv).(var_series)(end - env.grid.oc + env.idx.pi_ind,:) );
        pol2_stds2.RR_new_std(vv) = manual_std(Pol2_results(vv).(var_series)(env.grid.numstates - env.grid.os + env.idx.R_cb_ind,:) - Pol2_results(vv).(var_series)(end - env.grid.oc + env.idx.pi_ind,:) );
    end
end

%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% (3) Plotting Function Calls: Non-Disaggregated Plots
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% (3.1) Create frontier plots for only Taylor rule results
%%%%%%%%%%

one_frontier_plot(pol1_stds, phipi1_vals, string(fieldnames(pol1_stds)), "SS_Dev", pltroot);
close all;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% (3.2) Create phi-pi and frontier volatility plots for both policy rules (using
% MATLAB std() function)
%%%%%%%%%%
if ~strcmp(simul_str, 'shock_series_rand_new')
    frontier_plots(pol1_stds, pol2_stds, pol1, plt_labels ,phipi1_vals, baseline_phi_ind, string(fieldnames(pol1_stds)), "SS_Dev", pltroot, 'kappa');
    close all;
end


%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% (4) Loading simulation results: Disaggregated
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


if plot_supply_demand
    for vv=1:length(phipi1)
        pp=phipi1(vv);
        filepath_pol1_d = results_folder + "/precovid_ugap_v4_pol=" + pol1 + "_shocks=" + simul_str +"_ELB=" + simul_ELB+ "_demand_phipi=" + pp + "_kappa=" + kappa_val1 + ".mat";
        filepath_pol1_s = results_folder + "/precovid_ugap_v4_pol=" + pol1 + "_shocks=" + simul_str +"_ELB=" + simul_ELB+ "_supply_phipi=" + pp + "_kappa=" + kappa_val1 + ".mat";

        vers_pol1_d(vv) = load(filepath_pol1_d).results;
        vers_pol1_s(vv) = load(filepath_pol1_s).results;

        %%% Policy 1:

        pol1_d_means.C_10_mean(vv) = mean(vers_pol1_d(vv).(var_series)(end - env.grid.oc + env.idx.C_10_ind,:));
        pol1_d_means.C_Q1_mean(vv) = mean(vers_pol1_d(vv).(var_series)(end - env.grid.oc + env.idx.C_Q1_ind,:));
        pol1_d_means.C_Q2_mean(vv) = mean(vers_pol1_d(vv).(var_series)(end - env.grid.oc + env.idx.C_Q2_ind,:));
        pol1_d_means.C_Q3_mean(vv) = mean(vers_pol1_d(vv).(var_series)(end - env.grid.oc + env.idx.C_Q3_ind,:));
        pol1_d_means.C_Q4_mean(vv) = mean(vers_pol1_d(vv).(var_series)(end - env.grid.oc + env.idx.C_Q4_ind,:));
        pol1_d_means.C_Q5_mean(vv) = mean(vers_pol1_d(vv).(var_series)(end - env.grid.oc + env.idx.C_Q5_ind,:));
        pol1_d_means.C_90_mean(vv) = mean(vers_pol1_d(vv).(var_series)(end - env.grid.oc + env.idx.C_90_ind,:));

        pol1_d_means.C_I_10_mean(vv) = mean(vers_pol1_d(vv).(var_series)(end - env.grid.oc + env.idx.C_I_10_ind,:));
        pol1_d_means.C_I_Q1_mean(vv) = mean(vers_pol1_d(vv).(var_series)(end - env.grid.oc + env.idx.C_I_Q1_ind,:));
        pol1_d_means.C_I_Q2_mean(vv) = mean(vers_pol1_d(vv).(var_series)(end - env.grid.oc + env.idx.C_I_Q2_ind,:));
        pol1_d_means.C_I_Q3_mean(vv) = mean(vers_pol1_d(vv).(var_series)(end - env.grid.oc + env.idx.C_I_Q3_ind,:));
        pol1_d_means.C_I_Q4_mean(vv) = mean(vers_pol1_d(vv).(var_series)(end - env.grid.oc + env.idx.C_I_Q4_ind,:));
        pol1_d_means.C_I_Q5_mean(vv) = mean(vers_pol1_d(vv).(var_series)(end - env.grid.oc + env.idx.C_I_Q5_ind,:));
        pol1_d_means.C_I_90_mean(vv) = mean(vers_pol1_d(vv).(var_series)(end - env.grid.oc + env.idx.C_I_90_ind,:));
        pol1_d_means.C_mean(vv)         = mean(vers_pol1_d(vv).(var_series)(end - env.grid.oc + env.idx.C_ind,:) );
        pol1_d_means.Y_mean(vv)         = mean(vers_pol1_d(vv).(var_series)(end - env.grid.oc + env.idx.Y_ind,:) );
        pol1_d_means.pi_mean(vv)        = mean(vers_pol1_d(vv).(var_series)(end - env.grid.oc + env.idx.pi_ind,:) );
        pol1_d_means.unemp_mean(vv)     = mean(vers_pol1_d(vv).(var_series)(end - env.grid.oc + env.idx.unemp_ind,:));
        pol1_d_means.RR_mean(vv)        = mean(vers_pol1_d(vv).(var_series)(end - env.grid.oc + env.idx.RR_ind,:));
        pol1_d_means.Profit_mean(vv)    = mean(vers_pol1_d(vv).(var_series)(end - env.grid.oc + env.idx.Profit_ind,:));
        pol1_d_means.Profit_FI_mean(vv) = mean(vers_pol1_d(vv).(var_series)(end - env.grid.oc + env.idx.unemp_ind,:));

        %scalars:
        pol1_d_means.Q_mean(vv) = mean(vers_pol1_d(vv).(var_series)(env.grid.numstates - env.grid.os + env.idx.Q_ind,:));
        pol1_d_means.R_cb_mean(vv) = mean(vers_pol1_d(vv).(var_series)(env.grid.numstates - env.grid.os + env.idx.R_cb_ind,:));
        pol1_d_means.W_mean(vv) = mean(vers_pol1_d(vv).(var_series)(env.grid.numstates - env.grid.os + env.idx.w_ind,:));

        %Standard Deviations
        pol1_d_stds.C_10_std(vv) = std(vers_pol1_d(vv).(var_series)(end - env.grid.oc + env.idx.C_10_ind,:));
        pol1_d_stds.C_Q1_std(vv) = std(vers_pol1_d(vv).(var_series)(end - env.grid.oc + env.idx.C_Q1_ind,:));
        pol1_d_stds.C_Q2_std(vv) = std(vers_pol1_d(vv).(var_series)(end - env.grid.oc + env.idx.C_Q2_ind,:));
        pol1_d_stds.C_Q3_std(vv) = std(vers_pol1_d(vv).(var_series)(end - env.grid.oc + env.idx.C_Q3_ind,:));
        pol1_d_stds.C_Q4_std(vv) = std(vers_pol1_d(vv).(var_series)(end - env.grid.oc + env.idx.C_Q4_ind,:));
        pol1_d_stds.C_Q5_std(vv) = std(vers_pol1_d(vv).(var_series)(end - env.grid.oc + env.idx.C_Q5_ind,:));
        pol1_d_stds.C_90_std(vv) = std(vers_pol1_d(vv).(var_series)(end - env.grid.oc + env.idx.C_90_ind,:));
        pol1_d_stds.C_I_10_std(vv) = std(vers_pol1_d(vv).(var_series)(end - env.grid.oc + env.idx.C_I_10_ind,:));
        pol1_d_stds.C_I_Q1_std(vv) = std(vers_pol1_d(vv).(var_series)(end - env.grid.oc + env.idx.C_I_Q1_ind,:));
        pol1_d_stds.C_I_Q2_std(vv) = std(vers_pol1_d(vv).(var_series)(end - env.grid.oc + env.idx.C_I_Q2_ind,:));
        pol1_d_stds.C_I_Q3_std(vv) = std(vers_pol1_d(vv).(var_series)(end - env.grid.oc + env.idx.C_I_Q3_ind,:));
        pol1_d_stds.C_I_Q4_std(vv) = std(vers_pol1_d(vv).(var_series)(end - env.grid.oc + env.idx.C_I_Q4_ind,:));
        pol1_d_stds.C_I_Q5_std(vv) = std(vers_pol1_d(vv).(var_series)(end - env.grid.oc + env.idx.C_I_Q5_ind,:));
        pol1_d_stds.C_I_90_std(vv) = std(vers_pol1_d(vv).(var_series)(end - env.grid.oc + env.idx.C_I_90_ind,:));
        pol1_d_stds.C_std(vv)         = std(vers_pol1_d(vv).(var_series)(end - env.grid.oc + env.idx.C_ind,:) );
        pol1_d_stds.Y_std(vv)         = std(vers_pol1_d(vv).(var_series)(end - env.grid.oc + env.idx.Y_ind,:) );
        pol1_d_stds.pi_std(vv)        = std(vers_pol1_d(vv).(var_series)(end - env.grid.oc + env.idx.pi_ind,:) );
        pol1_d_stds.unemp_std(vv)     = std(vers_pol1_d(vv).(var_series)(end - env.grid.oc + env.idx.unemp_ind,:));
        pol1_d_stds.RR_std(vv)        = std(vers_pol1_d(vv).(var_series)(end - env.grid.oc + env.idx.RR_ind,:));
        pol1_d_stds.Profit_std(vv)    = std(vers_pol1_d(vv).(var_series)(end - env.grid.oc + env.idx.Profit_ind,:));
        pol1_d_stds.Profit_FI_std(vv) = std(vers_pol1_d(vv).(var_series)(end - env.grid.oc + env.idx.unemp_ind,:));

        %scalars:
        pol1_d_stds.Q_std(vv) = std(vers_pol1_d(vv).(var_series)(env.grid.numstates - env.grid.os + env.idx.Q_ind,:));
        pol1_d_stds.R_cb_std(vv) = std(vers_pol1_d(vv).(var_series)(env.grid.numstates - env.grid.os + env.idx.R_cb_ind,:));
        pol1_d_stds.W_std(vv) = std(vers_pol1_d(vv).(var_series)(env.grid.numstates - env.grid.os + env.idx.w_ind,:));


        % Manual Standard Deviations:
        pol1_d_stds2.C_10_std(vv) = manual_std(vers_pol1_d(vv).(var_series)(end - env.grid.oc + env.idx.C_10_ind,:));
        pol1_d_stds2.C_Q1_std(vv) = manual_std(vers_pol1_d(vv).(var_series)(end - env.grid.oc + env.idx.C_Q1_ind,:));
        pol1_d_stds2.C_Q2_std(vv) = manual_std(vers_pol1_d(vv).(var_series)(end - env.grid.oc + env.idx.C_Q2_ind,:));
        pol1_d_stds2.C_Q3_std(vv) = manual_std(vers_pol1_d(vv).(var_series)(end - env.grid.oc + env.idx.C_Q3_ind,:));
        pol1_d_stds2.C_Q4_std(vv) = manual_std(vers_pol1_d(vv).(var_series)(end - env.grid.oc + env.idx.C_Q4_ind,:));
        pol1_d_stds2.C_Q5_std(vv) = manual_std(vers_pol1_d(vv).(var_series)(end - env.grid.oc + env.idx.C_Q5_ind,:));
        pol1_d_stds2.C_90_std(vv) = manual_std(vers_pol1_d(vv).(var_series)(end - env.grid.oc + env.idx.C_90_ind,:));
        pol1_d_stds2.C_I_10_std(vv) = manual_std(vers_pol1_d(vv).(var_series)(end - env.grid.oc + env.idx.C_I_10_ind,:));
        pol1_d_stds2.C_I_Q1_std(vv) = manual_std(vers_pol1_d(vv).(var_series)(end - env.grid.oc + env.idx.C_I_Q1_ind,:));
        pol1_d_stds2.C_I_Q2_std(vv) = manual_std(vers_pol1_d(vv).(var_series)(end - env.grid.oc + env.idx.C_I_Q2_ind,:));
        pol1_d_stds2.C_I_Q3_std(vv) = manual_std(vers_pol1_d(vv).(var_series)(end - env.grid.oc + env.idx.C_I_Q3_ind,:));
        pol1_d_stds2.C_I_Q4_std(vv) = manual_std(vers_pol1_d(vv).(var_series)(end - env.grid.oc + env.idx.C_I_Q4_ind,:));
        pol1_d_stds2.C_I_Q5_std(vv) = manual_std(vers_pol1_d(vv).(var_series)(end - env.grid.oc + env.idx.C_I_Q5_ind,:));
        pol1_d_stds2.C_I_90_std(vv) = manual_std(vers_pol1_d(vv).(var_series)(end - env.grid.oc + env.idx.C_I_90_ind,:));
        pol1_d_stds2.C_std(vv)         = manual_std(vers_pol1_d(vv).(var_series)(end - env.grid.oc + env.idx.C_ind,:) );
        pol1_d_stds2.Y_std(vv)         = manual_std(vers_pol1_d(vv).(var_series)(end - env.grid.oc + env.idx.Y_ind,:) );
        pol1_d_stds2.pi_std(vv)        = manual_std(vers_pol1_d(vv).(var_series)(end - env.grid.oc + env.idx.pi_ind,:) );
        pol1_d_stds2.unemp_std(vv)     = manual_std(vers_pol1_d(vv).(var_series)(end - env.grid.oc + env.idx.unemp_ind,:));
        pol1_d_stds2.RR_std(vv)        = manual_std(vers_pol1_d(vv).(var_series)(end - env.grid.oc + env.idx.RR_ind,:));
        pol1_d_stds2.Profit_std(vv)    = manual_std(vers_pol1_d(vv).(var_series)(end - env.grid.oc + env.idx.Profit_ind,:));
        pol1_d_stds2.Profit_FI_std(vv) = manual_std(vers_pol1_d(vv).(var_series)(end - env.grid.oc + env.idx.unemp_ind,:));

        %scalars:
        pol1_d_stds2.Q_std(vv) = manual_std(vers_pol1_d(vv).(var_series)(env.grid.numstates - env.grid.os + env.idx.Q_ind,:));
        pol1_d_stds2.R_cb_std(vv) = manual_std(vers_pol1_d(vv).(var_series)(env.grid.numstates - env.grid.os + env.idx.R_cb_ind,:));
        pol1_d_stds2.W_std(vv) = manual_std(vers_pol1_d(vv).(var_series)(env.grid.numstates - env.grid.os + env.idx.w_ind,:));

        pol1_d_means.W_10_std(vv) = mean(vers_pol1_d(vv).(var_series)(end - env.grid.oc + env.idx.w_10_ind,:));
        pol1_d_means.W_Q1_std(vv) = mean(vers_pol1_d(vv).(var_series)(end - env.grid.oc + env.idx.w_Q1_ind,:));
        pol1_d_means.W_Q2_std(vv) = mean(vers_pol1_d(vv).(var_series)(end - env.grid.oc + env.idx.w_Q2_ind,:));
        pol1_d_means.W_Q3_std(vv) = mean(vers_pol1_d(vv).(var_series)(end - env.grid.oc + env.idx.w_Q3_ind,:));
        pol1_d_means.W_Q4_std(vv) = mean(vers_pol1_d(vv).(var_series)(end - env.grid.oc + env.idx.w_Q4_ind,:));
        pol1_d_means.W_Q5_std(vv) = mean(vers_pol1_d(vv).(var_series)(end - env.grid.oc + env.idx.w_Q5_ind,:));
        pol1_d_means.W_90_std(vv) = mean(vers_pol1_d(vv).(var_series)(end - env.grid.oc + env.idx.w_90_ind,:));
        pol1_d_stds.W_10_std(vv) = std(vers_pol1_d(vv).(var_series)(end - env.grid.oc + env.idx.w_10_ind,:));
        pol1_d_stds.W_Q1_std(vv) = std(vers_pol1_d(vv).(var_series)(end - env.grid.oc + env.idx.w_Q1_ind,:));
        pol1_d_stds.W_Q2_std(vv) = std(vers_pol1_d(vv).(var_series)(end - env.grid.oc + env.idx.w_Q2_ind,:));
        pol1_d_stds.W_Q3_std(vv) = std(vers_pol1_d(vv).(var_series)(end - env.grid.oc + env.idx.w_Q3_ind,:));
        pol1_d_stds.W_Q4_std(vv) = std(vers_pol1_d(vv).(var_series)(end - env.grid.oc + env.idx.w_Q4_ind,:));
        pol1_d_stds.W_Q5_std(vv) = std(vers_pol1_d(vv).(var_series)(end - env.grid.oc + env.idx.w_Q5_ind,:));
        pol1_d_stds.W_90_std(vv) = std(vers_pol1_d(vv).(var_series)(end - env.grid.oc + env.idx.w_90_ind,:));
        pol1_d_stds2.W_10_std(vv) = manual_std(vers_pol1_d(vv).(var_series)(end - env.grid.oc + env.idx.w_10_ind,:));
        pol1_d_stds2.W_Q1_std(vv) = manual_std(vers_pol1_d(vv).(var_series)(end - env.grid.oc + env.idx.w_Q1_ind,:));
        pol1_d_stds2.W_Q2_std(vv) = manual_std(vers_pol1_d(vv).(var_series)(end - env.grid.oc + env.idx.w_Q2_ind,:));
        pol1_d_stds2.W_Q3_std(vv) = manual_std(vers_pol1_d(vv).(var_series)(end - env.grid.oc + env.idx.w_Q3_ind,:));
        pol1_d_stds2.W_Q4_std(vv) = manual_std(vers_pol1_d(vv).(var_series)(end - env.grid.oc + env.idx.w_Q4_ind,:));
        pol1_d_stds2.W_Q5_std(vv) = manual_std(vers_pol1_d(vv).(var_series)(end - env.grid.oc + env.idx.w_Q5_ind,:));
        pol1_d_stds2.W_90_std(vv) = manual_std(vers_pol1_d(vv).(var_series)(end - env.grid.oc + env.idx.w_90_ind,:));

        pol1_d_means.R_a_mean(vv)     = mean(vers_pol1_d(vv).(var_series)(end - env.grid.oc + env.idx.r_a_ind,:));
        pol1_d_stds.R_a_std(vv)     = std(vers_pol1_d(vv).(var_series)(end - env.grid.oc + env.idx.r_a_ind,:));
        pol1_d_stds2.R_a_std(vv)     = manual_std(vers_pol1_d(vv).(var_series)(end - env.grid.oc + env.idx.r_a_ind,:));
        pol1_d_means.RR_new_mean(vv)        = mean(vers_pol1_d(vv).(var_series)(env.grid.numstates - env.grid.os + env.idx.R_cb_ind,:) - vers_pol1_d(vv).(var_series)(end - env.grid.oc + env.idx.pi_ind,:) );
        pol1_d_stds.RR_new_std(vv)        = std(vers_pol1_d(vv).(var_series)(env.grid.numstates - env.grid.os + env.idx.R_cb_ind,:) - vers_pol1_d(vv).(var_series)(end - env.grid.oc + env.idx.pi_ind,:) );
        pol1_d_stds2.RR_new_std(vv)        = manual_std(vers_pol1_d(vv).(var_series)(env.grid.numstates - env.grid.os + env.idx.R_cb_ind,:) - vers_pol1_d(vv).(var_series)(end - env.grid.oc + env.idx.pi_ind,:) );

        pol1_s_means.C_10_mean(vv) = mean(vers_pol1_s(vv).(var_series)(end - env.grid.oc + env.idx.C_10_ind,:));
        pol1_s_means.C_Q1_mean(vv) = mean(vers_pol1_s(vv).(var_series)(end - env.grid.oc + env.idx.C_Q1_ind,:));
        pol1_s_means.C_Q2_mean(vv) = mean(vers_pol1_s(vv).(var_series)(end - env.grid.oc + env.idx.C_Q2_ind,:));
        pol1_s_means.C_Q3_mean(vv) = mean(vers_pol1_s(vv).(var_series)(end - env.grid.oc + env.idx.C_Q3_ind,:));
        pol1_s_means.C_Q4_mean(vv) = mean(vers_pol1_s(vv).(var_series)(end - env.grid.oc + env.idx.C_Q4_ind,:));
        pol1_s_means.C_Q5_mean(vv) = mean(vers_pol1_s(vv).(var_series)(end - env.grid.oc + env.idx.C_Q5_ind,:));
        pol1_s_means.C_90_mean(vv) = mean(vers_pol1_s(vv).(var_series)(end - env.grid.oc + env.idx.C_90_ind,:));
        pol1_s_means.C_I_10_mean(vv) = mean(vers_pol1_s(vv).(var_series)(end - env.grid.oc + env.idx.C_I_10_ind,:));
        pol1_s_means.C_I_Q1_mean(vv) = mean(vers_pol1_s(vv).(var_series)(end - env.grid.oc + env.idx.C_I_Q1_ind,:));
        pol1_s_means.C_I_Q2_mean(vv) = mean(vers_pol1_s(vv).(var_series)(end - env.grid.oc + env.idx.C_I_Q2_ind,:));
        pol1_s_means.C_I_Q3_mean(vv) = mean(vers_pol1_s(vv).(var_series)(end - env.grid.oc + env.idx.C_I_Q3_ind,:));
        pol1_s_means.C_I_Q4_mean(vv) = mean(vers_pol1_s(vv).(var_series)(end - env.grid.oc + env.idx.C_I_Q4_ind,:));
        pol1_s_means.C_I_Q5_mean(vv) = mean(vers_pol1_s(vv).(var_series)(end - env.grid.oc + env.idx.C_I_Q5_ind,:));
        pol1_s_means.C_I_90_mean(vv) = mean(vers_pol1_s(vv).(var_series)(end - env.grid.oc + env.idx.C_I_90_ind,:));
        pol1_s_means.C_mean(vv)         = mean(vers_pol1_s(vv).(var_series)(end - env.grid.oc + env.idx.C_ind,:) );
        pol1_s_means.Y_mean(vv)         = mean(vers_pol1_s(vv).(var_series)(end - env.grid.oc + env.idx.Y_ind,:) );
        pol1_s_means.pi_mean(vv)        = mean(vers_pol1_s(vv).(var_series)(end - env.grid.oc + env.idx.pi_ind,:) );
        pol1_s_means.unemp_mean(vv)     = mean(vers_pol1_s(vv).(var_series)(end - env.grid.oc + env.idx.unemp_ind,:));
        pol1_s_means.RR_mean(vv)        = mean(vers_pol1_s(vv).(var_series)(end - env.grid.oc + env.idx.RR_ind,:));
        pol1_s_means.Profit_mean(vv)    = mean(vers_pol1_s(vv).(var_series)(end - env.grid.oc + env.idx.Profit_ind,:));
        pol1_s_means.Profit_FI_mean(vv) = mean(vers_pol1_s(vv).(var_series)(end - env.grid.oc + env.idx.unemp_ind,:));

        %scalars:
        pol1_s_means.Q_mean(vv) = mean(vers_pol1_s(vv).(var_series)(env.grid.numstates - env.grid.os + env.idx.Q_ind,:));
        pol1_s_means.R_cb_mean(vv) = mean(vers_pol1_s(vv).(var_series)(env.grid.numstates - env.grid.os + env.idx.R_cb_ind,:));
        pol1_s_means.W_mean(vv) = mean(vers_pol1_s(vv).(var_series)(env.grid.numstates - env.grid.os + env.idx.w_ind,:));


        %Standard Deviations
        pol1_s_stds.C_10_std(vv) = std(vers_pol1_s(vv).(var_series)(end - env.grid.oc + env.idx.C_10_ind,:));
        pol1_s_stds.C_Q1_std(vv) = std(vers_pol1_s(vv).(var_series)(end - env.grid.oc + env.idx.C_Q1_ind,:));
        pol1_s_stds.C_Q2_std(vv) = std(vers_pol1_s(vv).(var_series)(end - env.grid.oc + env.idx.C_Q2_ind,:));
        pol1_s_stds.C_Q3_std(vv) = std(vers_pol1_s(vv).(var_series)(end - env.grid.oc + env.idx.C_Q3_ind,:));
        pol1_s_stds.C_Q4_std(vv) = std(vers_pol1_s(vv).(var_series)(end - env.grid.oc + env.idx.C_Q4_ind,:));
        pol1_s_stds.C_Q5_std(vv) = std(vers_pol1_s(vv).(var_series)(end - env.grid.oc + env.idx.C_Q5_ind,:));
        pol1_s_stds.C_90_std(vv) = std(vers_pol1_s(vv).(var_series)(end - env.grid.oc + env.idx.C_90_ind,:));
        pol1_s_stds.C_I_10_std(vv) = std(vers_pol1_s(vv).(var_series)(end - env.grid.oc + env.idx.C_I_10_ind,:));
        pol1_s_stds.C_I_Q1_std(vv) = std(vers_pol1_s(vv).(var_series)(end - env.grid.oc + env.idx.C_I_Q1_ind,:));
        pol1_s_stds.C_I_Q2_std(vv) = std(vers_pol1_s(vv).(var_series)(end - env.grid.oc + env.idx.C_I_Q2_ind,:));
        pol1_s_stds.C_I_Q3_std(vv) = std(vers_pol1_s(vv).(var_series)(end - env.grid.oc + env.idx.C_I_Q3_ind,:));
        pol1_s_stds.C_I_Q4_std(vv) = std(vers_pol1_s(vv).(var_series)(end - env.grid.oc + env.idx.C_I_Q4_ind,:));
        pol1_s_stds.C_I_Q5_std(vv) = std(vers_pol1_s(vv).(var_series)(end - env.grid.oc + env.idx.C_I_Q5_ind,:));
        pol1_s_stds.C_I_90_std(vv) = std(vers_pol1_s(vv).(var_series)(end - env.grid.oc + env.idx.C_I_90_ind,:));
        pol1_s_stds.C_std(vv)         = std(vers_pol1_s(vv).(var_series)(end - env.grid.oc + env.idx.C_ind,:) );
        pol1_s_stds.Y_std(vv)         = std(vers_pol1_s(vv).(var_series)(end - env.grid.oc + env.idx.Y_ind,:) );
        pol1_s_stds.pi_std(vv)        = std(vers_pol1_s(vv).(var_series)(end - env.grid.oc + env.idx.pi_ind,:) );
        pol1_s_stds.unemp_std(vv)     = std(vers_pol1_s(vv).(var_series)(end - env.grid.oc + env.idx.unemp_ind,:));
        pol1_s_stds.RR_std(vv)        = std(vers_pol1_s(vv).(var_series)(end - env.grid.oc + env.idx.RR_ind,:));
        pol1_s_stds.Profit_std(vv)    = std(vers_pol1_s(vv).(var_series)(end - env.grid.oc + env.idx.Profit_ind,:));
        pol1_s_stds.Profit_FI_std(vv) = std(vers_pol1_s(vv).(var_series)(end - env.grid.oc + env.idx.unemp_ind,:));

        %scalars:
        pol1_s_stds.Q_std(vv) = std(vers_pol1_s(vv).(var_series)(env.grid.numstates - env.grid.os + env.idx.Q_ind,:));
        pol1_s_stds.R_cb_std(vv) = std(vers_pol1_s(vv).(var_series)(env.grid.numstates - env.grid.os + env.idx.R_cb_ind,:));
        pol1_s_stds.W_std(vv) = std(vers_pol1_s(vv).(var_series)(env.grid.numstates - env.grid.os + env.idx.w_ind,:));


        % Manual Standard Deviations:
        pol1_s_stds2.C_10_std(vv) = manual_std(vers_pol1_s(vv).(var_series)(end - env.grid.oc + env.idx.C_10_ind,:));
        pol1_s_stds2.C_Q1_std(vv) = manual_std(vers_pol1_s(vv).(var_series)(end - env.grid.oc + env.idx.C_Q1_ind,:));
        pol1_s_stds2.C_Q2_std(vv) = manual_std(vers_pol1_s(vv).(var_series)(end - env.grid.oc + env.idx.C_Q2_ind,:));
        pol1_s_stds2.C_Q3_std(vv) = manual_std(vers_pol1_s(vv).(var_series)(end - env.grid.oc + env.idx.C_Q3_ind,:));
        pol1_s_stds2.C_Q4_std(vv) = manual_std(vers_pol1_s(vv).(var_series)(end - env.grid.oc + env.idx.C_Q4_ind,:));
        pol1_s_stds2.C_Q5_std(vv) = manual_std(vers_pol1_s(vv).(var_series)(end - env.grid.oc + env.idx.C_Q5_ind,:));
        pol1_s_stds2.C_90_std(vv) = manual_std(vers_pol1_s(vv).(var_series)(end - env.grid.oc + env.idx.C_90_ind,:));
        pol1_s_stds2.C_I_10_std(vv) = manual_std(vers_pol1_s(vv).(var_series)(end - env.grid.oc + env.idx.C_I_10_ind,:));
        pol1_s_stds2.C_I_Q1_std(vv) = manual_std(vers_pol1_s(vv).(var_series)(end - env.grid.oc + env.idx.C_I_Q1_ind,:));
        pol1_s_stds2.C_I_Q2_std(vv) = manual_std(vers_pol1_s(vv).(var_series)(end - env.grid.oc + env.idx.C_I_Q2_ind,:));
        pol1_s_stds2.C_I_Q3_std(vv) = manual_std(vers_pol1_s(vv).(var_series)(end - env.grid.oc + env.idx.C_I_Q3_ind,:));
        pol1_s_stds2.C_I_Q4_std(vv) = manual_std(vers_pol1_s(vv).(var_series)(end - env.grid.oc + env.idx.C_I_Q4_ind,:));
        pol1_s_stds2.C_I_Q5_std(vv) = manual_std(vers_pol1_s(vv).(var_series)(end - env.grid.oc + env.idx.C_I_Q5_ind,:));
        pol1_s_stds2.C_I_90_std(vv) = manual_std(vers_pol1_s(vv).(var_series)(end - env.grid.oc + env.idx.C_I_90_ind,:));
        pol1_s_stds2.C_std(vv)         = manual_std(vers_pol1_s(vv).(var_series)(end - env.grid.oc + env.idx.C_ind,:) );
        pol1_s_stds2.Y_std(vv)         = manual_std(vers_pol1_s(vv).(var_series)(end - env.grid.oc + env.idx.Y_ind,:) );
        pol1_s_stds2.pi_std(vv)        = manual_std(vers_pol1_s(vv).(var_series)(end - env.grid.oc + env.idx.pi_ind,:) );
        pol1_s_stds2.unemp_std(vv)     = manual_std(vers_pol1_s(vv).(var_series)(end - env.grid.oc + env.idx.unemp_ind,:));
        pol1_s_stds2.RR_std(vv)        = manual_std(vers_pol1_s(vv).(var_series)(end - env.grid.oc + env.idx.RR_ind,:));
        pol1_s_stds2.Profit_std(vv)    = manual_std(vers_pol1_s(vv).(var_series)(end - env.grid.oc + env.idx.Profit_ind,:));
        pol1_s_stds2.Profit_FI_std(vv) = manual_std(vers_pol1_s(vv).(var_series)(end - env.grid.oc + env.idx.unemp_ind,:));

        %scalars:
        pol1_s_stds2.Q_std(vv) = manual_std(vers_pol1_s(vv).(var_series)(env.grid.numstates - env.grid.os + env.idx.Q_ind,:));
        pol1_s_stds2.R_cb_std(vv) = manual_std(vers_pol1_s(vv).(var_series)(env.grid.numstates - env.grid.os + env.idx.R_cb_ind,:));
        pol1_s_stds2.W_std(vv) = manual_std(vers_pol1_s(vv).(var_series)(env.grid.numstates - env.grid.os + env.idx.w_ind,:));

        pol1_s_means.W_10_std(vv) = mean(vers_pol1_s(vv).(var_series)(end - env.grid.oc + env.idx.w_10_ind,:));
        pol1_s_means.W_Q1_std(vv) = mean(vers_pol1_s(vv).(var_series)(end - env.grid.oc + env.idx.w_Q1_ind,:));
        pol1_s_means.W_Q2_std(vv) = mean(vers_pol1_s(vv).(var_series)(end - env.grid.oc + env.idx.w_Q2_ind,:));
        pol1_s_means.W_Q3_std(vv) = mean(vers_pol1_s(vv).(var_series)(end - env.grid.oc + env.idx.w_Q3_ind,:));
        pol1_s_means.W_Q4_std(vv) = mean(vers_pol1_s(vv).(var_series)(end - env.grid.oc + env.idx.w_Q4_ind,:));
        pol1_s_means.W_Q5_std(vv) = mean(vers_pol1_s(vv).(var_series)(end - env.grid.oc + env.idx.w_Q5_ind,:));
        pol1_s_means.W_90_std(vv) = mean(vers_pol1_s(vv).(var_series)(end - env.grid.oc + env.idx.w_90_ind,:));
        pol1_s_stds.W_10_std(vv) = std(vers_pol1_s(vv).(var_series)(end - env.grid.oc + env.idx.w_10_ind,:));
        pol1_s_stds.W_Q1_std(vv) = std(vers_pol1_s(vv).(var_series)(end - env.grid.oc + env.idx.w_Q1_ind,:));
        pol1_s_stds.W_Q2_std(vv) = std(vers_pol1_s(vv).(var_series)(end - env.grid.oc + env.idx.w_Q2_ind,:));
        pol1_s_stds.W_Q3_std(vv) = std(vers_pol1_s(vv).(var_series)(end - env.grid.oc + env.idx.w_Q3_ind,:));
        pol1_s_stds.W_Q4_std(vv) = std(vers_pol1_s(vv).(var_series)(end - env.grid.oc + env.idx.w_Q4_ind,:));
        pol1_s_stds.W_Q5_std(vv) = std(vers_pol1_s(vv).(var_series)(end - env.grid.oc + env.idx.w_Q5_ind,:));
        pol1_s_stds.W_90_std(vv) = std(vers_pol1_s(vv).(var_series)(end - env.grid.oc + env.idx.w_90_ind,:));
        pol1_s_stds2.W_10_std(vv) = manual_std(vers_pol1_s(vv).(var_series)(end - env.grid.oc + env.idx.w_10_ind,:));
        pol1_s_stds2.W_Q1_std(vv) = manual_std(vers_pol1_s(vv).(var_series)(end - env.grid.oc + env.idx.w_Q1_ind,:));
        pol1_s_stds2.W_Q2_std(vv) = manual_std(vers_pol1_s(vv).(var_series)(end - env.grid.oc + env.idx.w_Q2_ind,:));
        pol1_s_stds2.W_Q3_std(vv) = manual_std(vers_pol1_s(vv).(var_series)(end - env.grid.oc + env.idx.w_Q3_ind,:));
        pol1_s_stds2.W_Q4_std(vv) = manual_std(vers_pol1_s(vv).(var_series)(end - env.grid.oc + env.idx.w_Q4_ind,:));
        pol1_s_stds2.W_Q5_std(vv) = manual_std(vers_pol1_s(vv).(var_series)(end - env.grid.oc + env.idx.w_Q5_ind,:));
        pol1_s_stds2.W_90_std(vv) = manual_std(vers_pol1_s(vv).(var_series)(end - env.grid.oc + env.idx.w_90_ind,:));

        pol1_s_means.R_a_mean(vv)    = mean(vers_pol1_s(vv).(var_series)(end - env.grid.oc + env.idx.r_a_ind,:));
        pol1_s_stds.R_a_std(vv)      = std(vers_pol1_s(vv).(var_series)(end - env.grid.oc + env.idx.r_a_ind,:));
        pol1_s_stds2.R_a_std(vv)     = manual_std(vers_pol1_s(vv).(var_series)(end - env.grid.oc + env.idx.r_a_ind,:));
        pol1_s_means.RR_new_mean(vv) = mean(vers_pol1_s(vv).(var_series)(env.grid.numstates - env.grid.os + env.idx.R_cb_ind,:) - vers_pol1_s(vv).(var_series)(end - env.grid.oc + env.idx.pi_ind,:) );
        pol1_s_stds.RR_new_std(vv)   = std(vers_pol1_s(vv).(var_series)(env.grid.numstates - env.grid.os + env.idx.R_cb_ind,:) - vers_pol1_s(vv).(var_series)(end - env.grid.oc + env.idx.pi_ind,:) );
        pol1_s_stds2.RR_new_std(vv)  = manual_std(vers_pol1_s(vv).(var_series)(env.grid.numstates - env.grid.os + env.idx.R_cb_ind,:) - vers_pol1_s(vv).(var_series)(end - env.grid.oc + env.idx.pi_ind,:) );

    end

    %%

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % (4) Plotting Function Calls: Disaggregated Plots
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % (4.1) Create  plots for disaggregated policy-rule plots
    %%%%%%%%%%

    disagg_plots(pol1_stds, pol1_s_stds, pol1_d_stds, phipi1_vals, string(fieldnames(pol1_stds)), "SS_Dev", pltroot)
    close all;
    frontier_plots(pol1_d_stds, pol1_s_stds, pol1, ["Demand", "Supply"] ,phipi1_vals, baseline_phi_ind, string(fieldnames(pol1_stds)), "SS_Dev", pltroot, 'supplydemand');
    close all;


end
%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% (6) Plotting functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% (6.1)
% function one_frontier_plot()
%
% Arguments:
% (1) pol_std: Struct of volatilities for each variable by phi_pi value for
% specified policy rule
% (2) phi_pi_vals: Iterables for phi_pi values under specified policy rule
% (3) type: String indentifier for plots (MATLAB std() = "SS_Dev, manual =
% "SS_Dev_manual)
% (4) pltroot: File path to save plots
%
% Description: Produces frontier plot for specified policy rule
%%%%%%%%%%

function one_frontier_plot(pol_std, phi_pi_vals, fields, type, pltroot)

% Point label offset for plotting:
dx = 0.0001;
dy = 0.0001;

if strcmp(type, "Mean")
    pol_infl = getfield(pol_std, "pi_mean");
else
    pol_infl = getfield(pol_std, "pi_std");
end

% Normalize plot relative to baseline phi_pi value
pol_infl_norm = 100 * (pol_infl - pol_infl(2));


for f = 1:length(fields)

    char_f = char(fields(f));

    if strcmp(type, "Mean")
        para_name = string(char_f(1:end-5)); %removes '_mean' from string name
    else
        para_name = string(char_f(1:end-4)); %removes '_std' from string name
    end
    x_lab = "pi";

    pol_vals = getfield(pol_std, fields(f));

    pol_vals_norm = 100 * (pol_vals - pol_vals(2));


    %%%%%% Actual Plotting %%%%%%%%%
    figure('name', strcat(para_name,"_frontier"));

    hold on
    plot(pol_infl_norm, pol_vals_norm, "--o", 'Color', "blue", "LineWidth", 1.5);
    yline(0);
    xline(0);
    text(pol_infl_norm+dx, pol_vals_norm+dy, string(phi_pi_vals));
    hold off

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
        base_fl = pltroot.onefrontierplot_mean;
    else
        base_fl = pltroot.onefrontierplot;
    end


    if ~exist(base_fl, 'dir')
        mkdir(base_fl)
    end
    exportgraphics(fi, base_fl + type + '_phi_pi_vs_' + para_name + '.png', 'Resolution', 300);

end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% (6.2)
% function frontier_plots()
%
% Arguments:
% (1) pol1_std: Struct containing volatility calculations for first policy rule
% (2) pol2_std: Struct containing volatility calculations for second policy rule
% (3) pol1: Name of first policy rule
% (4) phi_pi_pol1 Iterables for phi_pi values the  policy
% rule
% (5) labels: Plot labels for comparative phi-pi plots
% (6) fields: Struct fieldnames for untransformed variables
% (7) type: String indentifier for plots (MATLAB std() = "SS_Dev, manual =
% "SS_Dev_manual)
% (8) pltroot: File path to save plots
%
% Description: Create frontier plots of policy rule 1 vs policy rule 2
% (normalized by reach rule's baseline phi-pi value)
%%%%%%%%%%

function frontier_plots(pol1_std, pol2_std, pol1, pol_names ,phi_pi_pol1, baseline_phi_ind, fields, type, pltroot, dirtree)

% Point label offset for plotting:
dx = 0.0001;
dy = 0.0001;

pol2 = "Taylor";

% Asserting that fields is of type string(array)
if strcmp(type, "Mean")
    pol2_infl = getfield(pol2_std, "pi_mean");
    pol1_infl = getfield(pol1_std, "pi_mean");
else
    pol2_infl = getfield(pol2_std, "pi_std");
    pol1_infl = getfield(pol1_std, "pi_std");
end

if strcmp(pol1, pol2)
    pol2_infl_norm = 100 * (pol2_infl - pol2_infl(baseline_phi_ind(2)));
    pol1_infl_norm = 100 * (pol1_infl - pol1_infl(baseline_phi_ind(1)));
else
    pol2_infl_norm = 100 * (pol2_infl - pol1_infl(baseline_phi_ind(1)));
    pol1_infl_norm = 100 * (pol1_infl - pol1_infl(baseline_phi_ind(1)));
end

for f = 1:length(fields)

    char_f = char(fields(f));
    para_name = string(char_f(1:end-4)); %removes '_std' from string name
    x_lab = "pi";


    %Variable relative to Taylor Baseline:
    pol2_vals = getfield(pol2_std, fields(f));
    pol1_vals = getfield(pol1_std, fields(f));

    if strcmp(pol1, pol2)
        pol2_vals_norm = 100 *(pol2_vals - pol2_vals(baseline_phi_ind(2)));
        pol1_vals_norm = 100 * (pol1_vals - pol1_vals(baseline_phi_ind(1)));
    else
        pol2_vals_norm = 100 *(pol2_vals - pol1_vals(baseline_phi_ind(1)));
        pol1_vals_norm = 100 * (pol1_vals - pol1_vals(baseline_phi_ind(1)));
    end

    %%%%%% Actual Plotting %%%%%%%%%
    figure('name', strcat(para_name,"_frontier"));

    hold on
    plot(pol1_infl_norm, pol1_vals_norm, "--o", 'Color', "blue", "LineWidth", 1.5);
    plot(pol2_infl_norm,  pol2_vals_norm, "--o", 'Color', "red", "LineWidth", 1.5);
    yline(0);
    xline(0);
    text(pol1_infl_norm+dx, pol1_vals_norm+dy, string(phi_pi_pol1));
    text(pol2_infl_norm+dx, pol2_vals_norm+dy, string(phi_pi_pol1));
    hold off

    legend(pol_names)
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
    elseif strcmp(dirtree, 'kappa')
        base_fl = pltroot.kappafrontierplots;
    elseif strcmp(dirtree, 'supplydemand')
        base_fl = pltroot.SDfrontierplots;
    else
        warn('Directory not defined!')
        assert(false);
    end
    if ~exist(base_fl, 'dir')
        mkdir(base_fl)
    end
    exportgraphics(fi, base_fl + type + '_phi_pi_vs_' + para_name + '.png', 'Resolution', 300);
end

end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% (6.3)
% function disagg_plots()
%
% Arguments:
% (1) all_std: Struct of volatilities for each variable by phi_pi value
% (under simulation where both supply and demand shocks are active)
% (2) supply_std: Struct of volatilities for each variable by phi_pi value
% (under simulation with isolated supply shocks)
% (3) demand_std: Struct of volatilities for each variable by phi_pi value
% (under simulation with isolated demand shocks)
% (4) phi_pi_s: Iterables for phi_pi values under specified policy rule
% (5) type: String indentifier for plots (MATLAB std() = "SS_Dev, manual =
% "SS_Dev_manual)
% (6) pltroot: File path to save plots
%
% Description: Produces phi-pi disaggregation plots
%%%%%%%%%%

function disagg_plots(all_std, supply_std, demand_std, phi_pis, fields, type, pltroot)
for f = 1:length(fields) %% Assuming that the two have same fields

    char_f = char(fields(f));
    para_name = string(char_f(1:end-4)); %removes '_std' from string name
    x_lab = "\phi_{\pi}";

    figure('name', strcat(para_name,"_phi_pi"));
    hold on
    plot(phi_pis, 100 * getfield(all_std, fields(f)), "--o", 'Color', "black", "LineWidth", 1.5)
    plot(phi_pis, 100 * getfield(supply_std, fields(f)), "--o", 'Color', "red", "LineWidth", 1.5)
    plot(phi_pis, 100 * getfield(demand_std, fields(f)), "--o", 'Color', "blue", "LineWidth", 1.5)
    hold off
    title('Deviation of ' + strrep(para_name, "_", "-") + " v.s " + x_lab);
    xlabel(x_lab);
    ylabel('Deviation of ' + strrep(para_name, "_", "-") + ' (in %)');

    l = legend('All Shocks', 'Supply Shocks', 'Demand Shocks');
    legend('boxoff')

    fi = gcf;
    base = pltroot.supply_demand_plots;

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

