 



    pol1_01_means.C_10_mean(vv) = mean(Pol1_results0_1(vv).(var_series)(end - env.grid.oc + env.idx.C_10_ind,:));
    pol1_01_means.C_Q1_mean(vv) = mean(Pol1_results0_1(vv).(var_series)(end - env.grid.oc + env.idx.C_Q1_ind,:));
    pol1_01_means.C_Q2_mean(vv) = mean(Pol1_results0_1(vv).(var_series)(end - env.grid.oc + env.idx.C_Q2_ind,:));
    pol1_01_means.C_Q3_mean(vv) = mean(Pol1_results0_1(vv).(var_series)(end - env.grid.oc + env.idx.C_Q3_ind,:));
    pol1_01_means.C_Q4_mean(vv) = mean(Pol1_results0_1(vv).(var_series)(end - env.grid.oc + env.idx.C_Q4_ind,:));
    pol1_01_means.C_Q5_mean(vv) = mean(Pol1_results0_1(vv).(var_series)(end - env.grid.oc + env.idx.C_Q5_ind,:));
    pol1_01_means.C_90_mean(vv) = mean(Pol1_results0_1(vv).(var_series)(end - env.grid.oc + env.idx.C_90_ind,:));


    pol1_01_means.W_10_mean(vv) = mean(Pol1_results0_1(vv).(var_series)(end - env.grid.oc + env.idx.w_10_ind,:));
    pol1_01_means.W_Q1_mean(vv) = mean(Pol1_results0_1(vv).(var_series)(end - env.grid.oc + env.idx.w_Q1_ind,:));
    pol1_01_means.W_Q2_mean(vv) = mean(Pol1_results0_1(vv).(var_series)(end - env.grid.oc + env.idx.w_Q2_ind,:));
    pol1_01_means.W_Q3_mean(vv) = mean(Pol1_results0_1(vv).(var_series)(end - env.grid.oc + env.idx.w_Q3_ind,:));
    pol1_01_means.W_Q4_mean(vv) = mean(Pol1_results0_1(vv).(var_series)(end - env.grid.oc + env.idx.w_Q4_ind,:));
    pol1_01_means.W_Q5_mean(vv) = mean(Pol1_results0_1(vv).(var_series)(end - env.grid.oc + env.idx.w_Q5_ind,:));
    pol1_01_means.W_90_mean(vv) = mean(Pol1_results0_1(vv).(var_series)(end - env.grid.oc + env.idx.w_90_ind,:));

    pol1_01_means.C_I_10_mean(vv) = mean(Pol1_results0_1(vv).(var_series)(end - env.grid.oc + env.idx.C_I_10_ind,:));
    pol1_01_means.C_I_Q1_mean(vv) = mean(Pol1_results0_1(vv).(var_series)(end - env.grid.oc + env.idx.C_I_Q1_ind,:));
    pol1_01_means.C_I_Q2_mean(vv) = mean(Pol1_results0_1(vv).(var_series)(end - env.grid.oc + env.idx.C_I_Q2_ind,:));
    pol1_01_means.C_I_Q3_mean(vv) = mean(Pol1_results0_1(vv).(var_series)(end - env.grid.oc + env.idx.C_I_Q3_ind,:));
    pol1_01_means.C_I_Q4_mean(vv) = mean(Pol1_results0_1(vv).(var_series)(end - env.grid.oc + env.idx.C_I_Q4_ind,:));
    pol1_01_means.C_I_Q5_mean(vv) = mean(Pol1_results0_1(vv).(var_series)(end - env.grid.oc + env.idx.C_I_Q5_ind,:));
    pol1_01_means.C_I_90_mean(vv) = mean(Pol1_results0_1(vv).(var_series)(end - env.grid.oc + env.idx.C_I_90_ind,:));
    pol1_01_means.C_mean(vv)         = mean(Pol1_results0_1(vv).(var_series)(end - env.grid.oc + env.idx.C_ind,:) );
    pol1_01_means.Y_mean(vv)         = mean(Pol1_results0_1(vv).(var_series)(end - env.grid.oc + env.idx.Y_ind,:) );
    pol1_01_means.pi_mean(vv)        = mean(Pol1_results0_1(vv).(var_series)(end - env.grid.oc + env.idx.pi_ind,:) );
    pol1_01_means.unemp_mean(vv)     = mean(Pol1_results0_1(vv).(var_series)(end - env.grid.oc + env.idx.unemp_ind,:));
    pol1_01_means.RR_mean(vv)        = mean(Pol1_results0_1(vv).(var_series)(end - env.grid.oc + env.idx.RR_ind,:));
    pol1_01_means.Profit_mean(vv)    = mean(Pol1_results0_1(vv).(var_series)(end - env.grid.oc + env.idx.Profit_ind,:));
    pol1_01_means.Profit_FI_mean(vv) = mean(Pol1_results0_1(vv).(var_series)(end - env.grid.oc + env.idx.unemp_ind,:));

    %scalars:
    pol1_01_means.Q_mean(vv) = mean(Pol1_results0_1(vv).(var_series)(env.grid.numstates - env.grid.os + env.idx.Q_ind,:));
    pol1_01_means.R_cb_mean(vv) = mean(Pol1_results0_1(vv).(var_series)(env.grid.numstates - env.grid.os + env.idx.R_cb_ind,:));
    pol1_01_means.W_mean(vv) = mean(Pol1_results0_1(vv).(var_series)(env.grid.numstates - env.grid.os + env.idx.w_ind,:));




    %Standard Deviations
    pol1_01_stds.C_10_std(vv) = std(Pol1_results0_1(vv).(var_series)(end - env.grid.oc + env.idx.C_10_ind,:));
    pol1_01_stds.C_Q1_std(vv) = std(Pol1_results0_1(vv).(var_series)(end - env.grid.oc + env.idx.C_Q1_ind,:));
    pol1_01_stds.C_Q2_std(vv) = std(Pol1_results0_1(vv).(var_series)(end - env.grid.oc + env.idx.C_Q2_ind,:));
    pol1_01_stds.C_Q3_std(vv) = std(Pol1_results0_1(vv).(var_series)(end - env.grid.oc + env.idx.C_Q3_ind,:));
    pol1_01_stds.C_Q4_std(vv) = std(Pol1_results0_1(vv).(var_series)(end - env.grid.oc + env.idx.C_Q4_ind,:));
    pol1_01_stds.C_Q5_std(vv) = std(Pol1_results0_1(vv).(var_series)(end - env.grid.oc + env.idx.C_Q5_ind,:));
    pol1_01_stds.C_90_std(vv) = std(Pol1_results0_1(vv).(var_series)(end - env.grid.oc + env.idx.C_90_ind,:));

    pol1_01_stds.W_10_std(vv) = std(Pol1_results0_1(vv).(var_series)(end - env.grid.oc + env.idx.w_10_ind,:));
    pol1_01_stds.W_Q1_std(vv) = std(Pol1_results0_1(vv).(var_series)(end - env.grid.oc + env.idx.w_Q1_ind,:));
    pol1_01_stds.W_Q2_std(vv) = std(Pol1_results0_1(vv).(var_series)(end - env.grid.oc + env.idx.w_Q2_ind,:));
    pol1_01_stds.W_Q3_std(vv) = std(Pol1_results0_1(vv).(var_series)(end - env.grid.oc + env.idx.w_Q3_ind,:));
    pol1_01_stds.W_Q4_std(vv) = std(Pol1_results0_1(vv).(var_series)(end - env.grid.oc + env.idx.w_Q4_ind,:));
    pol1_01_stds.W_Q5_std(vv) = std(Pol1_results0_1(vv).(var_series)(end - env.grid.oc + env.idx.w_Q5_ind,:));
    pol1_01_stds.W_90_std(vv) = std(Pol1_results0_1(vv).(var_series)(end - env.grid.oc + env.idx.w_90_ind,:));



    pol1_01_stds.C_I_10_std(vv) = std(Pol1_results0_1(vv).(var_series)(end - env.grid.oc + env.idx.C_I_10_ind,:));
    pol1_01_stds.C_I_Q1_std(vv) = std(Pol1_results0_1(vv).(var_series)(end - env.grid.oc + env.idx.C_I_Q1_ind,:));
    pol1_01_stds.C_I_Q2_std(vv) = std(Pol1_results0_1(vv).(var_series)(end - env.grid.oc + env.idx.C_I_Q2_ind,:));
    pol1_01_stds.C_I_Q3_std(vv) = std(Pol1_results0_1(vv).(var_series)(end - env.grid.oc + env.idx.C_I_Q3_ind,:));
    pol1_01_stds.C_I_Q4_std(vv) = std(Pol1_results0_1(vv).(var_series)(end - env.grid.oc + env.idx.C_I_Q4_ind,:));
    pol1_01_stds.C_I_Q5_std(vv) = std(Pol1_results0_1(vv).(var_series)(end - env.grid.oc + env.idx.C_I_Q5_ind,:));
    pol1_01_stds.C_I_90_std(vv) = std(Pol1_results0_1(vv).(var_series)(end - env.grid.oc + env.idx.C_I_90_ind,:));
    pol1_01_stds.C_std(vv)         = std(Pol1_results0_1(vv).(var_series)(end - env.grid.oc + env.idx.C_ind,:) );
    pol1_01_stds.Y_std(vv)         = std(Pol1_results0_1(vv).(var_series)(end - env.grid.oc + env.idx.Y_ind,:) );
    pol1_01_stds.pi_std(vv)        = std(Pol1_results0_1(vv).(var_series)(end - env.grid.oc + env.idx.pi_ind,:) );
    pol1_01_stds.unemp_std(vv)     = std(Pol1_results0_1(vv).(var_series)(end - env.grid.oc + env.idx.unemp_ind,:));
    pol1_01_stds.RR_std(vv)        = std(Pol1_results0_1(vv).(var_series)(end - env.grid.oc + env.idx.RR_ind,:));
    pol1_01_stds.Profit_std(vv)    = std(Pol1_results0_1(vv).(var_series)(end - env.grid.oc + env.idx.Profit_ind,:));
    pol1_01_stds.Profit_FI_std(vv) = std(Pol1_results0_1(vv).(var_series)(end - env.grid.oc + env.idx.unemp_ind,:));

    %scalars:
    pol1_01_stds.Q_std(vv) = std(Pol1_results0_1(vv).(var_series)(env.grid.numstates - env.grid.os + env.idx.Q_ind,:));
    pol1_01_stds.R_cb_std(vv) = std(Pol1_results0_1(vv).(var_series)(env.grid.numstates - env.grid.os + env.idx.R_cb_ind,:));
    pol1_01_stds.W_std(vv) = std(Pol1_results0_1(vv).(var_series)(env.grid.numstates - env.grid.os + env.idx.w_ind,:));


    %Manual Standard Deviations:
    pol1_01_stds2.C_10_std(vv) = manual_std(Pol1_results0_1(vv).(var_series)(end - env.grid.oc + env.idx.C_10_ind,:));
    pol1_01_stds2.C_Q1_std(vv) = manual_std(Pol1_results0_1(vv).(var_series)(end - env.grid.oc + env.idx.C_Q1_ind,:));
    pol1_01_stds2.C_Q2_std(vv) = manual_std(Pol1_results0_1(vv).(var_series)(end - env.grid.oc + env.idx.C_Q2_ind,:));
    pol1_01_stds2.C_Q3_std(vv) = manual_std(Pol1_results0_1(vv).(var_series)(end - env.grid.oc + env.idx.C_Q3_ind,:));
    pol1_01_stds2.C_Q4_std(vv) = manual_std(Pol1_results0_1(vv).(var_series)(end - env.grid.oc + env.idx.C_Q4_ind,:));
    pol1_01_stds2.C_Q5_std(vv) = manual_std(Pol1_results0_1(vv).(var_series)(end - env.grid.oc + env.idx.C_Q5_ind,:));
    pol1_01_stds2.C_90_std(vv) = manual_std(Pol1_results0_1(vv).(var_series)(end - env.grid.oc + env.idx.C_90_ind,:));


    pol1_01_stds2.W_10_std(vv) = manual_std(Pol1_results0_1(vv).(var_series)(end - env.grid.oc + env.idx.w_10_ind,:));
    pol1_01_stds2.W_Q1_std(vv) = manual_std(Pol1_results0_1(vv).(var_series)(end - env.grid.oc + env.idx.w_Q1_ind,:));
    pol1_01_stds2.W_Q2_std(vv) = manual_std(Pol1_results0_1(vv).(var_series)(end - env.grid.oc + env.idx.w_Q2_ind,:));
    pol1_01_stds2.W_Q3_std(vv) = manual_std(Pol1_results0_1(vv).(var_series)(end - env.grid.oc + env.idx.w_Q3_ind,:));
    pol1_01_stds2.W_Q4_std(vv) = manual_std(Pol1_results0_1(vv).(var_series)(end - env.grid.oc + env.idx.w_Q4_ind,:));
    pol1_01_stds2.W_Q5_std(vv) = manual_std(Pol1_results0_1(vv).(var_series)(end - env.grid.oc + env.idx.w_Q5_ind,:));
    pol1_01_stds2.W_90_std(vv) = manual_std(Pol1_results0_1(vv).(var_series)(end - env.grid.oc + env.idx.w_90_ind,:));

    pol1_01_stds2.C_I_10_std(vv) = manual_std(Pol1_results0_1(vv).(var_series)(end - env.grid.oc + env.idx.C_I_10_ind,:));
    pol1_01_stds2.C_I_Q1_std(vv) = manual_std(Pol1_results0_1(vv).(var_series)(end - env.grid.oc + env.idx.C_I_Q1_ind,:));
    pol1_01_stds2.C_I_Q2_std(vv) = manual_std(Pol1_results0_1(vv).(var_series)(end - env.grid.oc + env.idx.C_I_Q2_ind,:));
    pol1_01_stds2.C_I_Q3_std(vv) = manual_std(Pol1_results0_1(vv).(var_series)(end - env.grid.oc + env.idx.C_I_Q3_ind,:));
    pol1_01_stds2.C_I_Q4_std(vv) = manual_std(Pol1_results0_1(vv).(var_series)(end - env.grid.oc + env.idx.C_I_Q4_ind,:));
    pol1_01_stds2.C_I_Q5_std(vv) = manual_std(Pol1_results0_1(vv).(var_series)(end - env.grid.oc + env.idx.C_I_Q5_ind,:));
    pol1_01_stds2.C_I_90_std(vv) = manual_std(Pol1_results0_1(vv).(var_series)(end - env.grid.oc + env.idx.C_I_90_ind,:));
    pol1_01_stds2.C_std(vv)         = manual_std(Pol1_results0_1(vv).(var_series)(end - env.grid.oc + env.idx.C_ind,:) );
    pol1_01_stds2.Y_std(vv)         = manual_std(Pol1_results0_1(vv).(var_series)(end - env.grid.oc + env.idx.Y_ind,:) );
    pol1_01_stds2.pi_std(vv)        = manual_std(Pol1_results0_1(vv).(var_series)(end - env.grid.oc + env.idx.pi_ind,:) );
    pol1_01_stds2.unemp_std(vv)     = manual_std(Pol1_results0_1(vv).(var_series)(end - env.grid.oc + env.idx.unemp_ind,:));
    pol1_01_stds2.RR_std(vv)        = manual_std(Pol1_results0_1(vv).(var_series)(end - env.grid.oc + env.idx.RR_ind,:));
    pol1_01_stds2.Profit_std(vv)    = manual_std(Pol1_results0_1(vv).(var_series)(end - env.grid.oc + env.idx.Profit_ind,:));
    pol1_01_stds2.Profit_FI_std(vv) = manual_std(Pol1_results0_1(vv).(var_series)(end - env.grid.oc + env.idx.unemp_ind,:));

    %scalars:
    pol1_01_stds2.Q_std(vv) = manual_std(Pol1_results0_1(vv).(var_series)(env.grid.numstates - env.grid.os + env.idx.Q_ind,:));
    pol1_01_stds2.R_cb_std(vv) = manual_std(Pol1_results0_1(vv).(var_series)(env.grid.numstates - env.grid.os + env.idx.R_cb_ind,:));
    pol1_01_stds2.W_std(vv) = manual_std(Pol1_results0_1(vv).(var_series)(env.grid.numstates - env.grid.os + env.idx.w_ind,:));

    pol1_01_means.R_a_mean(vv)     = mean(Pol1_results0_1(vv).(var_series)(end - env.grid.oc + env.idx.r_a_ind,:));
    pol1_01_stds.R_a_std(vv)     = std(Pol1_results0_1(vv).(var_series)(end - env.grid.oc + env.idx.r_a_ind,:));
    pol1_01_stds2.R_a_std(vv)     = manual_std(Pol1_results0_1(vv).(var_series)(end - env.grid.oc + env.idx.r_a_ind,:));



    pol1_01_means.RR_new_mean(vv)        = mean(Pol1_results0_1(vv).(var_series)(env.grid.numstates - env.grid.os + env.idx.R_cb_ind,:) - Pol1_results0_1(vv).(var_series)(end - env.grid.oc + env.idx.pi_ind,:) );
    pol1_01_stds.RR_new_std(vv)        = std(Pol1_results0_1(vv).(var_series)(env.grid.numstates - env.grid.os + env.idx.R_cb_ind,:) - Pol1_results0_1(vv).(var_series)(end - env.grid.oc + env.idx.pi_ind,:) );
    pol1_01_stds2.RR_new_std(vv)        = manual_std(Pol1_results0_1(vv).(var_series)(env.grid.numstates - env.grid.os + env.idx.R_cb_ind,:) - Pol1_results0_1(vv).(var_series)(end - env.grid.oc + env.idx.pi_ind,:) );









%%%%%%%%%%%%%%%%%%%%%%%%

%var_series = 'X_t_series';

    pol1_03_means.C_10_mean(vv) = mean(Pol1_results0_3(vv).(var_series)(end - env.grid.oc + env.idx.C_10_ind,:));
    pol1_03_means.C_Q1_mean(vv) = mean(Pol1_results0_3(vv).(var_series)(end - env.grid.oc + env.idx.C_Q1_ind,:));
    pol1_03_means.C_Q2_mean(vv) = mean(Pol1_results0_3(vv).(var_series)(end - env.grid.oc + env.idx.C_Q2_ind,:));
    pol1_03_means.C_Q3_mean(vv) = mean(Pol1_results0_3(vv).(var_series)(end - env.grid.oc + env.idx.C_Q3_ind,:));
    pol1_03_means.C_Q4_mean(vv) = mean(Pol1_results0_3(vv).(var_series)(end - env.grid.oc + env.idx.C_Q4_ind,:));
    pol1_03_means.C_Q5_mean(vv) = mean(Pol1_results0_3(vv).(var_series)(end - env.grid.oc + env.idx.C_Q5_ind,:));
    pol1_03_means.C_90_mean(vv) = mean(Pol1_results0_3(vv).(var_series)(end - env.grid.oc + env.idx.C_90_ind,:));


    pol1_03_means.W_10_mean(vv) = mean(Pol1_results0_3(vv).(var_series)(end - env.grid.oc + env.idx.w_10_ind,:));
    pol1_03_means.W_Q1_mean(vv) = mean(Pol1_results0_3(vv).(var_series)(end - env.grid.oc + env.idx.w_Q1_ind,:));
    pol1_03_means.W_Q2_mean(vv) = mean(Pol1_results0_3(vv).(var_series)(end - env.grid.oc + env.idx.w_Q2_ind,:));
    pol1_03_means.W_Q3_mean(vv) = mean(Pol1_results0_3(vv).(var_series)(end - env.grid.oc + env.idx.w_Q3_ind,:));
    pol1_03_means.W_Q4_mean(vv) = mean(Pol1_results0_3(vv).(var_series)(end - env.grid.oc + env.idx.w_Q4_ind,:));
    pol1_03_means.W_Q5_mean(vv) = mean(Pol1_results0_3(vv).(var_series)(end - env.grid.oc + env.idx.w_Q5_ind,:));
    pol1_03_means.W_90_mean(vv) = mean(Pol1_results0_3(vv).(var_series)(end - env.grid.oc + env.idx.w_90_ind,:));

    pol1_03_means.C_I_10_mean(vv) = mean(Pol1_results0_3(vv).(var_series)(end - env.grid.oc + env.idx.C_I_10_ind,:));
    pol1_03_means.C_I_Q1_mean(vv) = mean(Pol1_results0_3(vv).(var_series)(end - env.grid.oc + env.idx.C_I_Q1_ind,:));
    pol1_03_means.C_I_Q2_mean(vv) = mean(Pol1_results0_3(vv).(var_series)(end - env.grid.oc + env.idx.C_I_Q2_ind,:));
    pol1_03_means.C_I_Q3_mean(vv) = mean(Pol1_results0_3(vv).(var_series)(end - env.grid.oc + env.idx.C_I_Q3_ind,:));
    pol1_03_means.C_I_Q4_mean(vv) = mean(Pol1_results0_3(vv).(var_series)(end - env.grid.oc + env.idx.C_I_Q4_ind,:));
    pol1_03_means.C_I_Q5_mean(vv) = mean(Pol1_results0_3(vv).(var_series)(end - env.grid.oc + env.idx.C_I_Q5_ind,:));
    pol1_03_means.C_I_90_mean(vv) = mean(Pol1_results0_3(vv).(var_series)(end - env.grid.oc + env.idx.C_I_90_ind,:));
    pol1_03_means.C_mean(vv)         = mean(Pol1_results0_3(vv).(var_series)(end - env.grid.oc + env.idx.C_ind,:) );
    pol1_03_means.Y_mean(vv)         = mean(Pol1_results0_3(vv).(var_series)(end - env.grid.oc + env.idx.Y_ind,:) );
    pol1_03_means.pi_mean(vv)        = mean(Pol1_results0_3(vv).(var_series)(end - env.grid.oc + env.idx.pi_ind,:) );
    pol1_03_means.unemp_mean(vv)     = mean(Pol1_results0_3(vv).(var_series)(end - env.grid.oc + env.idx.unemp_ind,:));
    pol1_03_means.RR_mean(vv)        = mean(Pol1_results0_3(vv).(var_series)(end - env.grid.oc + env.idx.RR_ind,:));
    pol1_03_means.Profit_mean(vv)    = mean(Pol1_results0_3(vv).(var_series)(end - env.grid.oc + env.idx.Profit_ind,:));
    pol1_03_means.Profit_FI_mean(vv) = mean(Pol1_results0_3(vv).(var_series)(end - env.grid.oc + env.idx.unemp_ind,:));

    %scalars:
    pol1_03_means.Q_mean(vv) = mean(Pol1_results0_3(vv).(var_series)(env.grid.numstates - env.grid.os + env.idx.Q_ind,:));
    pol1_03_means.R_cb_mean(vv) = mean(Pol1_results0_3(vv).(var_series)(env.grid.numstates - env.grid.os + env.idx.R_cb_ind,:));
    pol1_03_means.W_mean(vv) = mean(Pol1_results0_3(vv).(var_series)(env.grid.numstates - env.grid.os + env.idx.w_ind,:));




    %Standard Deviations
    pol1_03_stds.C_10_std(vv) = std(Pol1_results0_3(vv).(var_series)(end - env.grid.oc + env.idx.C_10_ind,:));
    pol1_03_stds.C_Q1_std(vv) = std(Pol1_results0_3(vv).(var_series)(end - env.grid.oc + env.idx.C_Q1_ind,:));
    pol1_03_stds.C_Q2_std(vv) = std(Pol1_results0_3(vv).(var_series)(end - env.grid.oc + env.idx.C_Q2_ind,:));
    pol1_03_stds.C_Q3_std(vv) = std(Pol1_results0_3(vv).(var_series)(end - env.grid.oc + env.idx.C_Q3_ind,:));
    pol1_03_stds.C_Q4_std(vv) = std(Pol1_results0_3(vv).(var_series)(end - env.grid.oc + env.idx.C_Q4_ind,:));
    pol1_03_stds.C_Q5_std(vv) = std(Pol1_results0_3(vv).(var_series)(end - env.grid.oc + env.idx.C_Q5_ind,:));
    pol1_03_stds.C_90_std(vv) = std(Pol1_results0_3(vv).(var_series)(end - env.grid.oc + env.idx.C_90_ind,:));

    pol1_03_stds.W_10_std(vv) = std(Pol1_results0_3(vv).(var_series)(end - env.grid.oc + env.idx.w_10_ind,:));
    pol1_03_stds.W_Q1_std(vv) = std(Pol1_results0_3(vv).(var_series)(end - env.grid.oc + env.idx.w_Q1_ind,:));
    pol1_03_stds.W_Q2_std(vv) = std(Pol1_results0_3(vv).(var_series)(end - env.grid.oc + env.idx.w_Q2_ind,:));
    pol1_03_stds.W_Q3_std(vv) = std(Pol1_results0_3(vv).(var_series)(end - env.grid.oc + env.idx.w_Q3_ind,:));
    pol1_03_stds.W_Q4_std(vv) = std(Pol1_results0_3(vv).(var_series)(end - env.grid.oc + env.idx.w_Q4_ind,:));
    pol1_03_stds.W_Q5_std(vv) = std(Pol1_results0_3(vv).(var_series)(end - env.grid.oc + env.idx.w_Q5_ind,:));
    pol1_03_stds.W_90_std(vv) = std(Pol1_results0_3(vv).(var_series)(end - env.grid.oc + env.idx.w_90_ind,:));



    pol1_03_stds.C_I_10_std(vv) = std(Pol1_results0_3(vv).(var_series)(end - env.grid.oc + env.idx.C_I_10_ind,:));
    pol1_03_stds.C_I_Q1_std(vv) = std(Pol1_results0_3(vv).(var_series)(end - env.grid.oc + env.idx.C_I_Q1_ind,:));
    pol1_03_stds.C_I_Q2_std(vv) = std(Pol1_results0_3(vv).(var_series)(end - env.grid.oc + env.idx.C_I_Q2_ind,:));
    pol1_03_stds.C_I_Q3_std(vv) = std(Pol1_results0_3(vv).(var_series)(end - env.grid.oc + env.idx.C_I_Q3_ind,:));
    pol1_03_stds.C_I_Q4_std(vv) = std(Pol1_results0_3(vv).(var_series)(end - env.grid.oc + env.idx.C_I_Q4_ind,:));
    pol1_03_stds.C_I_Q5_std(vv) = std(Pol1_results0_3(vv).(var_series)(end - env.grid.oc + env.idx.C_I_Q5_ind,:));
    pol1_03_stds.C_I_90_std(vv) = std(Pol1_results0_3(vv).(var_series)(end - env.grid.oc + env.idx.C_I_90_ind,:));
    pol1_03_stds.C_std(vv)         = std(Pol1_results0_3(vv).(var_series)(end - env.grid.oc + env.idx.C_ind,:) );
    pol1_03_stds.Y_std(vv)         = std(Pol1_results0_3(vv).(var_series)(end - env.grid.oc + env.idx.Y_ind,:) );
    pol1_03_stds.pi_std(vv)        = std(Pol1_results0_3(vv).(var_series)(end - env.grid.oc + env.idx.pi_ind,:) );
    pol1_03_stds.unemp_std(vv)     = std(Pol1_results0_3(vv).(var_series)(end - env.grid.oc + env.idx.unemp_ind,:));
    pol1_03_stds.RR_std(vv)        = std(Pol1_results0_3(vv).(var_series)(end - env.grid.oc + env.idx.RR_ind,:));
    pol1_03_stds.Profit_std(vv)    = std(Pol1_results0_3(vv).(var_series)(end - env.grid.oc + env.idx.Profit_ind,:));
    pol1_03_stds.Profit_FI_std(vv) = std(Pol1_results0_3(vv).(var_series)(end - env.grid.oc + env.idx.unemp_ind,:));

    %scalars:
    pol1_03_stds.Q_std(vv) = std(Pol1_results0_3(vv).(var_series)(env.grid.numstates - env.grid.os + env.idx.Q_ind,:));
    pol1_03_stds.R_cb_std(vv) = std(Pol1_results0_3(vv).(var_series)(env.grid.numstates - env.grid.os + env.idx.R_cb_ind,:));
    pol1_03_stds.W_std(vv) = std(Pol1_results0_3(vv).(var_series)(env.grid.numstates - env.grid.os + env.idx.w_ind,:));


    %Manual Standard Deviations:
    pol1_03_stds2.C_10_std(vv) = manual_std(Pol1_results0_3(vv).(var_series)(end - env.grid.oc + env.idx.C_10_ind,:));
    pol1_03_stds2.C_Q1_std(vv) = manual_std(Pol1_results0_3(vv).(var_series)(end - env.grid.oc + env.idx.C_Q1_ind,:));
    pol1_03_stds2.C_Q2_std(vv) = manual_std(Pol1_results0_3(vv).(var_series)(end - env.grid.oc + env.idx.C_Q2_ind,:));
    pol1_03_stds2.C_Q3_std(vv) = manual_std(Pol1_results0_3(vv).(var_series)(end - env.grid.oc + env.idx.C_Q3_ind,:));
    pol1_03_stds2.C_Q4_std(vv) = manual_std(Pol1_results0_3(vv).(var_series)(end - env.grid.oc + env.idx.C_Q4_ind,:));
    pol1_03_stds2.C_Q5_std(vv) = manual_std(Pol1_results0_3(vv).(var_series)(end - env.grid.oc + env.idx.C_Q5_ind,:));
    pol1_03_stds2.C_90_std(vv) = manual_std(Pol1_results0_3(vv).(var_series)(end - env.grid.oc + env.idx.C_90_ind,:));


    pol1_03_stds2.W_10_std(vv) = manual_std(Pol1_results0_3(vv).(var_series)(end - env.grid.oc + env.idx.w_10_ind,:));
    pol1_03_stds2.W_Q1_std(vv) = manual_std(Pol1_results0_3(vv).(var_series)(end - env.grid.oc + env.idx.w_Q1_ind,:));
    pol1_03_stds2.W_Q2_std(vv) = manual_std(Pol1_results0_3(vv).(var_series)(end - env.grid.oc + env.idx.w_Q2_ind,:));
    pol1_03_stds2.W_Q3_std(vv) = manual_std(Pol1_results0_3(vv).(var_series)(end - env.grid.oc + env.idx.w_Q3_ind,:));
    pol1_03_stds2.W_Q4_std(vv) = manual_std(Pol1_results0_3(vv).(var_series)(end - env.grid.oc + env.idx.w_Q4_ind,:));
    pol1_03_stds2.W_Q5_std(vv) = manual_std(Pol1_results0_3(vv).(var_series)(end - env.grid.oc + env.idx.w_Q5_ind,:));
    pol1_03_stds2.W_90_std(vv) = manual_std(Pol1_results0_3(vv).(var_series)(end - env.grid.oc + env.idx.w_90_ind,:));

    pol1_03_stds2.C_I_10_std(vv) = manual_std(Pol1_results0_3(vv).(var_series)(end - env.grid.oc + env.idx.C_I_10_ind,:));
    pol1_03_stds2.C_I_Q1_std(vv) = manual_std(Pol1_results0_3(vv).(var_series)(end - env.grid.oc + env.idx.C_I_Q1_ind,:));
    pol1_03_stds2.C_I_Q2_std(vv) = manual_std(Pol1_results0_3(vv).(var_series)(end - env.grid.oc + env.idx.C_I_Q2_ind,:));
    pol1_03_stds2.C_I_Q3_std(vv) = manual_std(Pol1_results0_3(vv).(var_series)(end - env.grid.oc + env.idx.C_I_Q3_ind,:));
    pol1_03_stds2.C_I_Q4_std(vv) = manual_std(Pol1_results0_3(vv).(var_series)(end - env.grid.oc + env.idx.C_I_Q4_ind,:));
    pol1_03_stds2.C_I_Q5_std(vv) = manual_std(Pol1_results0_3(vv).(var_series)(end - env.grid.oc + env.idx.C_I_Q5_ind,:));
    pol1_03_stds2.C_I_90_std(vv) = manual_std(Pol1_results0_3(vv).(var_series)(end - env.grid.oc + env.idx.C_I_90_ind,:));
    pol1_03_stds2.C_std(vv)         = manual_std(Pol1_results0_3(vv).(var_series)(end - env.grid.oc + env.idx.C_ind,:) );
    pol1_03_stds2.Y_std(vv)         = manual_std(Pol1_results0_3(vv).(var_series)(end - env.grid.oc + env.idx.Y_ind,:) );
    pol1_03_stds2.pi_std(vv)        = manual_std(Pol1_results0_3(vv).(var_series)(end - env.grid.oc + env.idx.pi_ind,:) );
    pol1_03_stds2.unemp_std(vv)     = manual_std(Pol1_results0_3(vv).(var_series)(end - env.grid.oc + env.idx.unemp_ind,:));
    pol1_03_stds2.RR_std(vv)        = manual_std(Pol1_results0_3(vv).(var_series)(end - env.grid.oc + env.idx.RR_ind,:));
    pol1_03_stds2.Profit_std(vv)    = manual_std(Pol1_results0_3(vv).(var_series)(end - env.grid.oc + env.idx.Profit_ind,:));
    pol1_03_stds2.Profit_FI_std(vv) = manual_std(Pol1_results0_3(vv).(var_series)(end - env.grid.oc + env.idx.unemp_ind,:));

    %scalars:
    pol1_03_stds2.Q_std(vv) = manual_std(Pol1_results0_3(vv).(var_series)(env.grid.numstates - env.grid.os + env.idx.Q_ind,:));
    pol1_03_stds2.R_cb_std(vv) = manual_std(Pol1_results0_3(vv).(var_series)(env.grid.numstates - env.grid.os + env.idx.R_cb_ind,:));
    pol1_03_stds2.W_std(vv) = manual_std(Pol1_results0_3(vv).(var_series)(env.grid.numstates - env.grid.os + env.idx.w_ind,:));

    pol1_03_means.R_a_mean(vv)     = mean(Pol1_results0_3(vv).(var_series)(end - env.grid.oc + env.idx.r_a_ind,:));
    pol1_03_stds.R_a_std(vv)     = std(Pol1_results0_3(vv).(var_series)(end - env.grid.oc + env.idx.r_a_ind,:));
    pol1_03_stds2.R_a_std(vv)     = manual_std(Pol1_results0_3(vv).(var_series)(end - env.grid.oc + env.idx.r_a_ind,:));



    pol1_03_means.RR_new_mean(vv)        = mean(Pol1_results0_3(vv).(var_series)(env.grid.numstates - env.grid.os + env.idx.R_cb_ind,:) - Pol1_results0_3(vv).(var_series)(end - env.grid.oc + env.idx.pi_ind,:) );
    pol1_03_stds.RR_new_std(vv)        = std(Pol1_results0_3(vv).(var_series)(env.grid.numstates - env.grid.os + env.idx.R_cb_ind,:) - Pol1_results0_3(vv).(var_series)(end - env.grid.oc + env.idx.pi_ind,:) );
    pol1_03_stds2.RR_new_std(vv)        = manual_std(Pol1_results0_3(vv).(var_series)(env.grid.numstates - env.grid.os + env.idx.R_cb_ind,:) - Pol1_results0_3(vv).(var_series)(end - env.grid.oc + env.idx.pi_ind,:) );





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%var_series = 'X_series';
 pol1_05_means.C_10_mean(vv) = mean(Pol1_results0_5(vv).(var_series)(end - env.grid.oc + env.idx.C_10_ind,:));
    pol1_05_means.C_Q1_mean(vv) = mean(Pol1_results0_5(vv).(var_series)(end - env.grid.oc + env.idx.C_Q1_ind,:));
    pol1_05_means.C_Q2_mean(vv) = mean(Pol1_results0_5(vv).(var_series)(end - env.grid.oc + env.idx.C_Q2_ind,:));
    pol1_05_means.C_Q3_mean(vv) = mean(Pol1_results0_5(vv).(var_series)(end - env.grid.oc + env.idx.C_Q3_ind,:));
    pol1_05_means.C_Q4_mean(vv) = mean(Pol1_results0_5(vv).(var_series)(end - env.grid.oc + env.idx.C_Q4_ind,:));
    pol1_05_means.C_Q5_mean(vv) = mean(Pol1_results0_5(vv).(var_series)(end - env.grid.oc + env.idx.C_Q5_ind,:));
    pol1_05_means.C_90_mean(vv) = mean(Pol1_results0_5(vv).(var_series)(end - env.grid.oc + env.idx.C_90_ind,:));


    pol1_05_means.W_10_mean(vv) = mean(Pol1_results0_5(vv).(var_series)(end - env.grid.oc + env.idx.w_10_ind,:));
    pol1_05_means.W_Q1_mean(vv) = mean(Pol1_results0_5(vv).(var_series)(end - env.grid.oc + env.idx.w_Q1_ind,:));
    pol1_05_means.W_Q2_mean(vv) = mean(Pol1_results0_5(vv).(var_series)(end - env.grid.oc + env.idx.w_Q2_ind,:));
    pol1_05_means.W_Q3_mean(vv) = mean(Pol1_results0_5(vv).(var_series)(end - env.grid.oc + env.idx.w_Q3_ind,:));
    pol1_05_means.W_Q4_mean(vv) = mean(Pol1_results0_5(vv).(var_series)(end - env.grid.oc + env.idx.w_Q4_ind,:));
    pol1_05_means.W_Q5_mean(vv) = mean(Pol1_results0_5(vv).(var_series)(end - env.grid.oc + env.idx.w_Q5_ind,:));
    pol1_05_means.W_90_mean(vv) = mean(Pol1_results0_5(vv).(var_series)(end - env.grid.oc + env.idx.w_90_ind,:));

    pol1_05_means.C_I_10_mean(vv) = mean(Pol1_results0_5(vv).(var_series)(end - env.grid.oc + env.idx.C_I_10_ind,:));
    pol1_05_means.C_I_Q1_mean(vv) = mean(Pol1_results0_5(vv).(var_series)(end - env.grid.oc + env.idx.C_I_Q1_ind,:));
    pol1_05_means.C_I_Q2_mean(vv) = mean(Pol1_results0_5(vv).(var_series)(end - env.grid.oc + env.idx.C_I_Q2_ind,:));
    pol1_05_means.C_I_Q3_mean(vv) = mean(Pol1_results0_5(vv).(var_series)(end - env.grid.oc + env.idx.C_I_Q3_ind,:));
    pol1_05_means.C_I_Q4_mean(vv) = mean(Pol1_results0_5(vv).(var_series)(end - env.grid.oc + env.idx.C_I_Q4_ind,:));
    pol1_05_means.C_I_Q5_mean(vv) = mean(Pol1_results0_5(vv).(var_series)(end - env.grid.oc + env.idx.C_I_Q5_ind,:));
    pol1_05_means.C_I_90_mean(vv) = mean(Pol1_results0_5(vv).(var_series)(end - env.grid.oc + env.idx.C_I_90_ind,:));
    pol1_05_means.C_mean(vv)         = mean(Pol1_results0_5(vv).(var_series)(end - env.grid.oc + env.idx.C_ind,:) );
    pol1_05_means.Y_mean(vv)         = mean(Pol1_results0_5(vv).(var_series)(end - env.grid.oc + env.idx.Y_ind,:) );
    pol1_05_means.pi_mean(vv)        = mean(Pol1_results0_5(vv).(var_series)(end - env.grid.oc + env.idx.pi_ind,:) );
    pol1_05_means.unemp_mean(vv)     = mean(Pol1_results0_5(vv).(var_series)(end - env.grid.oc + env.idx.unemp_ind,:));
    pol1_05_means.RR_mean(vv)        = mean(Pol1_results0_5(vv).(var_series)(end - env.grid.oc + env.idx.RR_ind,:));
    pol1_05_means.Profit_mean(vv)    = mean(Pol1_results0_5(vv).(var_series)(end - env.grid.oc + env.idx.Profit_ind,:));
    pol1_05_means.Profit_FI_mean(vv) = mean(Pol1_results0_5(vv).(var_series)(end - env.grid.oc + env.idx.unemp_ind,:));

    %scalars:
    pol1_05_means.Q_mean(vv) = mean(Pol1_results0_5(vv).(var_series)(env.grid.numstates - env.grid.os + env.idx.Q_ind,:));
    pol1_05_means.R_cb_mean(vv) = mean(Pol1_results0_5(vv).(var_series)(env.grid.numstates - env.grid.os + env.idx.R_cb_ind,:));
    pol1_05_means.W_mean(vv) = mean(Pol1_results0_5(vv).(var_series)(env.grid.numstates - env.grid.os + env.idx.w_ind,:));




    %Standard Deviations
    pol1_05_stds.C_10_std(vv) = std(Pol1_results0_5(vv).(var_series)(end - env.grid.oc + env.idx.C_10_ind,:));
    pol1_05_stds.C_Q1_std(vv) = std(Pol1_results0_5(vv).(var_series)(end - env.grid.oc + env.idx.C_Q1_ind,:));
    pol1_05_stds.C_Q2_std(vv) = std(Pol1_results0_5(vv).(var_series)(end - env.grid.oc + env.idx.C_Q2_ind,:));
    pol1_05_stds.C_Q3_std(vv) = std(Pol1_results0_5(vv).(var_series)(end - env.grid.oc + env.idx.C_Q3_ind,:));
    pol1_05_stds.C_Q4_std(vv) = std(Pol1_results0_5(vv).(var_series)(end - env.grid.oc + env.idx.C_Q4_ind,:));
    pol1_05_stds.C_Q5_std(vv) = std(Pol1_results0_5(vv).(var_series)(end - env.grid.oc + env.idx.C_Q5_ind,:));
    pol1_05_stds.C_90_std(vv) = std(Pol1_results0_5(vv).(var_series)(end - env.grid.oc + env.idx.C_90_ind,:));

    pol1_05_stds.W_10_std(vv) = std(Pol1_results0_5(vv).(var_series)(end - env.grid.oc + env.idx.w_10_ind,:));
    pol1_05_stds.W_Q1_std(vv) = std(Pol1_results0_5(vv).(var_series)(end - env.grid.oc + env.idx.w_Q1_ind,:));
    pol1_05_stds.W_Q2_std(vv) = std(Pol1_results0_5(vv).(var_series)(end - env.grid.oc + env.idx.w_Q2_ind,:));
    pol1_05_stds.W_Q3_std(vv) = std(Pol1_results0_5(vv).(var_series)(end - env.grid.oc + env.idx.w_Q3_ind,:));
    pol1_05_stds.W_Q4_std(vv) = std(Pol1_results0_5(vv).(var_series)(end - env.grid.oc + env.idx.w_Q4_ind,:));
    pol1_05_stds.W_Q5_std(vv) = std(Pol1_results0_5(vv).(var_series)(end - env.grid.oc + env.idx.w_Q5_ind,:));
    pol1_05_stds.W_90_std(vv) = std(Pol1_results0_5(vv).(var_series)(end - env.grid.oc + env.idx.w_90_ind,:));



    pol1_05_stds.C_I_10_std(vv) = std(Pol1_results0_5(vv).(var_series)(end - env.grid.oc + env.idx.C_I_10_ind,:));
    pol1_05_stds.C_I_Q1_std(vv) = std(Pol1_results0_5(vv).(var_series)(end - env.grid.oc + env.idx.C_I_Q1_ind,:));
    pol1_05_stds.C_I_Q2_std(vv) = std(Pol1_results0_5(vv).(var_series)(end - env.grid.oc + env.idx.C_I_Q2_ind,:));
    pol1_05_stds.C_I_Q3_std(vv) = std(Pol1_results0_5(vv).(var_series)(end - env.grid.oc + env.idx.C_I_Q3_ind,:));
    pol1_05_stds.C_I_Q4_std(vv) = std(Pol1_results0_5(vv).(var_series)(end - env.grid.oc + env.idx.C_I_Q4_ind,:));
    pol1_05_stds.C_I_Q5_std(vv) = std(Pol1_results0_5(vv).(var_series)(end - env.grid.oc + env.idx.C_I_Q5_ind,:));
    pol1_05_stds.C_I_90_std(vv) = std(Pol1_results0_5(vv).(var_series)(end - env.grid.oc + env.idx.C_I_90_ind,:));
    pol1_05_stds.C_std(vv)         = std(Pol1_results0_5(vv).(var_series)(end - env.grid.oc + env.idx.C_ind,:) );
    pol1_05_stds.Y_std(vv)         = std(Pol1_results0_5(vv).(var_series)(end - env.grid.oc + env.idx.Y_ind,:) );
    pol1_05_stds.pi_std(vv)        = std(Pol1_results0_5(vv).(var_series)(end - env.grid.oc + env.idx.pi_ind,:) );
    pol1_05_stds.unemp_std(vv)     = std(Pol1_results0_5(vv).(var_series)(end - env.grid.oc + env.idx.unemp_ind,:));
    pol1_05_stds.RR_std(vv)        = std(Pol1_results0_5(vv).(var_series)(end - env.grid.oc + env.idx.RR_ind,:));
    pol1_05_stds.Profit_std(vv)    = std(Pol1_results0_5(vv).(var_series)(end - env.grid.oc + env.idx.Profit_ind,:));
    pol1_05_stds.Profit_FI_std(vv) = std(Pol1_results0_5(vv).(var_series)(end - env.grid.oc + env.idx.unemp_ind,:));

    %scalars:
    pol1_05_stds.Q_std(vv) = std(Pol1_results0_5(vv).(var_series)(env.grid.numstates - env.grid.os + env.idx.Q_ind,:));
    pol1_05_stds.R_cb_std(vv) = std(Pol1_results0_5(vv).(var_series)(env.grid.numstates - env.grid.os + env.idx.R_cb_ind,:));
    pol1_05_stds.W_std(vv) = std(Pol1_results0_5(vv).(var_series)(env.grid.numstates - env.grid.os + env.idx.w_ind,:));


    %Manual Standard Deviations:
    pol1_05_stds2.C_10_std(vv) = manual_std(Pol1_results0_5(vv).(var_series)(end - env.grid.oc + env.idx.C_10_ind,:));
    pol1_05_stds2.C_Q1_std(vv) = manual_std(Pol1_results0_5(vv).(var_series)(end - env.grid.oc + env.idx.C_Q1_ind,:));
    pol1_05_stds2.C_Q2_std(vv) = manual_std(Pol1_results0_5(vv).(var_series)(end - env.grid.oc + env.idx.C_Q2_ind,:));
    pol1_05_stds2.C_Q3_std(vv) = manual_std(Pol1_results0_5(vv).(var_series)(end - env.grid.oc + env.idx.C_Q3_ind,:));
    pol1_05_stds2.C_Q4_std(vv) = manual_std(Pol1_results0_5(vv).(var_series)(end - env.grid.oc + env.idx.C_Q4_ind,:));
    pol1_05_stds2.C_Q5_std(vv) = manual_std(Pol1_results0_5(vv).(var_series)(end - env.grid.oc + env.idx.C_Q5_ind,:));
    pol1_05_stds2.C_90_std(vv) = manual_std(Pol1_results0_5(vv).(var_series)(end - env.grid.oc + env.idx.C_90_ind,:));


    pol1_05_stds2.W_10_std(vv) = manual_std(Pol1_results0_5(vv).(var_series)(end - env.grid.oc + env.idx.w_10_ind,:));
    pol1_05_stds2.W_Q1_std(vv) = manual_std(Pol1_results0_5(vv).(var_series)(end - env.grid.oc + env.idx.w_Q1_ind,:));
    pol1_05_stds2.W_Q2_std(vv) = manual_std(Pol1_results0_5(vv).(var_series)(end - env.grid.oc + env.idx.w_Q2_ind,:));
    pol1_05_stds2.W_Q3_std(vv) = manual_std(Pol1_results0_5(vv).(var_series)(end - env.grid.oc + env.idx.w_Q3_ind,:));
    pol1_05_stds2.W_Q4_std(vv) = manual_std(Pol1_results0_5(vv).(var_series)(end - env.grid.oc + env.idx.w_Q4_ind,:));
    pol1_05_stds2.W_Q5_std(vv) = manual_std(Pol1_results0_5(vv).(var_series)(end - env.grid.oc + env.idx.w_Q5_ind,:));
    pol1_05_stds2.W_90_std(vv) = manual_std(Pol1_results0_5(vv).(var_series)(end - env.grid.oc + env.idx.w_90_ind,:));

    pol1_05_stds2.C_I_10_std(vv) = manual_std(Pol1_results0_5(vv).(var_series)(end - env.grid.oc + env.idx.C_I_10_ind,:));
    pol1_05_stds2.C_I_Q1_std(vv) = manual_std(Pol1_results0_5(vv).(var_series)(end - env.grid.oc + env.idx.C_I_Q1_ind,:));
    pol1_05_stds2.C_I_Q2_std(vv) = manual_std(Pol1_results0_5(vv).(var_series)(end - env.grid.oc + env.idx.C_I_Q2_ind,:));
    pol1_05_stds2.C_I_Q3_std(vv) = manual_std(Pol1_results0_5(vv).(var_series)(end - env.grid.oc + env.idx.C_I_Q3_ind,:));
    pol1_05_stds2.C_I_Q4_std(vv) = manual_std(Pol1_results0_5(vv).(var_series)(end - env.grid.oc + env.idx.C_I_Q4_ind,:));
    pol1_05_stds2.C_I_Q5_std(vv) = manual_std(Pol1_results0_5(vv).(var_series)(end - env.grid.oc + env.idx.C_I_Q5_ind,:));
    pol1_05_stds2.C_I_90_std(vv) = manual_std(Pol1_results0_5(vv).(var_series)(end - env.grid.oc + env.idx.C_I_90_ind,:));
    pol1_05_stds2.C_std(vv)         = manual_std(Pol1_results0_5(vv).(var_series)(end - env.grid.oc + env.idx.C_ind,:) );
    pol1_05_stds2.Y_std(vv)         = manual_std(Pol1_results0_5(vv).(var_series)(end - env.grid.oc + env.idx.Y_ind,:) );
    pol1_05_stds2.pi_std(vv)        = manual_std(Pol1_results0_5(vv).(var_series)(end - env.grid.oc + env.idx.pi_ind,:) );
    pol1_05_stds2.unemp_std(vv)     = manual_std(Pol1_results0_5(vv).(var_series)(end - env.grid.oc + env.idx.unemp_ind,:));
    pol1_05_stds2.RR_std(vv)        = manual_std(Pol1_results0_5(vv).(var_series)(end - env.grid.oc + env.idx.RR_ind,:));
    pol1_05_stds2.Profit_std(vv)    = manual_std(Pol1_results0_5(vv).(var_series)(end - env.grid.oc + env.idx.Profit_ind,:));
    pol1_05_stds2.Profit_FI_std(vv) = manual_std(Pol1_results0_5(vv).(var_series)(end - env.grid.oc + env.idx.unemp_ind,:));

    %scalars:
    pol1_05_stds2.Q_std(vv) = manual_std(Pol1_results0_5(vv).(var_series)(env.grid.numstates - env.grid.os + env.idx.Q_ind,:));
    pol1_05_stds2.R_cb_std(vv) = manual_std(Pol1_results0_5(vv).(var_series)(env.grid.numstates - env.grid.os + env.idx.R_cb_ind,:));
    pol1_05_stds2.W_std(vv) = manual_std(Pol1_results0_5(vv).(var_series)(env.grid.numstates - env.grid.os + env.idx.w_ind,:));

    pol1_05_means.R_a_mean(vv)     = mean(Pol1_results0_5(vv).(var_series)(end - env.grid.oc + env.idx.r_a_ind,:));
    pol1_05_stds.R_a_std(vv)     = std(Pol1_results0_5(vv).(var_series)(end - env.grid.oc + env.idx.r_a_ind,:));
    pol1_05_stds2.R_a_std(vv)     = manual_std(Pol1_results0_5(vv).(var_series)(end - env.grid.oc + env.idx.r_a_ind,:));



    pol1_05_means.RR_new_mean(vv)        = mean(Pol1_results0_5(vv).(var_series)(env.grid.numstates - env.grid.os + env.idx.R_cb_ind,:) - Pol1_results0_5(vv).(var_series)(end - env.grid.oc + env.idx.pi_ind,:) );
    pol1_05_stds.RR_new_std(vv)        = std(Pol1_results0_5(vv).(var_series)(env.grid.numstates - env.grid.os + env.idx.R_cb_ind,:) - Pol1_results0_5(vv).(var_series)(end - env.grid.oc + env.idx.pi_ind,:) );
    pol1_05_stds2.RR_new_std(vv)        = manual_std(Pol1_results0_5(vv).(var_series)(env.grid.numstates - env.grid.os + env.idx.R_cb_ind,:) - Pol1_results0_5(vv).(var_series)(end - env.grid.oc + env.idx.pi_ind,:) );

















