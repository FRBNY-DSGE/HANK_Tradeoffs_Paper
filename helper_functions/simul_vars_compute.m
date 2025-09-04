% clear simul_mean_collet_base simul_std_collet_base simul_mean_collet_CF1 simul_std_collet_CF1
function [mean_simul_base,std_simul_base,std_simul_base_annual,std_simul_base_quarterly] = simul_vars_compute(IRFs_base,num_simul,SS_stats,param)

simul_mean_collet_base.Y         = zeros(num_simul,1);
simul_mean_collet_base.C         = zeros(num_simul,1);
simul_mean_collet_base.pi        = zeros(num_simul,1);
simul_mean_collet_base.R_cb      = zeros(num_simul,1);
simul_mean_collet_base.r_a       = zeros(num_simul,1);
simul_mean_collet_base.Q         = zeros(num_simul,1);
simul_mean_collet_base.Profit    = zeros(num_simul,1);
simul_mean_collet_base.Profit_FI = zeros(num_simul,1);
simul_mean_collet_base.unemp     = zeros(num_simul,1);
simul_mean_collet_base.W         = zeros(num_simul,1);
simul_mean_collet_base.I         = zeros(num_simul,1);
simul_mean_collet_base.GiniW     = zeros(num_simul,1);
simul_mean_collet_base.GiniI_1   = zeros(num_simul,1);
simul_mean_collet_base.GiniI_2   = zeros(num_simul,1);
simul_mean_collet_base.GiniC     = zeros(num_simul,1);

simul_mean_collet_base.R         = zeros(num_simul,1);
simul_mean_collet_base.RR        = zeros(num_simul,1);
simul_mean_collet_base.Ra        = zeros(num_simul,1);


simul_std_collet_base.Y         = zeros(num_simul,1);
simul_std_collet_base.C         = zeros(num_simul,1);
simul_std_collet_base.pi        = zeros(num_simul,1);
simul_std_collet_base.R_cb      = zeros(num_simul,1);
simul_std_collet_base.r_a       = zeros(num_simul,1);
simul_std_collet_base.Q         = zeros(num_simul,1);
simul_std_collet_base.Profit    = zeros(num_simul,1);
simul_std_collet_base.Profit_FI = zeros(num_simul,1);
simul_std_collet_base.unemp     = zeros(num_simul,1);
simul_std_collet_base.W         = zeros(num_simul,1);
simul_std_collet_base.I         = zeros(num_simul,1);
simul_std_collet_base.GiniW     = zeros(num_simul,1);
simul_std_collet_base.GiniI_1   = zeros(num_simul,1);
simul_std_collet_base.GiniI_2   = zeros(num_simul,1);
simul_std_collet_base.GiniC     = zeros(num_simul,1);

simul_std_collet_base.R         = zeros(num_simul,1);
simul_std_collet_base.RR        = zeros(num_simul,1);
simul_std_collet_base.Ra        = zeros(num_simul,1);

simul_mean_collet_base.C_1      = zeros(num_simul,1);
simul_mean_collet_base.C_10      = zeros(num_simul,1);
simul_mean_collet_base.C_Q1      = zeros(num_simul,1);
simul_mean_collet_base.C_Q2      = zeros(num_simul,1);
simul_mean_collet_base.C_Q3      = zeros(num_simul,1);
simul_mean_collet_base.C_Q4      = zeros(num_simul,1);
simul_mean_collet_base.C_Q5      = zeros(num_simul,1);
simul_mean_collet_base.C_90      = zeros(num_simul,1);
simul_mean_collet_base.C_99      = zeros(num_simul,1);

simul_mean_collet_base.C_I_1      = zeros(num_simul,1);
simul_mean_collet_base.C_I_10      = zeros(num_simul,1);
simul_mean_collet_base.C_I_Q1      = zeros(num_simul,1);
simul_mean_collet_base.C_I_Q2      = zeros(num_simul,1);
simul_mean_collet_base.C_I_Q3      = zeros(num_simul,1);
simul_mean_collet_base.C_I_Q4      = zeros(num_simul,1);
simul_mean_collet_base.C_I_Q5      = zeros(num_simul,1);
simul_mean_collet_base.C_I_90      = zeros(num_simul,1);
simul_mean_collet_base.C_I_99      = zeros(num_simul,1);

simul_mean_collet_base.I_1_1      = zeros(num_simul,1);
simul_mean_collet_base.I_1_10      = zeros(num_simul,1);
simul_mean_collet_base.I_1_Q1      = zeros(num_simul,1);
simul_mean_collet_base.I_1_Q2      = zeros(num_simul,1);
simul_mean_collet_base.I_1_Q3      = zeros(num_simul,1);
simul_mean_collet_base.I_1_Q4      = zeros(num_simul,1);
simul_mean_collet_base.I_1_Q5      = zeros(num_simul,1);
simul_mean_collet_base.I_1_90      = zeros(num_simul,1);
simul_mean_collet_base.I_1_99      = zeros(num_simul,1);

simul_mean_collet_base.I_I_1_1      = zeros(num_simul,1);
simul_mean_collet_base.I_I_1_10      = zeros(num_simul,1);
simul_mean_collet_base.I_I_1_Q1      = zeros(num_simul,1);
simul_mean_collet_base.I_I_1_Q2      = zeros(num_simul,1);
simul_mean_collet_base.I_I_1_Q3      = zeros(num_simul,1);
simul_mean_collet_base.I_I_1_Q4      = zeros(num_simul,1);
simul_mean_collet_base.I_I_1_Q5      = zeros(num_simul,1);
simul_mean_collet_base.I_I_1_90      = zeros(num_simul,1);
simul_mean_collet_base.I_I_1_99      = zeros(num_simul,1);

simul_std_collet_base.C_1      = zeros(num_simul,1);
simul_std_collet_base.C_10      = zeros(num_simul,1);
simul_std_collet_base.C_Q1      = zeros(num_simul,1);
simul_std_collet_base.C_Q2      = zeros(num_simul,1);
simul_std_collet_base.C_Q3      = zeros(num_simul,1);
simul_std_collet_base.C_Q4      = zeros(num_simul,1);
simul_std_collet_base.C_Q5      = zeros(num_simul,1);
simul_std_collet_base.C_90      = zeros(num_simul,1);
simul_std_collet_base.C_99      = zeros(num_simul,1);

simul_std_collet_base.C_I_1      = zeros(num_simul,1);
simul_std_collet_base.C_I_10      = zeros(num_simul,1);
simul_std_collet_base.C_I_Q1      = zeros(num_simul,1);
simul_std_collet_base.C_I_Q2      = zeros(num_simul,1);
simul_std_collet_base.C_I_Q3      = zeros(num_simul,1);
simul_std_collet_base.C_I_Q4      = zeros(num_simul,1);
simul_std_collet_base.C_I_Q5      = zeros(num_simul,1);
simul_std_collet_base.C_I_90      = zeros(num_simul,1);
simul_std_collet_base.C_I_99      = zeros(num_simul,1);

simul_std_collet_base.I_1_1      = zeros(num_simul,1);
simul_std_collet_base.I_1_10      = zeros(num_simul,1);
simul_std_collet_base.I_1_Q1      = zeros(num_simul,1);
simul_std_collet_base.I_1_Q2      = zeros(num_simul,1);
simul_std_collet_base.I_1_Q3      = zeros(num_simul,1);
simul_std_collet_base.I_1_Q4      = zeros(num_simul,1);
simul_std_collet_base.I_1_Q5      = zeros(num_simul,1);
simul_std_collet_base.I_1_90      = zeros(num_simul,1);
simul_std_collet_base.I_1_99      = zeros(num_simul,1);

simul_std_collet_base.I_I_1_1      = zeros(num_simul,1);
simul_std_collet_base.I_I_1_10      = zeros(num_simul,1);
simul_std_collet_base.I_I_1_Q1      = zeros(num_simul,1);
simul_std_collet_base.I_I_1_Q2      = zeros(num_simul,1);
simul_std_collet_base.I_I_1_Q3      = zeros(num_simul,1);
simul_std_collet_base.I_I_1_Q4      = zeros(num_simul,1);
simul_std_collet_base.I_I_1_Q5      = zeros(num_simul,1);
simul_std_collet_base.I_I_1_90      = zeros(num_simul,1);
simul_std_collet_base.I_I_1_99      = zeros(num_simul,1);

simul_mean_collet_base.K         = zeros(num_simul,1);
simul_std_collet_base.K          = zeros(num_simul,1);

simul_mean_collet_base.A_b       = zeros(num_simul,1);
simul_mean_collet_base.A_hh      = zeros(num_simul,1);
simul_std_collet_base.B_hh       = zeros(num_simul,1);

simul_std_collet_base.A_b       = zeros(num_simul,1);
simul_std_collet_base.A_hh      = zeros(num_simul,1);
simul_std_collet_base.B_hh      = zeros(num_simul,1);

for jj = 1:num_simul

	simul_mean_collet_base.Y(jj,1)         = mean(IRFs_base(jj).IRFs_Y_p);
	simul_mean_collet_base.C(jj,1)         = mean(IRFs_base(jj).IRFs_C_p);
	simul_mean_collet_base.pi(jj,1)        = mean(IRFs_base(jj).IRFs_pi_p);
	simul_mean_collet_base.R_cb(jj,1)      = mean(IRFs_base(jj).IRFs_R_cb_p);
	simul_mean_collet_base.r_a(jj,1)       = mean(IRFs_base(jj).IRFs_r_a_p);
	simul_mean_collet_base.Q(jj,1)         = mean(IRFs_base(jj).IRFs_Q_p);
	simul_mean_collet_base.Profit(jj,1)    = mean(IRFs_base(jj).IRFs_Profit_p);
	simul_mean_collet_base.Profit_FI(jj,1) = mean(IRFs_base(jj).IRFs_Profit_FI_p);
	simul_mean_collet_base.unemp(jj,1)     = mean(IRFs_base(jj).IRFs_unemp_p);
	simul_mean_collet_base.W(jj,1)         = mean(IRFs_base(jj).IRFs_w_p);
    simul_mean_collet_base.I(jj,1)         = mean(IRFs_base(jj).IRFs_I_p);
	simul_mean_collet_base.GiniW(jj,1)     = mean(IRFs_base(jj).IRFs_GiniW_p);
	simul_mean_collet_base.GiniI_1(jj,1)   = mean(IRFs_base(jj).IRFs_GiniI_1_p);
	simul_mean_collet_base.GiniI_2(jj,1)   = mean(IRFs_base(jj).IRFs_GiniI_2_p);
	simul_mean_collet_base.GiniC(jj,1)     = mean(IRFs_base(jj).IRFs_GiniC_p);

	simul_mean_collet_base.R(jj,1)         = mean(IRFs_base(jj).IRFs_RR_p);
	simul_mean_collet_base.RR(jj,1)        = mean(IRFs_base(jj).IRFs_R_cb_p-IRFs_base(jj).IRFs_pi_p);
	simul_mean_collet_base.Ra(jj,1)        = mean(IRFs_base(jj).IRFs_RRa_p);

	simul_mean_collet_base.K(jj,1)       = mean(IRFs_base(jj).IRFs_K_p); 

	simul_mean_collet_base.C_1(jj,1)      = mean(IRFs_base(jj).IRFs_C_1_p);
    simul_mean_collet_base.C_10(jj,1)      = mean(IRFs_base(jj).IRFs_C_10_p);
	simul_mean_collet_base.C_Q1(jj,1)      = mean(IRFs_base(jj).IRFs_C_Q1_p);
    simul_mean_collet_base.C_Q2(jj,1)      = mean(IRFs_base(jj).IRFs_C_Q2_p);
    simul_mean_collet_base.C_Q3(jj,1)      = mean(IRFs_base(jj).IRFs_C_Q3_p);
    simul_mean_collet_base.C_Q4(jj,1)      = mean(IRFs_base(jj).IRFs_C_Q4_p);
    simul_mean_collet_base.C_Q5(jj,1)      = mean(IRFs_base(jj).IRFs_C_Q5_p);
	simul_mean_collet_base.C_90(jj,1)      = mean(IRFs_base(jj).IRFs_C_90_p);
    simul_mean_collet_base.C_99(jj,1)      = mean(IRFs_base(jj).IRFs_C_99_p);

    simul_mean_collet_base.C_I_1(jj,1)      = mean(IRFs_base(jj).IRFs_C_I_1_p);
    simul_mean_collet_base.C_I_10(jj,1)      = mean(IRFs_base(jj).IRFs_C_I_10_p);
	simul_mean_collet_base.C_I_Q1(jj,1)      = mean(IRFs_base(jj).IRFs_C_I_Q1_p);
    simul_mean_collet_base.C_I_Q2(jj,1)      = mean(IRFs_base(jj).IRFs_C_I_Q2_p);
    simul_mean_collet_base.C_I_Q3(jj,1)      = mean(IRFs_base(jj).IRFs_C_I_Q3_p);
    simul_mean_collet_base.C_I_Q4(jj,1)      = mean(IRFs_base(jj).IRFs_C_I_Q4_p);
    simul_mean_collet_base.C_I_Q5(jj,1)      = mean(IRFs_base(jj).IRFs_C_I_Q5_p);
	simul_mean_collet_base.C_I_90(jj,1)      = mean(IRFs_base(jj).IRFs_C_I_90_p);
    simul_mean_collet_base.C_I_99(jj,1)      = mean(IRFs_base(jj).IRFs_C_I_99_p);


    simul_mean_collet_base.I_1_1(jj,1)      = mean(IRFs_base(jj).IRFs_I_1_1_p);
	simul_mean_collet_base.I_1_10(jj,1)      = mean(IRFs_base(jj).IRFs_I_1_10_p);
	simul_mean_collet_base.I_1_Q1(jj,1)      = mean(IRFs_base(jj).IRFs_I_1_Q1_p);
	simul_mean_collet_base.I_1_Q2(jj,1)      = mean(IRFs_base(jj).IRFs_I_1_Q2_p);
	simul_mean_collet_base.I_1_Q3(jj,1)      = mean(IRFs_base(jj).IRFs_I_1_Q3_p);
	simul_mean_collet_base.I_1_Q4(jj,1)      = mean(IRFs_base(jj).IRFs_I_1_Q4_p);
	simul_mean_collet_base.I_1_Q5(jj,1)      = mean(IRFs_base(jj).IRFs_I_1_Q5_p);
	simul_mean_collet_base.I_1_90(jj,1)      = mean(IRFs_base(jj).IRFs_I_1_90_p);
    simul_mean_collet_base.I_1_99(jj,1)      = mean(IRFs_base(jj).IRFs_I_1_99_p);

    simul_mean_collet_base.I_I_1_1(jj,1)      = mean(IRFs_base(jj).IRFs_I_I_1_1_p);
	simul_mean_collet_base.I_I_1_10(jj,1)      = mean(IRFs_base(jj).IRFs_I_I_1_10_p);
	simul_mean_collet_base.I_I_1_Q1(jj,1)      = mean(IRFs_base(jj).IRFs_I_I_1_Q1_p);
	simul_mean_collet_base.I_I_1_Q2(jj,1)      = mean(IRFs_base(jj).IRFs_I_I_1_Q2_p);
	simul_mean_collet_base.I_I_I_1_Q3(jj,1)      = mean(IRFs_base(jj).IRFs_I_I_1_Q3_p);
	simul_mean_collet_base.I_I_1_Q4(jj,1)      = mean(IRFs_base(jj).IRFs_I_I_1_Q4_p);
	simul_mean_collet_base.I_I_1_Q5(jj,1)      = mean(IRFs_base(jj).IRFs_I_I_1_Q5_p);
	simul_mean_collet_base.I_I_1_90(jj,1)      = mean(IRFs_base(jj).IRFs_I_I_1_90_p);
    simul_mean_collet_base.I_I_1_99(jj,1)      = mean(IRFs_base(jj).IRFs_I_I_1_99_p);

	simul_mean_collet_base.A_b   = mean(IRFs_base(jj).IRFs_A_b_p);
	simul_mean_collet_base.A_hh  = mean(IRFs_base(jj).IRFs_A_hh_p);
	simul_mean_collet_base.B_hh   = mean(IRFs_base(jj).IRFs_B_hh_p);

	simul_std_collet_base.Y(jj,1)         = std(IRFs_base(jj).IRFs_Y_p);
	simul_std_collet_base.C(jj,1)         = std(IRFs_base(jj).IRFs_C_p);
	simul_std_collet_base.pi(jj,1)        = std(IRFs_base(jj).pi_IRFs_p);
	simul_std_collet_base.R_cb(jj,1)      = std(IRFs_base(jj).R_cb_IRFs_p);
	simul_std_collet_base.r_a(jj,1)       = std(IRFs_base(jj).r_a_IRFs_p);
	simul_std_collet_base.Q(jj,1)         = std(IRFs_base(jj).IRFs_Q_p);
	simul_std_collet_base.Profit(jj,1)    = std(IRFs_base(jj).IRFs_Profit_p);
	simul_std_collet_base.Profit_FI(jj,1) = std(IRFs_base(jj).IRFs_Profit_FI_p);
	simul_std_collet_base.unemp(jj,1)     = std(IRFs_base(jj).unemp_IRFs_p);
	simul_std_collet_base.I(jj,1)         = std(IRFs_base(jj).IRFs_I_p);
    simul_std_collet_base.W(jj,1)         = std(IRFs_base(jj).IRFs_w_p);
	simul_std_collet_base.GiniW(jj,1)     = std(IRFs_base(jj).GiniW_IRFs_p);
	simul_std_collet_base.GiniI_1(jj,1)   = std(IRFs_base(jj).GiniI_1_IRFs_p);
	simul_std_collet_base.GiniI_2(jj,1)   = std(IRFs_base(jj).GiniI_2_IRFs_p);
	simul_std_collet_base.GiniC(jj,1)     = std(IRFs_base(jj).GiniC_IRFs_p);

	simul_std_collet_base.R(jj,1)         = std(IRFs_base(jj).RR_IRFs_p);
	simul_std_collet_base.RR(jj,1)        = std(IRFs_base(jj).R_cb_IRFs_p./IRFs_base(jj).pi_IRFs_p);
	simul_std_collet_base.Ra(jj,1)        = std(IRFs_base(jj).RRa_IRFs_p);

	simul_std_collet_base.K(jj,1)         = std(IRFs_base(jj).IRFs_K_p); 
	simul_std_collet_base.C_1(jj,1)      = std(IRFs_base(jj).IRFs_C_1_p);
    simul_std_collet_base.C_10(jj,1)      = std(IRFs_base(jj).IRFs_C_10_p);
	simul_std_collet_base.C_Q1(jj,1)      = std(IRFs_base(jj).IRFs_C_Q1_p);
    simul_std_collet_base.C_Q2(jj,1)      = std(IRFs_base(jj).IRFs_C_Q2_p);
    simul_std_collet_base.C_Q3(jj,1)      = std(IRFs_base(jj).IRFs_C_Q3_p);
    simul_std_collet_base.C_Q4(jj,1)      = std(IRFs_base(jj).IRFs_C_Q4_p);
    simul_std_collet_base.C_Q5(jj,1)      = std(IRFs_base(jj).IRFs_C_Q5_p);
	simul_std_collet_base.C_90(jj,1)      = std(IRFs_base(jj).IRFs_C_90_p);
	simul_std_collet_base.C_99(jj,1)      = std(IRFs_base(jj).IRFs_C_99_p);
	simul_std_collet_base.C_I_1(jj,1)      = std(IRFs_base(jj).IRFs_C_I_1_p);
    simul_std_collet_base.C_I_10(jj,1)      = std(IRFs_base(jj).IRFs_C_I_10_p);
	simul_std_collet_base.C_I_Q1(jj,1)      = std(IRFs_base(jj).IRFs_C_I_Q1_p);
    simul_std_collet_base.C_I_Q2(jj,1)      = std(IRFs_base(jj).IRFs_C_I_Q2_p);
    simul_std_collet_base.C_I_Q3(jj,1)      = std(IRFs_base(jj).IRFs_C_I_Q3_p);
    simul_std_collet_base.C_I_Q4(jj,1)      = std(IRFs_base(jj).IRFs_C_I_Q4_p);
    simul_std_collet_base.C_I_Q5(jj,1)      = std(IRFs_base(jj).IRFs_C_I_Q5_p);
	simul_std_collet_base.C_I_90(jj,1)      = std(IRFs_base(jj).IRFs_C_I_90_p);
	simul_std_collet_base.C_I_99(jj,1)      = std(IRFs_base(jj).IRFs_C_I_99_p);
	
    simul_std_collet_base.I_1_1(jj,1)    = std(IRFs_base(jj).IRFs_I_1_1_p);
    simul_std_collet_base.I_1_10(jj,1)    = std(IRFs_base(jj).IRFs_I_1_10_p);
	simul_std_collet_base.I_1_Q1(jj,1)    = std(IRFs_base(jj).IRFs_I_1_Q1_p);
	simul_std_collet_base.I_1_Q2(jj,1)    = std(IRFs_base(jj).IRFs_I_1_Q2_p);
	simul_std_collet_base.I_1_Q3(jj,1)    = std(IRFs_base(jj).IRFs_I_1_Q3_p);
	simul_std_collet_base.I_1_Q4(jj,1)    = std(IRFs_base(jj).IRFs_I_1_Q4_p);
	simul_std_collet_base.I_1_Q5(jj,1)    = std(IRFs_base(jj).IRFs_I_1_Q5_p);
	simul_std_collet_base.I_1_90(jj,1)    = std(IRFs_base(jj).IRFs_I_1_90_p);
	simul_std_collet_base.I_1_99(jj,1)    = std(IRFs_base(jj).IRFs_I_1_99_p);

    simul_std_collet_base.I_I_1_1(jj,1)    = std(IRFs_base(jj).IRFs_I_I_1_1_p);
    simul_std_collet_base.I_I_1_10(jj,1)    = std(IRFs_base(jj).IRFs_I_I_1_10_p);
	simul_std_collet_base.I_I_1_Q1(jj,1)    = std(IRFs_base(jj).IRFs_I_I_1_Q1_p);
	simul_std_collet_base.I_I_1_Q2(jj,1)    = std(IRFs_base(jj).IRFs_I_I_1_Q2_p);
	simul_std_collet_base.I_I_1_Q3(jj,1)    = std(IRFs_base(jj).IRFs_I_I_1_Q3_p);
	simul_std_collet_base.I_I_1_Q4(jj,1)    = std(IRFs_base(jj).IRFs_I_I_1_Q4_p);
	simul_std_collet_base.I_I_1_Q5(jj,1)    = std(IRFs_base(jj).IRFs_I_I_1_Q5_p);
	simul_std_collet_base.I_I_1_90(jj,1)    = std(IRFs_base(jj).IRFs_I_I_1_90_p);
	simul_std_collet_base.I_I_1_99(jj,1)    = std(IRFs_base(jj).IRFs_I_I_1_99_p);
	
    simul_std_collet_base.A_b             = std(IRFs_base(jj).IRFs_A_b_p);
	simul_std_collet_base.A_hh            = std(IRFs_base(jj).IRFs_A_hh_p);
	simul_std_collet_base.B_hh            = std(IRFs_base(jj).IRFs_B_hh_p);

	simul_std_collet_base_annual.Y(jj,1)         = std(IRFs_base(jj).Y_IRFs_p(5:end)./IRFs_base(jj).Y_IRFs_p(1:end-4)-1);
	simul_std_collet_base_annual.C(jj,1)         = std(IRFs_base(jj).C_IRFs_p(5:end)./IRFs_base(jj).C_IRFs_p(1:end-4)-1);
	simul_std_collet_base_annual.pi(jj,1)        = std(400*(IRFs_base(jj).pi_IRFs_p(5:end)-IRFs_base(jj).pi_IRFs_p(1:end-4)));
	simul_std_collet_base_annual.R_cb(jj,1)      = std(400*(IRFs_base(jj).R_cb_IRFs_p(5:end)-IRFs_base(jj).R_cb_IRFs_p(1:end-4)));
	simul_std_collet_base_annual.r_a(jj,1)       = std(400*(IRFs_base(jj).r_a_IRFs_p(5:end)-IRFs_base(jj).r_a_IRFs_p(1:end-4)));

	simul_std_collet_base_annual.Q(jj,1)         = std(IRFs_base(jj).Q_IRFs_p(5:end)./IRFs_base(jj).Q_IRFs_p(1:end-4)-1);
	simul_std_collet_base_annual.Profit(jj,1)    = std(IRFs_base(jj).Profit_IRFs_p(5:end)./IRFs_base(jj).Profit_IRFs_p(1:end-4)-1);
	simul_std_collet_base_annual.Profit_FI(jj,1) = std(IRFs_base(jj).Profit_FI_IRFs_p(5:end)./IRFs_base(jj).Profit_FI_IRFs_p(1:end-4)-1);
	simul_std_collet_base_annual.unemp(jj,1)     = std(IRFs_base(jj).unemp_IRFs_p(5:end)-IRFs_base(jj).unemp_IRFs_p(1:end-4));
	simul_std_collet_base_annual.W(jj,1)         = std(IRFs_base(jj).w_IRFs_p(5:end)./IRFs_base(jj).w_IRFs_p(1:end-4)-1);
    simul_std_collet_base_annual.I(jj,1)         = std(IRFs_base(jj).I_IRFs_p(5:end)./IRFs_base(jj).I_IRFs_p(1:end-4)-1);
	simul_std_collet_base_annual.GiniW(jj,1)     = std(IRFs_base(jj).GiniW_IRFs_p(5:end)-IRFs_base(jj).GiniW_IRFs_p(1:end-4));
	simul_std_collet_base_annual.GiniI_1(jj,1)   = std(IRFs_base(jj).GiniI_1_IRFs_p(5:end)-IRFs_base(jj).GiniI_1_IRFs_p(1:end-4));
	simul_std_collet_base_annual.GiniI_2(jj,1)   = std(IRFs_base(jj).GiniI_2_IRFs_p(5:end)-IRFs_base(jj).GiniI_2_IRFs_p(1:end-4));
	simul_std_collet_base_annual.GiniC(jj,1)     = std(IRFs_base(jj).GiniC_IRFs_p(5:end)-IRFs_base(jj).GiniC_IRFs_p(1:end-4));

	simul_std_collet_base_annual.R(jj,1)         = std(400*(IRFs_base(jj).RR_IRFs_p(5:end)-IRFs_base(jj).RR_IRFs_p(1:end-4)));
	simul_std_collet_base_annual.RR(jj,1)        = std(400*(IRFs_base(jj).R_cb_IRFs_p(5:end)./IRFs_base(jj).pi_IRFs_p(5:end)-IRFs_base(jj).R_cb_IRFs_p(1:end-4)./IRFs_base(jj).pi_IRFs_p(1:end-4)));
	simul_std_collet_base_annual.Ra(jj,1)       = std(400*(IRFs_base(jj).RRa_IRFs_p(5:end)-IRFs_base(jj).RRa_IRFs_p(1:end-4)));

	simul_std_collet_base_annual.K(jj,1)         = std(IRFs_base(jj).K_IRFs_p(5:end)./IRFs_base(jj).K_IRFs_p(1:end-4)-1);

    simul_std_collet_base_annual.C_1(jj,1)      = std(IRFs_base(jj).C_1_IRFs_p(5:end)./IRFs_base(jj).C_1_IRFs_p(1:end-4)-1);
	simul_std_collet_base_annual.C_10(jj,1)      = std(IRFs_base(jj).C_10_IRFs_p(5:end)./IRFs_base(jj).C_10_IRFs_p(1:end-4)-1);
	simul_std_collet_base_annual.C_Q1(jj,1)      = std(IRFs_base(jj).C_Q1_IRFs_p(5:end)./IRFs_base(jj).C_Q1_IRFs_p(1:end-4)-1);
	simul_std_collet_base_annual.C_Q2(jj,1)      = std(IRFs_base(jj).C_Q2_IRFs_p(5:end)./IRFs_base(jj).C_Q2_IRFs_p(1:end-4)-1);
	simul_std_collet_base_annual.C_Q3(jj,1)      = std(IRFs_base(jj).C_Q3_IRFs_p(5:end)./IRFs_base(jj).C_Q3_IRFs_p(1:end-4)-1);
	simul_std_collet_base_annual.C_Q4(jj,1)      = std(IRFs_base(jj).C_Q4_IRFs_p(5:end)./IRFs_base(jj).C_Q4_IRFs_p(1:end-4)-1);
	simul_std_collet_base_annual.C_Q5(jj,1)      = std(IRFs_base(jj).C_Q5_IRFs_p(5:end)./IRFs_base(jj).C_Q5_IRFs_p(1:end-4)-1);
	simul_std_collet_base_annual.C_90(jj,1)      = std(IRFs_base(jj).C_90_IRFs_p(5:end)./IRFs_base(jj).C_90_IRFs_p(1:end-4)-1);
    simul_std_collet_base_annual.C_99(jj,1)      = std(IRFs_base(jj).C_99_IRFs_p(5:end)./IRFs_base(jj).C_99_IRFs_p(1:end-4)-1);

    simul_std_collet_base_annual.C_I_1(jj,1)      = std(IRFs_base(jj).C_I_1_IRFs_p(5:end)./IRFs_base(jj).C_I_1_IRFs_p(1:end-4)-1);
	simul_std_collet_base_annual.C_I_10(jj,1)      = std(IRFs_base(jj).C_I_10_IRFs_p(5:end)./IRFs_base(jj).C_I_10_IRFs_p(1:end-4)-1);
	simul_std_collet_base_annual.C_I_Q1(jj,1)      = std(IRFs_base(jj).C_I_Q1_IRFs_p(5:end)./IRFs_base(jj).C_I_Q1_IRFs_p(1:end-4)-1);
	simul_std_collet_base_annual.C_I_Q2(jj,1)      = std(IRFs_base(jj).C_I_Q2_IRFs_p(5:end)./IRFs_base(jj).C_I_Q2_IRFs_p(1:end-4)-1);
	simul_std_collet_base_annual.C_I_Q3(jj,1)      = std(IRFs_base(jj).C_I_Q3_IRFs_p(5:end)./IRFs_base(jj).C_I_Q3_IRFs_p(1:end-4)-1);
	simul_std_collet_base_annual.C_I_Q4(jj,1)      = std(IRFs_base(jj).C_I_Q4_IRFs_p(5:end)./IRFs_base(jj).C_I_Q4_IRFs_p(1:end-4)-1);
	simul_std_collet_base_annual.C_I_Q5(jj,1)      = std(IRFs_base(jj).C_I_Q5_IRFs_p(5:end)./IRFs_base(jj).C_I_Q5_IRFs_p(1:end-4)-1);
	simul_std_collet_base_annual.C_I_90(jj,1)      = std(IRFs_base(jj).C_I_90_IRFs_p(5:end)./IRFs_base(jj).C_I_90_IRFs_p(1:end-4)-1);
    simul_std_collet_base_annual.C_I_99(jj,1)      = std(IRFs_base(jj).C_I_99_IRFs_p(5:end)./IRFs_base(jj).C_I_99_IRFs_p(1:end-4)-1);

    simul_std_collet_base_annual.I_1_1(jj,1)    = std(IRFs_base(jj).I_1_1_IRFs_p(5:end)./IRFs_base(jj).I_1_1_IRFs_p(1:end-4)-1);
	simul_std_collet_base_annual.I_1_10(jj,1)    = std(IRFs_base(jj).I_1_10_IRFs_p(5:end)./IRFs_base(jj).I_1_10_IRFs_p(1:end-4)-1);
	simul_std_collet_base_annual.I_1_Q1(jj,1)    = std(IRFs_base(jj).I_1_Q1_IRFs_p(5:end)./IRFs_base(jj).I_1_Q1_IRFs_p(1:end-4)-1);
	simul_std_collet_base_annual.I_1_Q2(jj,1)    = std(IRFs_base(jj).I_1_Q2_IRFs_p(5:end)./IRFs_base(jj).I_1_Q2_IRFs_p(1:end-4)-1);
	simul_std_collet_base_annual.I_1_Q3(jj,1)    = std(IRFs_base(jj).I_1_Q3_IRFs_p(5:end)./IRFs_base(jj).I_1_Q3_IRFs_p(1:end-4)-1);
	simul_std_collet_base_annual.I_1_Q4(jj,1)    = std(IRFs_base(jj).I_1_Q4_IRFs_p(5:end)./IRFs_base(jj).I_1_Q4_IRFs_p(1:end-4)-1);
	simul_std_collet_base_annual.I_1_Q5(jj,1)    = std(IRFs_base(jj).I_1_Q5_IRFs_p(5:end)./IRFs_base(jj).I_1_Q5_IRFs_p(1:end-4)-1);
	simul_std_collet_base_annual.I_1_90(jj,1)    = std(IRFs_base(jj).I_1_90_IRFs_p(5:end)./IRFs_base(jj).I_1_90_IRFs_p(1:end-4)-1);
    simul_std_collet_base_annual.I_1_99(jj,1)    = std(IRFs_base(jj).I_1_99_IRFs_p(5:end)./IRFs_base(jj).I_1_99_IRFs_p(1:end-4)-1);

    simul_std_collet_base_annual.I_I_1_1(jj,1)    = std(IRFs_base(jj).I_I_1_1_IRFs_p(5:end)./IRFs_base(jj).I_I_1_1_IRFs_p(1:end-4)-1);
	simul_std_collet_base_annual.I_I_1_10(jj,1)    = std(IRFs_base(jj).I_I_1_10_IRFs_p(5:end)./IRFs_base(jj).I_I_1_10_IRFs_p(1:end-4)-1);
	simul_std_collet_base_annual.I_I_1_Q1(jj,1)    = std(IRFs_base(jj).I_I_1_Q1_IRFs_p(5:end)./IRFs_base(jj).I_I_1_Q1_IRFs_p(1:end-4)-1);
	simul_std_collet_base_annual.I_I_1_Q2(jj,1)    = std(IRFs_base(jj).I_I_1_Q2_IRFs_p(5:end)./IRFs_base(jj).I_I_1_Q2_IRFs_p(1:end-4)-1);
	simul_std_collet_base_annual.I_I_1_Q3(jj,1)    = std(IRFs_base(jj).I_I_1_Q3_IRFs_p(5:end)./IRFs_base(jj).I_I_1_Q3_IRFs_p(1:end-4)-1);
	simul_std_collet_base_annual.I_I_1_Q4(jj,1)    = std(IRFs_base(jj).I_I_1_Q4_IRFs_p(5:end)./IRFs_base(jj).I_I_1_Q4_IRFs_p(1:end-4)-1);
	simul_std_collet_base_annual.I_I_1_Q5(jj,1)    = std(IRFs_base(jj).I_I_1_Q5_IRFs_p(5:end)./IRFs_base(jj).I_I_1_Q5_IRFs_p(1:end-4)-1);
	simul_std_collet_base_annual.I_I_1_90(jj,1)    = std(IRFs_base(jj).I_I_1_90_IRFs_p(5:end)./IRFs_base(jj).I_I_1_90_IRFs_p(1:end-4)-1);
    simul_std_collet_base_annual.I_I_1_99(jj,1)    = std(IRFs_base(jj).I_I_1_99_IRFs_p(5:end)./IRFs_base(jj).I_I_1_99_IRFs_p(1:end-4)-1);



	simul_std_collet_base_annual.A_b(jj,1)       = std(IRFs_base(jj).A_b_IRFs_p(5:end)./IRFs_base(jj).A_b_IRFs_p(1:end-4)-1);
	simul_std_collet_base_annual.A_hh(jj,1)      = std(IRFs_base(jj).A_hh_IRFs_p(5:end)./IRFs_base(jj).A_hh_IRFs_p(1:end-4)-1);
	simul_std_collet_base_annual.B_hh(jj,1)      = std(IRFs_base(jj).B_hh_IRFs_p(5:end)./IRFs_base(jj).B_hh_IRFs_p(1:end-4)-1);



	simul_std_collet_base_quarterly.Y(jj,1)         = std(IRFs_base(jj).Y_IRFs_p(2:end)./IRFs_base(jj).Y_IRFs_p(1:end-1)-1);
	simul_std_collet_base_quarterly.C(jj,1)         = std(IRFs_base(jj).C_IRFs_p(2:end)./IRFs_base(jj).C_IRFs_p(1:end-1)-1);
	simul_std_collet_base_quarterly.pi(jj,1)        = std(400*(IRFs_base(jj).pi_IRFs_p(2:end)-IRFs_base(jj).pi_IRFs_p(1:end-1)));
	simul_std_collet_base_quarterly.R_cb(jj,1)      = std(400*(IRFs_base(jj).R_cb_IRFs_p(2:end)-IRFs_base(jj).R_cb_IRFs_p(1:end-1)));
	simul_std_collet_base_quarterly.r_a(jj,1)       = std(400*(IRFs_base(jj).r_a_IRFs_p(2:end)-IRFs_base(jj).r_a_IRFs_p(1:end-1)));

	simul_std_collet_base_quarterly.Q(jj,1)         = std(IRFs_base(jj).Q_IRFs_p(2:end)./IRFs_base(jj).Q_IRFs_p(1:end-1)-1);
	simul_std_collet_base_quarterly.Profit(jj,1)    = std(IRFs_base(jj).Profit_IRFs_p(2:end)./IRFs_base(jj).Profit_IRFs_p(1:end-1)-1);
	simul_std_collet_base_quarterly.Profit_FI(jj,1) = std(IRFs_base(jj).Profit_FI_IRFs_p(2:end)./IRFs_base(jj).Profit_FI_IRFs_p(1:end-1)-1);
	simul_std_collet_base_quarterly.unemp(jj,1)     = std(IRFs_base(jj).unemp_IRFs_p(2:end)-IRFs_base(jj).unemp_IRFs_p(1:end-1));
	simul_std_collet_base_quarterly.W(jj,1)         = std(IRFs_base(jj).w_IRFs_p(2:end)./IRFs_base(jj).w_IRFs_p(1:end-1)-1);
    simul_std_collet_base_quarterly.I(jj,1)         = std(IRFs_base(jj).I_IRFs_p(2:end)./IRFs_base(jj).I_IRFs_p(1:end-1)-1);
	simul_std_collet_base_quarterly.GiniW(jj,1)     = std(IRFs_base(jj).GiniW_IRFs_p(2:end)-IRFs_base(jj).GiniW_IRFs_p(1:end-1));
	simul_std_collet_base_quarterly.GiniI_1(jj,1)   = std(IRFs_base(jj).GiniI_1_IRFs_p(2:end)-IRFs_base(jj).GiniI_1_IRFs_p(1:end-1));
	simul_std_collet_base_quarterly.GiniI_2(jj,1)   = std(IRFs_base(jj).GiniI_2_IRFs_p(2:end)-IRFs_base(jj).GiniI_2_IRFs_p(1:end-1));
	simul_std_collet_base_quarterly.GiniC(jj,1)     = std(IRFs_base(jj).GiniC_IRFs_p(2:end)-IRFs_base(jj).GiniC_IRFs_p(1:end-1));

	simul_std_collet_base_quarterly.R(jj,1)         = std(400*(IRFs_base(jj).RR_IRFs_p(2:end)-IRFs_base(jj).RR_IRFs_p(1:end-1)));
	simul_std_collet_base_quarterly.RR(jj,1)        = std(400*(IRFs_base(jj).R_cb_IRFs_p(2:end)./IRFs_base(jj).pi_IRFs_p(2:end)-IRFs_base(jj).R_cb_IRFs_p(1:end-1)./IRFs_base(jj).pi_IRFs_p(1:end-1)));
	simul_std_collet_base_quarterly.Ra(jj,1)        = std(400*(IRFs_base(jj).RRa_IRFs_p(2:end)-IRFs_base(jj).RRa_IRFs_p(1:end-1)));

	simul_std_collet_base_quarterly.K(jj,1)         = std(IRFs_base(jj).K_IRFs_p(2:end)./IRFs_base(jj).K_IRFs_p(1:end-1)-1);

	simul_std_collet_base_quarterly.C_1(jj,1)       = std(IRFs_base(jj).C_1_IRFs_p(2:end)./IRFs_base(jj).C_1_IRFs_p(1:end-1)-1);
    simul_std_collet_base_quarterly.C_10(jj,1)      = std(IRFs_base(jj).C_10_IRFs_p(2:end)./IRFs_base(jj).C_10_IRFs_p(1:end-1)-1);
	simul_std_collet_base_quarterly.C_Q1(jj,1)      = std(IRFs_base(jj).C_Q1_IRFs_p(2:end)./IRFs_base(jj).C_Q1_IRFs_p(1:end-1)-1);
	simul_std_collet_base_quarterly.C_Q2(jj,1)      = std(IRFs_base(jj).C_Q2_IRFs_p(2:end)./IRFs_base(jj).C_Q2_IRFs_p(1:end-1)-1);
	simul_std_collet_base_quarterly.C_Q3(jj,1)      = std(IRFs_base(jj).C_Q3_IRFs_p(2:end)./IRFs_base(jj).C_Q3_IRFs_p(1:end-1)-1);
	simul_std_collet_base_quarterly.C_Q4(jj,1)      = std(IRFs_base(jj).C_Q4_IRFs_p(2:end)./IRFs_base(jj).C_Q4_IRFs_p(1:end-1)-1);
	simul_std_collet_base_quarterly.C_Q5(jj,1)      = std(IRFs_base(jj).C_Q5_IRFs_p(2:end)./IRFs_base(jj).C_Q5_IRFs_p(1:end-1)-1);
	simul_std_collet_base_quarterly.C_90(jj,1)      = std(IRFs_base(jj).C_90_IRFs_p(2:end)./IRFs_base(jj).C_90_IRFs_p(1:end-1)-1);
    simul_std_collet_base_quarterly.C_99(jj,1)      = std(IRFs_base(jj).C_99_IRFs_p(2:end)./IRFs_base(jj).C_99_IRFs_p(1:end-1)-1);

    simul_std_collet_base_quarterly.C_I_1(jj,1)      = std(IRFs_base(jj).C_I_1_IRFs_p(2:end)./IRFs_base(jj).C_I_1_IRFs_p(1:end-1)-1);
    simul_std_collet_base_quarterly.C_I_10(jj,1)      = std(IRFs_base(jj).C_I_10_IRFs_p(2:end)./IRFs_base(jj).C_I_10_IRFs_p(1:end-1)-1);
	simul_std_collet_base_quarterly.C_I_Q1(jj,1)      = std(IRFs_base(jj).C_I_Q1_IRFs_p(2:end)./IRFs_base(jj).C_I_Q1_IRFs_p(1:end-1)-1);
	simul_std_collet_base_quarterly.C_I_Q2(jj,1)      = std(IRFs_base(jj).C_I_Q2_IRFs_p(2:end)./IRFs_base(jj).C_I_Q2_IRFs_p(1:end-1)-1);
	simul_std_collet_base_quarterly.C_I_Q3(jj,1)      = std(IRFs_base(jj).C_I_Q3_IRFs_p(2:end)./IRFs_base(jj).C_I_Q3_IRFs_p(1:end-1)-1);
	simul_std_collet_base_quarterly.C_I_Q4(jj,1)      = std(IRFs_base(jj).C_I_Q4_IRFs_p(2:end)./IRFs_base(jj).C_I_Q4_IRFs_p(1:end-1)-1);
	simul_std_collet_base_quarterly.C_I_Q5(jj,1)      = std(IRFs_base(jj).C_I_Q5_IRFs_p(2:end)./IRFs_base(jj).C_I_Q5_IRFs_p(1:end-1)-1);
	simul_std_collet_base_quarterly.C_I_90(jj,1)      = std(IRFs_base(jj).C_I_90_IRFs_p(2:end)./IRFs_base(jj).C_I_90_IRFs_p(1:end-1)-1);
    simul_std_collet_base_quarterly.C_I_99(jj,1)      = std(IRFs_base(jj).C_I_99_IRFs_p(2:end)./IRFs_base(jj).C_I_99_IRFs_p(1:end-1)-1);

    simul_std_collet_base_quarterly.I_1_1(jj,1)    = std(IRFs_base(jj).I_1_1_IRFs_p(2:end)./IRFs_base(jj).I_1_1_IRFs_p(1:end-1)-1);
	simul_std_collet_base_quarterly.I_1_10(jj,1)    = std(IRFs_base(jj).I_1_10_IRFs_p(2:end)./IRFs_base(jj).I_1_10_IRFs_p(1:end-1)-1);
	simul_std_collet_base_quarterly.I_1_Q1(jj,1)    = std(IRFs_base(jj).I_1_Q1_IRFs_p(2:end)./IRFs_base(jj).I_1_Q1_IRFs_p(1:end-1)-1);
	simul_std_collet_base_quarterly.I_1_Q2(jj,1)    = std(IRFs_base(jj).I_1_Q2_IRFs_p(2:end)./IRFs_base(jj).I_1_Q2_IRFs_p(1:end-1)-1);
	simul_std_collet_base_quarterly.I_1_Q3(jj,1)    = std(IRFs_base(jj).I_1_Q3_IRFs_p(2:end)./IRFs_base(jj).I_1_Q3_IRFs_p(1:end-1)-1);
	simul_std_collet_base_quarterly.I_1_Q4(jj,1)    = std(IRFs_base(jj).I_1_Q4_IRFs_p(2:end)./IRFs_base(jj).I_1_Q4_IRFs_p(1:end-1)-1);
	simul_std_collet_base_quarterly.I_1_Q5(jj,1)    = std(IRFs_base(jj).I_1_Q5_IRFs_p(2:end)./IRFs_base(jj).I_1_Q5_IRFs_p(1:end-1)-1);
	simul_std_collet_base_quarterly.I_1_90(jj,1)    = std(IRFs_base(jj).I_1_90_IRFs_p(2:end)./IRFs_base(jj).I_1_90_IRFs_p(1:end-1)-1);
    simul_std_collet_base_quarterly.I_1_99(jj,1)    = std(IRFs_base(jj).I_1_99_IRFs_p(2:end)./IRFs_base(jj).I_1_99_IRFs_p(1:end-1)-1);

    simul_std_collet_base_quarterly.I_I_1_1(jj,1)    = std(IRFs_base(jj).I_I_1_1_IRFs_p(2:end)./IRFs_base(jj).I_I_1_1_IRFs_p(1:end-1)-1);
	simul_std_collet_base_quarterly.I_I_1_10(jj,1)    = std(IRFs_base(jj).I_I_1_10_IRFs_p(2:end)./IRFs_base(jj).I_I_1_10_IRFs_p(1:end-1)-1);
	simul_std_collet_base_quarterly.I_I_1_Q1(jj,1)    = std(IRFs_base(jj).I_I_1_Q1_IRFs_p(2:end)./IRFs_base(jj).I_I_1_Q1_IRFs_p(1:end-1)-1);
	simul_std_collet_base_quarterly.I_I_1_Q2(jj,1)    = std(IRFs_base(jj).I_I_1_Q2_IRFs_p(2:end)./IRFs_base(jj).I_I_1_Q2_IRFs_p(1:end-1)-1);
	simul_std_collet_base_quarterly.I_I_1_Q3(jj,1)    = std(IRFs_base(jj).I_I_1_Q3_IRFs_p(2:end)./IRFs_base(jj).I_I_1_Q3_IRFs_p(1:end-1)-1);
	simul_std_collet_base_quarterly.I_I_1_Q4(jj,1)    = std(IRFs_base(jj).I_I_1_Q4_IRFs_p(2:end)./IRFs_base(jj).I_I_1_Q4_IRFs_p(1:end-1)-1);
	simul_std_collet_base_quarterly.I_I_1_Q5(jj,1)    = std(IRFs_base(jj).I_I_1_Q5_IRFs_p(2:end)./IRFs_base(jj).I_I_1_Q5_IRFs_p(1:end-1)-1);
	simul_std_collet_base_quarterly.I_I_1_90(jj,1)    = std(IRFs_base(jj).I_I_1_90_IRFs_p(2:end)./IRFs_base(jj).I_I_1_90_IRFs_p(1:end-1)-1);
    simul_std_collet_base_quarterly.I_I_1_99(jj,1)    = std(IRFs_base(jj).I_I_1_99_IRFs_p(2:end)./IRFs_base(jj).I_I_1_99_IRFs_p(1:end-1)-1);

	simul_std_collet_base_quarterly.A_b(jj,1)       = std(IRFs_base(jj).A_b_IRFs_p(2:end)./IRFs_base(jj).A_b_IRFs_p(1:end-1)-1);
	simul_std_collet_base_quarterly.A_hh(jj,1)      = std(IRFs_base(jj).A_hh_IRFs_p(2:end)./IRFs_base(jj).A_hh_IRFs_p(1:end-1)-1);
	simul_std_collet_base_quarterly.B_hh(jj,1)      = std(IRFs_base(jj).B_hh_IRFs_p(2:end)./IRFs_base(jj).B_hh_IRFs_p(1:end-1)-1);

end


mean_simul_base.Y         = exp(mean(simul_mean_collet_base.Y))*SS_stats.Y;
mean_simul_base.C         = exp(mean(simul_mean_collet_base.C))*SS_stats.C;
mean_simul_base.pi        = exp(mean(simul_mean_collet_base.pi))*param.pi_bar;
mean_simul_base.R_cb      = exp(mean(simul_mean_collet_base.R_cb))*param.R_cb;
mean_simul_base.R         = exp(mean(simul_mean_collet_base.R_cb))*param.R_cb/param.pi_bar;
mean_simul_base.Ra        = exp(mean(simul_mean_collet_base.R_cb))*SS_stats.RRa;
mean_simul_base.Q         = exp(mean(simul_mean_collet_base.Q));
mean_simul_base.Profit    = exp(mean(simul_mean_collet_base.Profit))*SS_stats.Profit;
mean_simul_base.Profit_FI = exp(mean(simul_mean_collet_base.Profit_FI))*SS_stats.Profit_FI;
mean_simul_base.unemp     = exp(mean(simul_mean_collet_base.unemp))*SS_stats.u;
mean_simul_base.W         = exp(mean(simul_mean_collet_base.W))*param.w_bar;
mean_simul_base.I         = exp(mean(simul_mean_collet_base.W))*SS_stats.I;
mean_simul_base.GiniW     = exp(mean(simul_mean_collet_base.GiniW))*SS_stats.GiniW;
mean_simul_base.GiniI_1   = exp(mean(simul_mean_collet_base.GiniI_1))*SS_stats.GiniIncome;
mean_simul_base.GiniI_2   = exp(mean(simul_mean_collet_base.GiniI_2))*SS_stats.GiniIncome;
mean_simul_base.GiniC     = exp(mean(simul_mean_collet_base.GiniC))*SS_stats.GiniCons;
mean_simul_base.K         = exp(mean(simul_mean_collet_base.K))*SS_stats.K;
mean_simul_base.C_1      = exp(mean(simul_mean_collet_base.C_1))*SS_stats.C_B1;
mean_simul_base.C_10      = exp(mean(simul_mean_collet_base.C_10))*SS_stats.C_B10;
mean_simul_base.C_Q1      = exp(mean(simul_mean_collet_base.C_Q1))*SS_stats.C_Q1;
mean_simul_base.C_Q2      = exp(mean(simul_mean_collet_base.C_Q2))*SS_stats.C_Q2;
mean_simul_base.C_Q3      = exp(mean(simul_mean_collet_base.C_Q3))*SS_stats.C_Q3;
mean_simul_base.C_Q4      = exp(mean(simul_mean_collet_base.C_Q4))*SS_stats.C_Q4;
mean_simul_base.C_Q5      = exp(mean(simul_mean_collet_base.C_Q5))*SS_stats.C_Q5;
mean_simul_base.C_90      = exp(mean(simul_mean_collet_base.C_90))*SS_stats.C_P90;
mean_simul_base.C_99      = exp(mean(simul_mean_collet_base.C_99))*SS_stats.C_P99;
mean_simul_base.C_I_1      = exp(mean(simul_mean_collet_base.C_I_1))*SS_stats.C_B1;
mean_simul_base.C_I_10      = exp(mean(simul_mean_collet_base.C_I_10))*SS_stats.C_I_B10;
mean_simul_base.C_I_Q1      = exp(mean(simul_mean_collet_base.C_I_Q1))*SS_stats.C_I_Q1;
mean_simul_base.C_I_Q2      = exp(mean(simul_mean_collet_base.C_I_Q2))*SS_stats.C_I_Q2;
mean_simul_base.C_I_Q3      = exp(mean(simul_mean_collet_base.C_I_Q3))*SS_stats.C_I_Q3;
mean_simul_base.C_I_Q4      = exp(mean(simul_mean_collet_base.C_I_Q4))*SS_stats.C_I_Q4;
mean_simul_base.C_I_Q5      = exp(mean(simul_mean_collet_base.C_I_Q5))*SS_stats.C_I_Q5;
mean_simul_base.C_I_90      = exp(mean(simul_mean_collet_base.C_I_90))*SS_stats.C_I_P90;
mean_simul_base.C_I_99      = exp(mean(simul_mean_collet_base.C_I_99))*SS_stats.C_I_P99;

mean_simul_base.I_1_1    = exp(mean(simul_mean_collet_base.I_1_1))*SS_stats.income_B1;
mean_simul_base.I_1_10    = exp(mean(simul_mean_collet_base.I_1_10))*SS_stats.income_B10;
mean_simul_base.I_1_Q1    = exp(mean(simul_mean_collet_base.I_1_Q1))*SS_stats.income_Q1;
mean_simul_base.I_1_Q2    = exp(mean(simul_mean_collet_base.I_1_Q2))*SS_stats.income_Q2;
mean_simul_base.I_1_Q3    = exp(mean(simul_mean_collet_base.I_1_Q3))*SS_stats.income_Q3;
mean_simul_base.I_1_Q4    = exp(mean(simul_mean_collet_base.I_1_Q4))*SS_stats.income_Q4;
mean_simul_base.I_1_Q5    = exp(mean(simul_mean_collet_base.I_1_Q5))*SS_stats.income_Q5;
mean_simul_base.I_1_90    = exp(mean(simul_mean_collet_base.I_1_90))*SS_stats.income_P90;
mean_simul_base.I_1_99    = exp(mean(simul_mean_collet_base.I_1_99))*SS_stats.income_P99;

mean_simul_base.I_I_1_1    = exp(mean(simul_mean_collet_base.I_I_1_1))*SS_stats.income_B1;
mean_simul_base.I_I_1_10    = exp(mean(simul_mean_collet_base.I_I_1_10))*SS_stats.income_B10;
mean_simul_base.I_I_1_Q1    = exp(mean(simul_mean_collet_base.I_I_1_Q1))*SS_stats.income_Q1;
mean_simul_base.I_I_1_Q2    = exp(mean(simul_mean_collet_base.I_I_1_Q2))*SS_stats.income_Q2;
mean_simul_base.I_I_1_Q3    = exp(mean(simul_mean_collet_base.I_I_1_Q3))*SS_stats.income_Q3;
mean_simul_base.I_I_1_Q4    = exp(mean(simul_mean_collet_base.I_I_1_Q4))*SS_stats.income_Q4;
mean_simul_base.I_I_1_Q5    = exp(mean(simul_mean_collet_base.I_I_1_Q5))*SS_stats.income_Q5;
mean_simul_base.I_I_1_90    = exp(mean(simul_mean_collet_base.I_I_1_90))*SS_stats.income_P90;
mean_simul_base.I_I_1_99    = exp(mean(simul_mean_collet_base.I_I_1_99))*SS_stats.income_P99;


mean_simul_base.A_b       = exp(mean(simul_mean_collet_base.A_b))*SS_stats.A_b;
mean_simul_base.A_hh      = exp(mean(simul_mean_collet_base.A_hh))*SS_stats.A_hh;
mean_simul_base.B_hh      = exp(mean(simul_mean_collet_base.B_hh))*SS_stats.B_hh;

std_simul_base.Y          = mean(simul_std_collet_base.Y);
std_simul_base.C          = mean(simul_std_collet_base.C);
std_simul_base.pi         = mean(simul_std_collet_base.pi);
std_simul_base.R_cb       = mean(simul_std_collet_base.R_cb);
std_simul_base.Q          = mean(simul_std_collet_base.Q);
std_simul_base.Profit     = mean(simul_std_collet_base.Profit);
std_simul_base.Profit_FI  = mean(simul_std_collet_base.Profit_FI);
std_simul_base.unemp      = mean(simul_std_collet_base.unemp);
std_simul_base.W          = mean(simul_std_collet_base.W);
std_simul_base.I          = mean(simul_std_collet_base.I);
std_simul_base.GiniW      = mean(simul_std_collet_base.GiniW);
std_simul_base.GiniI_1    = mean(simul_std_collet_base.GiniI_1);
std_simul_base.GiniI_2    = mean(simul_std_collet_base.GiniI_2);
std_simul_base.GiniC      = mean(simul_std_collet_base.GiniC);
std_simul_base.R          = mean(simul_std_collet_base.R);
std_simul_base.RR         = mean(simul_std_collet_base.RR);
std_simul_base.Ra         = mean(simul_std_collet_base.Ra);
std_simul_base.K          = mean(simul_std_collet_base.K);
std_simul_base.C_1       = mean(simul_std_collet_base.C_1);
std_simul_base.C_10       = mean(simul_std_collet_base.C_10);
std_simul_base.C_Q1       = mean(simul_std_collet_base.C_Q1);
std_simul_base.C_Q2       = mean(simul_std_collet_base.C_Q2);
std_simul_base.C_Q3       = mean(simul_std_collet_base.C_Q3);
std_simul_base.C_Q4       = mean(simul_std_collet_base.C_Q4);
std_simul_base.C_Q5       = mean(simul_std_collet_base.C_Q5);
std_simul_base.C_90       = mean(simul_std_collet_base.C_90);
std_simul_base.C_99       = mean(simul_std_collet_base.C_99);
std_simul_base.C_I_1       = mean(simul_std_collet_base.C_I_1);
std_simul_base.C_I_10       = mean(simul_std_collet_base.C_I_10);
std_simul_base.C_I_Q1       = mean(simul_std_collet_base.C_I_Q1);
std_simul_base.C_I_Q2       = mean(simul_std_collet_base.C_I_Q2);
std_simul_base.C_I_Q3       = mean(simul_std_collet_base.C_I_Q3);
std_simul_base.C_I_Q4       = mean(simul_std_collet_base.C_I_Q4);
std_simul_base.C_I_Q5       = mean(simul_std_collet_base.C_I_Q5);
std_simul_base.C_I_90       = mean(simul_std_collet_base.C_I_90);
std_simul_base.C_I_99       = mean(simul_std_collet_base.C_I_99);

std_simul_base.I_1_1     = mean(simul_std_collet_base.I_1_1);
std_simul_base.I_1_10     = mean(simul_std_collet_base.I_1_10);
std_simul_base.I_1_Q1     = mean(simul_std_collet_base.I_1_Q1);
std_simul_base.I_1_Q2     = mean(simul_std_collet_base.I_1_Q2);
std_simul_base.I_1_Q3     = mean(simul_std_collet_base.I_1_Q3);
std_simul_base.I_1_Q4     = mean(simul_std_collet_base.I_1_Q4);
std_simul_base.I_1_Q5     = mean(simul_std_collet_base.I_1_Q5);
std_simul_base.I_1_90     = mean(simul_std_collet_base.I_1_90);
std_simul_base.I_1_99     = mean(simul_std_collet_base.I_1_99);

std_simul_base.I_I_1_1     = mean(simul_std_collet_base.I_I_1_1);
std_simul_base.I_I_1_10     = mean(simul_std_collet_base.I_I_1_10);
std_simul_base.I_I_1_Q1     = mean(simul_std_collet_base.I_I_1_Q1);
std_simul_base.I_I_1_Q2     = mean(simul_std_collet_base.I_I_1_Q2);
std_simul_base.I_I_1_Q3     = mean(simul_std_collet_base.I_I_1_Q3);
std_simul_base.I_I_1_Q4     = mean(simul_std_collet_base.I_I_1_Q4);
std_simul_base.I_I_1_Q5     = mean(simul_std_collet_base.I_I_1_Q5);
std_simul_base.I_I_1_90     = mean(simul_std_collet_base.I_I_1_90);
std_simul_base.I_I_1_99     = mean(simul_std_collet_base.I_I_1_99);


std_simul_base.A_b        = mean(simul_std_collet_base.A_b);
std_simul_base.A_hh       = mean(simul_std_collet_base.A_hh);
std_simul_base.B_hh       = mean(simul_std_collet_base.B_hh);

std_simul_base_annual.Y          = mean(simul_std_collet_base_annual.Y);
std_simul_base_annual.C          = mean(simul_std_collet_base_annual.C);
std_simul_base_annual.pi         = mean(simul_std_collet_base_annual.pi);
std_simul_base_annual.R_cb       = mean(simul_std_collet_base_annual.R_cb);
std_simul_base_annual.Q          = mean(simul_std_collet_base_annual.Q);
std_simul_base_annual.Profit     = mean(simul_std_collet_base_annual.Profit);
std_simul_base_annual.Profit_FI  = mean(simul_std_collet_base_annual.Profit_FI);
std_simul_base_annual.unemp      = mean(simul_std_collet_base_annual.unemp);
std_simul_base_annual.W          = mean(simul_std_collet_base_annual.W);
std_simul_base_annual.I          = mean(simul_std_collet_base_annual.I);
std_simul_base_annual.GiniW      = mean(simul_std_collet_base_annual.GiniW);
std_simul_base_annual.GiniI_1    = mean(simul_std_collet_base_annual.GiniI_1);
std_simul_base_annual.GiniI_2    = mean(simul_std_collet_base_annual.GiniI_2);
std_simul_base_annual.GiniC      = mean(simul_std_collet_base_annual.GiniC);
std_simul_base_annual.R          = mean(simul_std_collet_base_annual.R);
std_simul_base_annual.RR         = mean(simul_std_collet_base_annual.RR);
std_simul_base_annual.Ra         = mean(simul_std_collet_base_annual.Ra);
std_simul_base_annual.K          = mean(simul_std_collet_base_annual.K);
std_simul_base_annual.C_1       = mean(simul_std_collet_base_annual.C_1);
std_simul_base_annual.C_10       = mean(simul_std_collet_base_annual.C_10);
std_simul_base_annual.C_Q1       = mean(simul_std_collet_base_annual.C_Q1);
std_simul_base_annual.C_Q2       = mean(simul_std_collet_base_annual.C_Q2);
std_simul_base_annual.C_Q3       = mean(simul_std_collet_base_annual.C_Q3);
std_simul_base_annual.C_Q4       = mean(simul_std_collet_base_annual.C_Q4);
std_simul_base_annual.C_Q5       = mean(simul_std_collet_base_annual.C_Q5);
std_simul_base_annual.C_90       = mean(simul_std_collet_base_annual.C_90);
std_simul_base_annual.C_99       = mean(simul_std_collet_base_annual.C_99);
std_simul_base_annual.C_I_1       = mean(simul_std_collet_base_annual.C_I_1);
std_simul_base_annual.C_I_10       = mean(simul_std_collet_base_annual.C_I_10);
std_simul_base_annual.C_I_Q1       = mean(simul_std_collet_base_annual.C_I_Q1);
std_simul_base_annual.C_I_Q2       = mean(simul_std_collet_base_annual.C_I_Q2);
std_simul_base_annual.C_I_Q3       = mean(simul_std_collet_base_annual.C_I_Q3);
std_simul_base_annual.C_I_Q4       = mean(simul_std_collet_base_annual.C_I_Q4);
std_simul_base_annual.C_I_Q5       = mean(simul_std_collet_base_annual.C_I_Q5);
std_simul_base_annual.C_I_90       = mean(simul_std_collet_base_annual.C_I_90);
std_simul_base_annual.C_I_99       = mean(simul_std_collet_base_annual.C_I_99);

std_simul_base_annual.I_1_1     = mean(simul_std_collet_base_annual.I_1_1);
std_simul_base_annual.I_1_10     = mean(simul_std_collet_base_annual.I_1_10);
std_simul_base_annual.I_1_Q1     = mean(simul_std_collet_base_annual.I_1_Q1);
std_simul_base_annual.I_1_Q2     = mean(simul_std_collet_base_annual.I_1_Q2);
std_simul_base_annual.I_1_Q3     = mean(simul_std_collet_base_annual.I_1_Q3);
std_simul_base_annual.I_1_Q4     = mean(simul_std_collet_base_annual.I_1_Q4);
std_simul_base_annual.I_1_Q5     = mean(simul_std_collet_base_annual.I_1_Q5);
std_simul_base_annual.I_1_90     = mean(simul_std_collet_base_annual.I_1_90);
std_simul_base_annual.I_1_99     = mean(simul_std_collet_base_annual.I_1_99);

std_simul_base_annual.I_I_1_1     = mean(simul_std_collet_base_annual.I_I_1_1);
std_simul_base_annual.I_I_1_10     = mean(simul_std_collet_base_annual.I_I_1_10);
std_simul_base_annual.I_I_1_Q1     = mean(simul_std_collet_base_annual.I_I_1_Q1);
std_simul_base_annual.I_I_1_Q2     = mean(simul_std_collet_base_annual.I_I_1_Q2);
std_simul_base_annual.I_I_1_Q3     = mean(simul_std_collet_base_annual.I_I_1_Q3);
std_simul_base_annual.I_I_1_Q4     = mean(simul_std_collet_base_annual.I_I_1_Q4);
std_simul_base_annual.I_I_1_Q5     = mean(simul_std_collet_base_annual.I_I_1_Q5);
std_simul_base_annual.I_I_1_90     = mean(simul_std_collet_base_annual.I_I_1_90);
std_simul_base_annual.I_I_1_99     = mean(simul_std_collet_base_annual.I_I_1_99);

std_simul_base_annual.A_b        = mean(simul_std_collet_base_annual.A_b);
std_simul_base_annual.A_hh       = mean(simul_std_collet_base_annual.A_hh);
std_simul_base_annual.B_hh       = mean(simul_std_collet_base_annual.B_hh);

std_simul_base_quarterly.Y          = mean(simul_std_collet_base_quarterly.Y);
std_simul_base_quarterly.C          = mean(simul_std_collet_base_quarterly.C);
std_simul_base_quarterly.pi         = mean(simul_std_collet_base_quarterly.pi);
std_simul_base_quarterly.R_cb       = mean(simul_std_collet_base_quarterly.R_cb);
std_simul_base_quarterly.Q          = mean(simul_std_collet_base_quarterly.Q);
std_simul_base_quarterly.Profit     = mean(simul_std_collet_base_quarterly.Profit);
std_simul_base_quarterly.Profit_FI  = mean(simul_std_collet_base_quarterly.Profit_FI);
std_simul_base_quarterly.unemp      = mean(simul_std_collet_base_quarterly.unemp);
std_simul_base_quarterly.W          = mean(simul_std_collet_base_quarterly.W);
std_simul_base_quarterly.I          = mean(simul_std_collet_base_quarterly.I);
std_simul_base_quarterly.GiniW      = mean(simul_std_collet_base_quarterly.GiniW);
std_simul_base_quarterly.GiniI_1    = mean(simul_std_collet_base_quarterly.GiniI_1);
std_simul_base_quarterly.GiniI_2    = mean(simul_std_collet_base_quarterly.GiniI_2);
std_simul_base_quarterly.GiniC      = mean(simul_std_collet_base_quarterly.GiniC);
std_simul_base_quarterly.R          = mean(simul_std_collet_base_quarterly.R);
std_simul_base_quarterly.RR         = mean(simul_std_collet_base_quarterly.RR);
std_simul_base_quarterly.Ra         = mean(simul_std_collet_base_quarterly.Ra);
std_simul_base_quarterly.K          = mean(simul_std_collet_base_quarterly.K);
std_simul_base_quarterly.C_1       = mean(simul_std_collet_base_quarterly.C_1);
std_simul_base_quarterly.C_10       = mean(simul_std_collet_base_quarterly.C_10);
std_simul_base_quarterly.C_Q1       = mean(simul_std_collet_base_quarterly.C_Q1);
std_simul_base_quarterly.C_Q2       = mean(simul_std_collet_base_quarterly.C_Q2);
std_simul_base_quarterly.C_Q3       = mean(simul_std_collet_base_quarterly.C_Q3);
std_simul_base_quarterly.C_Q4       = mean(simul_std_collet_base_quarterly.C_Q4);
std_simul_base_quarterly.C_Q5       = mean(simul_std_collet_base_quarterly.C_Q5);
std_simul_base_quarterly.C_90       = mean(simul_std_collet_base_quarterly.C_90);
std_simul_base_quarterly.C_99       = mean(simul_std_collet_base_quarterly.C_99);
std_simul_base_quarterly.C_I_1       = mean(simul_std_collet_base_quarterly.C_I_1);
std_simul_base_quarterly.C_I_10       = mean(simul_std_collet_base_quarterly.C_I_10);
std_simul_base_quarterly.C_I_Q1       = mean(simul_std_collet_base_quarterly.C_I_Q1);
std_simul_base_quarterly.C_I_Q2       = mean(simul_std_collet_base_quarterly.C_I_Q2);
std_simul_base_quarterly.C_I_Q3       = mean(simul_std_collet_base_quarterly.C_I_Q3);
std_simul_base_quarterly.C_I_Q4       = mean(simul_std_collet_base_quarterly.C_I_Q4);
std_simul_base_quarterly.C_I_Q5       = mean(simul_std_collet_base_quarterly.C_I_Q5);
std_simul_base_quarterly.C_I_90       = mean(simul_std_collet_base_quarterly.C_I_90);
std_simul_base_quarterly.C_I_99       = mean(simul_std_collet_base_quarterly.C_I_99);

std_simul_base_quarterly.I_1_1     = mean(simul_std_collet_base_quarterly.I_1_1);
std_simul_base_quarterly.I_1_10     = mean(simul_std_collet_base_quarterly.I_1_10);
std_simul_base_quarterly.I_1_Q1     = mean(simul_std_collet_base_quarterly.I_1_Q1);
std_simul_base_quarterly.I_1_Q2     = mean(simul_std_collet_base_quarterly.I_1_Q2);
std_simul_base_quarterly.I_1_Q3     = mean(simul_std_collet_base_quarterly.I_1_Q3);
std_simul_base_quarterly.I_1_Q4     = mean(simul_std_collet_base_quarterly.I_1_Q4);
std_simul_base_quarterly.I_1_Q5     = mean(simul_std_collet_base_quarterly.I_1_Q5);
std_simul_base_quarterly.I_1_90     = mean(simul_std_collet_base_quarterly.I_1_90);
std_simul_base_quarterly.I_1_99     = mean(simul_std_collet_base_quarterly.I_1_99);

std_simul_base_quarterly.I_I_1_1     = mean(simul_std_collet_base_quarterly.I_I_1_1);
std_simul_base_quarterly.I_I_1_10     = mean(simul_std_collet_base_quarterly.I_I_1_10);
std_simul_base_quarterly.I_I_1_Q1     = mean(simul_std_collet_base_quarterly.I_I_1_Q1);
std_simul_base_quarterly.I_I_1_Q2     = mean(simul_std_collet_base_quarterly.I_I_1_Q2);
std_simul_base_quarterly.I_I_1_Q3     = mean(simul_std_collet_base_quarterly.I_I_1_Q3);
std_simul_base_quarterly.I_I_1_Q4     = mean(simul_std_collet_base_quarterly.I_I_1_Q4);
std_simul_base_quarterly.I_I_1_Q5     = mean(simul_std_collet_base_quarterly.I_I_1_Q5);
std_simul_base_quarterly.I_I_1_90     = mean(simul_std_collet_base_quarterly.I_I_1_90);
std_simul_base_quarterly.I_I_1_99     = mean(simul_std_collet_base_quarterly.I_I_1_99);

std_simul_base_quarterly.A_b        = mean(simul_std_collet_base_quarterly.A_b);
std_simul_base_quarterly.A_hh       = mean(simul_std_collet_base_quarterly.A_hh);
std_simul_base_quarterly.B_hh       = mean(simul_std_collet_base_quarterly.B_hh);

%out.mean = mean_simul_base;
%out.std  = std_simul_base;