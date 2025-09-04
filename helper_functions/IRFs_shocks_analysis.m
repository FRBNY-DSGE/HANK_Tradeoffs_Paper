grid.maxlag = 1001;

x0 = zeros(grid.numstates,1);

if lllllll <= 9

    x0(end-9+lllllll) = shock_sets(lllllll);

elseif lllllll == 10

    x0(end-3) = 0.2;

end

clear MX IRF_state_sparse x

MX               = [eye(length(x0));gx];
IRF_state_sparse = [];
x                = x0;

for t=1:grid.maxlag
    IRF_state_sparse(:,t)=(MX*x)';
    x=hx*x;

end

IRF_distr=Gamma_state*IRF_state_sparse(1:grid.numstates-grid.os,1:grid.maxlag);

Statenext_IRFs   = repmat(Xss(end-grid.os+1:end),[1 grid.maxlag]) + IRF_state_sparse(grid.numstates-grid.os+1:grid.numstates,:);
Control_IRFs     = repmat(Yss(end-grid.oc+1:end),[1 grid.maxlag]) + Gamma_control(end-grid.oc+1:end,:)*IRF_state_sparse(grid.numstates+1:end,:);

Statenext_IRFs   = exp(Statenext_IRFs);
Control_IRFs     = exp(Control_IRFs);

% preparation

dist_a  = diff([0;marginal_a']);
dist_b  = diff([0;marginal_b]);
dist_se = diff([0;marginal_se]);

IRF_distr_full = repmat([dist_b;dist_a;dist_se],[1 grid.maxlag]) + IRF_distr;

IRF_distr_full_b  = IRF_distr_full(1:grid.nb,:);
IRF_distr_full_a  = IRF_distr_full(grid.nb+(1:grid.na),:);
IRF_distr_full_se = IRF_distr_full(grid.nb+grid.na+(1:grid.nse),:);

indexes_analysis;

Rcb_IRFs_p      = Statenext_IRFs(R_cb_ind,1:end);
W_IRFs_p        = Statenext_IRFs(w_ind,1:end);
Q_IRFs_p        = Statenext_IRFs(Q_ind,1:end);
A_g_IRFs_p      = Statenext_IRFs(A_g_ind,1:end);
B_F_IRFs_p      = Statenext_IRFs(B_F_ind,1:end);
eta_IRFs_p      = Statenext_IRFs(eta_ind,1:end);
Lambda_IRFs_p   = Control_IRFs(Lambda_ind,:);
pi_IRFs_p       = Control_IRFs(pi_ind,:);
V_IRFs_p        = Control_IRFs(V_ind,:);
L_IRFs_p        = Control_IRFs(L_ind,:);
J_IRFs_p        = Control_IRFs(J_ind,:);
h_IRFs_p        = Control_IRFs(h_ind,:);
v_IRFs_p        = Control_IRFs(v_ind,:);
Y_IRFs_p        = Control_IRFs(Y_ind,:);
Profit_IRFs_p   = Control_IRFs(Profit_ind,:);
r_k_IRFs_p      = Control_IRFs(r_k_ind,:);
r_a_IRFs_p      = Control_IRFs(r_a_ind,:);
T_IRFs_p        = Control_IRFs(T_ind,:);
B_IRFs_p        = Control_IRFs(B_ind,:);
K_IRFs_p        = Control_IRFs(K_ind,:);
N_IRFs_p        = Control_IRFs(N_ind,:);
LT_IRFs_p       = Control_IRFs(LT_ind,:);
G_IRFs_p        = Control_IRFs(G_ind,:);
C_IRFs_p        = Control_IRFs(C_ind,:);
MC_IRFs_p       = Control_IRFs(MC_ind,:);
RRa_IRFs_p      = Control_IRFs(RRa_ind,:);
RR_IRFs_p       = Control_IRFs(RR_ind,:);
I_IRFs_p        = Control_IRFs(I_ind,:);
A_hh_IRFs_p     = Control_IRFs(A_hh_ind,:);
unemp_IRFs_p    = Control_IRFs(unemp_ind,:);
B_hh_IRFs_p     = Control_IRFs(B_hh_ind,:);
MRS_IRFs_p      = Control_IRFs(MRS_ind,:);

w2_IRFs_p         = Control_IRFs(w2_ind,:);
Profit2_IRFs_p    = Control_IRFs(Profit2_ind,:);
Profit3_IRFs_p    = Control_IRFs(Profit3_ind,:);
Profit_FI2_IRFs_p = Control_IRFs(Profit_FI2_ind,:);
r_a2_IRFs_p       = Control_IRFs(r_a2_ind,:);


U_IRFs_p = 1-N_IRFs_p-param.Eshare;
M_IRFs_p = (U_IRFs_p + param.lambda*N_IRFs_p).*(V_IRFs_p)./(((U_IRFs_p + param.lambda*N_IRFs_p).^param.alpha + (V_IRFs_p).^param.alpha).^(1/param.alpha));
f_IRFs_p = M_IRFs_p./(U_IRFs_p+param.lambda*N_IRFs_p);
nn_IRFs_p     = ((1-param.tau_w).*W_IRFs_p/param.psi).^(1/param.xi);
Ntilde_IRFs_p = (1-param.lambda).*N_IRFs_p + M_IRFs_p;
IRFs_M_p = 100*(M_IRFs_p/SS_stats.M-1);
IRFs_f_p = 100*(f_IRFs_p/SS_stats.f-1);
IRFs_Y_p = 100*(Y_IRFs_p/SS_stats.Y-1);
IRFs_C_p = 100*(C_IRFs_p/SS_stats.C-1);


IRFs_Rcb_p       = 100*IRF_state_sparse(grid.numstates-grid.os+R_cb_ind,:);
IRFs_W_p         = 100*IRF_state_sparse(grid.numstates-grid.os+w_ind,:);
IRFs_Q_p         = 100*IRF_state_sparse(grid.numstates-grid.os+Q_ind,:);
IRFs_A_g_p       = 100*IRF_state_sparse(grid.numstates-grid.os+A_g_ind,:);
IRFs_Lambda_p    = 100*IRF_state_sparse(end-grid.oc+Lambda_ind,:);
IRFs_pi_p        = 100*IRF_state_sparse(end-grid.oc+pi_ind,:);
IRFs_V_p         = 100*IRF_state_sparse(end-grid.oc+V_ind,:);
IRFs_J_p         = 100*IRF_state_sparse(end-grid.oc+J_ind,:);
IRFs_h_p         = 100*IRF_state_sparse(end-grid.oc+h_ind,:);
IRFs_v_p         = 100*IRF_state_sparse(end-grid.oc+v_ind,:);
IRFs_Profit_p    = 100*IRF_state_sparse(end-grid.oc+Profit_ind,:);
IRFs_r_k_p       = 100*IRF_state_sparse(end-grid.oc+r_k_ind,:);
IRFs_r_a_p       = 100*IRF_state_sparse(end-grid.oc+r_a_ind,:);
IRFs_T_p         = 100*IRF_state_sparse(end-grid.oc+T_ind,:);
IRFs_B_p         = 100*IRF_state_sparse(end-grid.oc+B_ind,:);
IRFs_K_p         = 100*IRF_state_sparse(end-grid.oc+K_ind,:);
IRFs_N_p         = 100*IRF_state_sparse(end-grid.oc+N_ind,:);
IRFs_LT_p        = 100*IRF_state_sparse(end-grid.oc+LT_ind,:);
IRFs_G_p         = 100*IRF_state_sparse(end-grid.oc+G_ind,:);
IRFs_C_hh_p      = 100*IRF_state_sparse(end-grid.oc+C_ind,:);
IRFs_MC_p        = 100*IRF_state_sparse(end-grid.oc+MC_ind,:);
IRFs_RRa_p       = 100*IRF_state_sparse(end-grid.oc+RRa_ind,:);
IRFs_RR_p        = 100*IRF_state_sparse(end-grid.oc+RR_ind,:);
IRFs_I_p         = 100*IRF_state_sparse(end-grid.oc+I_ind,:);
IRFs_A_hh_p      = 100*IRF_state_sparse(end-grid.oc+A_hh_ind,:);
IRFs_unemp_p     = 100*(unemp_IRFs_p-SS_stats.u);
IRFs_B_hh_p      = 100*IRF_state_sparse(end-grid.oc+B_hh_ind,:);
IRFs_MRS_p       = 100*IRF_state_sparse(end-grid.oc+MRS_ind,:);

IRFs_w2_p           = 100*(w2_IRFs_p/param.w_bar-1);
IRFs_Profit2_p      = 100*(Profit2_IRFs_p/SS_stats.Profit-1);
IRFs_Profit3_p      = 100*(Profit3_IRFs_p/SS_stats.Profit-1);
IRFs_Profit_FI2_p   = 100*(Profit_FI2_IRFs_p/SS_stats.Profit_FI-1);
IRFs_r_a2_p         = 400*(r_a2_IRFs_p-SS_stats.r_a);


A_b_IRFs_p     = Statenext_IRFs(A_b_ind,1:end);
B_b_IRFs_p     = Statenext_IRFs(B_b_ind,1:end);
lev_IRFs_p     = Statenext_IRFs(lev_ind,1:end);
NW_b_IRFs_p    = Statenext_IRFs(NW_b_ind,1:end);
zz_IRFs_p        = Control_IRFs(zz_ind,:);
xx_IRFs_p        = Control_IRFs(xx_ind,:);
Profit_FI_IRFs_p = Control_IRFs(Profit_FI_ind,:);
vv_IRFs_p        = Control_IRFs(vv_ind,:);
ee_IRFs_p        = Control_IRFs(ee_ind,:);
C_b_IRFs_p       = Control_IRFs(C_b_ind,:);



IRFs_A_b_p       = 100*IRF_state_sparse(grid.numstates-grid.os+A_b_ind,:);
IRFs_B_b_p       = 100*IRF_state_sparse(grid.numstates-grid.os+B_b_ind,:);
IRFs_lev_p       = 100*IRF_state_sparse(grid.numstates-grid.os+lev_ind,:);
IRFs_NW_b_p      = 100*IRF_state_sparse(grid.numstates-grid.os+NW_b_ind,:);
IRFs_zz_p        = 100*IRF_state_sparse(end-grid.oc+zz_ind,:);
IRFs_xx_p        = 100*IRF_state_sparse(end-grid.oc+xx_ind,:);
IRFs_Profit_FI_p = 100*IRF_state_sparse(end-grid.oc+Profit_FI_ind,:);
IRFs_vv_p        = 100*IRF_state_sparse(end-grid.oc+vv_ind,:);
IRFs_ee_p        = 100*IRF_state_sparse(end-grid.oc+ee_ind,:);
IRFs_C_b_p       = 100*IRF_state_sparse(end-grid.oc+C_b_ind,:);

TC_ss = param.w_bar*SS_stats.L+param.delta_0*param.v^param.delta_1*grid.K+param.iota*SS_stats.V+param.fix;

Amarkup = (1-TC_ss./SS_stats.Y);

Markup_IRFs_p    = 1-(W_IRFs_p.*L_IRFs_p+param.delta_0*v_IRFs_p.^param.delta_1.*K_IRFs_p+param.iota*V_IRFs_p+param.fix)./Y_IRFs_p;
IRFs_Markups_p   = 100*(Markup_IRFs_p-Amarkup);
IRFs_Profit_NF_p = 100*((Profit_IRFs_p-Profit_FI_IRFs_p+param.fix2)/SS_stats.Profit-1);
IRFs_Profit_F_p  = 100*((SS_stats.Profit-param.fix2+Profit_FI_IRFs_p)/SS_stats.Profit-1);

AC_IRFs_p = (param.delta_0*v_IRFs_p(1:end-1).^param.delta_1.*K_IRFs_p(1:end-1) + W_IRFs_p(2:end).*L_IRFs_p(1:end-1) + param.fix + param.iota*V_IRFs_p(1:end-1))./Y_IRFs_p(1:end-1);
IRFs_AC_p = 100*(AC_IRFs_p/SS_stats.AvgC-1);

IRFs_Profit_int_p = 100*(((1-MC_IRFs_p).*Y_IRFs_p-param.fix)/SS_stats.Profit_int-1);

IRFs.AC_IRFs_p = AC_IRFs_p;
IRFs.IRFs_AC_p = IRFs_AC_p;
IRFs.IRFs_Profit_int_p = IRFs_Profit_int_p;


IRFs.Rcb_IRFs_p = Rcb_IRFs_p;
IRFs.W_IRFs_p = W_IRFs_p;
IRFs.Q_IRFs_p = Q_IRFs_p;
IRFs.A_g_IRFs_p = A_g_IRFs_p;
IRFs.Lambda_IRFs_p = Lambda_IRFs_p;
IRFs.pi_IRFs_p = pi_IRFs_p;
IRFs.V_IRFs_p = V_IRFs_p;
IRFs.J_IRFs_p = J_IRFs_p;
IRFs.h_IRFs_p = h_IRFs_p;
IRFs.v_IRFs_p = v_IRFs_p;
IRFs.Y_IRFs_p = Y_IRFs_p;
IRFs.Profit_IRFs_p = Profit_IRFs_p;
IRFs.r_k_IRFs_p = r_k_IRFs_p;
IRFs.r_a_IRFs_p = r_a_IRFs_p;
IRFs.T_IRFs_p = T_IRFs_p;
IRFs.B_IRFs_p = B_IRFs_p;
IRFs.K_IRFs_p = K_IRFs_p;
IRFs.N_IRFs_p = N_IRFs_p;
IRFs.LT_IRFs_p = LT_IRFs_p;
IRFs.G_IRFs_p = G_IRFs_p;
IRFs.C_IRFs_p = C_IRFs_p;
IRFs.MC_IRFs_p = MC_IRFs_p;
IRFs.RRa_IRFs_p = RRa_IRFs_p;
IRFs.RR_IRFs_p = RR_IRFs_p;
IRFs.I_IRFs_p = I_IRFs_p;
IRFs.A_hh_IRFs_p = A_hh_IRFs_p;
IRFs.unemp_IRFs_p = unemp_IRFs_p;
IRFs.B_hh_IRFs_p = B_hh_IRFs_p;
IRFs.MRS_IRFs_p = MRS_IRFs_p;

IRFs.Markup_IRFs_p = Markup_IRFs_p;
IRFs.IRFs_Markups_p = IRFs_Markups_p;
IRFs.IRFs_Profit_NF_p = IRFs_Profit_NF_p;
IRFs.IRFs_Profit_F_p = IRFs_Profit_F_p;

IRFs.U_IRFs_p = U_IRFs_p;
IRFs.M_IRFs_p = M_IRFs_p;
IRFs.f_IRFs_p = f_IRFs_p;
IRFs.nn_IRFs_p = nn_IRFs_p;
IRFs.Ntilde_IRFs_p = Ntilde_IRFs_p;
IRFs.IRFs_M_p = IRFs_M_p;
IRFs.IRFs_f_p = IRFs_f_p;
IRFs.IRFs_Y_p = IRFs_Y_p;

IRFs.IRFs_Rcb_p = IRFs_Rcb_p;
IRFs.IRFs_W_p = IRFs_W_p;
IRFs.IRFs_Q_p = IRFs_Q_p;
IRFs.IRFs_A_g_p = IRFs_A_g_p;
IRFs.IRFs_Lambda_p = IRFs_Lambda_p;
IRFs.IRFs_pi_p = IRFs_pi_p;
IRFs.IRFs_V_p = IRFs_V_p;
IRFs.IRFs_J_p = IRFs_J_p;
IRFs.IRFs_h_p = IRFs_h_p;
IRFs.IRFs_v_p = IRFs_v_p;
IRFs.IRFs_Profit_p = IRFs_Profit_p;
IRFs.IRFs_r_k_p = IRFs_r_k_p;
IRFs.IRFs_r_a_p = IRFs_r_a_p;
IRFs.IRFs_T_p = IRFs_T_p;
IRFs.IRFs_B_p = IRFs_B_p;
IRFs.IRFs_K_p = IRFs_K_p;
IRFs.IRFs_N_p = IRFs_N_p;
IRFs.IRFs_LT_p = IRFs_LT_p;
IRFs.IRFs_G_p = IRFs_G_p;
IRFs.IRFs_C_hh_p = IRFs_C_hh_p;
IRFs.IRFs_MC_p = IRFs_MC_p;
plot_irfIRFs.IRFs_RRa_p = IRFs_RRa_p;
IRFs.IRFs_RR_p = IRFs_RR_p;
IRFs.IRFs_I_p = IRFs_I_p;
IRFs.IRFs_A_hh_p = IRFs_A_hh_p;
IRFs.IRFs_unemp_p = IRFs_unemp_p;
IRFs.IRFs_B_hh_p = IRFs_B_hh_p;
IRFs.IRFs_MRS_p = IRFs_MRS_p;

IRFs.A_b_IRFs_p       = A_b_IRFs_p;
IRFs.B_b_IRFs_p       = B_b_IRFs_p;
IRFs.lev_IRFs_p       = lev_IRFs_p;
IRFs.NW_b_IRFs_p      = NW_b_IRFs_p;
IRFs.zz_IRFs_p        = zz_IRFs_p;
IRFs.xx_IRFs_p        = xx_IRFs_p;
IRFs.Profit_FI_IRFs_p = Profit_FI_IRFs_p;
IRFs.vv_IRFs_p        = vv_IRFs_p;
IRFs.ee_IRFs_p        = ee_IRFs_p;
IRFs.C_b_IRFs_p       = C_b_IRFs_p;

IRFs.IRFs_A_b_p       =IRFs_A_b_p;
IRFs.IRFs_B_b_p       =IRFs_B_b_p;
IRFs.IRFs_lev_p       =IRFs_lev_p;
IRFs.IRFs_NW_b_p      =IRFs_NW_b_p;
IRFs.IRFs_zz_p        =IRFs_zz_p;
IRFs.IRFs_xx_p        =IRFs_xx_p;
IRFs.IRFs_Profit_FI_p =IRFs_Profit_FI_p;
IRFs.IRFs_vv_p        =IRFs_vv_p;
IRFs.IRFs_ee_p        =IRFs_ee_p;
IRFs.IRFs_C_b_p       =IRFs_C_b_p;

B_FI_IRFs_p      = Q_IRFs_p.*A_b_IRFs_p-NW_b_IRFs_p;
B_cp_IRFs_p      = Q_IRFs_p.*(1+param.tau_cp).*(A_g_IRFs_p-1);
IRFs_B_FI_p      = 100*(B_FI_IRFs_p/SS_stats.B_FI-1);
B_gov_IRFs_p     = B_b_IRFs_p(2:end) + B_hh_IRFs_p(1:end-1)/(1-param.death_rate) + SS_stats.B_F - B_FI_IRFs_p(2:end);
B_gov_ncp_IRFs_p = B_gov_IRFs_p - B_cp_IRFs_p(2:end);
IRFs_B_gov_p     = 100*(B_gov_IRFs_p/SS_stats.B_gov-1);
IRFs_B_gov_ncp_p = 100*(B_gov_ncp_IRFs_p/SS_stats.B_gov_ncp-1);

IRFs.B_FI_IRFs_p  = B_FI_IRFs_p;
IRFs.B_cp_IRFs_p  = B_cp_IRFs_p;
IRFs.IRFs_B_FI_p  = IRFs_B_FI_p;
IRFs.B_gov_IRFs_p  = B_gov_IRFs_p;
IRFs.B_gov_ncp_IRFs_p  = B_gov_ncp_IRFs_p;
IRFs.IRFs_B_gov_p  = IRFs_B_gov_p;
IRFs.IRFs_B_gov_ncp_p  = IRFs_B_gov_ncp_p;

IRFs.w2_IRFs_p         = w2_IRFs_p;
IRFs.Profit2_IRFs_p    = Profit2_IRFs_p;
IRFs.Profit3_IRFs_p    = Profit3_IRFs_p;
IRFs.Profit_FI2_IRFs_p = Profit_FI2_IRFs_p;
IRFs.r_a2_IRFs_p       = r_a2_IRFs_p;

IRFs.IRFs_w2_p         = IRFs_w2_p;
IRFs.IRFs_Profit2_p    = IRFs_Profit2_p;
IRFs.IRFs_Profit3_p    = IRFs_Profit3_p;
IRFs.IRFs_Profit_FI2_p = IRFs_Profit_FI2_p;
IRFs.IRFs_r_a2_p       = IRFs_r_a2_p;

%%Adding B10 Consumption IRF
C_10_IRFs_p = Control_IRFs(C_10_ind, :);
IRFs.IRFs_C_10_p = 100 *((C_10_IRFs_p/SS_stats.C_B10) -1);
IRFs.C_10_IRFs_p = C_10_IRFs_p;

C_Q1_IRFs_p = Control_IRFs(C_Q1_ind, :);
C_Q5_IRFs_p = Control_IRFs(C_Q5_ind, :);

IRFs.IRFs_C_Q1_p = 100 * ((C_Q1_IRFs_p/SS_stats.C_Q1) -1);
IRFs.C_Q1_IRFs_p = C_Q1_IRFs_p;
IRFs.IRFs_C_Q5_p = 100 * ((C_Q5_IRFs_p/SS_stats.C_Q5) -1);
IRFs.C_Q5_IRFs_p = C_Q5_IRFs_p;


%% Q3 Consumption IRF
C_Q3_IRFs_p = Control_IRFs(C_Q3_ind, :);
IRFs.IRFs_C_Q3_p = 100 * ((C_Q3_IRFs_p/SS_stats.C_Q3) -1);
IRFs.C_Q3_IRFs_p = C_Q3_IRFs_p;


%P10 Consumption IRF
C_90_IRFs_p = Control_IRFs(C_90_ind, :);
IRFs.IRFs_C_90_p = 100 * ((C_90_IRFs_p/SS_stats.C_P90) -1);
IRFs.C_90_IRFs_p = C_90_IRFs_p;



I_90to10_IRFs_p = Control_IRFs(I_90to10_ind,:);
I_50to10_IRFs_p = Control_IRFs(I_50to10_ind,:);
I_90to50_IRFs_p = Control_IRFs(I_90to50_ind,:);

I_1_01_IRFs_p = Control_IRFs(I_1_01_ind,:);
I_1_1_IRFs_p = Control_IRFs(I_1_1_ind,:);
I_1_10_IRFs_p = Control_IRFs(I_1_10_ind,:);
I_1_Q1_IRFs_p = Control_IRFs(I_1_Q1_ind,:);
I_1_Q2_IRFs_p = Control_IRFs(I_1_Q2_ind,:);
I_1_Q3_IRFs_p = Control_IRFs(I_1_Q3_ind,:);
I_1_Q4_IRFs_p = Control_IRFs(I_1_Q4_ind,:);
I_1_Q5_IRFs_p = Control_IRFs(I_1_Q5_ind,:);
I_1_90_IRFs_p = Control_IRFs(I_1_90_ind,:);
I_1_99_IRFs_p = Control_IRFs(I_1_99_ind,:);
I_1_999_IRFs_p = Control_IRFs(I_1_999_ind,:);

I_2_01_IRFs_p = Control_IRFs(I_2_01_ind,:);
I_2_1_IRFs_p = Control_IRFs(I_2_1_ind,:);
I_2_10_IRFs_p = Control_IRFs(I_2_10_ind,:);
I_2_Q1_IRFs_p = Control_IRFs(I_2_Q1_ind,:);
I_2_Q2_IRFs_p = Control_IRFs(I_2_Q2_ind,:);
I_2_Q3_IRFs_p = Control_IRFs(I_2_Q3_ind,:);
I_2_Q4_IRFs_p = Control_IRFs(I_2_Q4_ind,:);
I_2_Q5_IRFs_p = Control_IRFs(I_2_Q5_ind,:);
I_2_90_IRFs_p = Control_IRFs(I_2_90_ind,:);
I_2_99_IRFs_p = Control_IRFs(I_2_99_ind,:);
I_2_999_IRFs_p = Control_IRFs(I_2_999_ind,:);


W_90to50_IRFs_p = Control_IRFs(W_90to50_ind, :);
Inc_I_1_T10toB10_IRFs_p =Control_IRFs(Inc_I_1_T10toB10_ind, :);
Inc_I_1_T10toQ3_IRFs_p = Control_IRFs(Inc_I_1_T10toQ3_ind, :);
Inc_I_1_Q3toB10_IRFs_p = Control_IRFs(Inc_I_1_Q3toB10_ind, :);
Inc_I_2_T10toB10_IRFs_p =Control_IRFs(Inc_I_2_T10toB10_ind, :);
Inc_I_2_T10toQ3_IRFs_p = Control_IRFs(Inc_I_2_T10toQ3_ind, :);
Inc_I_2_Q3toB10_IRFs_p = Control_IRFs(Inc_I_2_Q3toB10_ind, :);
Inc_1_T10toB10_IRFs_p =Control_IRFs(Inc_1_T10toB10_ind, :);
Inc_1_T10toQ3_IRFs_p = Control_IRFs(Inc_1_T10toQ3_ind, :);
Inc_1_Q3toB10_IRFs_p = Control_IRFs(Inc_1_Q3toB10_ind, :);
Inc_2_T10toB10_IRFs_p =Control_IRFs(Inc_2_T10toB10_ind, :);
Inc_2_T10toQ3_IRFs_p = Control_IRFs(Inc_2_T10toQ3_ind, :);
Inc_2_Q3toB10_IRFs_p = Control_IRFs(Inc_2_Q3toB10_ind, :);
Wlth_T10toQ3_IRFs_p = Control_IRFs(Wlth_T10toQ3_ind, :);

C_T10toB10_IRFs_p =Control_IRFs(C_T10toB10_ind, :);
C_T10toQ3_IRFs_p = Control_IRFs(C_T10toQ3_ind, :);
C_Q3toB10_IRFs_p = Control_IRFs(C_Q3toB10_ind, :);

C_I_T10toB10_IRFs_p =Control_IRFs(C_I_T10toB10_ind, :);
C_I_T10toQ3_IRFs_p = Control_IRFs(C_I_T10toQ3_ind, :);
C_I_Q3toB10_IRFs_p = Control_IRFs(C_I_Q3toB10_ind, :);


GiniW_IRFs_p    = Control_IRFs(GiniW_ind,:);
GiniI1_IRFs_p    = Control_IRFs(GiniI_1_ind,:);
GiniI2_IRFs_p = Control_IRFs(GiniI_2_ind,:);
GiniC_IRFs_p    = Control_IRFs(GiniC_ind,:);
IRFs_GiniW_p    = GiniW_IRFs_p - SS_stats.GiniW;
IRFs_GiniI1_p    = GiniI1_IRFs_p - SS_stats.GiniIncome;
IRFs_GiniI2_p    = GiniI2_IRFs_p - SS_stats.GiniIncome;
IRFs_GiniC_p    = GiniC_IRFs_p - SS_stats.GiniCons;

IRFs.GiniW_IRFs_p = GiniW_IRFs_p;
IRFs.GiniI1_IRFs_p = GiniI1_IRFs_p;
IRFs.GiniI2_IRFs_p = GiniI2_IRFs_p;
IRFs.GiniC_IRFs_p = GiniC_IRFs_p;
IRFs.IRFs_GiniW_p = IRFs_GiniW_p;
IRFs.IRFs_GiniI1_p = IRFs_GiniI1_p;
IRFs.IRFs_GiniI2_p = IRFs_GiniI2_p;
IRFs.IRFs_GiniC_p = IRFs_GiniC_p;

IRFs_I_90to10_p = I_90to10_IRFs_p - (SS_stats.income90/SS_stats.income10);
IRFs_I_50to10_p = I_50to10_IRFs_p - (SS_stats.income50/SS_stats.income10);
IRFs_I_90to50_p = I_90to50_IRFs_p - (SS_stats.income90/SS_stats.income50);


IRFs_W_90to50_p = W_90to50_IRFs_p - (SS_stats.w90/SS_stats.w50);
IRFs_Inc_I_1_T10toB10_p =Inc_I_1_T10toB10_IRFs_p - (SS_stats.income_I_P90/SS_stats.income_I_B10);
IRFs_Inc_I_1_T10toQ3_p = Inc_I_1_T10toQ3_IRFs_p - (SS_stats.income_I_P90/SS_stats.income_I_Q3);
IRFs_Inc_I_1_Q3toB10_p = Inc_I_1_Q3toB10_IRFs_p - (SS_stats.income_I_Q3/SS_stats.income_I_B10);
IRFs_Inc_I_2_T10toB10_p =Inc_I_2_T10toB10_IRFs_p - (SS_stats.income_I_P90/SS_stats.income_I_B10);
IRFs_Inc_I_2_T10toQ3_p = Inc_I_2_T10toQ3_IRFs_p - (SS_stats.income_I_P90/SS_stats.income_I_Q3);
IRFs_Inc_I_2_Q3toB10_p = Inc_I_2_Q3toB10_IRFs_p - (SS_stats.income_I_Q3/SS_stats.income_I_B10);
IRFs_Inc_1_T10toB10_p =Inc_1_T10toB10_IRFs_p - (SS_stats.income_P90/SS_stats.income_B10);
IRFs_Inc_1_T10toQ3_p = Inc_1_T10toQ3_IRFs_p - (SS_stats.income_P90/SS_stats.income_Q3);
IRFs_Inc_1_Q3toB10_p = Inc_1_Q3toB10_IRFs_p - (SS_stats.income_Q3/SS_stats.income_B10);
IRFs_Inc_2_T10toB10_p =Inc_2_T10toB10_IRFs_p - (SS_stats.income_P90/SS_stats.income_B10);
IRFs_Inc_2_T10toQ3_p = Inc_2_T10toQ3_IRFs_p - (SS_stats.income_P90/SS_stats.income_Q3);
IRFs_Inc_2_Q3toB10_p = Inc_2_Q3toB10_IRFs_p - (SS_stats.income_Q3/SS_stats.income_B10);
IRFs_Wlth_T10toQ3_p = Wlth_T10toQ3_IRFs_p - (SS_stats.w_P90/SS_stats.w_Q3);




IRFs_C_T10toB10_p = C_T10toB10_IRFs_p - (SS_stats.C_P90/SS_stats.C_B10);
IRFs_C_T10toQ3_p = C_T10toQ3_IRFs_p - (SS_stats.C_P90/SS_stats.C_Q3);
IRFs_C_Q3toB10_p = C_Q3toB10_IRFs_p - (SS_stats.C_Q3/SS_stats.C_B10);
