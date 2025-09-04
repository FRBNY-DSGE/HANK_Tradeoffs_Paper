function out = compute_Jacob_GR_Taylor_v1(param,grid,SS_stats,idx,StateSS,ControlSS)

% indexes

R_cb_ind        = idx.R_cb_ind        ;
w_ind           = idx.w_ind           ;
A_b_ind         = idx.A_b_ind         ;
B_b_ind         = idx.B_b_ind         ;
A_g_ind         = idx.A_g_ind         ;
Q_ind           = idx.Q_ind           ;
lev_ind         = idx.lev_ind         ;
NW_b_ind        = idx.NW_b_ind        ;
R_tilde_ind     = idx.R_tilde_ind     ;
x_cb_ind        = idx.x_cb_ind        ;
pastpi_ind      = idx.pastpi_ind      ;
pastY_ind       = idx.pastY_ind       ;
pastC_ind       = idx.pastC_ind       ;
pastI_ind       = idx.pastI_ind       ;
pastPROFIT_ind  = idx.pastPROFIT_ind  ;
pastunemp_ind   = idx.pastunemp_ind   ;
pastG_ind       = idx.pastG_ind       ;
pastLT_ind      = idx.pastLT_ind      ;
R_star_ind      = idx.R_star_ind      ;
B_F_ind         = idx.B_F_ind         ;
Z_ind           = idx.Z_ind           ;
PSI_RP_ind      = idx.PSI_RP_ind      ;
eta_ind         = idx.eta_ind         ;
D_ind           = idx.D_ind           ;
GG_ind          = idx.GG_ind          ;
iota_ind        = idx.iota_ind        ;
BB_ind          = idx.BB_ind          ;
PSI_W_ind       = idx.PSI_W_ind       ;
MP_ind          = idx.MP_ind          ;
p_mark_ind      = idx.p_mark_ind      ;
pgap_ind        = idx.pgap_ind        ;
ygap_ind        = idx.ygap_ind        ;
eps_QE_ind      = idx.eps_QE_ind      ;
eps_B_F_ind     = idx.eps_B_F_ind     ;
eps_RP_ind      = idx.eps_RP_ind      ;
eps_BB_ind      = idx.eps_BB_ind      ;
eps_Z_ind       = idx.eps_Z_ind       ;
eps_G_ind       = idx.eps_G_ind       ;
eps_D_ind       = idx.eps_D_ind       ;
eps_R_ind       = idx.eps_R_ind       ;
eps_iota_ind    = idx.eps_iota_ind    ;
eps_eta_ind     = idx.eps_eta_ind     ;
eps_w_ind       = idx.eps_w_ind       ;

MRS_ind         = idx.MRS_ind;
A_hh_ind        = idx.A_hh_ind;
B_hh_ind        = idx.B_hh_ind;
C_ind           = idx.C_ind;
N_ind           = idx.N_ind;
L_ind           = idx.L_ind;
UB_ind          = idx.UB_ind;
K_ind           = idx.K_ind;
B_ind           = idx.B_ind;
B_gov_ncp_ind   = idx.B_gov_ncp_ind;
T_ind           = idx.T_ind;
LT_ind          = idx.LT_ind;
G_ind           = idx.G_ind;
Lambda_ind      = idx.Lambda_ind;
pi_ind          = idx.pi_ind;
V_ind           = idx.V_ind;
J_ind           = idx.J_ind;
h_ind           = idx.h_ind;
v_ind           = idx.v_ind;
Y_ind           = idx.Y_ind;
Profit_ind      = idx.Profit_ind;
r_k_ind         = idx.r_k_ind;
r_a_ind         = idx.r_a_ind;
MC_ind          = idx.MC_ind;
unemp_ind       = idx.unemp_ind;
nn_ind          = idx.nn_ind;
M_ind           = idx.M_ind;
f_ind           = idx.f_ind;
zz_ind          = idx.zz_ind;
xx_ind          = idx.xx_ind;
vv_ind          = idx.vv_ind;
ee_ind          = idx.ee_ind;
C_b_ind         = idx.C_b_ind;
Profit_FI_ind   = idx.Profit_FI_ind;
RRa_ind         = idx.RRa_ind;
RR_ind          = idx.RR_ind;
I_ind           = idx.I_ind;
x_k_ind         = idx.x_k_ind;
A_g_obs_ind     = idx.A_g_obs_ind;
Y_obs_ind       = idx.Y_obs_ind;
C_obs_ind       = idx.C_obs_ind;
I_obs_ind       = idx.I_obs_ind;
w_obs_ind       = idx.w_obs_ind;
PROFIT_obs_ind  = idx.PROFIT_obs_ind;
unemp_obs_ind   = idx.unemp_obs_ind;
inf_obs_ind     = idx.inf_obs_ind;
R_obs_ind       = idx.R_obs_ind;
pastw_ind       = idx.pastw_ind;
G_obs_ind       = idx.G_obs_ind;
pastYY_ind      = idx.pastYY_ind;
pastCC_ind      = idx.pastCC_ind;
pastII_ind      = idx.pastII_ind;
pastPPROFIT_ind = idx.pastPPROFIT_ind;
pastuu_ind      = idx.pastuu_ind;
pastGG_ind      = idx.pastGG_ind;
pastA_g_ind     = idx.pastA_g_ind;
l_lambda_ind    = idx.l_lambda_ind;
pastQ_ind       = idx.pastQ_ind;
x_I_ind         = idx.x_I_ind; 
eta2_ind        = idx.eta2_ind; 
iota2_ind       = idx.iota2_ind; 
pastLT2_ind     = idx.pastLT2_ind; 
pastG2_ind      = idx.pastG2_ind;
LT_obs_ind      = idx.LT_obs_ind;
B_gov_ncp2_ind  = idx.B_gov_ncp2_ind;

w2_ind          = idx.w2_ind;
Profit2_ind     = idx.Profit2_ind;
Profit3_ind     = idx.Profit3_ind;
Profit_FI2_ind  = idx.Profit_FI2_ind;
r_a2_ind        = idx.r_a2_ind;

perturb_size  = 1e-5;

State_perturb   = zeros(grid.os,1) + perturb_size;
Control_perturb = zeros(grid.oc,1) + perturb_size;

State_aux       = exp( StateSS(end-grid.os+1:end)   + State_perturb );
Control_aux     = exp( ControlSS(end-grid.oc+1:end) + Control_perturb );

State_ss        = exp( StateSS(end-grid.os+1:end)    );
Control_ss      = exp( ControlSS(end-grid.oc+1:end)  );

R_cb_aux    = State_aux(R_cb_ind);
w_aux       = State_aux(w_ind);
A_b_aux     = State_aux(A_b_ind);
B_b_aux     = State_aux(B_b_ind);
A_g_aux     = State_aux(A_g_ind);
Q_aux       = State_aux(Q_ind);
lev_aux     = State_aux(lev_ind);
NW_b_aux    = State_aux(NW_b_ind);
R_tilde_aux = State_aux(R_tilde_ind);
x_cb_aux    = State_aux(x_cb_ind);
pastpi_aux  = State_aux(pastpi_ind);
pastY_aux   = State_aux(pastY_ind);
pastC_aux   = State_aux(pastC_ind);
pastI_aux   = State_aux(pastI_ind);
pastPROFIT_aux   = State_aux(pastPROFIT_ind);
pastunemp_aux   = State_aux(pastunemp_ind);
pastG_aux   = State_aux(pastG_ind);
pastLT_aux   = State_aux(pastLT_ind);
R_star_aux  = State_aux(R_star_ind);
B_F_aux  = State_aux(B_F_ind);
Z_aux       = State_aux(Z_ind);
PSI_RP_aux  = State_aux(PSI_RP_ind);
eta_aux     = State_aux(eta_ind);
D_aux       = State_aux(D_ind);
GG_aux      = State_aux(GG_ind);
iota_aux    = State_aux(iota_ind);
BB_aux      = State_aux(BB_ind);
PSI_W_aux   = State_aux(PSI_W_ind); 
MP_aux      = State_aux(MP_ind); 
p_mark_aux  = State_aux(p_mark_ind); 

pgap_aux    = State_aux(pgap_ind);
ygap_aux    = State_aux(ygap_ind);

eps_QE_aux     = perturb_size;
eps_B_F_aux     = perturb_size;
eps_RP_aux     = perturb_size;
eps_eta_aux    = perturb_size;
eps_D_aux      = perturb_size;
eps_R_aux      = perturb_size;
eps_Z_aux      = perturb_size;
eps_G_aux      = perturb_size;
eps_iota_aux   = perturb_size;
eps_BB_aux     = perturb_size;
eps_w_aux      = perturb_size;

MRS_aux       = Control_aux(MRS_ind);
A_hh_aux      = Control_aux(A_hh_ind);
B_hh_aux      = Control_aux(B_hh_ind);
C_aux         = Control_aux(C_ind);
N_aux         = Control_aux(N_ind);
L_aux         = Control_aux(L_ind);
UB_aux        = Control_aux(UB_ind);
K_aux         = Control_aux(K_ind);
B_aux         = Control_aux(B_ind);
B_gov_ncp_aux = Control_aux(B_gov_ncp_ind);
T_aux         = Control_aux(T_ind);
LT_aux        = Control_aux(LT_ind);
G_aux         = Control_aux(G_ind);
Lambda_aux    = Control_aux(Lambda_ind);
pi_aux        = Control_aux(pi_ind);
V_aux         = Control_aux(V_ind);
J_aux         = Control_aux(J_ind);
h_aux         = Control_aux(h_ind);
v_aux         = Control_aux(v_ind);
Y_aux         = Control_aux(Y_ind);
Profit_aux    = Control_aux(Profit_ind);
r_k_aux       = Control_aux(r_k_ind);
r_a_aux       = Control_aux(r_a_ind);
MC_aux        = Control_aux(MC_ind);
unemp_aux     = Control_aux(unemp_ind);
nn_aux        = Control_aux(nn_ind);
M_aux         = Control_aux(M_ind);
f_aux         = Control_aux(f_ind);
zz_aux        = Control_aux(zz_ind);
xx_aux        = Control_aux(xx_ind);
vv_aux        = Control_aux(vv_ind);
ee_aux        = Control_aux(ee_ind);
C_b_aux       = Control_aux(C_b_ind);
Profit_FI_aux = Control_aux(Profit_FI_ind);
RRa_aux       = Control_aux(RRa_ind);
RR_aux        = Control_aux(RR_ind);
I_aux         = Control_aux(I_ind);
x_k_aux       = Control_aux(x_k_ind);
A_g_obs_aux   = Control_aux(A_g_obs_ind);
Y_obs_aux     = Control_aux(Y_obs_ind);
C_obs_aux     = Control_aux(C_obs_ind);
I_obs_aux     = Control_aux(I_obs_ind);
w_obs_aux     = Control_aux(w_obs_ind);
PROFIT_obs_aux     = Control_aux(PROFIT_obs_ind);
unemp_obs_aux     = Control_aux(unemp_obs_ind);
inf_obs_aux     = Control_aux(inf_obs_ind);
R_obs_aux     = Control_aux(R_obs_ind);
pastw_aux     = Control_aux(pastw_ind);
G_obs_aux     = Control_aux(G_obs_ind);

pastYY_aux = Control_aux(pastYY_ind);
pastCC_aux = Control_aux(pastCC_ind);
pastII_aux = Control_aux(pastII_ind);
pastPPROFIT_aux = Control_aux(pastPPROFIT_ind);
pastuu_aux = Control_aux(pastuu_ind);
pastGG_aux = Control_aux(pastGG_ind);
pastA_g_aux = Control_aux(pastA_g_ind);

l_lambda_aux = Control_aux(l_lambda_ind);
pastQ_aux    = Control_aux(pastQ_ind);
x_I_aux      = Control_aux(x_I_ind); 
eta2_aux     = Control_aux(eta2_ind); 
iota2_aux    = Control_aux(iota2_ind);
pastLT2_aux  = Control_aux(pastLT2_ind); 
pastG2_aux   = Control_aux(pastG2_ind);
LT_obs_aux   = Control_aux(LT_obs_ind);
B_gov_ncp2_aux = Control_aux(B_gov_ncp2_ind);

w2_aux          = Control_aux(w2_ind);
Profit2_aux     = Control_aux(Profit2_ind);
Profit3_aux     = Control_aux(Profit3_ind);
Profit_FI2_aux  = Control_aux(Profit_FI2_ind);
r_a2_aux        = Control_aux(r_a2_ind);

tau_w_aux = param.tau_w;
tau_a_aux = param.tau_a;

R_cb_ss    = State_ss(R_cb_ind);
w_ss       = State_ss(w_ind);
A_b_ss     = State_ss(A_b_ind);
B_b_ss     = State_ss(B_b_ind);
A_g_ss     = State_ss(A_g_ind);
Q_ss       = State_ss(Q_ind);
lev_ss     = State_ss(lev_ind);
NW_b_ss    = State_ss(NW_b_ind);
R_tilde_ss = State_ss(R_tilde_ind);
x_cb_ss    = State_ss(x_cb_ind);
pastpi_ss  = State_ss(pastpi_ind);
pastY_ss   = State_ss(pastY_ind);
pastC_ss   = State_ss(pastC_ind);
pastI_ss   = State_ss(pastI_ind);
pastPROFIT_ss   = State_ss(pastPROFIT_ind);
pastunemp_ss   = State_ss(pastunemp_ind);
pastG_ss   = State_ss(pastG_ind);
pastLT_ss  = State_ss(pastLT_ind);
R_star_ss  = State_ss(R_star_ind);
B_F_ss     = State_ss(B_F_ind);


Z_ss       = State_ss(Z_ind);
PSI_RP_ss  = State_ss(PSI_RP_ind);
eta_ss     = State_ss(eta_ind);
D_ss       = State_ss(D_ind);
GG_ss      = State_ss(GG_ind);
iota_ss    = State_ss(iota_ind);
BB_ss      = State_ss(BB_ind);
PSI_W_ss   = State_ss(PSI_W_ind); 
MP_ss      = State_ss(MP_ind); 
p_mark_ss  = State_ss(p_mark_ind); 

pgap_ss    = State_ss(pgap_ind);
ygap_ss    = State_ss(ygap_ind);


eps_RP_ss   = 0;
eps_B_F_ss  = 0;
eps_eta_ss  = 0;
eps_D_ss    = 0;
eps_R_ss    = 0;
eps_Z_ss    = 0;
eps_G_ss    = 0;
eps_iota_ss = 0;
eps_QE_ss   = 0;
eps_w_ss    = 0;

eps_varphi_ss = 0;
eps_BB_ss     = 0;

MRS_ss       = Control_ss(MRS_ind);
A_hh_ss      = Control_ss(A_hh_ind);
B_hh_ss      = Control_ss(B_hh_ind);
C_ss         = Control_ss(C_ind);
N_ss         = Control_ss(N_ind);
L_ss         = Control_ss(L_ind);
UB_ss        = Control_ss(UB_ind);
K_ss         = Control_ss(K_ind);
B_ss         = Control_ss(B_ind);
B_gov_ncp_ss = Control_ss(B_gov_ncp_ind);
T_ss         = Control_ss(T_ind);
LT_ss        = Control_ss(LT_ind);
G_ss         = Control_ss(G_ind);
Lambda_ss    = Control_ss(Lambda_ind);
pi_ss        = Control_ss(pi_ind);
V_ss         = Control_ss(V_ind);
J_ss         = Control_ss(J_ind);
h_ss         = Control_ss(h_ind);
v_ss         = Control_ss(v_ind);
Y_ss         = Control_ss(Y_ind);
Profit_ss    = Control_ss(Profit_ind);
r_k_ss       = Control_ss(r_k_ind);
r_a_ss       = Control_ss(r_a_ind);
MC_ss        = Control_ss(MC_ind);
unemp_ss     = Control_ss(unemp_ind);
nn_ss        = Control_ss(nn_ind);
M_ss         = Control_ss(M_ind);
f_ss         = Control_ss(f_ind);
zz_ss        = Control_ss(zz_ind);
xx_ss        = Control_ss(xx_ind);
vv_ss        = Control_ss(vv_ind);
ee_ss        = Control_ss(ee_ind);
C_b_ss       = Control_ss(C_b_ind);
Profit_FI_ss = Control_ss(Profit_FI_ind);
RRa_ss       = Control_ss(RRa_ind);
RR_ss        = Control_ss(RR_ind);
I_ss         = Control_ss(I_ind);
x_k_ss       = Control_ss(x_k_ind);
A_g_obs_ss   = Control_ss(A_g_obs_ind);
Y_obs_ss     = Control_ss(Y_obs_ind);
C_obs_ss     = Control_ss(C_obs_ind);
I_obs_ss     = Control_ss(I_obs_ind);
w_obs_ss     = Control_ss(w_obs_ind);
PROFIT_obs_ss= Control_ss(PROFIT_obs_ind);
unemp_obs_ss = Control_ss(unemp_obs_ind);
inf_obs_ss   = Control_ss(inf_obs_ind);
R_obs_ss     = Control_ss(R_obs_ind);
pastw_ss     = Control_ss(pastw_ind);
G_obs_ss     = Control_ss(G_obs_ind);

pastYY_ss = Control_ss(pastYY_ind);
pastCC_ss = Control_ss(pastCC_ind);
pastII_ss = Control_ss(pastII_ind);
pastPPROFIT_ss = Control_ss(pastPPROFIT_ind);
pastuu_ss = Control_ss(pastuu_ind);
pastGG_ss = Control_ss(pastGG_ind);
pastA_g_ss = Control_ss(pastA_g_ind);

l_lambda_ss = Control_ss(l_lambda_ind);
pastQ_ss    = Control_ss(pastQ_ind);
x_I_ss      = Control_ss(x_I_ind); 
eta2_ss     = Control_ss(eta2_ind); 
iota2_ss    = Control_ss(iota2_ind); 
pastLT2_ss  = Control_ss(pastLT2_ind); 

pastG2_ss   = Control_ss(pastG2_ind);
LT_obs_ss   = Control_ss(LT_obs_ind);

B_gov_ncp2_ss = Control_ss(B_gov_ncp2_ind);

w2_ss          = Control_ss(w2_ind);
Profit2_ss     = Control_ss(Profit2_ind);
Profit3_ss     = Control_ss(Profit3_ind);
Profit_FI2_ss  = Control_ss(Profit_FI2_ind);
r_a2_ss        = Control_ss(r_a2_ind);


tau_w_ss = param.tau_w;
tau_a_ss = param.tau_a;

F21_aux = zeros(grid.os,grid.numstates);
F22_aux = zeros(grid.os,grid.numcontrols);
F23_aux = zeros(grid.os,grid.numstates);
F24_aux = zeros(grid.os,grid.numcontrols);

F41_aux = zeros(grid.oc,grid.numstates);
F42_aux = zeros(grid.oc,grid.numcontrols);
F43_aux = zeros(grid.oc,grid.numstates);
F44_aux = zeros(grid.oc,grid.numcontrols);



% 1) State variables

% log(R_cb) = log(param.R_cb) + (1-param.rho_R)*(param.phi_pi*log(PI/param.pi_cb)-param.phi_y*(unemp-SS_stats.u)) + param.rho_R*(log(R_cbminus/param.R_cb)) + eps_R

% log(R_cb_aux) - ( log(param.R_cb) + (1-param.rho_R)*(param.phi_pi*log(pi_ss/param.pi_cb)-param.phi_y*(unemp_ss-SS_stats.u)) + param.rho_R*(log(R_cb_ss/param.R_cb)) + eps_R_ss )

% switch(param.regime)
% 	case('QE')
% 		log(R_cb)  = log(param.R_cb);
% 	case('Partial QE')
% 		log(R_cb)  = log(param.R_cb) + 1*log(PI/param.pi_cb);
% end

% (1) R_{t+1}

F21_aux(R_cb_ind,end-grid.os+R_cb_ind)  = (log(R_cb_aux)-log(R_cb_ss))/perturb_size;
F23_aux(R_cb_ind,end-grid.os+eps_R_ind) = -(eps_R_aux-eps_R_ss)/perturb_size;

F21_aux(R_star_ind,end-grid.os+R_star_ind) = (log(R_star_aux)-log(R_star_ss))/perturb_size;
F23_aux(R_star_ind,end-grid.os+R_star_ind) = (-param.rho_R*(log(R_star_aux/param.R_cb)-log(R_star_ss/param.R_cb)))./perturb_size;
F23_aux(R_star_ind,end-grid.os+eps_R_ind)  = -(eps_R_aux-eps_R_ss)/perturb_size;
F24_aux(R_star_ind,end-grid.oc+pi_ind)     = - ( (1-param.rho_R)*(param.phi_pi*(log(pi_aux/param.pi_cb) - log(pi_ss/param.pi_cb))) )/perturb_size;
F24_aux(R_star_ind,end-grid.oc+unemp_ind)  = - (1-param.rho_R)*(-param.phi_u*((unemp_aux-unemp_ss)))./perturb_size;

% (2) W_t

% log(W) = log(param.w_bar) + param.rho_w * log(Wminus/param.w_bar) + param.rho_w*(param.d*log(param.pi_bar/PI)+(1-param.d)*log(pastpiminus/PI))+ (1-param.rho_w)*param.eps_w*log(H/SS_stats.r_l) + eps_w;

F21_aux(w_ind,end-grid.os+w_ind)        = (log(w_aux)-log(w_ss))/perturb_size;
F21_aux(w_ind,end-grid.os+PSI_W_ind)    = - (1-param.rho_w)*param.eps_w*(log(1/PSI_W_aux)-log(1/PSI_W_ss))/perturb_size;
F23_aux(w_ind,end-grid.os+w_ind)        = - param.rho_w * (log(w_aux/param.w_bar) - log(w_ss/param.w_bar))/perturb_size;
F23_aux(w_ind,end-grid.os+pastpi_ind)   = - param.rho_w*(1-param.d)*(log(pi_aux/param.pi_bar)-log(pi_ss/param.pi_bar))/perturb_size;
F24_aux(w_ind,end-grid.oc+pi_ind)       = - param.rho_w*(param.d*log(param.pi_bar/pi_aux)+(1-param.d)*log(param.pi_bar/pi_aux))/perturb_size;
F24_aux(w_ind,end-grid.oc+h_ind)        = - (1-param.rho_w)*param.eps_w*(log(h_aux/SS_stats.r_l)-log(h_ss/SS_stats.r_l))/perturb_size;

% (3) A_t+1^b

% log(A_bnext) = log( lev*NW_b/Q );

F21_aux(A_b_ind,end-grid.os+A_b_ind)    = (log(A_b_aux)-log(A_b_ss))/perturb_size;
F21_aux(A_b_ind,end-grid.os+Q_ind)      = -(log(lev_ss*NW_b_ss/Q_aux)-log(lev_ss*NW_b_ss/Q_ss))/perturb_size;
F21_aux(A_b_ind,end-grid.os+lev_ind)    = -(log(lev_aux*NW_b_ss/Q_ss)-log(lev_ss*NW_b_ss/Q_ss))/perturb_size;
F21_aux(A_b_ind,end-grid.os+NW_b_ind)   = -(log(lev_ss*NW_b_aux/Q_ss)-log(lev_ss*NW_b_ss/Q_ss))/perturb_size;


% (4) B_t+1^b

% log(B_bnext) = log( param.MMF_ratio_1*PROFIT + param.MMF_ratio_2*Profit_FI + R_cbminus/PI*B_b - C_b );


F21_aux(B_b_ind,end-grid.os+B_b_ind)       = (log(B_b_aux)-log(B_b_ss))/perturb_size;

% log(B_b) - log( param.MMF_cont_ratio*T + R_cbminus/PI*B_b - C_b );

% log( T + B_gov_ncpnext - B_gov_ncp*R_cbminus/PI + (Q+R_A-R_cbminus/PI*Qminus)*A_g - UB - LT - param.tau_cp*Q*A_gnext - G + param.b_a_aux2*R_cbminus/PI*B_b - C_b );

F23_aux(B_b_ind,end-grid.os+B_b_ind)       = -( log( param.MMF_cont_ratio*T_ss  + param.b_a_aux2*R_cb_ss/pi_ss*B_b_aux - C_b_ss )  - log( param.MMF_cont_ratio*T_ss + param.b_a_aux2*R_cb_ss/pi_ss*B_b_ss - C_b_ss ))/perturb_size;
F24_aux(B_b_ind,end-grid.oc+C_b_ind)       = -( log( param.MMF_cont_ratio*T_ss  + param.b_a_aux2*R_cb_ss/pi_ss*B_b_ss  - C_b_aux ) - log( param.MMF_cont_ratio*T_ss + param.b_a_aux2*R_cb_ss/pi_ss*B_b_ss - C_b_ss ))/perturb_size;

F24_aux(B_b_ind,end-grid.oc+T_ind)         = -( log( T_aux + B_gov_ncp_ss - B_gov_ncp_ss*R_cb_ss/pi_ss + (Q_ss+r_a_ss-R_cb_ss/pi_ss*Q_ss)*(A_g_ss) - UB_ss - LT_ss - param.tau_cp*Q_ss*(A_g_ss) - G_ss + param.b_a_aux2*R_cb_ss/pi_ss*B_b_ss - C_b_ss )  - log( T_ss + B_gov_ncp_ss - B_gov_ncp_ss*R_cb_ss/pi_ss + (Q_ss+r_a_ss-R_cb_ss/pi_ss*Q_ss)*(A_g_ss) - UB_ss - LT_ss - param.tau_cp*Q_ss*(A_g_ss) - G_ss + param.b_a_aux2*R_cb_ss/pi_ss*B_b_ss - C_b_ss ))/perturb_size;
F23_aux(B_b_ind,end-grid.os+R_cb_ind)      = -( log( T_ss + B_gov_ncp_ss  - B_gov_ncp_ss*R_cb_aux/pi_ss + (Q_ss+r_a_ss-R_cb_aux/pi_ss*Q_ss)*(A_g_ss) - UB_ss - LT_ss - param.tau_cp*Q_ss*(A_g_ss) - G_ss + param.b_a_aux2*R_cb_aux/pi_ss*B_b_ss - C_b_ss )  - log( T_ss + B_gov_ncp_ss - B_gov_ncp_ss*R_cb_ss/pi_ss + (Q_ss+r_a_ss-R_cb_ss/pi_ss*Q_ss)*(A_g_ss) - UB_ss - LT_ss - param.tau_cp*Q_ss*(A_g_ss) - G_ss + param.b_a_aux2*R_cb_ss/pi_ss*B_b_ss - C_b_ss ))/perturb_size;
F24_aux(B_b_ind,end-grid.oc+pi_ind)        = -( log( T_ss + B_gov_ncp_ss  - B_gov_ncp_ss*R_cb_ss/pi_aux + (Q_ss+r_a_ss-R_cb_ss/pi_aux*Q_ss)*(A_g_ss) - UB_ss - LT_ss - param.tau_cp*Q_ss*(A_g_ss) - G_ss + param.b_a_aux2*R_cb_ss/pi_aux*B_b_ss - C_b_ss )  - log( T_ss + B_gov_ncp_ss - B_gov_ncp_ss*R_cb_ss/pi_ss + (Q_ss+r_a_ss-R_cb_ss/pi_ss*Q_ss)*(A_g_ss) - UB_ss - LT_ss - param.tau_cp*Q_ss*(A_g_ss) - G_ss + param.b_a_aux2*R_cb_ss/pi_ss*B_b_ss - C_b_ss ))/perturb_size;
F22_aux(B_b_ind,end-grid.oc+B_gov_ncp_ind) = -( log( T_ss + B_gov_ncp_aux - B_gov_ncp_ss*R_cb_ss/pi_ss + (Q_ss+r_a_ss-R_cb_ss/pi_ss*Q_ss)*(A_g_ss) - UB_ss - LT_ss - param.tau_cp*Q_ss*(A_g_ss) - G_ss + param.b_a_aux2*R_cb_ss/pi_ss*B_b_ss - C_b_ss )  - log( T_ss + B_gov_ncp_ss - B_gov_ncp_ss*R_cb_ss/pi_ss + (Q_ss+r_a_ss-R_cb_ss/pi_ss*Q_ss)*(A_g_ss) - UB_ss - LT_ss - param.tau_cp*Q_ss*(A_g_ss) - G_ss + param.b_a_aux2*R_cb_ss/pi_ss*B_b_ss - C_b_ss ))/perturb_size;
F24_aux(B_b_ind,end-grid.oc+B_gov_ncp_ind) = -( log( T_ss + B_gov_ncp_ss - B_gov_ncp_aux*R_cb_ss/pi_ss + (Q_ss+r_a_ss-R_cb_ss/pi_ss*Q_ss)*(A_g_ss) - UB_ss - LT_ss - param.tau_cp*Q_ss*(A_g_ss) - G_ss + param.b_a_aux2*R_cb_ss/pi_ss*B_b_ss - C_b_ss )  - log( T_ss + B_gov_ncp_ss - B_gov_ncp_ss*R_cb_ss/pi_ss + (Q_ss+r_a_ss-R_cb_ss/pi_ss*Q_ss)*(A_g_ss) - UB_ss - LT_ss - param.tau_cp*Q_ss*(A_g_ss) - G_ss + param.b_a_aux2*R_cb_ss/pi_ss*B_b_ss - C_b_ss ))/perturb_size;
F21_aux(B_b_ind,end-grid.os+Q_ind)         = -( log( T_ss + B_gov_ncp_ss - B_gov_ncp_ss*R_cb_ss/pi_ss + (Q_aux+r_a_ss-R_cb_ss/pi_ss*Q_ss)*(A_g_ss) - UB_ss - LT_ss - param.tau_cp*Q_aux*(A_g_ss) - G_ss + param.b_a_aux2*R_cb_ss/pi_ss*B_b_ss - C_b_ss )  - log( T_ss + B_gov_ncp_ss - B_gov_ncp_ss*R_cb_ss/pi_ss + (Q_ss+r_a_ss-R_cb_ss/pi_ss*Q_ss)*(A_g_ss) - UB_ss - LT_ss - param.tau_cp*Q_ss*(A_g_ss) - G_ss + param.b_a_aux2*R_cb_ss/pi_ss*B_b_ss - C_b_ss ))/perturb_size;
F23_aux(B_b_ind,end-grid.os+Q_ind)         = -( log( T_ss + B_gov_ncp_ss - B_gov_ncp_ss*R_cb_ss/pi_ss + (Q_ss+r_a_ss-R_cb_ss/pi_ss*Q_aux)*(A_g_ss) - UB_ss - LT_ss - param.tau_cp*Q_ss*(A_g_ss) - G_ss + param.b_a_aux2*R_cb_ss/pi_ss*B_b_ss - C_b_ss )  - log( T_ss + B_gov_ncp_ss - B_gov_ncp_ss*R_cb_ss/pi_ss + (Q_ss+r_a_ss-R_cb_ss/pi_ss*Q_ss)*(A_g_ss) - UB_ss - LT_ss - param.tau_cp*Q_ss*(A_g_ss) - G_ss + param.b_a_aux2*R_cb_ss/pi_ss*B_b_ss - C_b_ss ))/perturb_size;
F24_aux(B_b_ind,end-grid.oc+r_a_ind)       = -( log( T_ss + B_gov_ncp_ss - B_gov_ncp_ss*R_cb_ss/pi_ss + (Q_ss+r_a_aux-R_cb_ss/pi_ss*Q_ss)*(A_g_ss) - UB_ss - LT_ss - param.tau_cp*Q_ss*(A_g_ss) - G_ss + param.b_a_aux2*R_cb_ss/pi_ss*B_b_ss - C_b_ss )  - log( T_ss + B_gov_ncp_ss - B_gov_ncp_ss*R_cb_ss/pi_ss + (Q_ss+r_a_ss-R_cb_ss/pi_ss*Q_ss)*(A_g_ss) - UB_ss - LT_ss - param.tau_cp*Q_ss*(A_g_ss) - G_ss + param.b_a_aux2*R_cb_ss/pi_ss*B_b_ss - C_b_ss ))/perturb_size;
F23_aux(B_b_ind,end-grid.os+A_g_ind)       = -( log( T_ss + B_gov_ncp_ss - B_gov_ncp_ss*R_cb_ss/pi_ss + (Q_ss+r_a_ss-R_cb_ss/pi_ss*Q_ss)*(A_g_aux) - UB_ss - LT_ss - param.tau_cp*Q_ss*(A_g_ss) - G_ss + param.b_a_aux2*R_cb_ss/pi_ss*B_b_ss - C_b_ss )  - log( T_ss + B_gov_ncp_ss - B_gov_ncp_ss*R_cb_ss/pi_ss + (Q_ss+r_a_ss-R_cb_ss/pi_ss*Q_ss)*(A_g_ss) - UB_ss - LT_ss - param.tau_cp*Q_ss*(A_g_ss) - G_ss + param.b_a_aux2*R_cb_ss/pi_ss*B_b_ss - C_b_ss ))/perturb_size;
F21_aux(B_b_ind,end-grid.os+A_g_ind)       = -( log( T_ss + B_gov_ncp_ss - B_gov_ncp_ss*R_cb_ss/pi_ss + (Q_ss+r_a_ss-R_cb_ss/pi_ss*Q_ss)*(A_g_ss) - UB_ss - LT_ss - param.tau_cp*Q_ss*(A_g_aux) - G_ss + param.b_a_aux2*R_cb_ss/pi_ss*B_b_ss - C_b_ss )  - log( T_ss + B_gov_ncp_ss - B_gov_ncp_ss*R_cb_ss/pi_ss + (Q_ss+r_a_ss-R_cb_ss/pi_ss*Q_ss)*(A_g_ss) - UB_ss - LT_ss - param.tau_cp*Q_ss*(A_g_ss) - G_ss + param.b_a_aux2*R_cb_ss/pi_ss*B_b_ss - C_b_ss ))/perturb_size;
F24_aux(B_b_ind,end-grid.oc+UB_ind)        = -( log( T_ss + B_gov_ncp_ss - B_gov_ncp_ss*R_cb_ss/pi_ss + (Q_ss+r_a_ss-R_cb_ss/pi_ss*Q_ss)*(A_g_ss) - UB_aux - LT_ss - param.tau_cp*Q_ss*(A_g_ss) - G_ss + param.b_a_aux2*R_cb_ss/pi_ss*B_b_ss - C_b_ss )  - log( T_ss + B_gov_ncp_ss - B_gov_ncp_ss*R_cb_ss/pi_ss + (Q_ss+r_a_ss-R_cb_ss/pi_ss*Q_ss)*(A_g_ss) - UB_ss - LT_ss - param.tau_cp*Q_ss*(A_g_ss) - G_ss + param.b_a_aux2*R_cb_ss/pi_ss*B_b_ss - C_b_ss ))/perturb_size;
F24_aux(B_b_ind,end-grid.oc+LT_ind)        = -( log( T_ss + B_gov_ncp_ss - B_gov_ncp_ss*R_cb_ss/pi_ss + (Q_ss+r_a_ss-R_cb_ss/pi_ss*Q_ss)*(A_g_ss) - UB_ss - LT_aux - param.tau_cp*Q_ss*(A_g_ss) - G_ss + param.b_a_aux2*R_cb_ss/pi_ss*B_b_ss - C_b_ss )  - log( T_ss + B_gov_ncp_ss - B_gov_ncp_ss*R_cb_ss/pi_ss + (Q_ss+r_a_ss-R_cb_ss/pi_ss*Q_ss)*(A_g_ss) - UB_ss - LT_ss - param.tau_cp*Q_ss*(A_g_ss) - G_ss + param.b_a_aux2*R_cb_ss/pi_ss*B_b_ss - C_b_ss ))/perturb_size;
F24_aux(B_b_ind,end-grid.oc+G_ind)         = -( log( T_ss + B_gov_ncp_ss - B_gov_ncp_ss*R_cb_ss/pi_ss + (Q_ss+r_a_ss-R_cb_ss/pi_ss*Q_ss)*(A_g_ss) - UB_ss - LT_ss - param.tau_cp*Q_ss*(A_g_ss) - G_aux + param.b_a_aux2*R_cb_ss/pi_ss*B_b_ss - C_b_ss )  - log( T_ss + B_gov_ncp_ss - B_gov_ncp_ss*R_cb_ss/pi_ss + (Q_ss+r_a_ss-R_cb_ss/pi_ss*Q_ss)*(A_g_ss) - UB_ss - LT_ss - param.tau_cp*Q_ss*(A_g_ss) - G_ss + param.b_a_aux2*R_cb_ss/pi_ss*B_b_ss - C_b_ss ))/perturb_size;

% (5) A_t+1^g

% log(A_gauxnext) = log( param.rho_passive_QE*A_g + (1-param.rho_passive_QE)*SS_stats.A_g + 1 ) ;
% log(A_gauxnext) = log( param.rho_passive_QE*A_g + (1-param.rho_passive_QE)*SS_stats.A_g + eps_QE + 1 ) ;
% log(A_gauxnext) = log( x_cb*Qminus*A_g/Q + eps_R+ 1 ) ;
% 
% log(A_gauxnext) - log( ((param.rho_passive_QE + eps_QE)*Qminus*A_g + (1-param.rho_passive_QE)*SS_stats.A_g)/Q + 1 ) = 0 

F21_aux(A_g_ind,end-grid.os+A_g_ind)       = (log(A_g_aux)-log(A_g_ss))/perturb_size;

switch(param.QE)
	case('QE')
        F21_aux(A_g_ind,end-grid.os+x_cb_ind)   = - (log( x_cb_aux*Q_ss *(A_g_ss) /Q_ss  )-log( x_cb_ss*Q_ss*(A_g_ss)/Q_ss  ))/perturb_size;
	case('No QE')
		F23_aux(A_g_ind,end-grid.os+A_g_ind)      = -( log( ((param.rho_passive_QE)*(A_g_aux) + (1-param.rho_passive_QE)*SS_stats.A_g) + eps_QE_ss*Y_ss ) - log( ((param.rho_passive_QE)*(A_g_ss) + (1-param.rho_passive_QE)*SS_stats.A_g) )  )/perturb_size;
end

% (6) Q

% log(Q) = log( 1+param.phi*(x_k-1) - LAMBDA/SS_stats.Lambda*((param.phi/2*(x_knext^2-1)) + Q2next - 1 ) );
% log( SS_stats.Lambda/SS_stats.Lambda*((param.phi/2*((1+1e-5)^2-1^2)) )


% log(Q) - log( 1+gamma_Q*param.phi*(x_k-1) - LAMBDA/SS_stats.Lambda*gamma_Q2next*(param.phi/2*(x_knext^2-1)) );
% log(Q) - log( eta2*(1+param.phi*log(x_k)+param.phi/2*(log(x_k))^2 ) - eta2next*LAMBDA/SS_stats.Lambda*param.phi*log(x_knext)*x_knext );

F21_aux(Q_ind,end-grid.os+Q_ind)           = (log(Q_aux)-log(Q_ss))/perturb_size;

F22_aux(Q_ind,end-grid.oc+x_k_ind)         = - ( log( eta2_ss*(1+param.phi*log(x_k_ss)+param.phi/2*(log(x_k_ss))^2 ) - eta2_ss*Lambda_ss/SS_stats.Lambda*param.phi*log(x_k_aux)*x_k_aux ) - log( eta2_ss*(1+param.phi*log(x_k_ss)+param.phi/2*(log(x_k_ss))^2 ) - eta2_ss*Lambda_ss/SS_stats.Lambda*param.phi*log(x_k_ss)*x_k_ss ) )/perturb_size;
F22_aux(Q_ind,end-grid.oc+Lambda_ind)      = - ( log( eta2_ss*(1+param.phi*log(x_k_ss)+param.phi/2*(log(x_k_ss))^2 ) - eta2_ss*Lambda_aux/SS_stats.Lambda*param.phi*log(x_k_ss)*x_k_ss ) - log( eta2_ss*(1+param.phi*log(x_k_ss)+param.phi/2*(log(x_k_ss))^2 ) - eta2_ss*Lambda_ss/SS_stats.Lambda*param.phi*log(x_k_ss)*x_k_ss ) )/perturb_size;
F24_aux(Q_ind,end-grid.oc+x_k_ind)         = - ( log( eta2_ss*(1+param.phi*log(x_k_aux)+param.phi/2*(log(x_k_aux))^2 ) - eta2_ss*Lambda_ss/SS_stats.Lambda*param.phi*log(x_k_ss)*x_k_ss ) - log( eta2_ss*(1+param.phi*log(x_k_ss)+param.phi/2*(log(x_k_ss))^2 ) - eta2_ss*Lambda_ss/SS_stats.Lambda*param.phi*log(x_k_ss)*x_k_ss ) )/perturb_size;
F24_aux(Q_ind,end-grid.oc+iota2_ind)       = - ( log( iota2_aux*(1+param.phi*log(x_k_ss)+param.phi/2*(log(x_k_ss))^2 ) - iota2_ss*Lambda_ss/SS_stats.Lambda*param.phi*log(x_k_ss)*x_k_ss ) - log( iota2_ss*(1+param.phi*log(x_k_ss)+param.phi/2*(log(x_k_ss))^2 ) - iota2_ss*Lambda_ss/SS_stats.Lambda*param.phi*log(x_k_ss)*x_k_ss ) )/perturb_size;
F22_aux(Q_ind,end-grid.oc+iota2_ind)       = - ( log( iota2_ss*(1+param.phi*log(x_k_ss)+param.phi/2*(log(x_k_ss))^2 ) - iota2_aux*Lambda_ss/SS_stats.Lambda*param.phi*log(x_k_ss)*x_k_ss ) - log( iota2_ss*(1+param.phi*log(x_k_ss)+param.phi/2*(log(x_k_ss))^2 ) - iota2_ss*Lambda_ss/SS_stats.Lambda*param.phi*log(x_k_ss)*x_k_ss ) )/perturb_size;

% (7) Leverage (THETA)

% log(lev) = log( ee/(param.DELTA-vv) )

F21_aux(lev_ind,end-grid.os+lev_ind)       = (log(lev_aux)-log(lev_ss))/perturb_size;
F24_aux(lev_ind,end-grid.oc+vv_ind)        = - (log( ee_ss/(param.DELTA-vv_aux))-log( ee_ss/(param.DELTA-vv_ss) ))/perturb_size;
F24_aux(lev_ind,end-grid.oc+ee_ind)        = - (log( ee_aux/(param.DELTA-vv_ss))-log( ee_ss/(param.DELTA-vv_ss) ))/perturb_size;

% (8) Net worth

% log(NW_b) = log( (param.theta_b*((RRa-R_cbminus/PI)*levminus+R_cbminus/PI)*NW_bminus + param.omega*Qminus*A_b) )

F21_aux(NW_b_ind,end-grid.os+NW_b_ind)     = (log(NW_b_aux)-log(NW_b_ss))/perturb_size;
F23_aux(NW_b_ind,end-grid.os+R_cb_ind)     = - ( log( (param.theta_b*((RRa_ss -param.b_a_aux2*R_cb_aux/pi_ss)*lev_ss +param.b_a_aux2*R_cb_aux/pi_ss) *NW_b_ss  + param.omega*Q_ss *A_b_ss) ) - log( (param.theta_b*((RRa_ss-param.b_a_aux2*R_cb_ss/pi_ss)*lev_ss+param.b_a_aux2*R_cb_ss/pi_ss)*NW_b_ss + param.omega*Q_ss*A_b_ss) )   )/perturb_size;
F23_aux(NW_b_ind,end-grid.os+A_b_ind)      = - ( log( (param.theta_b*((RRa_ss -param.b_a_aux2*R_cb_ss/pi_ss) *lev_ss +param.b_a_aux2*R_cb_ss /pi_ss) *NW_b_ss  + param.omega*Q_ss *A_b_aux)) - log( (param.theta_b*((RRa_ss-param.b_a_aux2*R_cb_ss/pi_ss)*lev_ss+param.b_a_aux2*R_cb_ss/pi_ss)*NW_b_ss + param.omega*Q_ss*A_b_ss) )   )/perturb_size;
F23_aux(NW_b_ind,end-grid.os+Q_ind)        = - ( log( (param.theta_b*((RRa_ss -param.b_a_aux2*R_cb_ss/pi_ss) *lev_ss +param.b_a_aux2*R_cb_ss /pi_ss) *NW_b_ss  + param.omega*Q_aux*A_b_ss) ) - log( (param.theta_b*((RRa_ss-param.b_a_aux2*R_cb_ss/pi_ss)*lev_ss+param.b_a_aux2*R_cb_ss/pi_ss)*NW_b_ss + param.omega*Q_ss*A_b_ss) )   )/perturb_size;
F23_aux(NW_b_ind,end-grid.os+NW_b_ind)     = - ( log( (param.theta_b*((RRa_ss -param.b_a_aux2*R_cb_ss/pi_ss) *lev_ss +param.b_a_aux2*R_cb_ss /pi_ss) *NW_b_aux + param.omega*Q_ss *A_b_ss) ) - log( (param.theta_b*((RRa_ss-param.b_a_aux2*R_cb_ss/pi_ss)*lev_ss+param.b_a_aux2*R_cb_ss/pi_ss)*NW_b_ss + param.omega*Q_ss*A_b_ss) )   )/perturb_size;
F23_aux(NW_b_ind,end-grid.os+lev_ind)      = - ( log( (param.theta_b*((RRa_ss -param.b_a_aux2*R_cb_ss/pi_ss) *lev_aux+param.b_a_aux2*R_cb_ss /pi_ss) *NW_b_ss  + param.omega*Q_ss *A_b_ss) ) - log( (param.theta_b*((RRa_ss-param.b_a_aux2*R_cb_ss/pi_ss)*lev_ss+param.b_a_aux2*R_cb_ss/pi_ss)*NW_b_ss + param.omega*Q_ss*A_b_ss) )   )/perturb_size;
F24_aux(NW_b_ind,end-grid.oc+pi_ind)       = - ( log( (param.theta_b*((RRa_ss -param.b_a_aux2*R_cb_ss/pi_aux)*lev_ss +param.b_a_aux2*R_cb_ss /pi_aux)*NW_b_ss  + param.omega*Q_ss *A_b_ss) ) - log( (param.theta_b*((RRa_ss-param.b_a_aux2*R_cb_ss/pi_ss)*lev_ss+param.b_a_aux2*R_cb_ss/pi_ss)*NW_b_ss + param.omega*Q_ss*A_b_ss) )   )/perturb_size;
F24_aux(NW_b_ind,end-grid.oc+RRa_ind)      = - ( log( (param.theta_b*((RRa_aux-param.b_a_aux2*R_cb_ss/pi_ss) *lev_ss +param.b_a_aux2*R_cb_ss /pi_ss) *NW_b_ss  + param.omega*Q_ss *A_b_ss) ) - log( (param.theta_b*((RRa_ss-param.b_a_aux2*R_cb_ss/pi_ss)*lev_ss+param.b_a_aux2*R_cb_ss/pi_ss)*NW_b_ss + param.omega*Q_ss*A_b_ss) )   )/perturb_size;


% (9) Inventory

% log(inven) =  log( stock - sale )

% R_tilde

% log( R_tilde ) = log( R_cb*D)

F21_aux(R_tilde_ind,end-grid.os+R_tilde_ind) = (log(R_tilde_aux)-log(R_tilde_ss))/perturb_size;
F21_aux(R_tilde_ind,end-grid.os+R_cb_ind)    = - (log( R_cb_aux)-log(R_cb_ss))/perturb_size;

% (10) X_cb

% log(x_cb) = log( 1 +param.rho_X_QE*(x_cbminus-1) + (1-param.rho_X_QE)*(-param.phi_pi_QE*log(PI/param.pi_cb)+param.phi_u_QE*(unemp-SS_stats.u))) 

F21_aux(x_cb_ind,end-grid.os+x_cb_ind)     = (log(x_cb_aux)-log(x_cb_ss))/perturb_size;
F23_aux(x_cb_ind,end-grid.os+x_cb_ind)     = - (log( 1 +param.rho_X_QE*(x_cb_aux-1) + (1-param.rho_X_QE)*(-param.phi_pi_QE*log(pi_ss/param.pi_cb) +param.phi_u_QE*(unemp_ss-SS_stats.u))) -  log( 1 +param.rho_X_QE*(x_cb_ss-1) + (1-param.rho_X_QE)*(-param.phi_pi_QE*log(pi_ss/param.pi_cb)+param.phi_u_QE*(unemp_ss-SS_stats.u))) )/perturb_size;
F24_aux(x_cb_ind,end-grid.oc+pi_ind)       = - (log( 1 +param.rho_X_QE*(x_cb_ss-1)  + (1-param.rho_X_QE)*(-param.phi_pi_QE*log(pi_aux/param.pi_cb)+param.phi_u_QE*(unemp_ss-SS_stats.u))) -  log( 1 +param.rho_X_QE*(x_cb_ss-1) + (1-param.rho_X_QE)*(-param.phi_pi_QE*log(pi_ss/param.pi_cb)+param.phi_u_QE*(unemp_ss-SS_stats.u))) )/perturb_size;
F24_aux(x_cb_ind,end-grid.oc+unemp_ind)    = - (log( 1 +param.rho_X_QE*(x_cb_ss-1)  + (1-param.rho_X_QE)*(-param.phi_pi_QE*log(pi_ss/param.pi_cb) +param.phi_u_QE*(unemp_aux-SS_stats.u))) -  log( 1 +param.rho_X_QE*(x_cb_ss-1) + (1-param.rho_X_QE)*(-param.phi_pi_QE*log(pi_ss/param.pi_cb)+param.phi_u_QE*(unemp_ss-SS_stats.u))) )/perturb_size;
F21_aux(x_cb_ind,end-grid.os+PSI_RP_ind)   = - (log( 1 +param.rho_X_QE*(x_cb_ss-1)  + (1-param.rho_X_QE)*(-param.phi_pi_QE*log(pi_ss/param.pi_cb) +param.phi_u_QE*(unemp_ss-SS_stats.u))  + log(PSI_RP_aux) ) -  log( 1 +param.rho_X_QE*(x_cb_ss-1) + (1-param.rho_X_QE)*(-param.phi_pi_QE*log(pi_ss/param.pi_cb)+param.phi_u_QE*(unemp_ss-SS_stats.u))) )/perturb_size;


% switch(param.regime)
% 	case('QE')
% 		F23_aux(x_cb_ind,end-grid.os+eps_R_ind)   = - (log( 1 +param.rho_X_QE*(x_cb_ss-1)  + (1-param.rho_X_QE)*(-param.phi_pi_QE*log(pi_ss/param.pi_cb) +param.phi_u_QE*(unemp_ss-SS_stats.u))  + eps_R_aux ) -  log( 1 +param.rho_X_QE*(x_cb_ss-1) + (1-param.rho_X_QE)*(-param.phi_pi_QE*log(pi_ss/param.pi_cb)+param.phi_u_QE*(unemp_ss-SS_stats.u))) )/perturb_size;
% 	case('Partial QE')
% 		F23_aux(x_cb_ind,end-grid.os+eps_R_ind)   = - (log( 1 +param.rho_X_QE*(x_cb_ss-1)  + (1-param.rho_X_QE)*(-param.phi_pi_QE*log(pi_ss/param.pi_cb) +param.phi_u_QE*(unemp_ss-SS_stats.u))  + eps_R_aux ) -  log( 1 +param.rho_X_QE*(x_cb_ss-1) + (1-param.rho_X_QE)*(-param.phi_pi_QE*log(pi_ss/param.pi_cb)+param.phi_u_QE*(unemp_ss-SS_stats.u))) )/perturb_size;		
% end



% (11) past s2s

% log(pasts2s) = log(s2s)

% (12) past pi

% log(pastpi) = log(pi)

F21_aux(pastpi_ind,end-grid.os+pastpi_ind) = (log(pastpi_aux)-log(pastpi_ss))/perturb_size;
F24_aux(pastpi_ind,end-grid.oc+pi_ind)     = -(log(pi_aux)-log(pi_ss))/perturb_size;

F21_aux(pastY_ind,end-grid.os+pastY_ind) = (log(pastY_aux)-log(pastY_ss))/perturb_size;
F24_aux(pastY_ind,end-grid.oc+Y_ind)     = -(log(Y_aux)-log(Y_ss))/perturb_size;

F21_aux(pastC_ind,end-grid.os+pastC_ind) = (log(pastC_aux)-log(pastC_ss))/perturb_size;
F24_aux(pastC_ind,end-grid.oc+C_ind)     = -(log(C_aux)-log(C_ss))/perturb_size;

F21_aux(pastI_ind,end-grid.os+pastI_ind) = (log(pastI_aux)-log(pastI_ss))/perturb_size;
F24_aux(pastI_ind,end-grid.oc+I_ind)     = -(log(I_aux)-log(I_ss))/perturb_size;

F21_aux(pastPROFIT_ind,end-grid.os+pastPROFIT_ind) = (log(pastPROFIT_aux)-log(pastPROFIT_ss))/perturb_size;
F24_aux(pastPROFIT_ind,end-grid.oc+Profit_ind)     = -(log(Profit_aux+Profit_FI_ss)-log(Profit_ss+Profit_FI_ss))/perturb_size;

F21_aux(pastunemp_ind,end-grid.os+pastunemp_ind) = (log(pastunemp_aux)-log(pastunemp_ss))/perturb_size;
F24_aux(pastunemp_ind,end-grid.oc+unemp_ind)     = -(log(unemp_aux)-log(unemp_ss))/perturb_size;

F21_aux(pastG_ind,end-grid.os+pastG_ind) = (log(pastG_aux)-log(pastG_ss))/perturb_size;
F24_aux(pastG_ind,end-grid.oc+G_ind)     = -(log(G_aux)-log(G_ss))/perturb_size;

F21_aux(pastLT_ind,end-grid.os+pastLT_ind) = (log(pastLT_aux)-log(pastLT_ss))/perturb_size;
F24_aux(pastLT_ind,end-grid.oc+LT_ind)     = -(log(LT_aux)-log(LT_ss))/perturb_size;




% (14)~(22) shock processes

F21_aux(B_F_ind,end-grid.os+B_F_ind)         = 1;
F23_aux(B_F_ind,end-grid.os+B_F_ind)         = - param.rho_B_F*(log(B_F_aux)-log(B_F_ss))/perturb_size;
F23_aux(B_F_ind,end-grid.os+eps_B_F_ind)     = - (eps_B_F_aux-eps_B_F_ss)/perturb_size;

F21_aux(Z_ind,end-grid.os+Z_ind)           = 1;
F23_aux(Z_ind,end-grid.os+Z_ind)           = - param.rho_Z*(log(Z_aux)-log(Z_ss))/perturb_size;
F23_aux(Z_ind,end-grid.os+eps_Z_ind)       = - (eps_Z_aux-eps_Z_ss)/perturb_size;

F21_aux(PSI_RP_ind,end-grid.os+PSI_RP_ind)        = 1;
F23_aux(PSI_RP_ind,end-grid.os+PSI_RP_ind)        = - param.rho_PSI_RP*(log(PSI_RP_aux)-log(PSI_RP_ss))/perturb_size;
F23_aux(PSI_RP_ind,end-grid.os+eps_RP_ind)        = - (eps_RP_aux-eps_RP_ss)/perturb_size;


F21_aux(eta_ind,end-grid.os+eta_ind)           = 1;
F21_aux(eta_ind,end-grid.os+p_mark_ind)        = - ( log(p_mark_aux/(p_mark_aux-1)) - log(p_mark_ss/(p_mark_ss-1)) )/perturb_size;

F21_aux(D_ind,end-grid.os+D_ind)           = 1;
F23_aux(D_ind,end-grid.os+D_ind)           = - param.rho_D*(log(D_aux)-log(D_ss))/perturb_size;
F23_aux(D_ind,end-grid.os+eps_D_ind)        = - (eps_D_aux-eps_D_ss)/perturb_size;


F21_aux(GG_ind,end-grid.os+GG_ind)           = 1;

switch(param.adjust)
	case('G')
		F23_aux(GG_ind,end-grid.os+GG_ind)           = - param.rho_G *(log(GG_aux)-log(GG_ss))/perturb_size;
	case('LT')
		F23_aux(GG_ind,end-grid.os+GG_ind)           = - param.rho_G *(log(GG_aux)-log(GG_ss))/perturb_size;
end

F23_aux(GG_ind,end-grid.os+eps_G_ind)        = - (eps_G_aux-eps_G_ss)/perturb_size;



F21_aux(iota_ind,end-grid.os+iota_ind)      = 1;
F23_aux(iota_ind,end-grid.os+iota_ind)      = - param.rho_iota*(log(iota_aux)-log(iota_ss))/perturb_size;
F23_aux(iota_ind,end-grid.os+eps_iota_ind)  = - (eps_iota_aux-eps_iota_ss)/perturb_size;

F21_aux(BB_ind,end-grid.os+BB_ind)          = 1;
F23_aux(BB_ind,end-grid.os+BB_ind)          = - param.rho_BB*(log(BB_aux)-log(BB_ss))/perturb_size;
F23_aux(BB_ind,end-grid.os+eps_BB_ind)      = - (eps_BB_aux-eps_BB_ss)/perturb_size;

F21_aux(PSI_W_ind,end-grid.os+PSI_W_ind)    = 1;
F23_aux(PSI_W_ind,end-grid.os+PSI_W_ind)    = - param.rho_PSI_W*(log(PSI_W_aux)-log(PSI_W_ss))/perturb_size;
F23_aux(PSI_W_ind,end-grid.os+eps_w_ind)    = - (eps_w_aux-eps_w_ss)/perturb_size;

F21_aux(MP_ind,end-grid.os+MP_ind)          = 1;
F23_aux(MP_ind,end-grid.os+MP_ind)          = - param.rho_MP*(log(MP_aux)-log(MP_ss))/perturb_size;
F23_aux(MP_ind,end-grid.os+eps_R_ind)       = - (eps_R_aux-eps_R_ss)/perturb_size;

F21_aux(p_mark_ind,end-grid.os+p_mark_ind)      = 1;
F23_aux(p_mark_ind,end-grid.os+p_mark_ind)      = - param.rho_eta*(log(p_mark_aux)-log(p_mark_ss))/perturb_size;
F23_aux(p_mark_ind,end-grid.os+eps_eta_ind)     = - (eps_eta_aux-eps_eta_ss)/perturb_size;

F21_aux(pgap_ind,end-grid.os+pgap_ind)      = (log(pgap_aux)-log(pgap_ss))/perturb_size;
F23_aux(pgap_ind,end-grid.os+pgap_ind)      = - param.rho_pgap*(log(pgap_aux)-log(pgap_ss))/perturb_size;
F24_aux(pgap_ind,end-grid.oc+pi_ind)        = - (log(pi_aux)-log(pi_ss))/perturb_size;

F21_aux(ygap_ind,end-grid.os+ygap_ind)      = (log(ygap_aux)-log(ygap_ss))/perturb_size;
F23_aux(ygap_ind,end-grid.os+ygap_ind)      = - param.rho_ygap*(log(ygap_aux)-log(ygap_ss))/perturb_size;
F24_aux(ygap_ind,end-grid.oc+unemp_ind)     = - (SS_stats.u-unemp_aux)/perturb_size;


% (23) Shocks

F21_aux(eps_QE_ind   ,end-grid.os+eps_QE_ind)      = 1;
F21_aux(eps_B_F_ind  ,end-grid.os+eps_B_F_ind)     = 1;
F21_aux(eps_BB_ind   ,end-grid.os+eps_BB_ind)      = 1;
F21_aux(eps_RP_ind   ,end-grid.os+eps_RP_ind)      = 1;
F21_aux(eps_eta_ind  ,end-grid.os+eps_eta_ind)     = 1;
F21_aux(eps_D_ind    ,end-grid.os+eps_D_ind)       = 1;
F21_aux(eps_R_ind    ,end-grid.os+eps_R_ind)       = 1;
F21_aux(eps_Z_ind    ,end-grid.os+eps_Z_ind)       = 1;
F21_aux(eps_G_ind    ,end-grid.os+eps_G_ind)       = 1;
F21_aux(eps_iota_ind ,end-grid.os+eps_iota_ind) = 1;
F21_aux(eps_w_ind    ,end-grid.os+eps_w_ind)       = 1;


F43_aux(K_ind,end-grid.os+A_b_ind)  = -(A_b_aux-A_b_ss)/perturb_size;
F43_aux(K_ind,end-grid.os+A_g_ind)  = -(A_g_aux-A_g_ss)/perturb_size;
F44_aux(K_ind,end-grid.oc+A_hh_ind) = -(A_hh_aux-A_hh_ss)/perturb_size;
F44_aux(K_ind,end-grid.oc+K_ind)    = (K_aux-K_ss)/perturb_size;

% (2) B_t

% B_gov_ncp = B_hh + B_b  + SS_stats.B_F


F43_aux(B_ind,end-grid.os+B_b_ind)  = -(B_b_aux-B_b_ss)/perturb_size;
F44_aux(B_ind,end-grid.oc+B_hh_ind) = -(B_hh_aux-B_hh_ss)/perturb_size;
F44_aux(B_ind,end-grid.oc+B_ind)    = (B_aux-B_ss)/perturb_size;

% (3) B_gov_ncp

% B_gov_ncp = B_b  + B_hh + SS_stats.B_F - (Qminus*A_b- NW_bminus) - Qminus *(1+param.tau_cp)* A_g;


F43_aux(B_gov_ncp_ind,end-grid.os+A_b_ind)       = - ( (B_b_ss + B_hh_ss + SS_stats.B_F  - (Q_ss*A_b_aux- NW_b_ss) - Q_ss *(1+param.tau_cp)* (A_g_ss)) - ( (B_b_ss + B_hh_ss + SS_stats.B_F - (Q_ss*A_b_ss- NW_b_ss) - Q_ss *(1+param.tau_cp)* (A_g_ss))  )  )/perturb_size; 
F43_aux(B_gov_ncp_ind,end-grid.os+B_b_ind)       = - ( (B_b_aux+ B_hh_ss + SS_stats.B_F - (Q_ss*A_b_ss- NW_b_ss)  - Q_ss *(1+param.tau_cp)* (A_g_ss)) - ( (B_b_ss + B_hh_ss + SS_stats.B_F - (Q_ss*A_b_ss- NW_b_ss) - Q_ss *(1+param.tau_cp)* (A_g_ss))  )  )/perturb_size;
F43_aux(B_gov_ncp_ind,end-grid.os+A_g_ind)       = - ( (B_b_ss + B_hh_ss + SS_stats.B_F  - (Q_ss*A_b_ss- NW_b_ss)  - Q_ss *(1+param.tau_cp)* (A_g_aux)) - ( (B_b_ss + B_hh_ss + SS_stats.B_F - (Q_ss*A_b_ss- NW_b_ss) - Q_ss *(1+param.tau_cp)* (A_g_ss))  )  )/perturb_size;
F43_aux(B_gov_ncp_ind,end-grid.os+Q_ind)         = - ( (B_b_ss + B_hh_ss + SS_stats.B_F  - (Q_aux*A_b_ss- NW_b_ss) - Q_aux*(1+param.tau_cp)* (A_g_ss)) - ( (B_b_ss + B_hh_ss + SS_stats.B_F - (Q_ss*A_b_ss- NW_b_ss) - Q_ss *(1+param.tau_cp)* (A_g_ss))  )  )/perturb_size;
F43_aux(B_gov_ncp_ind,end-grid.os+NW_b_ind)      = - ( (B_b_ss + B_hh_ss + SS_stats.B_F  - (Q_ss*A_b_ss- NW_b_aux) - Q_ss *(1+param.tau_cp)* (A_g_ss)) - ( (B_b_ss + B_hh_ss + SS_stats.B_F - (Q_ss*A_b_ss- NW_b_ss) - Q_ss *(1+param.tau_cp)* (A_g_ss))  )  )/perturb_size;
F44_aux(B_gov_ncp_ind,end-grid.oc+B_hh_ind)      = - (B_hh_aux-B_hh_ss)/perturb_size;
F44_aux(B_gov_ncp_ind,end-grid.oc+B_gov_ncp_ind) = (B_gov_ncp_aux-B_gov_ncp_ss)/perturb_size;

% (4) Tax

% T = tau_w.*(W.*L + UB) + tau_a.*(1-param.MMF_ratio_1)*PROFIT + tau_a.*param.Eratio_2*Profit_FI;
% T - (tau_w.*(W.*L + UB) + tau_a.*(PROFIT + Profit_FI - param.fix2))


F41_aux(T_ind,end-grid.os+w_ind)         = - ((tau_w_ss.*(w_aux.*L_ss + UB_ss) + tau_a_ss.*Profit_ss + tau_a_ss.*(Profit_FI_ss-param.fix2))-(tau_w_ss.*(w_ss.*L_ss + UB_ss) + tau_a_ss.*Profit_ss + tau_a_ss.*(Profit_FI_ss-param.fix2)))/perturb_size;
F44_aux(T_ind,end-grid.oc+T_ind)         = (T_aux-T_ss)/perturb_size;
F44_aux(T_ind,end-grid.oc+L_ind)         = - ((tau_w_ss.*(w_ss.*L_aux + UB_ss) + tau_a_ss.*Profit_ss + tau_a_ss.*(Profit_FI_ss-param.fix2))-(tau_w_ss.*(w_ss.*L_ss + UB_ss) + tau_a_ss.*Profit_ss + tau_a_ss.*(Profit_FI_ss-param.fix2)))/perturb_size;
F44_aux(T_ind,end-grid.oc+UB_ind)        = - ((tau_w_ss.*(w_ss.*L_ss + UB_aux) + tau_a_ss.*Profit_ss + tau_a_ss.*(Profit_FI_ss-param.fix2))-(tau_w_ss.*(w_ss.*L_ss + UB_ss) + tau_a_ss.*Profit_ss + tau_a_ss.*(Profit_FI_ss-param.fix2)))/perturb_size;
F44_aux(T_ind,end-grid.oc+Profit_ind)    = - ((tau_w_ss.*(w_ss.*L_ss + UB_ss) + tau_a_ss.*Profit_aux + tau_a_ss.*(Profit_FI_ss-param.fix2))-(tau_w_ss.*(w_ss.*L_ss + UB_ss) + tau_a_ss.*Profit_ss + tau_a_ss.*(Profit_FI_ss-param.fix2)))/perturb_size;

% (5-6) LT and G

switch(param.adjust)
	case('G')

		% LT = GG - (C_b-param.fix2);
		% G  = T + B_gov_ncpnext - B_gov_ncp*R_cbminus/PI + (Q+R_A-R_cbminus/PI*Qminus)*A_g - UB - LT - param.tau_cp*Q*A_gnext;        		           


		F41_aux(LT_ind,end-grid.os+GG_ind)  = -(-1/GG_aux+1/GG_ss)*Y_ss/perturb_size;
		F44_aux(LT_ind,end-grid.oc+LT_ind)  =  (LT_aux-LT_ss) /perturb_size;

		F44_aux(G_ind,end-grid.oc+G_ind)         = (G_aux-G_ss)/perturb_size;
		F44_aux(G_ind,end-grid.oc+C_b_ind) = -(C_b_aux-C_b_ss)/perturb_size;

    case('LT')

    	% G  = SS_stats.G;
		% LT = T + B_gov_ncpnext - B_gov_ncp*R_cbminus/PI + (Q+R_A-R_cbminus/PI*Qminus)*A_g - UB - G - param.tau_cp*Q*A_gnext;        		           

    	F41_aux(G_ind,end-grid.os+GG_ind)  = -(-1/GG_aux+1/GG_ss)*Y_ss/perturb_size;
    	F44_aux(G_ind,end-grid.oc+G_ind)   = (G_aux-G_ss)/perturb_size;

		F44_aux(LT_ind,end-grid.oc+LT_ind)        = (LT_aux-LT_ss)/perturb_size;
end

% (7) Lambda

% Lambda = MRS*param.Lambda_aux/param.Lambda_b_aux

F44_aux(Lambda_ind,end-grid.oc+Lambda_ind)        = (Lambda_aux - Lambda_ss)/perturb_size;
F44_aux(Lambda_ind,end-grid.oc+MRS_ind)           = -(MRS_aux*param.Lambda_aux/param.Lambda_b_aux-MRS_ss*param.Lambda_aux/param.Lambda_b_aux)/perturb_size;

% (8) pi

% pi = param.pi_cb*exp((param.rho_B/(param.rho_B+param.gamma_pi))*log(B_gov_ncp*R_cbminus/(SS_stats.B_gov_ncp*param.R_cb))-(param.gamma_T/(param.rho_B+param.gamma_pi))*log(T/SS_stats.T)-1/(param.rho_B+param.gamma_pi)*log(B_gov_ncpnext/SS_stats.B_gov_ncp))

F44_aux(pi_ind,end-grid.oc+pi_ind)                = (pi_aux-pi_ss)/perturb_size;

F42_aux(pi_ind,end-grid.oc+B_gov_ncp_ind)         = - (param.pi_cb*exp((param.rho_B/(param.rho_B+param.gamma_pi))*log(B_gov_ncp_ss*R_cb_ss/(SS_stats.B_gov_ncp*param.R_cb))-(param.gamma_T/(param.rho_B+param.gamma_pi))*log(T_ss/SS_stats.T)-1/(param.rho_B+param.gamma_pi)*log(B_gov_ncp_aux/SS_stats.B_gov_ncp)) - ...
													param.pi_cb*exp((param.rho_B/(param.rho_B+param.gamma_pi))*log(B_gov_ncp_ss*R_cb_ss/(SS_stats.B_gov_ncp*param.R_cb))-(param.gamma_T/(param.rho_B+param.gamma_pi))*log(T_ss/SS_stats.T)-1/(param.rho_B+param.gamma_pi)*log(B_gov_ncp_ss/SS_stats.B_gov_ncp)))/perturb_size;
F43_aux(pi_ind,end-grid.os+R_cb_ind)              = - (param.pi_cb*exp((param.rho_B/(param.rho_B+param.gamma_pi))*log(B_gov_ncp_ss*R_cb_aux/(SS_stats.B_gov_ncp*param.R_cb))-(param.gamma_T/(param.rho_B+param.gamma_pi))*log(T_ss/SS_stats.T)-1/(param.rho_B+param.gamma_pi)*log(B_gov_ncp_ss/SS_stats.B_gov_ncp)) - ...
													param.pi_cb*exp((param.rho_B/(param.rho_B+param.gamma_pi))*log(B_gov_ncp_ss*R_cb_ss/(SS_stats.B_gov_ncp*param.R_cb))-(param.gamma_T/(param.rho_B+param.gamma_pi))*log(T_ss/SS_stats.T)-1/(param.rho_B+param.gamma_pi)*log(B_gov_ncp_ss/SS_stats.B_gov_ncp)))/perturb_size;
F44_aux(pi_ind,end-grid.oc+B_gov_ncp_ind)         = - (param.pi_cb*exp((param.rho_B/(param.rho_B+param.gamma_pi))*log(B_gov_ncp_aux*R_cb_ss/(SS_stats.B_gov_ncp*param.R_cb))-(param.gamma_T/(param.rho_B+param.gamma_pi))*log(T_ss/SS_stats.T)-1/(param.rho_B+param.gamma_pi)*log(B_gov_ncp_ss/SS_stats.B_gov_ncp)) - ...
													param.pi_cb*exp((param.rho_B/(param.rho_B+param.gamma_pi))*log(B_gov_ncp_ss*R_cb_ss/(SS_stats.B_gov_ncp*param.R_cb))-(param.gamma_T/(param.rho_B+param.gamma_pi))*log(T_ss/SS_stats.T)-1/(param.rho_B+param.gamma_pi)*log(B_gov_ncp_ss/SS_stats.B_gov_ncp)))/perturb_size;
F44_aux(pi_ind,end-grid.oc+T_ind)                 = - (param.pi_cb*exp((param.rho_B/(param.rho_B+param.gamma_pi))*log(B_gov_ncp_ss*R_cb_ss/(SS_stats.B_gov_ncp*param.R_cb))-(param.gamma_T/(param.rho_B+param.gamma_pi))*log(T_aux/SS_stats.T)-1/(param.rho_B+param.gamma_pi)*log(B_gov_ncp_ss/SS_stats.B_gov_ncp)) - ...
													param.pi_cb*exp((param.rho_B/(param.rho_B+param.gamma_pi))*log(B_gov_ncp_ss*R_cb_ss/(SS_stats.B_gov_ncp*param.R_cb))-(param.gamma_T/(param.rho_B+param.gamma_pi))*log(T_ss/SS_stats.T)-1/(param.rho_B+param.gamma_pi)*log(B_gov_ncp2_ss/SS_stats.B_gov_ncp2)))/perturb_size;


% param.pi_cb*exp(param.gamma_B/param.gamma_pi*log(B_gov_ncp_ss/B_gov_ncp_ss)+1/param.gamma_pi*log(B_gov_ncp_ss/B_gov_ncp_ss)-param.gamma_Y/param.gamma_pi*log(Y_ss/Y_ss));


% (9) Vacancies

% V = M/param.iota*(J'*grid.s_dist)

F44_aux(V_ind,end-grid.oc+V_ind)                  = (V_aux-V_ss)/perturb_size;
F44_aux(V_ind,end-grid.oc+M_ind)                  = - ( M_aux/param.iota*(J_ss'*grid.s_dist) - M_ss/param.iota*(J_ss'*grid.s_dist)  )/perturb_size;

for iii = 1:grid.ns

	J_auxs{iii}      = J_ss;
	J_auxs{iii}(iii) = J_aux(iii);

	F44_aux(V_ind,end-grid.oc+J_ind(iii))               = - (M_ss/param.iota*(J_auxs{iii}'*grid.s_dist)-M_ss/param.iota*(J_ss'*grid.s_dist))/perturb_size;

end


% (10-14) J

% J = (H-W-param.fix_L).*(grid.s')*nn + LAMBDA*(1-param.death_rate)*(1-param.lambda)*(1-param.in)*param.P_SS*Jnext

for iii = 1:grid.ns


	F41_aux(J_ind(iii),end-grid.os+w_ind)      = (w_aux-w_ss)*grid.s(iii)*nn_ss/perturb_size;
	F44_aux(J_ind(iii),end-grid.oc+J_ind(iii))   = (J_auxs{iii}(iii)-J_ss(iii))/perturb_size;
	F44_aux(J_ind(iii),end-grid.oc+h_ind)        = - (h_aux-h_ss)*grid.s(iii)*nn_ss/perturb_size;
	F44_aux(J_ind(iii),end-grid.oc+nn_ind)       = - (h_ss-param.fix_L-w_ss).*grid.s(iii)*(nn_aux-nn_ss)/perturb_size;
	F44_aux(J_ind(iii),end-grid.oc+Lambda_ind)   = - (Lambda_aux*(1-param.death_rate)*(1-param.lambda)*(1-param.in)*param.P_SS(iii,:)*J_ss - Lambda_ss*(1-param.death_rate)*(1-param.lambda)*(1-param.in)*param.P_SS(iii,:)*J_ss )/perturb_size;
	F42_aux(J_ind(iii),end-grid.oc+l_lambda_ind) = - (Lambda_ss*(1-param.death_rate)*(1-l_lambda_aux)*(1-param.in)*param.P_SS(iii,:)*J_ss - Lambda_ss*(1-param.death_rate)*(1-param.lambda)*(1-param.in)*param.P_SS(iii,:)*J_ss )/perturb_size;

	for jjj  = 1:grid.ns

		F42_aux(J_ind(iii),end-grid.oc+J_ind(jjj)) = - (Lambda_ss*(1-param.death_rate)*(1-param.lambda)*(1-param.in)*param.P_SS(iii,:)*J_auxs{jjj} - Lambda_ss*(1-param.death_rate)*(1-param.lambda)*(1-param.in)*param.P_SS(iii,:)*J_ss )/perturb_size;

	end

end

% (15) r^l

% r_l = Z*MC*(1-param.theta)*(v*K).^param.theta.*(L).^(-param.theta)

% Z*MC*(v*K).^param.theta.*(L).^(param.theta_2-1)*(param.theta_2+param.alpha_lk*log(v*K/SS_stats.v/grid.K)+param.alpha_ll*log(L/SS_stats.L)); 


F41_aux(h_ind,end-grid.os+Z_ind)       = - ( Z_aux*MC_ss*param.theta_2*(v_ss*K_ss).^param.theta.*(L_ss).^(param.theta_2-1) - Z_ss*MC_ss*param.theta_2*(v_ss*K_ss).^param.theta.*(L_ss).^(param.theta_2-1)   )/perturb_size;
F44_aux(h_ind,end-grid.oc+h_ind)        = (h_aux-h_ss)/perturb_size;
F44_aux(h_ind,end-grid.oc+MC_ind)       = - ( Z_ss*MC_aux*param.theta_2*(v_ss *K_ss) .^param.theta.*(L_ss) .^(param.theta_2-1) - Z_ss*MC_ss*param.theta_2*(v_ss*K_ss).^param.theta.*(L_ss).^(param.theta_2-1)   )/perturb_size;
F44_aux(h_ind,end-grid.oc+v_ind)        = - ( Z_ss*MC_ss *(v_aux*K_ss) .^param.theta.*(L_ss) .^(param.theta_2-1)*(param.theta_2+param.alpha_lk*log(v_aux*K_ss/SS_stats.v/grid.K)+param.alpha_ll*log(L_ss /SS_stats.L)) - Z_ss*MC_ss*param.theta_2*(v_ss*K_ss).^param.theta.*(L_ss).^(param.theta_2-1)   )/perturb_size;
F44_aux(h_ind,end-grid.oc+K_ind)        = - ( Z_ss*MC_ss *(v_ss *K_aux).^param.theta.*(L_ss) .^(param.theta_2-1)*(param.theta_2+param.alpha_lk*log(v_ss*K_aux/SS_stats.v/grid.K)+param.alpha_ll*log(L_ss /SS_stats.L)) - Z_ss*MC_ss*param.theta_2*(v_ss*K_ss).^param.theta.*(L_ss).^(param.theta_2-1)   )/perturb_size;
F44_aux(h_ind,end-grid.oc+L_ind)        = - ( Z_ss*MC_ss *(v_ss *K_ss) .^param.theta.*(L_aux).^(param.theta_2-1)*(param.theta_2+param.alpha_lk*log(v_ss *K_ss/SS_stats.v/grid.K)+param.alpha_ll*log(L_aux/SS_stats.L)) - Z_ss*MC_ss*param.theta_2*(v_ss*K_ss).^param.theta.*(L_ss).^(param.theta_2-1)   )/perturb_size;

% (16) v

% v = (R_K/(param.delta_0*param.delta_1))^(1/(param.delta_1-1))

% v = SS_stats.v + (R_K-param.delta_1)/param.delta_2


F44_aux(v_ind,end-grid.oc+v_ind)        = (v_aux-v_ss)/perturb_size;
F44_aux(v_ind,end-grid.oc+r_k_ind)      = - ((r_k_aux/(param.delta_0*param.delta_1))^(1/(param.delta_1-1)) - (r_k_ss/(param.delta_0*param.delta_1))^(1/(param.delta_1-1)) )/perturb_size;

% (17) Y

% Y = Z*(v*K).^(param.theta).*(L).^(1-param.theta)


F41_aux(Y_ind,end-grid.os+Z_ind)       = - (Z_aux*(v_ss*K_ss).^(param.theta).*(L_ss).^param.theta_2 - Z_ss*(v_ss*K_ss).^(param.theta).*(L_ss).^param.theta_2)/perturb_size;
F44_aux(Y_ind,end-grid.oc+Y_ind)        = (Y_aux-Y_ss)/perturb_size;
F44_aux(Y_ind,end-grid.oc+v_ind)        = - (Z_ss*(v_aux*K_ss).^(param.theta).*(L_ss).^param.theta_2 - Z_ss*(v_ss*K_ss).^(param.theta).*(L_ss).^param.theta_2)/perturb_size;
F44_aux(Y_ind,end-grid.oc+K_ind)        = - (Z_ss*(v_ss*K_aux).^(param.theta).*(L_ss).^param.theta_2 - Z_ss*(v_ss*K_ss).^(param.theta).*(L_ss).^param.theta_2)/perturb_size;
F44_aux(Y_ind,end-grid.oc+L_ind)        = - (Z_ss*(v_ss*K_ss).^(param.theta).*(L_aux).^param.theta_2 - Z_ss*(v_ss*K_ss).^(param.theta).*(L_ss).^param.theta_2)/perturb_size;

% (18) Profit

% Profit = (sale *(1-param.eta/(2*param.kappa)*(log(PI)-(1-param.GAMMA)*log(param.pi_cb)-param.GAMMA*log(pastpiminus)).^2)+ Y.*(-MC) - param.fix  ...
% 		                - param.phi_x/2*(log(s2s)-log(pasts2sminus))^2 * sale + (H-W-param.fix_L).*L - iota.*V)  ...
% 		                + (R_K*v - param.delta_0.*v.^param.delta_1).*K +  Q*(Knext-K)-(Knext-K) - param.phi/2*(Knext/K-1)^2*K  ; 		 		 

% ((sale *(1-param.eta/(2*param.kappa)*(log(PI)-(1-param.GAMMA)*log(param.pi_cb)-param.GAMMA*log(pastpiminus)).^2)+ Y.*(-MC) - param.fix  ...
% 		                - param.phi_x/2*(log(s2s)-log(pasts2sminus))^2 * sale + (H-W-param.fix_L).*L - param.iota.*V)  ...
% 		                + (R_K*v - param.delta_0.*v.^param.delta_1).*K +  Q*(Knext-K)-(Knext-K) - param.phi/2*(Knext/K-1)^2*K )



F41_aux(Profit_ind,end-grid.os+eta_ind)       = - ( ((Y_ss *(1-eta_aux/(2*param.kappa)*(log(pi_ss)-(1-param.GAMMA)*log(param.pi_cb)-param.GAMMA*log(pastpi_ss)).^2)+ Y_ss.*(-MC_ss) - param.fix  ...
		                						    + (h_ss-w_ss).*L_ss - param.iota.*V_ss)  ...
		                							+ (r_k_ss*v_ss - param.delta_0*v_ss^param.delta_1).*K_ss +  Q_ss*(K_ss-K_ss)-(K_ss-K_ss) - param.phi/2*(K_ss/K_ss-1)^2*K_ss ) - ...
													((Y_ss *(1-param.eta/(2*param.kappa)*(log(pi_ss)-(1-param.GAMMA)*log(param.pi_cb)-param.GAMMA*log(pastpi_ss)).^2)+ Y_ss.*(-MC_ss) - param.fix  ...
		                						    + (h_ss-w_ss).*L_ss - param.iota.*V_ss)  ...
		                							+ (r_k_ss*v_ss - param.delta_0*v_ss^param.delta_1).*K_ss +  Q_ss*(K_ss-K_ss)-(K_ss-K_ss) - param.phi/2*(K_ss/K_ss-1)^2*K_ss ) )/perturb_size;

F41_aux(Profit_ind,end-grid.os+B_F_ind)       = - ( ((Y_ss *(1-param.eta/(2*param.kappa)*(log(pi_ss)-(1-param.GAMMA)*log(param.pi_cb)-param.GAMMA*log(pastpi_ss)).^2)+ Y_ss.*(-MC_ss) - (param.fix_ratio + log(B_F_aux))*Y_ss   ...
		                						    + (h_ss-w_ss).*L_ss - param.iota.*V_ss)  ...
		                							+ (r_k_ss*v_ss - param.delta_0*v_ss^param.delta_1).*K_ss +  Q_ss*(K_ss-K_ss)-(K_ss-K_ss) - param.phi/2*(K_ss/K_ss-1)^2*K_ss ) - ...
													((Y_ss *(1-param.eta/(2*param.kappa)*(log(pi_ss)-(1-param.GAMMA)*log(param.pi_cb)-param.GAMMA*log(pastpi_ss)).^2)+ Y_ss.*(-MC_ss) - param.fix  ...
		                						    + (h_ss-w_ss).*L_ss - param.iota.*V_ss)  ...
		                							+ (r_k_ss*v_ss - param.delta_0*v_ss^param.delta_1).*K_ss +  Q_ss*(K_ss-K_ss)-(K_ss-K_ss) - param.phi/2*(K_ss/K_ss-1)^2*K_ss ) )/perturb_size;
F41_aux(Profit_ind,end-grid.os+w_ind)         = - ( ((Y_ss *(1-param.eta/(2*param.kappa)*(log(pi_ss)-(1-param.GAMMA)*log(param.pi_cb)-param.GAMMA*log(pastpi_ss)).^2)+ Y_ss.*(-MC_ss) - param.fix  ...
		                						    + (h_ss-w_aux).*L_ss - param.iota.*V_ss)  ...
		                							+ (r_k_ss*v_ss - param.delta_0*v_ss^param.delta_1).*K_ss +  Q_ss*(K_ss-K_ss)-(K_ss-K_ss) - param.phi/2*(K_ss/K_ss-1)^2*K_ss ) - ...
													((Y_ss *(1-param.eta/(2*param.kappa)*(log(pi_ss)-(1-param.GAMMA)*log(param.pi_cb)-param.GAMMA*log(pastpi_ss)).^2)+ Y_ss.*(-MC_ss) - param.fix  ...
		                						    + (h_ss-w_ss).*L_ss - param.iota.*V_ss)  ...
		                							+ (r_k_ss*v_ss - param.delta_0*v_ss^param.delta_1).*K_ss +  Q_ss*(K_ss-K_ss)-(K_ss-K_ss) - param.phi/2*(K_ss/K_ss-1)^2*K_ss ) )/perturb_size;
F41_aux(Profit_ind,end-grid.os+Q_ind)         = - ( ((Y_ss *(1-param.eta/(2*param.kappa)*(log(pi_ss)-(1-param.GAMMA)*log(param.pi_cb)-param.GAMMA*log(pastpi_ss)).^2)+ Y_ss.*(-MC_ss) - param.fix  ...
		                						    + (h_ss-w_ss).*L_ss - param.iota.*V_ss)  ...
		                							+ (r_k_ss*v_ss - param.delta_0*v_ss^param.delta_1).*K_ss +  Q_aux*(K_ss-K_ss)-(K_ss-K_ss) - param.phi/2*(K_ss/K_ss-1)^2*K_ss ) - ...
													((Y_ss *(1-param.eta/(2*param.kappa)*(log(pi_ss)-(1-param.GAMMA)*log(param.pi_cb)-param.GAMMA*log(pastpi_ss)).^2)+ Y_ss.*(-MC_ss) - param.fix  ...
		                						    + (h_ss-w_ss).*L_ss - param.iota.*V_ss)  ...
		                							+ (r_k_ss*v_ss - param.delta_0*v_ss^param.delta_1).*K_ss +  Q_ss*(K_ss-K_ss)-(K_ss-K_ss) - param.phi/2*(K_ss/K_ss-1)^2*K_ss ) )/perturb_size;
F42_aux(Profit_ind,end-grid.oc+K_ind)         = - ( ((Y_ss *(1-param.eta/(2*param.kappa)*(log(pi_ss)-(1-param.GAMMA)*log(param.pi_cb)-param.GAMMA*log(pastpi_ss)).^2)+ Y_ss.*(-MC_ss) - param.fix  ...
		                						    + (h_ss-w_ss).*L_ss - param.iota.*V_ss)  ...
		                							+ (r_k_ss*v_ss - param.delta_0*v_ss^param.delta_1).*K_ss +  Q_ss*(K_aux-K_ss)-(K_aux-K_ss) - param.phi/2*(K_aux/K_ss-1)^2*K_ss ) - ...
													((Y_ss *(1-param.eta/(2*param.kappa)*(log(pi_ss)-(1-param.GAMMA)*log(param.pi_cb)-param.GAMMA*log(pastpi_ss)).^2)+ Y_ss.*(-MC_ss) - param.fix  ...
		                						    + (h_ss-w_ss).*L_ss - param.iota.*V_ss)  ...
		                							+ (r_k_ss*v_ss - param.delta_0*v_ss^param.delta_1).*K_ss +  Q_ss*(K_ss-K_ss)-(K_ss-K_ss) - param.phi/2*(K_ss/K_ss-1)^2*K_ss ) )/perturb_size;
F43_aux(Profit_ind,end-grid.os+pastpi_ind)    = - ( ((Y_ss *(1-param.eta/(2*param.kappa)*(log(pi_ss)-(1-param.GAMMA)*log(param.pi_cb)-param.GAMMA*log(pastpi_aux)).^2)+ Y_ss.*(-MC_ss) - param.fix  ...
		                						    + (h_ss-w_ss).*L_ss - param.iota.*V_ss)  ...
		                							+ (r_k_ss*v_ss - param.delta_0*v_ss^param.delta_1).*K_ss +  Q_ss*(K_ss-K_ss)-(K_ss-K_ss) - param.phi/2*(K_ss/K_ss-1)^2*K_ss ) - ...
													((Y_ss *(1-param.eta/(2*param.kappa)*(log(pi_ss)-(1-param.GAMMA)*log(param.pi_cb)-param.GAMMA*log(pastpi_ss)).^2)+ Y_ss.*(-MC_ss) - param.fix  ...
		                						    + (h_ss-w_ss).*L_ss - param.iota.*V_ss)  ...
		                							+ (r_k_ss*v_ss - param.delta_0*v_ss^param.delta_1).*K_ss +  Q_ss*(K_ss-K_ss)-(K_ss-K_ss) - param.phi/2*(K_ss/K_ss-1)^2*K_ss ) )/perturb_size;
F44_aux(Profit_ind,end-grid.oc+Profit_ind)    = (Profit_aux-Profit_ss)/perturb_size;
F44_aux(Profit_ind,end-grid.oc+pi_ind)        = - ( ((Y_ss *(1-param.eta/(2*param.kappa)*(log(pi_aux)-(1-param.GAMMA)*log(param.pi_cb)-param.GAMMA*log(pastpi_ss)).^2)+ Y_ss.*(-MC_ss) - param.fix  ...
		                						    + (h_ss-w_ss).*L_ss - param.iota.*V_ss)  ...
		                							+ (r_k_ss*v_ss - param.delta_0*v_ss^param.delta_1).*K_ss +  Q_ss*(K_ss-K_ss)-(K_ss-K_ss) - param.phi/2*(K_ss/K_ss-1)^2*K_ss ) - ...
													((Y_ss *(1-param.eta/(2*param.kappa)*(log(pi_ss)-(1-param.GAMMA)*log(param.pi_cb)-param.GAMMA*log(pastpi_ss)).^2)+ Y_ss.*(-MC_ss) - param.fix  ...
		                						    + (h_ss-w_ss).*L_ss - param.iota.*V_ss)  ...
		                							+ (r_k_ss*v_ss - param.delta_0*v_ss^param.delta_1).*K_ss +  Q_ss*(K_ss-K_ss)-(K_ss-K_ss) - param.phi/2*(K_ss/K_ss-1)^2*K_ss ) )/perturb_size;
F44_aux(Profit_ind,end-grid.oc+MC_ind)        = - ( ((Y_ss *(1-param.eta/(2*param.kappa)*(log(pi_ss)-(1-param.GAMMA)*log(param.pi_cb)-param.GAMMA*log(pastpi_ss)).^2)+ Y_ss.*(-MC_aux) - param.fix  ...
		                						    + (h_ss-w_ss).*L_ss - param.iota.*V_ss)  ...
		                							+ (r_k_ss*v_ss - param.delta_0*v_ss^param.delta_1).*K_ss +  Q_ss*(K_ss-K_ss)-(K_ss-K_ss) - param.phi/2*(K_ss/K_ss-1)^2*K_ss ) - ...
													((Y_ss *(1-param.eta/(2*param.kappa)*(log(pi_ss)-(1-param.GAMMA)*log(param.pi_cb)-param.GAMMA*log(pastpi_ss)).^2)+ Y_ss.*(-MC_ss) - param.fix  ...
		                						    + (h_ss-w_ss).*L_ss - param.iota.*V_ss)  ...
		                							+ (r_k_ss*v_ss - param.delta_0*v_ss^param.delta_1).*K_ss +  Q_ss*(K_ss-K_ss)-(K_ss-K_ss) - param.phi/2*(K_ss/K_ss-1)^2*K_ss ) )/perturb_size;
F44_aux(Profit_ind,end-grid.oc+Y_ind)         = - ( ((Y_aux *(1-param.eta/(2*param.kappa)*(log(pi_ss)-(1-param.GAMMA)*log(param.pi_cb)-param.GAMMA*log(pastpi_ss)).^2)+ Y_aux.*(-MC_ss) - param.fix  ...
		                						    + (h_ss-w_ss).*L_ss - param.iota.*V_ss)  ...
		                							+ (r_k_ss*v_ss - param.delta_0*v_ss^param.delta_1).*K_ss +  Q_ss*(K_ss-K_ss)-(K_ss-K_ss) - param.phi/2*(K_ss/K_ss-1)^2*K_ss ) - ...
													((Y_ss *(1-param.eta/(2*param.kappa)*(log(pi_ss)-(1-param.GAMMA)*log(param.pi_cb)-param.GAMMA*log(pastpi_ss)).^2)+ Y_ss.*(-MC_ss) - param.fix  ...
		                						    + (h_ss-w_ss).*L_ss - param.iota.*V_ss)  ...
		                							+ (r_k_ss*v_ss - param.delta_0*v_ss^param.delta_1).*K_ss +  Q_ss*(K_ss-K_ss)-(K_ss-K_ss) - param.phi/2*(K_ss/K_ss-1)^2*K_ss ) )/perturb_size;
F44_aux(Profit_ind,end-grid.oc+h_ind)         = - ( ((Y_ss *(1-param.eta/(2*param.kappa)*(log(pi_ss)-(1-param.GAMMA)*log(param.pi_cb)-param.GAMMA*log(pastpi_ss)).^2)+ Y_ss.*(-MC_ss) - param.fix  ...
		                						    + (h_aux-w_ss).*L_ss - param.iota.*V_ss)  ...
		                							+ (r_k_ss*v_ss - param.delta_0*v_ss^param.delta_1).*K_ss +  Q_ss*(K_ss-K_ss)-(K_ss-K_ss) - param.phi/2*(K_ss/K_ss-1)^2*K_ss ) - ...
													((Y_ss *(1-param.eta/(2*param.kappa)*(log(pi_ss)-(1-param.GAMMA)*log(param.pi_cb)-param.GAMMA*log(pastpi_ss)).^2)+ Y_ss.*(-MC_ss) - param.fix  ...
		                						    + (h_ss-w_ss).*L_ss - param.iota.*V_ss)  ...
		                							+ (r_k_ss*v_ss - param.delta_0*v_ss^param.delta_1).*K_ss +  Q_ss*(K_ss-K_ss)-(K_ss-K_ss) - param.phi/2*(K_ss/K_ss-1)^2*K_ss ) )/perturb_size;
F44_aux(Profit_ind,end-grid.oc+L_ind)         = - ( ((Y_ss *(1-param.eta/(2*param.kappa)*(log(pi_ss)-(1-param.GAMMA)*log(param.pi_cb)-param.GAMMA*log(pastpi_ss)).^2)+ Y_ss.*(-MC_ss) - param.fix  ...
		                						    + (h_ss-param.fix_L-w_ss).*L_aux - param.iota.*V_ss)  ...
		                							+ (r_k_ss*v_ss - param.delta_0*v_ss^param.delta_1).*K_ss +  Q_ss*(K_ss-K_ss)-(K_ss-K_ss) - param.phi/2*(K_ss/K_ss-1)^2*K_ss ) - ...
													((Y_ss *(1-param.eta/(2*param.kappa)*(log(pi_ss)-(1-param.GAMMA)*log(param.pi_cb)-param.GAMMA*log(pastpi_ss)).^2)+ Y_ss.*(-MC_ss) - param.fix  ...
		                						    + (h_ss-param.fix_L-w_ss).*L_ss - param.iota.*V_ss)  ...
		                							+ (r_k_ss*v_ss - param.delta_0*v_ss^param.delta_1).*K_ss +  Q_ss*(K_ss-K_ss)-(K_ss-K_ss) - param.phi/2*(K_ss/K_ss-1)^2*K_ss ) )/perturb_size;
F44_aux(Profit_ind,end-grid.oc+V_ind)         = - ( ((Y_ss *(1-param.eta/(2*param.kappa)*(log(pi_ss)-(1-param.GAMMA)*log(param.pi_cb)-param.GAMMA*log(pastpi_ss)).^2)+ Y_ss.*(-MC_ss) - param.fix  ...
		                						    + (h_ss-w_ss).*L_ss - param.iota.*V_aux)  ...
		                							+ (r_k_ss*v_ss - param.delta_0*v_ss^param.delta_1).*K_ss +  Q_ss*(K_ss-K_ss)-(K_ss-K_ss) - param.phi/2*(K_ss/K_ss-1)^2*K_ss ) - ...
													((Y_ss *(1-param.eta/(2*param.kappa)*(log(pi_ss)-(1-param.GAMMA)*log(param.pi_cb)-param.GAMMA*log(pastpi_ss)).^2)+ Y_ss.*(-MC_ss) - param.fix  ...
		                						    + (h_ss-w_ss).*L_ss - param.iota.*V_ss)  ...
		                							+ (r_k_ss*v_ss - param.delta_0*v_ss^param.delta_1).*K_ss +  Q_ss*(K_ss-K_ss)-(K_ss-K_ss) - param.phi/2*(K_ss/K_ss-1)^2*K_ss ) )/perturb_size;
F44_aux(Profit_ind,end-grid.oc+v_ind)         = - ( ((Y_ss *(1-param.eta/(2*param.kappa)*(log(pi_ss)-(1-param.GAMMA)*log(param.pi_cb)-param.GAMMA*log(pastpi_ss)).^2)+ Y_ss.*(-MC_ss) - param.fix  ...
		                						    + (h_ss-w_ss).*L_ss - param.iota.*V_ss)  ...
		                							+ (r_k_ss*v_aux - param.delta_0*v_aux^param.delta_1).*K_ss +  Q_ss*(K_ss-K_ss)-(K_ss-K_ss) - param.phi/2*(K_ss/K_ss-1)^2*K_ss ) - ...
													((Y_ss *(1-param.eta/(2*param.kappa)*(log(pi_ss)-(1-param.GAMMA)*log(param.pi_cb)-param.GAMMA*log(pastpi_ss)).^2)+ Y_ss.*(-MC_ss) - param.fix  ...
		                						    + (h_ss-w_ss).*L_ss - param.iota.*V_ss)  ...
		                							+ (r_k_ss*v_ss - param.delta_0*v_ss^param.delta_1).*K_ss +  Q_ss*(K_ss-K_ss)-(K_ss-K_ss) - param.phi/2*(K_ss/K_ss-1)^2*K_ss ) )/perturb_size;
F44_aux(Profit_ind,end-grid.oc+r_k_ind)       = - ( ((Y_ss *(1-param.eta/(2*param.kappa)*(log(pi_ss)-(1-param.GAMMA)*log(param.pi_cb)-param.GAMMA*log(pastpi_ss)).^2)+ Y_ss.*(-MC_ss) - param.fix  ...
		                						    + (h_ss-w_ss).*L_ss - param.iota.*V_ss)  ...
		                							+ (r_k_aux*v_ss - param.delta_0*v_ss^param.delta_1).*K_ss +  Q_ss*(K_ss-K_ss)-(K_ss-K_ss) - param.phi/2*(K_ss/K_ss-1)^2*K_ss ) - ...
													((Y_ss *(1-param.eta/(2*param.kappa)*(log(pi_ss)-(1-param.GAMMA)*log(param.pi_cb)-param.GAMMA*log(pastpi_ss)).^2)+ Y_ss.*(-MC_ss) - param.fix  ...
		                						    + (h_ss-w_ss).*L_ss - param.iota.*V_ss)  ...
		                							+ (r_k_ss*v_ss - param.delta_0*v_ss^param.delta_1).*K_ss +  Q_ss*(K_ss-K_ss)-(K_ss-K_ss) - param.phi/2*(K_ss/K_ss-1)^2*K_ss ) )/perturb_size;
F44_aux(Profit_ind,end-grid.oc+K_ind)         = - ( ((Y_ss *(1-param.eta/(2*param.kappa)*(log(pi_ss)-(1-param.GAMMA)*log(param.pi_cb)-param.GAMMA*log(pastpi_ss)).^2)+ Y_ss.*(-MC_ss) - param.fix  ...
		                						    + (h_ss-w_ss).*L_ss - param.iota.*V_ss)  ...
		                							+ (r_k_ss*v_ss - param.delta_0*v_ss^param.delta_1).*K_aux +  Q_ss*(K_ss-K_aux)-(K_ss-K_aux) - param.phi/2*(K_ss/K_aux-1)^2*K_aux ) - ...
													((Y_ss *(1-param.eta/(2*param.kappa)*(log(pi_ss)-(1-param.GAMMA)*log(param.pi_cb)-param.GAMMA*log(pastpi_ss)).^2)+ Y_ss.*(-MC_ss) - param.fix  ...
		                						    + (h_ss-w_ss).*L_ss - param.iota.*V_ss)  ...
		                							+ (r_k_ss*v_ss - param.delta_0*v_ss^param.delta_1).*K_ss +  Q_ss*(K_ss-K_ss)-(K_ss-K_ss) - param.phi/2*(K_ss/K_ss-1)^2*K_ss ) )/perturb_size;

% (19) r^k

% r_k = Z*MC*param.theta*(v*K).^(param.theta-1).*(L).^(1-param.theta)

F41_aux(r_k_ind,end-grid.os+Z_ind)            = -( Z_aux*MC_ss*param.theta*(v_ss*K_ss).^(param.theta-1).*(L_ss).^param.theta_2 - Z_ss*MC_ss*param.theta*(v_ss*K_ss).^(param.theta-1).*(L_ss).^param.theta_2)/perturb_size;
F44_aux(r_k_ind,end-grid.oc+r_k_ind)          = (r_k_aux-r_k_ss)/perturb_size;
F44_aux(r_k_ind,end-grid.oc+MC_ind)           = -( Z_ss*MC_aux*param.theta*(v_ss*K_ss).^(param.theta-1).*(L_ss).^param.theta_2 - Z_ss*MC_ss*param.theta*(v_ss*K_ss).^(param.theta-1).*(L_ss).^param.theta_2)/perturb_size;
F44_aux(r_k_ind,end-grid.oc+v_ind)            = -( Z_ss*MC_ss*(v_aux*K_ss).^(param.theta-1).*(L_ss) .^param.theta_2*(param.theta+param.alpha_lk*log(L_ss/SS_stats.L)+param.alpha_kk*log(v_aux*K_ss/SS_stats.v/grid.K)) - Z_ss*MC_ss*param.theta*(v_ss*K_ss).^(param.theta-1).*(L_ss).^param.theta_2)/perturb_size;
F44_aux(r_k_ind,end-grid.oc+K_ind)            = -( Z_ss*MC_ss*(v_ss*K_aux).^(param.theta-1).*(L_ss) .^param.theta_2*(param.theta+param.alpha_lk*log(L_ss/SS_stats.L)+param.alpha_kk*log(v_ss*K_aux/SS_stats.v/grid.K)) - Z_ss*MC_ss*param.theta*(v_ss*K_ss).^(param.theta-1).*(L_ss).^param.theta_2)/perturb_size;
F44_aux(r_k_ind,end-grid.oc+L_ind)            = -( Z_ss*MC_ss*(v_ss*K_ss) .^(param.theta-1).*(L_aux).^param.theta_2*(param.theta+param.alpha_lk*log(L_aux/SS_stats.L)+param.alpha_kk*log(v_ss*K_ss/SS_stats.v/grid.K)) - Z_ss*MC_ss*param.theta*(v_ss*K_ss).^(param.theta-1).*(L_ss).^param.theta_2)/perturb_size;

% (20) r^a

% r_a = (1-tau_a)*(1-param.Eratio-param.b_share)*PROFIT/K

F44_aux(r_a_ind,end-grid.oc+r_a_ind)          = (r_a_aux-r_a_ss)/perturb_size;
F44_aux(r_a_ind,end-grid.oc+Profit_ind)       = - ((1-tau_a_ss)*(1-param.Eratio-param.b_share)*Profit_aux/K_ss - (1-tau_a_ss)*(1-param.Eratio-param.b_share)*Profit_ss/K_ss )/perturb_size;
F44_aux(r_a_ind,end-grid.oc+K_ind)            = - ((1-tau_a_ss)*(1-param.Eratio-param.b_share)*Profit_ss/K_aux - (1-tau_a_ss)*(1-param.Eratio-param.b_share)*Profit_ss/K_ss )/perturb_size;

% (21) MC

% 1-1/eta + (log(PI/param.pi_bar)-LAMBDA*eta2next/eta2*Ynext/Y*log(PInext/param.pi_bar))/param.kappa

% RHS(nx+MC_ind) 	   = 1-1/eta + (log(PI/(pastpiminus^param.GAMMA*param.pi_bar^param.pi_bar^(1-param.GAMMA)))-LAMBDA*eta2next/eta2*Ynext/Y*log(PInext/(PI^param.GAMMA*param.pi_bar^(1-param.GAMMA)) ))/param.kappa;

F44_aux(MC_ind,end-grid.oc+MC_ind)            = (MC_aux-MC_ss)/perturb_size;
F41_aux(MC_ind,end-grid.os+eta_ind)           = - ( (1-1/eta_aux)-(1-1/eta_ss) )/perturb_size;
F44_aux(MC_ind,end-grid.oc+pi_ind)            = - ( ( 1-1/eta_ss + (log(pi_aux/(pastpi_ss ^param.GAMMA*param.pi_bar^(1-param.GAMMA)))-Lambda_ss*eta2_ss /eta2_ss *Y_ss/Y_ss *log(pi_ss /(pi_aux^param.GAMMA*param.pi_bar^(1-param.GAMMA)) ))/param.kappa  ) - ...
													( 1-1/eta_ss + (log(pi_ss /(pastpi_ss ^param.GAMMA*param.pi_bar^(1-param.GAMMA)))-Lambda_ss*eta2_ss /eta2_ss *Y_ss/Y_ss *log(pi_ss /(pi_ss ^param.GAMMA*param.pi_bar^(1-param.GAMMA)) ))/param.kappa  ) )/perturb_size;
F42_aux(MC_ind,end-grid.oc+pi_ind)            = - ( ( 1-1/eta_ss + (log(pi_ss /(pastpi_ss ^param.GAMMA*param.pi_bar^(1-param.GAMMA)))-Lambda_ss*eta2_ss /eta2_ss *Y_ss/Y_ss *log(pi_aux/(pi_ss ^param.GAMMA*param.pi_bar^(1-param.GAMMA)) ))/param.kappa  ) - ...
													( 1-1/eta_ss + (log(pi_ss /(pastpi_ss ^param.GAMMA*param.pi_bar^(1-param.GAMMA)))-Lambda_ss*eta2_ss /eta2_ss *Y_ss/Y_ss *log(pi_ss /(pi_ss ^param.GAMMA*param.pi_bar^(1-param.GAMMA)) ))/param.kappa  ) )/perturb_size;
F43_aux(MC_ind,end-grid.os+pastpi_ind)        = - ( ( 1-1/eta_ss + (log(pi_ss /(pastpi_aux^param.GAMMA*param.pi_bar^(1-param.GAMMA)))-Lambda_ss*eta2_ss /eta2_ss *Y_ss/Y_ss *log(pi_ss /(pi_ss ^param.GAMMA*param.pi_bar^(1-param.GAMMA)) ))/param.kappa  ) - ...
													( 1-1/eta_ss + (log(pi_ss /(pastpi_ss ^param.GAMMA*param.pi_bar^(1-param.GAMMA)))-Lambda_ss*eta2_ss /eta2_ss *Y_ss/Y_ss *log(pi_ss /(pi_ss ^param.GAMMA*param.pi_bar^(1-param.GAMMA)) ))/param.kappa  ) )/perturb_size;
F44_aux(MC_ind,end-grid.oc+Y_ind)             = - ( ( 1-1/eta_ss + (log(pi_ss /(pastpi_ss ^param.GAMMA*param.pi_bar^(1-param.GAMMA)))-Lambda_ss*eta2_ss /eta2_ss *Y_ss/Y_aux*log(pi_ss /(pi_ss ^param.GAMMA*param.pi_bar^(1-param.GAMMA)) ))/param.kappa  ) - ...
													( 1-1/eta_ss + (log(pi_ss /(pastpi_ss ^param.GAMMA*param.pi_bar^(1-param.GAMMA)))-Lambda_ss*eta2_ss /eta2_ss *Y_ss/Y_ss *log(pi_ss /(pi_ss ^param.GAMMA*param.pi_bar^(1-param.GAMMA)) ))/param.kappa  ) )/perturb_size;
F42_aux(MC_ind,end-grid.oc+Y_ind)             = - ( ( 1-1/eta_ss + (log(pi_ss /(pastpi_ss ^param.GAMMA*param.pi_bar^(1-param.GAMMA)))-Lambda_ss*eta2_ss /eta2_ss *Y_aux/Y_ss*log(pi_ss /(pi_ss ^param.GAMMA*param.pi_bar^(1-param.GAMMA)) ))/param.kappa  ) - ...
													( 1-1/eta_ss + (log(pi_ss /(pastpi_ss ^param.GAMMA*param.pi_bar^(1-param.GAMMA)))-Lambda_ss*eta2_ss /eta2_ss *Y_ss/Y_ss *log(pi_ss /(pi_ss ^param.GAMMA*param.pi_bar^(1-param.GAMMA)) ))/param.kappa  ) )/perturb_size;
F42_aux(MC_ind,end-grid.oc+eta2_ind)          = - ( ( 1-1/eta_ss + (log(pi_ss /(pastpi_ss ^param.GAMMA*param.pi_bar^(1-param.GAMMA)))-Lambda_ss*eta2_aux/eta2_ss *Y_ss/Y_ss *log(pi_ss /(pi_ss ^param.GAMMA*param.pi_bar^(1-param.GAMMA)) ))/param.kappa  ) - ...
													( 1-1/eta_ss + (log(pi_ss /(pastpi_ss ^param.GAMMA*param.pi_bar^(1-param.GAMMA)))-Lambda_ss*eta2_ss /eta2_ss *Y_ss/Y_ss *log(pi_ss /(pi_ss ^param.GAMMA*param.pi_bar^(1-param.GAMMA)) ))/param.kappa  ) )/perturb_size;
F44_aux(MC_ind,end-grid.oc+eta2_ind)          = - ( ( 1-1/eta_ss + (log(pi_ss /(pastpi_ss ^param.GAMMA*param.pi_bar^(1-param.GAMMA)))-Lambda_ss*eta2_ss /eta2_aux*Y_ss/Y_ss *log(pi_ss /(pi_ss ^param.GAMMA*param.pi_bar^(1-param.GAMMA)) ))/param.kappa  ) - ...
													( 1-1/eta_ss + (log(pi_ss /(pastpi_ss ^param.GAMMA*param.pi_bar^(1-param.GAMMA)))-Lambda_ss*eta2_ss /eta2_ss *Y_ss/Y_ss *log(pi_ss /(pi_ss ^param.GAMMA*param.pi_bar^(1-param.GAMMA)) ))/param.kappa  ) )/perturb_size;

% % (22) sale

% % sale = stock*s2s

% % (23) s2s : sales-to-stock ratio

% % A_aux    = 0;
% % B_aux    = param.kappa*((1-param.delta_i)*SS_stats.Lambda*SS_stats.mc - (param.eta-1)/param.eta - param.phi_x*log(param.s2s) - param.phi_x*SS_stats.Lambda*log(param.s2s));
% % C_aux    = param.kappa*param.phi_x*(1+SS_stats.Lambda);

% % A_aux = LAMBDA*salenext/sale*(log(PInext)-log(PI^param.GAMMA*param.pi_bar^(1-param.GAMMA)));
% % B_aux = param.kappa*( (1-param.delta_i)*LAMBDA*MCnext - param.phi_x*log(pasts2sminus) - param.phi_x*LAMBDA*salenext/sale*log(s2snext) - (param.eta-1)/param.eta  );

% % s2s = (exp( (log(PI)-log(pastpiminus^param.GAMMA*param.pi_bar^(1-param.GAMMA))-LAMBDA*salenext/sale*(log(PInext)-log(PI^param.GAMMA*param.pi_bar^(1-param.GAMMA)))- ...
% 			% param.kappa*( (1-param.delta_i)*LAMBDA*MCnext - param.phi_x*log(pasts2sminus) - param.phi_x*LAMBDA*salenext/sale*log(s2snext) - (param.eta-1)/param.eta  ))/( param.kappa*param.phi_x*(1+LAMBDA*salenext/sale) ) ));


% % (24) stock

% % stock = Y + (1-param.delta_i)*invenminus

% (25) Unemployment rate

% unemp = (1-((1-param.lambda)*N + M)-param.Eshare)/(1-param.Eshare)


F44_aux(unemp_ind,end-grid.oc+unemp_ind)      = (unemp_aux-unemp_ss)/perturb_size;
F44_aux(unemp_ind,end-grid.oc+N_ind)          = - ((1-((1-param.lambda)*N_aux + M_ss)-param.Eshare)/(1-param.Eshare) - (1-((1-param.lambda)*N_ss + M_ss)-param.Eshare)/(1-param.Eshare))/perturb_size;
F44_aux(unemp_ind,end-grid.oc+M_ind)          = - ((1-((1-param.lambda)*N_ss + M_aux)-param.Eshare)/(1-param.Eshare) - (1-((1-param.lambda)*N_ss + M_ss)-param.Eshare)/(1-param.Eshare))/perturb_size;
F44_aux(unemp_ind,end-grid.oc+l_lambda_ind)   = - ((1-((1-l_lambda_aux)*N_ss + M_ss)-param.Eshare)/(1-param.Eshare) - (1-((1-param.lambda)*N_ss + M_ss)-param.Eshare)/(1-param.Eshare))/perturb_size;

% (26) n

% nn = ((1-tau_w)*W/param.psi)^(1/param.xi)

F41_aux(nn_ind,end-grid.os+w_ind)             = - ( ((1-tau_w_ss)*w_aux/param.psi)^(1/param.xi)- ((1-tau_w_ss)*w_ss/param.psi)^(1/param.xi))/perturb_size;
F44_aux(nn_ind,end-grid.oc+nn_ind)            = (nn_aux-nn_ss)/perturb_size;

% (27) M

% M = ((1-N-Eshare)+param.lambda.*(N)).*V/((((1-N-Eshare)+param.lambda.*(N))^param.alpha+V^param.alpha)^(1/param.alpha))

% aux_M_1 = 1-SS_stats.N-param.Eshare + param.lambda*SS_stats.N;

F44_aux(M_ind,end-grid.oc+M_ind)    = (M_aux-M_ss)/perturb_size;
F44_aux(M_ind,end-grid.oc+N_ind)    = - ( ((1-N_aux-param.Eshare)+param.lambda.*(N_aux)).*V_ss/((((1-N_aux-param.Eshare)+param.lambda.*(N_aux))^param.alpha+V_ss^param.alpha)^(1/param.alpha)) - ...
										  ((1-N_ss-param.Eshare)+param.lambda.*(N_ss)).*V_ss/((((1-N_ss-param.Eshare)+param.lambda.*(N_ss))^param.alpha+V_ss^param.alpha)^(1/param.alpha)) )/perturb_size;
F44_aux(M_ind,end-grid.oc+V_ind)    = - ( ((1-N_ss-param.Eshare)+param.lambda.*(N_ss)).*V_aux/((((1-N_ss-param.Eshare)+param.lambda.*(N_ss))^param.alpha+V_aux^param.alpha)^(1/param.alpha)) - ...
										  ((1-N_ss-param.Eshare)+param.lambda.*(N_ss)).*V_ss/((((1-N_ss-param.Eshare)+param.lambda.*(N_ss))^param.alpha+V_ss^param.alpha)^(1/param.alpha)) )/perturb_size;
F44_aux(M_ind,end-grid.oc+l_lambda_ind)    = - ( ((1-N_ss-param.Eshare)+l_lambda_aux.*(N_ss)).*V_ss/((((1-N_ss-param.Eshare)+l_lambda_aux.*(N_ss))^param.alpha+V_ss^param.alpha)^(1/param.alpha)) - ...
										  ((1-N_ss-param.Eshare)+param.lambda.*(N_ss)).*V_ss/((((1-N_ss-param.Eshare)+param.lambda.*(N_ss))^param.alpha+V_ss^param.alpha)^(1/param.alpha)) )/perturb_size;

% (24) f

% f = M/((1-N-Eshare)+param.lambda*(N));

F44_aux(f_ind,end-grid.oc+f_ind)         = (f_aux-f_ss)/perturb_size;
F44_aux(f_ind,end-grid.oc+M_ind)         = - (M_aux-M_ss)/((1-N_ss-param.Eshare)+param.lambda*(N_ss))/perturb_size;
F44_aux(f_ind,end-grid.oc+N_ind)         = - ( M_ss/((1-N_aux-param.Eshare)+param.lambda*(N_aux)) - M_ss/((1-N_ss-param.Eshare)+param.lambda*(N_ss))   )/perturb_size;
F44_aux(f_ind,end-grid.oc+l_lambda_ind)  = - ( M_ss/((1-N_ss-param.Eshare)+l_lambda_aux*(N_ss)) - M_ss/((1-N_ss-param.Eshare)+param.lambda*(N_ss))   )/perturb_size;

% (25) zz

% zz = ((RRa-R_cbminus/PI)*levminus + R_cbminus/PI)

F44_aux(zz_ind,end-grid.oc+zz_ind)            = (zz_aux-zz_ss)/perturb_size;
F44_aux(zz_ind,end-grid.oc+RRa_ind)           = - (((RRa_aux-param.b_a_aux2*R_cb_ss /pi_ss) *lev_ss  + param.b_a_aux2*R_cb_ss /pi_ss)  - ((RRa_ss-R_cb_ss/pi_ss)*lev_ss + R_cb_ss/pi_ss)  )/perturb_size;
F43_aux(zz_ind,end-grid.os+R_cb_ind)          = - (((RRa_ss- param.b_a_aux2*R_cb_aux/pi_ss) *lev_ss  + param.b_a_aux2*R_cb_aux/pi_ss)  - ((RRa_ss-R_cb_ss/pi_ss)*lev_ss + R_cb_ss/pi_ss)  )/perturb_size;
F43_aux(zz_ind,end-grid.os+lev_ind)           = - (((RRa_ss- param.b_a_aux2*R_cb_ss /pi_ss) *lev_aux + param.b_a_aux2*R_cb_ss /pi_ss)  - ((RRa_ss-R_cb_ss/pi_ss)*lev_ss + R_cb_ss/pi_ss)  )/perturb_size;
F44_aux(zz_ind,end-grid.oc+pi_ind)            = - (((RRa_ss- param.b_a_aux2*R_cb_ss /pi_aux)*lev_ss  + param.b_a_aux2*R_cb_ss /pi_aux) - ((RRa_ss-R_cb_ss/pi_ss)*lev_ss + R_cb_ss/pi_ss)  )/perturb_size;

% (26) xx

% xx = (lev/levminus)*zz

F44_aux(xx_ind,end-grid.oc+xx_ind)            = (xx_aux-xx_ss)/perturb_size;
F44_aux(xx_ind,end-grid.oc+zz_ind)            = - ((lev_ss/lev_ss)*zz_aux-(lev_ss/lev_ss)*zz_ss)/perturb_size;
F41_aux(xx_ind,end-grid.os+lev_ind)           = - ((lev_aux/lev_ss)*zz_ss-(lev_ss/lev_ss)*zz_ss)/perturb_size;
F43_aux(xx_ind,end-grid.os+lev_ind)           = - ((lev_ss/lev_aux)*zz_ss-(lev_ss/lev_ss)*zz_ss)/perturb_size;

% (27) vv

% vv = ((1-param.theta_b)*PSI*LAMBDA/SS_stats.Lambda*(RRanext-R_cb/PInext)+param.theta_b*PSI*LAMBDA/SS_stats.Lambda*xxnext*vvnext)

F41_aux(vv_ind,end-grid.os+R_cb_ind)   = - ( ((1-param.theta_b)*Lambda_ss       *(RRa_ss -param.b_a_aux2*R_cb_aux/pi_ss)+param.theta_b*Lambda_ss       *xx_ss*vv_ss)  - ((1-param.theta_b)*Lambda_ss*(RRa_ss-R_cb_ss/pi_ss)+param.theta_b*Lambda_ss*xx_ss*vv_ss)  )/perturb_size;
F41_aux(vv_ind,end-grid.os+BB_ind)     = - ( ((1-param.theta_b)*BB_aux*Lambda_ss*(RRa_ss -param.b_a_aux2*R_cb_ss/pi_ss) +param.theta_b*BB_aux*Lambda_ss*xx_ss*vv_ss)  - ((1-param.theta_b)*Lambda_ss*(RRa_ss-R_cb_ss/pi_ss)+param.theta_b*Lambda_ss*xx_ss*vv_ss)  )/perturb_size;
F42_aux(vv_ind,end-grid.oc+RRa_ind)    = - ( ((1-param.theta_b)*Lambda_ss       *(RRa_aux-param.b_a_aux2*R_cb_ss/pi_ss) +param.theta_b*Lambda_ss       *xx_ss*vv_ss)  - ((1-param.theta_b)*Lambda_ss*(RRa_ss-R_cb_ss/pi_ss)+param.theta_b*Lambda_ss*xx_ss*vv_ss)  )/perturb_size;
F42_aux(vv_ind,end-grid.oc+pi_ind)     = - ( ((1-param.theta_b)*Lambda_ss       *(RRa_ss -param.b_a_aux2*R_cb_ss/pi_aux)+param.theta_b*Lambda_ss       *xx_ss*vv_ss)  - ((1-param.theta_b)*Lambda_ss*(RRa_ss-R_cb_ss/pi_ss)+param.theta_b*Lambda_ss*xx_ss*vv_ss)  )/perturb_size;
F42_aux(vv_ind,end-grid.oc+xx_ind)     = - ( ((1-param.theta_b)*Lambda_ss       *(RRa_ss -param.b_a_aux2*R_cb_ss/pi_ss) +param.theta_b*Lambda_ss       *xx_aux*vv_ss) - ((1-param.theta_b)*Lambda_ss*(RRa_ss-R_cb_ss/pi_ss)+param.theta_b*Lambda_ss*xx_ss*vv_ss)  )/perturb_size;
F42_aux(vv_ind,end-grid.oc+vv_ind)     = - ( ((1-param.theta_b)*Lambda_ss       *(RRa_ss -param.b_a_aux2*R_cb_ss/pi_ss) +param.theta_b*Lambda_ss       *xx_ss*vv_aux) - ((1-param.theta_b)*Lambda_ss*(RRa_ss-R_cb_ss/pi_ss)+param.theta_b*Lambda_ss*xx_ss*vv_ss)  )/perturb_size;
F44_aux(vv_ind,end-grid.oc+Lambda_ind) = - ( ((1-param.theta_b)*Lambda_aux      *(RRa_ss -param.b_a_aux2*R_cb_ss/pi_ss) +param.theta_b*Lambda_aux      *xx_ss*vv_ss)  - ((1-param.theta_b)*Lambda_ss*(RRa_ss-R_cb_ss/pi_ss)+param.theta_b*Lambda_ss*xx_ss*vv_ss)  )/perturb_size;
F44_aux(vv_ind,end-grid.oc+vv_ind)     = (vv_aux-vv_ss)/perturb_size;


% (28) ee

% ee = ((1-param.theta_b)*PSI_RP*LAMBDA/SS_stats.Lambda*R_cb/PInext+param.theta_b*PSI_RP*LAMBDA/SS_stats.Lambda*zznext*eenext)

F41_aux(ee_ind,end-grid.os+R_cb_ind)          = -((1-param.theta_b)*Lambda_ss*(R_cb_aux-param.b_a_aux2*R_cb_ss)/pi_ss)/perturb_size;
F41_aux(ee_ind,end-grid.os+BB_ind)            = -(((1-param.theta_b)*BB_aux*Lambda_ss*param.b_a_aux2*R_cb_ss/pi_ss+param.theta_b*BB_aux*Lambda_ss*zz_ss*ee_ss)-((1-param.theta_b)*Lambda_ss*R_cb_ss/pi_ss+param.theta_b*Lambda_ss*zz_ss*ee_ss))/perturb_size;
F44_aux(ee_ind,end-grid.oc+ee_ind)            = (ee_aux-ee_ss)/perturb_size;
F44_aux(ee_ind,end-grid.oc+Lambda_ind)        = -(((1-param.theta_b)*Lambda_aux*param.b_a_aux2*R_cb_ss/pi_ss +param.theta_b*Lambda_aux*zz_ss *ee_ss) -((1-param.theta_b)*Lambda_ss*R_cb_ss/pi_ss+param.theta_b*Lambda_ss*zz_ss*ee_ss))/perturb_size;
F42_aux(ee_ind,end-grid.oc+pi_ind)            = -(((1-param.theta_b)*Lambda_ss *param.b_a_aux2*R_cb_ss/pi_aux+param.theta_b*Lambda_ss *zz_ss *ee_ss) -((1-param.theta_b)*Lambda_ss*R_cb_ss/pi_ss+param.theta_b*Lambda_ss*zz_ss*ee_ss))/perturb_size;
F42_aux(ee_ind,end-grid.oc+zz_ind)            = -(((1-param.theta_b)*Lambda_ss *param.b_a_aux2*R_cb_ss/pi_ss +param.theta_b*Lambda_ss *zz_aux*ee_ss) -((1-param.theta_b)*Lambda_ss*R_cb_ss/pi_ss+param.theta_b*Lambda_ss*zz_ss*ee_ss))/perturb_size;
F42_aux(ee_ind,end-grid.oc+ee_ind)            = -(((1-param.theta_b)*Lambda_ss *param.b_a_aux2*R_cb_ss/pi_ss +param.theta_b*Lambda_ss *zz_ss *ee_aux)-((1-param.theta_b)*Lambda_ss*R_cb_ss/pi_ss+param.theta_b*Lambda_ss*zz_ss*ee_ss))/perturb_size;

% (29) C_b

% C_b = ((1*param.beta_b*(1-param.death_rate)*R_cb/PInext)/(1+param.phi_B*(B_bnext/B_b-1))*C_bnext^(-param.sigma2))^(-1/param.sigma2);

F44_aux(C_b_ind,end-grid.oc+C_b_ind)          = (C_b_aux-C_b_ss)/perturb_size;
F41_aux(C_b_ind,end-grid.os+R_cb_ind)         = - ( ((1*param.beta_b*R_cb_aux/pi_ss) /(1+param.phi_B*(B_b_ss/B_b_ss-1))*C_b_ss^(-param.sigma2))^(-1/param.sigma2) -  ((1*param.beta_b*R_cb_ss/pi_ss)/(1+param.phi_B*(B_b_ss/B_b_ss-1))*C_b_ss^(-param.sigma2))^(-1/param.sigma2) )/perturb_size;
F42_aux(C_b_ind,end-grid.oc+pi_ind)           = - ( ((1*param.beta_b*R_cb_ss /pi_aux)/(1+param.phi_B*(B_b_ss/B_b_ss-1))*C_b_ss^(-param.sigma2))^(-1/param.sigma2) -  ((1*param.beta_b*R_cb_ss/pi_ss)/(1+param.phi_B*(B_b_ss/B_b_ss-1))*C_b_ss^(-param.sigma2))^(-1/param.sigma2) )/perturb_size;
F42_aux(C_b_ind,end-grid.oc+C_b_ind)          = - ( ((1*param.beta_b*R_cb_ss /pi_ss) /(1+param.phi_B*(B_b_ss/B_b_ss-1))*C_b_aux^(-param.sigma2))^(-1/param.sigma2) - ((1*param.beta_b*R_cb_ss/pi_ss)/(1+param.phi_B*(B_b_ss/B_b_ss-1))*C_b_ss^(-param.sigma2))^(-1/param.sigma2) )/perturb_size;

% (30) Profit_FI

% Profit_FI = ((1-param.theta_b)*((RRa - R_cbminus/PI)*levminus+R_cbminus/PI)*NW_bminus - param.omega*Qminus*A_b);

F44_aux(Profit_FI_ind,end-grid.oc+Profit_FI_ind) = (Profit_FI_aux-Profit_FI_ss)/perturb_size;
F44_aux(Profit_FI_ind,end-grid.oc+RRa_ind)       = - ( ((1-param.theta_b)*((RRa_aux - param.b_a_aux2*R_cb_ss/pi_ss) *lev_ss +param.b_a_aux2*R_cb_ss/pi_ss) *NW_b_ss  - param.omega*Q_ss*A_b_ss)  - ((1-param.theta_b)*((RRa_ss - param.b_a_aux2*R_cb_ss/pi_ss)*lev_ss+R_cb_ss/pi_ss)*NW_b_ss - param.omega*Q_ss*A_b_ss)     )/perturb_size;
F44_aux(Profit_FI_ind,end-grid.oc+pi_ind)        = - ( ((1-param.theta_b)*((RRa_ss  - param.b_a_aux2*R_cb_ss/pi_aux)*lev_ss +param.b_a_aux2*R_cb_ss/pi_aux)*NW_b_ss  - param.omega*Q_ss*A_b_ss)  - ((1-param.theta_b)*((RRa_ss - param.b_a_aux2*R_cb_ss/pi_ss)*lev_ss+R_cb_ss/pi_ss)*NW_b_ss - param.omega*Q_ss*A_b_ss)     )/perturb_size;
F43_aux(Profit_FI_ind,end-grid.os+lev_ind)       = - ( ((1-param.theta_b)*((RRa_ss  - param.b_a_aux2*R_cb_ss/pi_ss) *lev_aux+param.b_a_aux2*R_cb_ss/pi_ss) *NW_b_ss  - param.omega*Q_ss*A_b_ss)  - ((1-param.theta_b)*((RRa_ss - param.b_a_aux2*R_cb_ss/pi_ss)*lev_ss+R_cb_ss/pi_ss)*NW_b_ss - param.omega*Q_ss*A_b_ss)     )/perturb_size;
F43_aux(Profit_FI_ind,end-grid.os+R_cb_ind)      = - ( ((1-param.theta_b)*((RRa_ss  - param.b_a_aux2*R_cb_aux/pi_ss)*lev_ss +param.b_a_aux2*R_cb_aux/pi_ss)*NW_b_ss  - param.omega*Q_ss*A_b_ss)  - ((1-param.theta_b)*((RRa_ss - param.b_a_aux2*R_cb_ss/pi_ss)*lev_ss+R_cb_ss/pi_ss)*NW_b_ss - param.omega*Q_ss*A_b_ss)     )/perturb_size;
F43_aux(Profit_FI_ind,end-grid.os+NW_b_ind)      = - ( ((1-param.theta_b)*((RRa_ss  - param.b_a_aux2*R_cb_ss/pi_ss) *lev_ss +param.b_a_aux2*R_cb_ss/pi_ss) *NW_b_aux - param.omega*Q_ss*A_b_ss)  - ((1-param.theta_b)*((RRa_ss - param.b_a_aux2*R_cb_ss/pi_ss)*lev_ss+R_cb_ss/pi_ss)*NW_b_ss - param.omega*Q_ss*A_b_ss)     )/perturb_size;
F43_aux(Profit_FI_ind,end-grid.os+Q_ind)         = - ( ((1-param.theta_b)*((RRa_ss  - param.b_a_aux2*R_cb_ss/pi_ss) *lev_ss +param.b_a_aux2*R_cb_ss/pi_ss) *NW_b_ss  - param.omega*Q_aux*A_b_ss) - ((1-param.theta_b)*((RRa_ss - param.b_a_aux2*R_cb_ss/pi_ss)*lev_ss+R_cb_ss/pi_ss)*NW_b_ss - param.omega*Q_ss*A_b_ss)     )/perturb_size;
F43_aux(Profit_FI_ind,end-grid.os+A_b_ind)       = - ( ((1-param.theta_b)*((RRa_ss  - param.b_a_aux2*R_cb_ss/pi_ss) *lev_ss +param.b_a_aux2*R_cb_ss/pi_ss) *NW_b_ss  - param.omega*Q_ss*A_b_aux) - ((1-param.theta_b)*((RRa_ss - param.b_a_aux2*R_cb_ss/pi_ss)*lev_ss+R_cb_ss/pi_ss)*NW_b_ss - param.omega*Q_ss*A_b_ss)     )/perturb_size;

% (31) RRa

% RRa = (Q+R_A)/Qminus

F41_aux(RRa_ind,end-grid.os+Q_ind)               = - ((Q_aux+r_a_ss)/Q_ss-(Q_ss+r_a_ss)/Q_ss)/perturb_size;
F43_aux(RRa_ind,end-grid.os+Q_ind)               = - ((Q_ss+r_a_ss)/Q_aux-(Q_ss+r_a_ss)/Q_ss)/perturb_size;
F44_aux(RRa_ind,end-grid.oc+r_a_ind)             = - ((Q_ss+r_a_aux)/Q_ss-(Q_ss+r_a_ss)/Q_ss)/perturb_size;
F44_aux(RRa_ind,end-grid.oc+RRa_ind)             = (RRa_aux-RRa_ss)/perturb_size;

% (32) RR

% RR = R_cbminus/PI

F43_aux(RR_ind,end-grid.os+R_cb_ind) = - (R_cb_aux/pi_ss-R_cb_ss/pi_ss)/perturb_size;
F44_aux(RR_ind,end-grid.oc+RR_ind)   = (RR_aux-RR_ss)/perturb_size;
F44_aux(RR_ind,end-grid.oc+pi_ind)   = - (R_cb_ss/pi_aux-R_cb_ss/pi_ss)/perturb_size;

% (33) I

% I = A_hhnext + A_gnext + A_bnext + SS_stats.A_F - (1-param.delta_0*v^param.delta_1)*K
% I = A_hhnext + A_gnext + A_bnext + SS_stats.A_F - (1-(param.delta_0+param.delta_1*(v-SS_stats.v)+param.delta_2/2*(v-SS_stats.v)^2))*K


F41_aux(I_ind,end-grid.os+iota_ind) = - (iota_aux*K_ss-K_ss)/perturb_size; 
F41_aux(I_ind,end-grid.os+A_b_ind)  = - (A_b_aux-A_b_ss)/perturb_size; 
F41_aux(I_ind,end-grid.os+A_g_ind)  = - (A_g_aux-A_g_ss)/perturb_size; 
F42_aux(I_ind,end-grid.oc+A_hh_ind) = - (A_hh_aux-A_hh_ss)/perturb_size; 
F44_aux(I_ind,end-grid.oc+I_ind)    = (I_aux-I_ss)/perturb_size;
F44_aux(I_ind,end-grid.oc+K_ind)    = (1-param.delta_0*v_ss^param.delta_1)*(K_aux-K_ss)/perturb_size;
F44_aux(I_ind,end-grid.oc+v_ind)    = - param.delta_0*(v_aux^param.delta_1-v_ss^param.delta_1)*K_ss/perturb_size;
F44_aux(I_ind,end-grid.oc+x_k_ind)  = - (K_ss*(1+param.phi/2*(log(x_k_aux))^2)-K_ss)/perturb_size; 

% (34) x_k

% x_k = Knext/K

F42_aux(x_k_ind,end-grid.oc+K_ind)    = - (K_aux/K_ss-K_ss/K_ss)/perturb_size;            
F44_aux(x_k_ind,end-grid.oc+x_k_ind)  = (x_k_aux-x_k_ss)/perturb_size;
F44_aux(x_k_ind,end-grid.oc+K_ind)    = - (K_ss/K_aux-K_ss/K_ss)/perturb_size;                          


% (39) A_g_obs

% A_g_obs - (A_gauxnext-1)/Y

F41_aux(A_g_obs_ind,end-grid.os+A_g_ind)         = - ( (A_g_aux) - (A_g_ss)  )/perturb_size;
F44_aux(A_g_obs_ind,end-grid.oc+A_g_obs_ind)     = (A_g_obs_aux - A_g_obs_ss)/perturb_size;

% % (41) Y_obs

% % Y_obs = Y/pastYminus;

F44_aux(Y_obs_ind,end-grid.oc+Y_ind)         = - (Y_aux-Y_ss)/perturb_size;
F44_aux(Y_obs_ind,end-grid.oc+Y_obs_ind)     = (Y_obs_aux-Y_obs_ss)/perturb_size;

% % (42) C_obs

% % C_obs = C/pastCminus;

F44_aux(C_obs_ind,end-grid.oc+C_ind)         = - (C_aux-C_ss)/perturb_size;
F44_aux(C_obs_ind,end-grid.oc+C_obs_ind)     = (C_obs_aux-C_obs_ss)/perturb_size;

% % (43) I_obs

% % I_obs = I/pastIminus;

F44_aux(I_obs_ind,end-grid.oc+I_ind)         = - (I_aux-I_ss)/perturb_size;
F44_aux(I_obs_ind,end-grid.oc+I_obs_ind)     = (I_obs_aux-I_obs_ss)/perturb_size;

% % (44) w_obs

% % w_obs = W/Wminus;

F41_aux(w_obs_ind,end-grid.os+w_ind)         = - (w_aux -  w_ss )/perturb_size;
F44_aux(w_obs_ind,end-grid.oc+w_obs_ind)     = (w_obs_aux-w_obs_ss)/perturb_size;

F44_aux(G_obs_ind,end-grid.oc+G_obs_ind)     = (G_obs_aux-G_obs_ss)/perturb_size;
F44_aux(G_obs_ind,end-grid.oc+G_ind)         = - (G_aux-G_ss)/perturb_size;

% % (45) PI_obs

% % PI_obs = PROFIT/pastPROFITminus;

F44_aux(PROFIT_obs_ind,end-grid.oc+Profit_ind)      = - ( Profit_aux - Profit_ss )/perturb_size;
F44_aux(PROFIT_obs_ind,end-grid.oc+PROFIT_obs_ind)  = ( PROFIT_obs_aux - PROFIT_obs_ss )/perturb_size;

F44_aux(unemp_obs_ind,end-grid.oc+unemp_ind)      = - ( 100*unemp_aux - 100*unemp_ss )/perturb_size;
F44_aux(unemp_obs_ind,end-grid.oc+unemp_obs_ind)  = ( unemp_obs_aux - unemp_obs_ss )/perturb_size;

F44_aux(inf_obs_ind,end-grid.oc+pi_ind)       = - ( pi_aux - pi_ss )/perturb_size;
F44_aux(inf_obs_ind,end-grid.oc+inf_obs_ind)  = ( inf_obs_aux - inf_obs_ss )/perturb_size;

F41_aux(R_obs_ind,end-grid.os+R_cb_ind)   = - ( R_cb_aux - R_cb_ss )/perturb_size;
F44_aux(R_obs_ind,end-grid.oc+R_obs_ind)  = ( R_obs_aux - R_obs_ss )/perturb_size;

F43_aux(pastw_ind,end-grid.os+w_ind)      = - ( w_aux - w_ss )/perturb_size;
F44_aux(pastw_ind,end-grid.oc+pastw_ind)  = ( pastw_aux - pastw_ss )/perturb_size;


F43_aux(pastYY_ind,end-grid.os+pastY_ind)   = - ( pastY_aux - pastY_ss )/perturb_size;
F44_aux(pastYY_ind,end-grid.oc+pastYY_ind)  = ( pastYY_aux - pastYY_ss )/perturb_size;

F43_aux(pastCC_ind,end-grid.os+pastC_ind)   = - ( pastC_aux - pastC_ss )/perturb_size;
F44_aux(pastCC_ind,end-grid.oc+pastCC_ind)  = ( pastCC_aux - pastCC_ss )/perturb_size;

F43_aux(pastII_ind,end-grid.os+pastI_ind)   = - ( pastI_aux - pastI_ss )/perturb_size;
F44_aux(pastII_ind,end-grid.oc+pastII_ind)  = ( pastII_aux - pastII_ss )/perturb_size;

F43_aux(pastPPROFIT_ind,end-grid.os+pastPROFIT_ind)   = - ( pastPROFIT_aux - pastPROFIT_ss )/perturb_size;
F44_aux(pastPPROFIT_ind,end-grid.oc+pastPPROFIT_ind)  = ( pastPPROFIT_aux - pastPPROFIT_ss )/perturb_size;

F43_aux(pastuu_ind,end-grid.os+pastunemp_ind) = - ( pastunemp_aux - pastunemp_ss )/perturb_size;
F44_aux(pastuu_ind,end-grid.oc+pastuu_ind)    = ( pastuu_aux - pastuu_ss )/perturb_size;

switch(param.adjust)
	case('G')
		F43_aux(pastGG_ind,end-grid.os+pastLT_ind) = - ( pastLT_aux - pastLT_ss )/perturb_size;
	case('LT')
		F43_aux(pastGG_ind,end-grid.os+pastG_ind)  = - ( pastG_aux - pastG_ss )/perturb_size;
end

F44_aux(pastGG_ind,end-grid.oc+pastGG_ind) = ( pastGG_aux - pastGG_ss )/perturb_size;

F43_aux(pastA_g_ind,end-grid.os+A_g_ind)     = - ( A_g_aux - A_g_ss )/perturb_size;
F44_aux(pastA_g_ind,end-grid.oc+pastA_g_ind) = ( pastA_g_aux - pastA_g_ss )/perturb_size;

F44_aux(l_lambda_ind,end-grid.oc+l_lambda_ind) = ( l_lambda_aux - l_lambda_ss )/perturb_size;

F43_aux(pastQ_ind,end-grid.os+Q_ind)     = - ( Q_aux - Q_ss )/perturb_size;
F44_aux(pastQ_ind,end-grid.oc+pastQ_ind) = ( pastQ_aux - pastQ_ss )/perturb_size;

F43_aux(x_I_ind,end-grid.os+pastI_ind) = - ( I_ss/pastI_aux - I_ss/pastI_ss )/perturb_size;
F44_aux(x_I_ind,end-grid.oc+I_ind)     = - ( I_aux/pastI_ss - I_ss/pastI_ss )/perturb_size;
F44_aux(x_I_ind,end-grid.oc+x_I_ind)   = ( x_I_aux - x_I_ss )/perturb_size;

F41_aux(eta2_ind,end-grid.os+eta_ind)   = - ( eta_aux - eta_ss )/perturb_size;
F44_aux(eta2_ind,end-grid.oc+eta2_ind)  = ( eta2_aux - eta2_ss )/perturb_size;

F41_aux(iota2_ind,end-grid.os+iota_ind)  = - ( iota_aux - iota_ss )/perturb_size;
F44_aux(iota2_ind,end-grid.oc+iota2_ind) = ( iota2_aux - iota2_ss )/perturb_size;

F44_aux(pastLT2_ind,end-grid.oc+pastLT2_ind)  = ( pastLT2_aux - pastLT2_ss )/perturb_size;
F43_aux(pastLT2_ind,end-grid.os+pastLT_ind)   = - ( pastLT_aux - pastLT_ss )/perturb_size;

F44_aux(pastG2_ind,end-grid.oc+pastG2_ind)  = ( pastG2_aux - pastG2_ss )/perturb_size;
F43_aux(pastG2_ind,end-grid.os+pastG_ind)   = - ( pastG_aux - pastG_ss )/perturb_size;

F44_aux(LT_obs_ind,end-grid.oc+LT_obs_ind)  = ( LT_obs_aux - LT_obs_ss )/perturb_size;
F44_aux(LT_obs_ind,end-grid.oc+LT_ind)      = - ( LT_aux - LT_ss )/perturb_size;

F43_aux(B_gov_ncp2_ind,end-grid.os+A_b_ind)       = - ( (B_b_ss + B_hh_ss + SS_stats.B_F - (Q_ss*A_b_aux- NW_b_ss) - Q_ss *(1+param.tau_cp)* (A_g_ss)) - ( (B_b_ss + B_hh_ss + SS_stats.B_F - (Q_ss*A_b_ss- NW_b_ss) - Q_ss *(1+param.tau_cp)* (A_g_ss))  )  )/perturb_size; 
F43_aux(B_gov_ncp2_ind,end-grid.os+B_b_ind)       = - ( (B_b_aux + B_hh_ss + SS_stats.B_F - (Q_ss*A_b_ss- NW_b_ss) - Q_ss *(1+param.tau_cp)* (A_g_ss)) - ( (B_b_ss + B_hh_ss + SS_stats.B_F - (Q_ss*A_b_ss- NW_b_ss) - Q_ss *(1+param.tau_cp)* (A_g_ss))  )  )/perturb_size;
F43_aux(B_gov_ncp2_ind,end-grid.os+A_g_ind)       = - ( (B_b_ss + B_hh_ss + SS_stats.B_F - (Q_ss*A_b_ss- NW_b_ss) - Q_ss *(1+param.tau_cp)* (A_g_aux)) - ( (B_b_ss + B_hh_ss + SS_stats.B_F - (Q_ss*A_b_ss- NW_b_ss) - Q_ss *(1+param.tau_cp)* (A_g_ss))  )  )/perturb_size;
F43_aux(B_gov_ncp2_ind,end-grid.os+Q_ind)         = - ( (B_b_ss + B_hh_ss + SS_stats.B_F - (Q_aux*A_b_ss- NW_b_ss) - Q_aux *(1+param.tau_cp)* (A_g_ss)) - ( (B_b_ss + B_hh_ss + SS_stats.B_F - (Q_ss*A_b_ss- NW_b_ss) - Q_ss *(1+param.tau_cp)* (A_g_ss))  )  )/perturb_size;
F43_aux(B_gov_ncp2_ind,end-grid.os+NW_b_ind)      = - ( (B_b_ss + B_hh_ss + SS_stats.B_F - (Q_ss*A_b_ss- NW_b_aux) - Q_ss *(1+param.tau_cp)* (A_g_ss)) - ( (B_b_ss + B_hh_ss + SS_stats.B_F - (Q_ss*A_b_ss- NW_b_ss) - Q_ss *(1+param.tau_cp)* (A_g_ss))  )  )/perturb_size;
F44_aux(B_gov_ncp2_ind,end-grid.oc+B_hh_ind)      = - (B_hh_aux-B_hh_ss)/perturb_size;
F44_aux(B_gov_ncp2_ind,end-grid.oc+B_gov_ncp2_ind) = (B_gov_ncp2_aux-B_gov_ncp2_ss)/perturb_size;



F44_aux(w2_ind,end-grid.oc+w2_ind)        = (w2_aux-w2_ss)/perturb_size;
F43_aux(w2_ind,end-grid.os+pastpi_ind)    = - (  param.w_bar*((param.pi_bar/pi_ss )^param.d*(pastpi_aux/pi_ss)^(1-param.d))^param.rho_w - param.w_bar*((param.pi_bar/pi_ss)^param.d*(pastpi_ss/pi_ss)^(1-param.d))^param.rho_w )/perturb_size;
F44_aux(w2_ind,end-grid.oc+pi_ind)        = - (  param.w_bar*((param.pi_bar/pi_aux)^param.d*(pastpi_ss/pi_aux)^(1-param.d))^param.rho_w - param.w_bar*((param.pi_bar/pi_ss)^param.d*(pastpi_ss/pi_ss)^(1-param.d))^param.rho_w )/perturb_size;

F44_aux(Profit2_ind,end-grid.oc+Profit2_ind)   = (Profit2_aux-Profit2_ss)/perturb_size;
F44_aux(Profit2_ind,end-grid.oc+w2_ind)        = - ( ((Y_ss *(1-param.eta/(2*param.kappa)*(log(pi_ss)-(1-param.GAMMA)*log(param.pi_cb)-param.GAMMA*log(pastpi_ss)).^2)+ Y_ss.*(-MC_ss) - param.fix  ...
		                						    + (h_ss-w2_aux).*L_ss - param.iota.*V_ss)  ...
		                							+ (r_k_ss*v_ss - param.delta_0*v_ss^param.delta_1).*K_ss +  Q_ss*(K_ss-K_ss)-(K_ss-K_ss) - param.phi/2*(K_ss/K_ss-1)^2*K_ss ) - ...
													((Y_ss *(1-param.eta/(2*param.kappa)*(log(pi_ss)-(1-param.GAMMA)*log(param.pi_cb)-param.GAMMA*log(pastpi_ss)).^2)+ Y_ss.*(-MC_ss) - param.fix  ...
		                						    + (h_ss-w2_ss).*L_ss - param.iota.*V_ss)  ...
		                							+ (r_k_ss*v_ss - param.delta_0*v_ss^param.delta_1).*K_ss +  Q_ss*(K_ss-K_ss)-(K_ss-K_ss) - param.phi/2*(K_ss/K_ss-1)^2*K_ss ) )/perturb_size;
F43_aux(Profit2_ind,end-grid.os+pastpi_ind)    = - ( ((Y_ss *(1-param.eta/(2*param.kappa)*(log(pi_ss)-(1-param.GAMMA)*log(param.pi_cb)-param.GAMMA*log(pastpi_aux)).^2)+ Y_ss.*(-MC_ss) - param.fix  ...
		                						    + (h_ss-w_ss).*L_ss - param.iota.*V_ss)  ...
		                							+ (r_k_ss*v_ss - param.delta_0*v_ss^param.delta_1).*K_ss +  Q_ss*(K_ss-K_ss)-(K_ss-K_ss) - param.phi/2*(K_ss/K_ss-1)^2*K_ss ) - ...
													((Y_ss *(1-param.eta/(2*param.kappa)*(log(pi_ss)-(1-param.GAMMA)*log(param.pi_cb)-param.GAMMA*log(pastpi_ss)).^2)+ Y_ss.*(-MC_ss) - param.fix  ...
		                						    + (h_ss-w_ss).*L_ss - param.iota.*V_ss)  ...
		                							+ (r_k_ss*v_ss - param.delta_0*v_ss^param.delta_1).*K_ss +  Q_ss*(K_ss-K_ss)-(K_ss-K_ss) - param.phi/2*(K_ss/K_ss-1)^2*K_ss ) )/perturb_size;
F44_aux(Profit2_ind,end-grid.oc+pi_ind)        = - ( ((Y_ss *(1-param.eta/(2*param.kappa)*(log(pi_aux)-(1-param.GAMMA)*log(param.pi_cb)-param.GAMMA*log(pastpi_ss)).^2)+ Y_ss.*(-MC_ss) - param.fix  ...
		                						    + (h_ss-w_ss).*L_ss - param.iota.*V_ss)  ...
		                							+ (r_k_ss*v_ss - param.delta_0*v_ss^param.delta_1).*K_ss +  Q_ss*(K_ss-K_ss)-(K_ss-K_ss) - param.phi/2*(K_ss/K_ss-1)^2*K_ss ) - ...
													((Y_ss *(1-param.eta/(2*param.kappa)*(log(pi_ss)-(1-param.GAMMA)*log(param.pi_cb)-param.GAMMA*log(pastpi_ss)).^2)+ Y_ss.*(-MC_ss) - param.fix  ...
		                						    + (h_ss-w_ss).*L_ss - param.iota.*V_ss)  ...
		                							+ (r_k_ss*v_ss - param.delta_0*v_ss^param.delta_1).*K_ss +  Q_ss*(K_ss-K_ss)-(K_ss-K_ss) - param.phi/2*(K_ss/K_ss-1)^2*K_ss ) )/perturb_size;

F44_aux(Profit3_ind,end-grid.oc+Profit3_ind)   = (Profit3_aux-Profit3_ss)/perturb_size;

F43_aux(Profit3_ind,end-grid.os+pastpi_ind)    = - ( ((Y_ss *(1-param.eta/(2*param.kappa)*(log(pi_ss)-(1-param.GAMMA)*log(param.pi_cb)-param.GAMMA*log(pastpi_aux)).^2)+ Y_ss.*(-MC_ss) - param.fix  ...
		                						    + (h_ss-w_ss).*L_ss - param.iota.*V_ss)  ...
		                							+ (r_k_ss*v_ss - param.delta_0*v_ss^param.delta_1).*K_ss +  Q_ss*(K_ss-K_ss)-(K_ss-K_ss) - param.phi/2*(K_ss/K_ss-1)^2*K_ss ) - ...
													((Y_ss *(1-param.eta/(2*param.kappa)*(log(pi_ss)-(1-param.GAMMA)*log(param.pi_cb)-param.GAMMA*log(pastpi_ss)).^2)+ Y_ss.*(-MC_ss) - param.fix  ...
		                						    + (h_ss-w_ss).*L_ss - param.iota.*V_ss)  ...
		                							+ (r_k_ss*v_ss - param.delta_0*v_ss^param.delta_1).*K_ss +  Q_ss*(K_ss-K_ss)-(K_ss-K_ss) - param.phi/2*(K_ss/K_ss-1)^2*K_ss ) )/perturb_size;
F44_aux(Profit3_ind,end-grid.oc+pi_ind)        = - ( ((Y_ss *(1-param.eta/(2*param.kappa)*(log(pi_aux)-(1-param.GAMMA)*log(param.pi_cb)-param.GAMMA*log(pastpi_ss)).^2)+ Y_ss.*(-MC_ss) - param.fix  ...
		                						    + (h_ss-w_ss).*L_ss - param.iota.*V_ss)  ...
		                							+ (r_k_ss*v_ss - param.delta_0*v_ss^param.delta_1).*K_ss +  Q_ss*(K_ss-K_ss)-(K_ss-K_ss) - param.phi/2*(K_ss/K_ss-1)^2*K_ss ) - ...
													((Y_ss *(1-param.eta/(2*param.kappa)*(log(pi_ss)-(1-param.GAMMA)*log(param.pi_cb)-param.GAMMA*log(pastpi_ss)).^2)+ Y_ss.*(-MC_ss) - param.fix  ...
		                						    + (h_ss-w_ss).*L_ss - param.iota.*V_ss)  ...
		                							+ (r_k_ss*v_ss - param.delta_0*v_ss^param.delta_1).*K_ss +  Q_ss*(K_ss-K_ss)-(K_ss-K_ss) - param.phi/2*(K_ss/K_ss-1)^2*K_ss ) )/perturb_size;

F44_aux(Profit_FI2_ind,end-grid.oc+Profit_FI2_ind) = (Profit_FI2_aux-Profit_FI2_ss)/perturb_size;
F44_aux(Profit_FI2_ind,end-grid.oc+pi_ind)        = - ( ((1-param.theta_b)*((RRa_ss  - param.b_a_aux2*R_cb_ss/pi_aux)*lev_ss +param.b_a_aux2*R_cb_ss/pi_aux)*NW_b_ss  - param.omega*Q_ss*A_b_ss)  - ((1-param.theta_b)*((RRa_ss - param.b_a_aux2*R_cb_ss/pi_ss)*lev_ss+R_cb_ss/pi_ss)*NW_b_ss - param.omega*Q_ss*A_b_ss)     )/perturb_size;

F44_aux(r_a2_ind,end-grid.oc+r_a2_ind)             = (r_a2_aux-r_a2_ss)/perturb_size;
F44_aux(r_a2_ind,end-grid.oc+Profit2_ind)          = - ((1-tau_a_ss)*(1-param.Eratio-param.b_share)*Profit2_aux/K_ss - (1-tau_a_ss)*(1-param.Eratio-param.b_share)*Profit2_ss/K_ss )/perturb_size;

out.F21_aux = F21_aux;
out.F22_aux = F22_aux;
out.F23_aux = F23_aux;
out.F24_aux = F24_aux;
out.F41_aux = F41_aux;
out.F42_aux = F42_aux;
out.F43_aux = F43_aux;
out.F44_aux = F44_aux;