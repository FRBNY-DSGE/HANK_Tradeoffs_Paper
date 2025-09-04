function [Difference,LHS,RHS,MUnext,c_a_star,b_a_star,a_a_star,c_n_star,b_n_star,AC,AProb,P_transition_exp,P_transition_dist] = F_sys_ref_anal_ELB(State,Stateminus,...
		Controlnext_sparse,Control_sparse,StateSS,...
		ControlSS,Gamma_state,Gamma_control,InvGamma,...
		param,grid,SS_stats,Copula,P_SE)

% System of equations written in Schmitt-Groh?-Uribe generic form with states and controls
% STATE: Vector of state variables t+1 (only marginal distributions for histogram)
% STATEMINUS: Vector of state variables t (only marginal distributions for histogram)
% CONTROL: Vector of state variables t+1 (only coefficients of sparse polynomial)
% CONTROLMINUS: Vector of state variables t (only coefficients of sparse polynomial)
% STATESS and CONTROLSS: Value of the state and control variables in steady
% state. For the Value functions these are at full grids.
% GAMMA_STATE: Mapping such that perturbationof marginals are still
% distributions (sum to 1).
% PARAM: Model and numerical parameters (structure)
% GRID: Liquid, illiquid and productivity grid
% SS_stats: Stores targets for government policy
% COPULA: Interpolant that allows to map marginals back to full-grid
% distribuitions
% P: steady state transition matrix
% aggshock: sets whether the Aggregate shock is Z or uncertainty
%
% =========================================================================
% Part of the Matlab code to aPSI_RPompany the paper
% 'Solving heterogeneous agent models in discrete time with
% many idiosyncratic states by perturbation methods', CEPR Discussion Paper
% DP13071
% by Christian Bayer and Ralph Luetticke
% http://www.erc-hump-uni-bonn.net/
% =========================================================================

%% Initializations

util     = @(c)  (c.^(1-param.sigma))./(1-param.sigma);
mutil    = @(c)  (1./(c.^param.sigma));
invutil  = @(u)  (((1-param.sigma).*u).^(1/(1-param.sigma)));
invmutil = @(mu) ((1./mu).^(1/param.sigma));

% Number of states, controls
nx   = grid.numstates; % Number of states
ny   = length(ControlSS); % number of Controls
NxNx = nx - grid.os; % Number of states without aggregate

NN   = grid.nb*grid.na*grid.nse; % Number of points in the full grid
Ny   = 3*NN;

%  Initialize LHS and RHS

LHS = zeros(nx+ny,1);
RHS = zeros(nx+ny,1);  

%  Indices for LHS/RHS

%  Indices for Controlsx

VALUE_ind   = 1:NN;
mutil_c_ind = NN + (1:NN);
Va_ind      = 2*NN + (1:NN);

% Summary variables

Profit_aux_ind       = Ny             + 1; % 1
GiniW_ind            = Profit_aux_ind + 1; % 2
GiniI_1_ind          = GiniW_ind      + 1; % 3
GiniI_2_ind          = GiniI_1_ind    + 1; % 4
GiniC_ind            = GiniI_2_ind    + 1; % 5
w_01_ind             = GiniC_ind      + 1; % 6
w_1_ind              = w_01_ind       + 1; % 7
w_10_ind             = w_1_ind        + 1; % 8
w_Q1_ind             = w_10_ind       + 1; % 9
w_Q2_ind             = w_Q1_ind       + 1; % 10
w_Q3_ind             = w_Q2_ind       + 1; % 11
w_Q4_ind             = w_Q3_ind       + 1; % 12
w_Q5_ind             = w_Q4_ind       + 1; % 13
w_90_ind             = w_Q5_ind       + 1; % 14
w_99_ind             = w_90_ind       + 1; % 15
w_999_ind            = w_99_ind       + 1; % 16
I_1_01_ind           = w_999_ind      + 1; % 17
I_1_1_ind            = I_1_01_ind     + 1; % 18
I_1_10_ind           = I_1_1_ind      + 1; % 19
I_1_Q1_ind           = I_1_10_ind     + 1; % 20
I_1_Q2_ind           = I_1_Q1_ind     + 1; % 21
I_1_Q3_ind           = I_1_Q2_ind     + 1; % 22
I_1_Q4_ind           = I_1_Q3_ind     + 1; % 23
I_1_Q5_ind           = I_1_Q4_ind     + 1; % 24
I_1_90_ind           = I_1_Q5_ind     + 1; % 25
I_1_99_ind           = I_1_90_ind     + 1; % 26
I_1_999_ind          = I_1_99_ind     + 1; % 27
I_2_01_ind           = I_1_999_ind    + 1; % 28
I_2_1_ind            = I_2_01_ind     + 1; % 29
I_2_10_ind           = I_2_1_ind      + 1; % 30
I_2_Q1_ind           = I_2_10_ind     + 1; % 31
I_2_Q2_ind           = I_2_Q1_ind     + 1; % 32
I_2_Q3_ind           = I_2_Q2_ind     + 1; % 33
I_2_Q4_ind           = I_2_Q3_ind     + 1; % 34
I_2_Q5_ind           = I_2_Q4_ind     + 1; % 35
I_2_90_ind           = I_2_Q5_ind     + 1; % 36
I_2_99_ind           = I_2_90_ind     + 1; % 37
I_2_999_ind          = I_2_99_ind     + 1; % 38
C_01_ind             = I_2_999_ind    + 1; % 39
C_1_ind              = C_01_ind       + 1; % 40
C_10_ind             = C_1_ind        + 1; % 41
C_Q1_ind             = C_10_ind       + 1; % 42
C_Q2_ind             = C_Q1_ind       + 1; % 43
C_Q3_ind             = C_Q2_ind       + 1; % 44
C_Q4_ind             = C_Q3_ind       + 1; % 45
C_Q5_ind             = C_Q4_ind       + 1; % 46
C_90_ind             = C_Q5_ind       + 1; % 47
C_99_ind             = C_90_ind       + 1; % 48
C_999_ind            = C_99_ind       + 1; % 49
w_I_01_ind           = C_999_ind     + 1; % 50
w_I_1_ind            = w_I_01_ind    + 1; % 51
w_I_10_ind           = w_I_1_ind     + 1; % 52
w_I_Q1_ind           = w_I_10_ind    + 1; % 53
w_I_Q2_ind           = w_I_Q1_ind    + 1; % 54
w_I_Q3_ind           = w_I_Q2_ind    + 1; % 55
w_I_Q4_ind           = w_I_Q3_ind    + 1; % 56
w_I_Q5_ind           = w_I_Q4_ind    + 1; % 57
w_I_90_ind           = w_I_Q5_ind    + 1; % 58
w_I_99_ind           = w_I_90_ind    + 1; % 59
w_I_999_ind          = w_I_99_ind    + 1; % 60
I_I_1_01_ind         = w_I_999_ind      + 1; % 61
I_I_1_1_ind          = I_I_1_01_ind       + 1; % 62
I_I_1_10_ind         = I_I_1_1_ind        + 1; % 63
I_I_1_Q1_ind         = I_I_1_10_ind       + 1; % 64
I_I_1_Q2_ind         = I_I_1_Q1_ind       + 1; % 65
I_I_1_Q3_ind         = I_I_1_Q2_ind       + 1; % 66
I_I_1_Q4_ind         = I_I_1_Q3_ind       + 1; % 67
I_I_1_Q5_ind         = I_I_1_Q4_ind       + 1; % 68
I_I_1_90_ind         = I_I_1_Q5_ind       + 1; % 69
I_I_1_99_ind         = I_I_1_90_ind       + 1; % 70
I_I_1_999_ind        = I_I_1_99_ind       + 1; % 71
I_I_2_01_ind         = I_I_1_999_ind    + 1; % 72
I_I_2_1_ind          = I_I_2_01_ind     + 1; % 73
I_I_2_10_ind         = I_I_2_1_ind      + 1; % 74
I_I_2_Q1_ind         = I_I_2_10_ind     + 1; % 75
I_I_2_Q2_ind         = I_I_2_Q1_ind     + 1; % 76
I_I_2_Q3_ind         = I_I_2_Q2_ind     + 1; % 77
I_I_2_Q4_ind         = I_I_2_Q3_ind     + 1; % 78
I_I_2_Q5_ind         = I_I_2_Q4_ind     + 1; % 79
I_I_2_90_ind         = I_I_2_Q5_ind     + 1; % 80
I_I_2_99_ind         = I_I_2_90_ind     + 1; % 81
I_I_2_999_ind        = I_I_2_99_ind     + 1; % 82
C_I_01_ind           = I_I_2_999_ind    + 1; % 83
C_I_1_ind            = C_I_01_ind       + 1; % 84
C_I_10_ind           = C_I_1_ind        + 1; % 85
C_I_Q1_ind           = C_I_10_ind       + 1; % 86
C_I_Q2_ind           = C_I_Q1_ind       + 1; % 87
C_I_Q3_ind           = C_I_Q2_ind       + 1; % 88
C_I_Q4_ind           = C_I_Q3_ind       + 1; % 89
C_I_Q5_ind           = C_I_Q4_ind       + 1; % 90
C_I_90_ind           = C_I_Q5_ind       + 1; % 91
C_I_99_ind           = C_I_90_ind       + 1; % 92
C_I_999_ind          = C_I_99_ind       + 1; % 93
I_90to10_ind         = C_I_999_ind    + 1; % 94
I_50to10_ind         = I_90to10_ind   + 1; % 95
I_90to50_ind         = I_50to10_ind   + 1; % 96
W_90to50_ind         = I_90to50_ind       + 1; % 97   %Differentiation between I/W and Inc/Wlth being between threshold vals (income50/w50 vs avg income_P90/B_10)
Inc_1_T10toB10_ind   = W_90to50_ind       + 1; % 98
Inc_1_T10toQ3_ind    = Inc_1_T10toB10_ind + 1; % 99
Inc_1_Q3toB10_ind    = Inc_1_T10toQ3_ind  + 1; % 100
Inc_2_T10toB10_ind   = Inc_1_Q3toB10_ind  + 1; % 101
Inc_2_T10toQ3_ind    = Inc_2_T10toB10_ind + 1; % 102
Inc_2_Q3toB10_ind    = Inc_2_T10toQ3_ind  + 1; % 103
Wlth_T10toQ3_ind     = Inc_2_Q3toB10_ind  + 1; % 104
C_T10toB10_ind       = Wlth_T10toQ3_ind   + 1; % 105
C_T10toQ3_ind        = C_T10toB10_ind     + 1; % 106
C_Q3toB10_ind        = C_T10toQ3_ind      + 1; % 107
Inc_I_1_T10toB10_ind = C_Q3toB10_ind        + 1; % 108
Inc_I_1_T10toQ3_ind  = Inc_I_1_T10toB10_ind + 1; % 109
Inc_I_1_Q3toB10_ind  = Inc_I_1_T10toQ3_ind  + 1; % 110
Inc_I_2_T10toB10_ind = Inc_I_1_Q3toB10_ind  + 1; % 111
Inc_I_2_T10toQ3_ind  = Inc_I_2_T10toB10_ind + 1; % 112
Inc_I_2_Q3toB10_ind  = Inc_I_2_T10toQ3_ind  + 1; % 113
Wlth_I_T10toQ3_ind   = Inc_I_2_Q3toB10_ind  + 1; % 114
C_I_T10toB10_ind     = Wlth_I_T10toQ3_ind   + 1; % 115
C_I_T10toQ3_ind      = C_I_T10toB10_ind     + 1; % 116
C_I_Q3toB10_ind      = C_I_T10toQ3_ind      + 1; % 117


MRS_ind       = C_I_Q3toB10_ind  + 1;  % 1
A_hh_ind      = MRS_ind       + 1;     % 2
B_hh_ind      = A_hh_ind      + 1;     % 3
C_ind         = B_hh_ind      + 1;     % 4
N_ind         = C_ind         + 1;     % 5
L_ind         = N_ind         + 1;     % 6
UB_ind        = L_ind         + 1;     % 7

% Other control variables

K_ind         = UB_ind        + 1;  % 1
B_ind         = K_ind         + 1;  %2 
B_gov_ncp_ind = B_ind         + 1;  %3
T_ind         = B_gov_ncp_ind + 1;  %4
LT_ind 	      = T_ind         + 1;  %5
G_ind         = LT_ind 	      + 1;  %6
Lambda_ind    = G_ind         + 1;  %7
pi_ind        = Lambda_ind    + 1;  %8
V_ind         = pi_ind        + 1;  %9
J_ind         = V_ind + (1:grid.ns)'; %10,11,12,13,14
h_ind         = J_ind(end)+1;         %15
v_ind         = h_ind         + 1;  %16
Y_ind         = v_ind         + 1;  %17
Profit_ind    = Y_ind         + 1;  %18
r_k_ind       = Profit_ind    + 1;  %19
r_a_ind       = r_k_ind       + 1;  %20
MC_ind        = r_a_ind       + 1;  %21
unemp_ind     = MC_ind        + 1;  %22
nn_ind        = unemp_ind     + 1;  %23
M_ind         = nn_ind        + 1;  %24
f_ind         = M_ind         + 1;  %25
zz_ind        = f_ind         + 1;  %26
xx_ind        = zz_ind        + 1;  %27
vv_ind        = xx_ind        + 1;  %28
ee_ind        = vv_ind        + 1;  %29
C_b_ind       = ee_ind        + 1;  %30
Profit_FI_ind = C_b_ind       + 1;  %31
RRa_ind       = Profit_FI_ind + 1;  %32
RR_ind        = RRa_ind       + 1;  %33
I_ind         = RR_ind        + 1;  %34
x_k_ind       = I_ind         + 1;  %35

% B_F_ind       = x_k_ind       + 1;

% Observables

A_g_obs_ind    = x_k_ind        + 1; % 36
Y_obs_ind      = A_g_obs_ind    + 1; % 37
C_obs_ind      = Y_obs_ind      + 1; % 38
I_obs_ind      = C_obs_ind      + 1; % 39
w_obs_ind      = I_obs_ind      + 1; % 40
PROFIT_obs_ind = w_obs_ind      + 1; % 41
unemp_obs_ind  = PROFIT_obs_ind + 1; % 42
inf_obs_ind    = unemp_obs_ind  + 1;
R_obs_ind      = inf_obs_ind    + 1;
pastw_ind      = R_obs_ind      + 1;
G_obs_ind      = pastw_ind      + 1;

pastYY_ind        = G_obs_ind       + 1;
pastCC_ind        = pastYY_ind      + 1; 
pastII_ind        = pastCC_ind      + 1; 
pastPPROFIT_ind   = pastII_ind      + 1;
pastuu_ind        = pastPPROFIT_ind + 1;
pastGG_ind        = pastuu_ind      + 1;
pastA_g_ind       = pastGG_ind      + 1;
l_lambda_ind      = pastA_g_ind     + 1;
pastQ_ind         = l_lambda_ind    + 1;
x_I_ind           = pastQ_ind       + 1;
eta2_ind          = x_I_ind         + 1;
iota2_ind         = eta2_ind        + 1;
pastLT2_ind       = iota2_ind       + 1;
pastG2_ind        = pastLT2_ind     + 1;
LT_obs_ind        = pastG2_ind      + 1;
B_gov_ncp2_ind    = LT_obs_ind      + 1;

w2_ind            = B_gov_ncp2_ind  + 1;
Profit2_ind       = w2_ind          + 1;
Profit3_ind       = Profit2_ind     + 1;
Profit_FI2_ind    = Profit3_ind     + 1;
r_a2_ind          = Profit_FI2_ind  + 1;


%  Indices for States

marginal_b_ind  = 1:grid.nb-1;
marginal_a_ind  = (grid.nb-1+(1:(grid.na-1)));
marginal_se_ind = (grid.nb+grid.na-2 + (1:(grid.nse-1)));

R_cb_ind        = NxNx+1;
w_ind           = NxNx+2;
A_b_ind         = NxNx+3;
B_b_ind         = NxNx+4;
A_g_ind         = NxNx+5;
Q_ind           = NxNx+6;
lev_ind         = NxNx+7;
NW_b_ind        = NxNx+8;
R_tilde_ind     = NxNx+9;
x_cb_ind        = NxNx+10;   
pastpi_ind      = NxNx+11;

pastY_ind       = NxNx+12;
pastC_ind       = NxNx+13;
pastI_ind       = NxNx+14;
pastPROFIT_ind  = NxNx+15;
pastunemp_ind   = NxNx+16;
pastG_ind       = NxNx+17;
pastLT_ind      = NxNx+18;
R_star_ind      = NxNx+19;

B_F_ind         = NxNx+20;
Z_ind           = NxNx+21; 
PSI_RP_ind      = NxNx+22; 
eta_ind         = NxNx+23; 
D_ind           = NxNx+24; 
GG_ind          = NxNx+25; 
iota_ind        = NxNx+26; 
BB_ind          = NxNx+27; 
PSI_W_ind       = NxNx+28;
MP_ind          = NxNx+29;
p_mark_ind      = NxNx+30;

pgap_ind        = NxNx+31;
ygap_ind        = NxNx+32;

eps_QE_ind      = NxNx+33; % 119
eps_RP_ind      = NxNx+34;  % Networth
eps_B_F_ind     = NxNx+35;
eps_BB_ind      = NxNx+36;  % Undertermined
eps_Z_ind       = NxNx+37;  % 1) TFP shock
eps_G_ind       = NxNx+38;  % 2) G shock
eps_D_ind       = NxNx+39;  % 3) Liquidity preference
eps_R_ind       = NxNx+40;  % 4) MP shock
eps_iota_ind    = NxNx+41;  % 5) MEI (or IST) shock
eps_eta_ind     = NxNx+42;  % 6) Price mark-up shock
eps_w_ind       = NxNx+43;  % 7) Wage (mark-up)  sjpcl

NxNx_aux = 43;


		


%% Control Variables (For value function and derivatives, DCT is used)

Controlnext  = ControlSS .* (1+Gamma_control*(Controlnext_sparse));
Control      = ControlSS .* (1+Gamma_control*(Control_sparse));

% Controlnext    = zeros(length(ControlSS));
% Control        = zeros(length(ControlSS));

Controlnext(end-grid.oc+1:end) = ControlSS(end-grid.oc+1:end) + Gamma_control(end-grid.oc+1:end,:)*(Controlnext_sparse);
Control(end-grid.oc+1:end)     = ControlSS(end-grid.oc+1:end) + Gamma_control(end-grid.oc+1:end,:)*(Control_sparse);

%% State Variables

%  Size of Distribution of nb+na+nse. Note that NxNx = nb+na+nse-3

Distribution       = StateSS(1:end-grid.os) + Gamma_state * State(1:NxNx);
Distributionminus  = StateSS(1:end-grid.os) + Gamma_state * Stateminus(1:NxNx);

%  Fiscal and Monetary Policy Shock

R_cb       		  = StateSS(end-(NxNx_aux-R_cb_ind+NxNx)) + (State(end-(NxNx_aux-R_cb_ind+NxNx)));
R_cbminus  		  = StateSS(end-(NxNx_aux-R_cb_ind+NxNx)) + (Stateminus(end-(NxNx_aux-R_cb_ind+NxNx)));
		
Wminus     		  = StateSS(end-(NxNx_aux-w_ind+NxNx))    + (Stateminus(end-(NxNx_aux-w_ind+NxNx)));
W          		  = StateSS(end-(NxNx_aux-w_ind+NxNx))    + (State(end-(NxNx_aux-w_ind+NxNx)));

A_bnext    		  = StateSS(end-(NxNx_aux-A_b_ind+NxNx)) + (State(end-(NxNx_aux-A_b_ind+NxNx)));
A_b        		  = StateSS(end-(NxNx_aux-A_b_ind+NxNx)) + (Stateminus(end-(NxNx_aux-A_b_ind+NxNx)));
		
B_bnext           = StateSS(end-(NxNx_aux-B_b_ind+NxNx)) + (State(end-(NxNx_aux-B_b_ind+NxNx)));
B_b               = StateSS(end-(NxNx_aux-B_b_ind+NxNx)) + (Stateminus(end-(NxNx_aux-B_b_ind+NxNx)));
		
A_gauxnext        = StateSS(end-(NxNx_aux-A_g_ind+NxNx)) + (State(end-(NxNx_aux-A_g_ind+NxNx)));
A_gaux            = StateSS(end-(NxNx_aux-A_g_ind+NxNx)) + (Stateminus(end-(NxNx_aux-A_g_ind+NxNx)));

Q                 = StateSS(end-(NxNx_aux-Q_ind+NxNx)) + (State(end-(NxNx_aux-Q_ind+NxNx)));
Qminus            = StateSS(end-(NxNx_aux-Q_ind+NxNx)) + (Stateminus(end-(NxNx_aux-Q_ind+NxNx)));
		
lev               = StateSS(end-(NxNx_aux-lev_ind+NxNx)) + (State(end-(NxNx_aux-lev_ind+NxNx)));
levminus          = StateSS(end-(NxNx_aux-lev_ind+NxNx)) + (Stateminus(end-(NxNx_aux-lev_ind+NxNx)));
		
NW_b              = StateSS(end-(NxNx_aux-NW_b_ind+NxNx)) + (State(end-(NxNx_aux-NW_b_ind+NxNx)));
NW_bminus         = StateSS(end-(NxNx_aux-NW_b_ind+NxNx)) + (Stateminus(end-(NxNx_aux-NW_b_ind+NxNx)));
		

R_tilde           = StateSS(end-(NxNx_aux-R_tilde_ind+NxNx)) + (State(end-(NxNx_aux-R_tilde_ind+NxNx)));
R_tildeminus      = StateSS(end-(NxNx_aux-R_tilde_ind+NxNx)) + (Stateminus(end-(NxNx_aux-R_tilde_ind+NxNx)));

x_cb              = StateSS(end-(NxNx_aux-x_cb_ind+NxNx)) + (State(end-(NxNx_aux-x_cb_ind+NxNx)));
x_cbminus         = StateSS(end-(NxNx_aux-x_cb_ind+NxNx)) + (Stateminus(end-(NxNx_aux-x_cb_ind+NxNx)));


pastpi            = StateSS(end-(NxNx_aux-pastpi_ind+NxNx)) + (State(end-(NxNx_aux-pastpi_ind+NxNx)));
pastpiminus       = StateSS(end-(NxNx_aux-pastpi_ind+NxNx)) + (Stateminus(end-(NxNx_aux-pastpi_ind+NxNx)));

pastY             = StateSS(end-(NxNx_aux-pastY_ind+NxNx)) + (State(end-(NxNx_aux-pastY_ind+NxNx)));
pastYminus        = StateSS(end-(NxNx_aux-pastY_ind+NxNx)) + (Stateminus(end-(NxNx_aux-pastY_ind+NxNx)));

pastC             = StateSS(end-(NxNx_aux-pastC_ind+NxNx)) + (State(end-(NxNx_aux-pastC_ind+NxNx)));
pastCminus        = StateSS(end-(NxNx_aux-pastC_ind+NxNx)) + (Stateminus(end-(NxNx_aux-pastC_ind+NxNx)));

pastI             = StateSS(end-(NxNx_aux-pastI_ind+NxNx)) + (State(end-(NxNx_aux-pastI_ind+NxNx)));
pastIminus        = StateSS(end-(NxNx_aux-pastI_ind+NxNx)) + (Stateminus(end-(NxNx_aux-pastI_ind+NxNx)));

pastPROFIT        = StateSS(end-(NxNx_aux-pastPROFIT_ind+NxNx)) + (State(end-(NxNx_aux-pastPROFIT_ind+NxNx)));
pastPROFITminus   = StateSS(end-(NxNx_aux-pastPROFIT_ind+NxNx)) + (Stateminus(end-(NxNx_aux-pastPROFIT_ind+NxNx)));

pastunemp         = StateSS(end-(NxNx_aux-pastunemp_ind+NxNx)) + (State(end-(NxNx_aux-pastunemp_ind+NxNx)));
pastunempminus    = StateSS(end-(NxNx_aux-pastunemp_ind+NxNx)) + (Stateminus(end-(NxNx_aux-pastunemp_ind+NxNx)));

pastG             = StateSS(end-(NxNx_aux-pastG_ind+NxNx)) + (State(end-(NxNx_aux-pastG_ind+NxNx)));
pastGminus        = StateSS(end-(NxNx_aux-pastG_ind+NxNx)) + (Stateminus(end-(NxNx_aux-pastG_ind+NxNx)));

pastLT             = StateSS(end-(NxNx_aux-pastLT_ind+NxNx)) + (State(end-(NxNx_aux-pastLT_ind+NxNx)));
pastLTminus        = StateSS(end-(NxNx_aux-pastLT_ind+NxNx)) + (Stateminus(end-(NxNx_aux-pastLT_ind+NxNx)));

R_star            = StateSS(end-(NxNx_aux-R_star_ind+NxNx)) + (State(end-(NxNx_aux-R_star_ind+NxNx)));
R_starminus       = StateSS(end-(NxNx_aux-R_star_ind+NxNx)) + (Stateminus(end-(NxNx_aux-R_star_ind+NxNx)));

B_Fnext           = StateSS(end-(NxNx_aux-B_F_ind+NxNx)) + (State(end-(NxNx_aux-B_F_ind+NxNx)));
B_F               = StateSS(end-(NxNx_aux-B_F_ind+NxNx)) + (Stateminus(end-(NxNx_aux-B_F_ind+NxNx)));

Z                 = StateSS(end-(NxNx_aux-Z_ind+NxNx)) + (State(end-(NxNx_aux-Z_ind+NxNx)));
Zminus            = StateSS(end-(NxNx_aux-Z_ind+NxNx)) + (Stateminus(end-(NxNx_aux-Z_ind+NxNx)));


PSI_RP            = StateSS(end-(NxNx_aux-PSI_RP_ind+NxNx)) + (State(end-(NxNx_aux-PSI_RP_ind+NxNx)));
PSI_RPminus       = StateSS(end-(NxNx_aux-PSI_RP_ind+NxNx)) + (Stateminus(end-(NxNx_aux-PSI_RP_ind+NxNx)));

eta               = StateSS(end-(NxNx_aux-eta_ind+NxNx)) + (State(end-(NxNx_aux-eta_ind+NxNx)));
etaminus          = StateSS(end-(NxNx_aux-eta_ind+NxNx)) + (Stateminus(end-(NxNx_aux-eta_ind+NxNx)));

D                 = StateSS(end-(NxNx_aux-D_ind+NxNx)) + (State(end-(NxNx_aux-D_ind+NxNx)));
Dminus            = StateSS(end-(NxNx_aux-D_ind+NxNx)) + (Stateminus(end-(NxNx_aux-D_ind+NxNx)));

GG                = StateSS(end-(NxNx_aux-GG_ind+NxNx)) + (State(end-(NxNx_aux-GG_ind+NxNx)));
GGminus           = StateSS(end-(NxNx_aux-GG_ind+NxNx)) + (Stateminus(end-(NxNx_aux-GG_ind+NxNx)));
		
iota              = StateSS(end-(NxNx_aux-iota_ind+NxNx)) + (State(end-(NxNx_aux-iota_ind+NxNx)));
iotaminus         = StateSS(end-(NxNx_aux-iota_ind+NxNx)) + (Stateminus(end-(NxNx_aux-iota_ind+NxNx)));

BB              = StateSS(end-(NxNx_aux-BB_ind+NxNx)) + (State(end-(NxNx_aux-BB_ind+NxNx)));
BBminus         = StateSS(end-(NxNx_aux-BB_ind+NxNx)) + (Stateminus(end-(NxNx_aux-BB_ind+NxNx)));

PSI_W              = StateSS(end-(NxNx_aux-PSI_W_ind+NxNx)) + (State(end-(NxNx_aux-PSI_W_ind+NxNx)));
PSI_Wminus         = StateSS(end-(NxNx_aux-PSI_W_ind+NxNx)) + (Stateminus(end-(NxNx_aux-PSI_W_ind+NxNx)));

MP              = StateSS(end-(NxNx_aux-MP_ind+NxNx)) + (State(end-(NxNx_aux-MP_ind+NxNx)));
MPminus         = StateSS(end-(NxNx_aux-MP_ind+NxNx)) + (Stateminus(end-(NxNx_aux-MP_ind+NxNx)));

p_mark            = StateSS(end-(NxNx_aux-p_mark_ind+NxNx)) + (State(end-(NxNx_aux-p_mark_ind+NxNx)));
p_markminus       = StateSS(end-(NxNx_aux-p_mark_ind+NxNx)) + (Stateminus(end-(NxNx_aux-p_mark_ind+NxNx)));

pgap            = StateSS(end-(NxNx_aux-pgap_ind+NxNx)) + (State(end-(NxNx_aux-pgap_ind+NxNx)));
pgapminus       = StateSS(end-(NxNx_aux-pgap_ind+NxNx)) + (Stateminus(end-(NxNx_aux-pgap_ind+NxNx)));

ygap            = StateSS(end-(NxNx_aux-ygap_ind+NxNx)) + (State(end-(NxNx_aux-ygap_ind+NxNx)));
ygapminus       = StateSS(end-(NxNx_aux-ygap_ind+NxNx)) + (Stateminus(end-(NxNx_aux-ygap_ind+NxNx)));
		
%  Shock

eps_B_Fnext        = StateSS(end-(NxNx_aux-eps_B_F_ind+NxNx)) + (State(end-(NxNx_aux-eps_B_F_ind+NxNx)));
eps_B_F            = StateSS(end-(NxNx_aux-eps_B_F_ind+NxNx)) + (Stateminus(end-(NxNx_aux-eps_B_F_ind+NxNx)));

eps_RPnext        = StateSS(end-(NxNx_aux-eps_RP_ind+NxNx)) + (State(end-(NxNx_aux-eps_RP_ind+NxNx)));
eps_RP            = StateSS(end-(NxNx_aux-eps_RP_ind+NxNx)) + (Stateminus(end-(NxNx_aux-eps_RP_ind+NxNx)));

eps_etanext       = StateSS(end-(NxNx_aux-eps_eta_ind+NxNx)) + (State(end-(NxNx_aux-eps_eta_ind+NxNx)));
eps_eta           = StateSS(end-(NxNx_aux-eps_eta_ind+NxNx)) + (Stateminus(end-(NxNx_aux-eps_eta_ind+NxNx)));

eps_Dnext         = StateSS(end-(NxNx_aux-eps_D_ind+NxNx)) + (State(end-(NxNx_aux-eps_D_ind+NxNx)));
eps_D             = StateSS(end-(NxNx_aux-eps_D_ind+NxNx)) + (Stateminus(end-(NxNx_aux-eps_D_ind+NxNx)));
		
eps_Rnext         = StateSS(end-(NxNx_aux-eps_R_ind+NxNx)) + (State(end-(NxNx_aux-eps_R_ind+NxNx)));
eps_R             = StateSS(end-(NxNx_aux-eps_R_ind+NxNx)) + (Stateminus(end-(NxNx_aux-eps_R_ind+NxNx)));
		
eps_Znext         = StateSS(end-(NxNx_aux-eps_Z_ind+NxNx)) + (State(end-(NxNx_aux-eps_Z_ind+NxNx)));
eps_Z             = StateSS(end-(NxNx_aux-eps_Z_ind+NxNx)) + (Stateminus(end-(NxNx_aux-eps_Z_ind+NxNx)));
		
eps_Gnext         = StateSS(end-(NxNx_aux-eps_G_ind+NxNx)) + (State(end-(NxNx_aux-eps_G_ind+NxNx)));
eps_G             = StateSS(end-(NxNx_aux-eps_G_ind+NxNx)) + (Stateminus(end-(NxNx_aux-eps_G_ind+NxNx)));

eps_iotanext      = StateSS(end-(NxNx_aux-eps_iota_ind+NxNx)) + (State(end-(NxNx_aux-eps_iota_ind+NxNx)));
eps_iota          = StateSS(end-(NxNx_aux-eps_iota_ind+NxNx)) + (Stateminus(end-(NxNx_aux-eps_iota_ind+NxNx)));

eps_QEnext        = StateSS(end-(NxNx_aux-eps_QE_ind+NxNx)) + (State(end-(NxNx_aux-eps_QE_ind+NxNx)));
eps_QE            = StateSS(end-(NxNx_aux-eps_QE_ind+NxNx)) + (Stateminus(end-(NxNx_aux-eps_QE_ind+NxNx)));

eps_wnext         = StateSS(end-(NxNx_aux-eps_w_ind+NxNx)) + (State(end-(NxNx_aux-eps_w_ind+NxNx)));
eps_w             = StateSS(end-(NxNx_aux-eps_w_ind+NxNx)) + (Stateminus(end-(NxNx_aux-eps_w_ind+NxNx)));

eps_BBnext        = StateSS(end-(NxNx_aux-eps_BB_ind+NxNx)) + (State(end-(NxNx_aux-eps_BB_ind+NxNx)));
eps_BB            = StateSS(end-(NxNx_aux-eps_BB_ind+NxNx)) + (Stateminus(end-(NxNx_aux-eps_BB_ind+NxNx)));

		

%% Split the Control vector into items with names

%  Controls

VALUEnext     = util(Controlnext(VALUE_ind));
VALUE         = util(Control(VALUE_ind));
Vanext        = mutil(Controlnext(Va_ind));
Va            = mutil(Control(Va_ind));
mutil_cnext   = mutil(Controlnext(mutil_c_ind));
mutil_c       = mutil(Control(mutil_c_ind));

%  Aggregate Controls (t+1);

A_hhnext      = exp(Controlnext(A_hh_ind));
B_hhnext      = exp(Controlnext(B_hh_ind));
Cnext         = exp(Controlnext(C_ind));
Knext         = exp(Controlnext(K_ind));
Bnext         = exp(Controlnext(B_ind));
Nnext         = exp(Controlnext(N_ind));

LAMBDAnext    = exp(Controlnext(Lambda_ind));
PInext        = exp(Controlnext(pi_ind));
Vnext         = exp(Controlnext(V_ind));
Jnext         = exp(Controlnext(J_ind));
Hnext         = exp(Controlnext(h_ind));
vnext         = exp(Controlnext(v_ind));
Ynext         = exp(Controlnext(Y_ind));
R_Knext       = exp(Controlnext(r_k_ind));
R_Anext       = exp(Controlnext(r_a_ind));
MCnext        = exp(Controlnext(MC_ind));
RRanext       = exp(Controlnext(RRa_ind));
RRnext        = exp(Controlnext(RR_ind));
Inext         = exp(Controlnext(I_ind));
unempnext     = exp(Controlnext(unemp_ind));
x_knext       = exp(Controlnext(x_k_ind));
nnnext        = exp(Controlnext(nn_ind));
Mnext         = exp(Controlnext(M_ind));
fnext         = exp(Controlnext(f_ind));
B_gov_ncpnext = exp(Controlnext(B_gov_ncp_ind));
UBnext        = exp(Controlnext(UB_ind));
l_lambdanext  = exp(Controlnext(l_lambda_ind));
x_Inext       = exp(Controlnext(x_I_ind));
eta2next      = exp(Controlnext(eta2_ind));
iota2next     = exp(Controlnext(iota2_ind));
B_gov_ncp2next= exp(Controlnext(B_gov_ncp2_ind));


%  Aggregate Controls (t);

Profit_aux = exp(Control(Profit_aux_ind));
GiniW = exp(Control(GiniW_ind));
GiniI_1 = exp(Control(GiniI_1_ind));
GiniI_2 = exp(Control(GiniI_2_ind));
GiniC = exp(Control(GiniC_ind));
w_01 = exp(Control(w_01_ind));
w_1 = exp(Control(w_1_ind));
w_10 = exp(Control(w_10_ind));
w_Q1 = exp(Control(w_Q1_ind));
w_Q2 = exp(Control(w_Q2_ind));
w_Q3 = exp(Control(w_Q3_ind));
w_Q4 = exp(Control(w_Q4_ind));
w_Q5 = exp(Control(w_Q5_ind));
w_90 = exp(Control(w_90_ind));
w_99 = exp(Control(w_99_ind));
w_999 = exp(Control(w_999_ind));
I_1_01 = exp(Control(I_1_01_ind));
I_1_1 = exp(Control(I_1_1_ind));
I_1_10 = exp(Control(I_1_10_ind));
I_1_Q1 = exp(Control(I_1_Q1_ind));
I_1_Q2 = exp(Control(I_1_Q2_ind));
I_1_Q3 = exp(Control(I_1_Q3_ind));
I_1_Q4 = exp(Control(I_1_Q4_ind));
I_1_Q5 = exp(Control(I_1_Q5_ind));
I_1_90 = exp(Control(I_1_90_ind));
I_1_99 = exp(Control(I_1_99_ind));
I_1_999 = exp(Control(I_1_999_ind));
I_2_01 = exp(Control(I_2_01_ind));
I_2_1 = exp(Control(I_2_1_ind));
I_2_10 = exp(Control(I_2_10_ind));
I_2_Q1 = exp(Control(I_2_Q1_ind));
I_2_Q2 = exp(Control(I_2_Q2_ind));
I_2_Q3 = exp(Control(I_2_Q3_ind));
I_2_Q4 = exp(Control(I_2_Q4_ind));
I_2_Q5 = exp(Control(I_2_Q5_ind));
I_2_90 = exp(Control(I_2_90_ind));
I_2_99 = exp(Control(I_2_99_ind));
I_2_999 = exp(Control(I_2_999_ind));
C_01 = exp(Control(C_01_ind));
C_1 = exp(Control(C_1_ind));
C_10 = exp(Control(C_10_ind));
C_Q1 = exp(Control(C_Q1_ind));
C_Q2 = exp(Control(C_Q2_ind));
C_Q3 = exp(Control(C_Q3_ind));
C_Q4 = exp(Control(C_Q4_ind));
C_Q5 = exp(Control(C_Q5_ind));
C_90 = exp(Control(C_90_ind));
C_99 = exp(Control(C_99_ind));
C_999 = exp(Control(C_999_ind));
w_I_01 = exp(Control(w_I_01_ind));
w_I_1 = exp(Control(w_I_1_ind));
w_I_10 = exp(Control(w_I_10_ind));
w_I_Q1 = exp(Control(w_I_Q1_ind));
w_I_Q2 = exp(Control(w_I_Q2_ind));
w_I_Q3 = exp(Control(w_I_Q3_ind));
w_I_Q4 = exp(Control(w_I_Q4_ind));
w_I_Q5 = exp(Control(w_I_Q5_ind));
w_I_90 = exp(Control(w_I_90_ind));
w_I_99 = exp(Control(w_I_99_ind));
w_I_999 = exp(Control(w_I_999_ind));
I_I_1_01 = exp(Control(I_I_1_01_ind));
I_I_1_1 = exp(Control(I_I_1_1_ind));
I_I_1_10 = exp(Control(I_I_1_10_ind));
I_I_1_Q1 = exp(Control(I_I_1_Q1_ind));
I_I_1_Q2 = exp(Control(I_I_1_Q2_ind));
I_I_1_Q3 = exp(Control(I_I_1_Q3_ind));
I_I_1_Q4 = exp(Control(I_I_1_Q4_ind));
I_I_1_Q5 = exp(Control(I_I_1_Q5_ind));
I_I_1_90 = exp(Control(I_I_1_90_ind));
I_I_1_99 = exp(Control(I_I_1_99_ind));
I_I_1_999 = exp(Control(I_I_1_999_ind));
I_I_2_01 = exp(Control(I_I_2_01_ind));
I_I_2_1 = exp(Control(I_I_2_1_ind));
I_I_2_10 = exp(Control(I_I_2_10_ind));
I_I_2_Q1 = exp(Control(I_I_2_Q1_ind));
I_I_2_Q2 = exp(Control(I_I_2_Q2_ind));
I_I_2_Q3 = exp(Control(I_I_2_Q3_ind));
I_I_2_Q4 = exp(Control(I_I_2_Q4_ind));
I_I_2_Q5 = exp(Control(I_I_2_Q5_ind));
I_I_2_90 = exp(Control(I_I_2_90_ind));
I_I_2_99 = exp(Control(I_I_2_99_ind));
I_I_2_999 = exp(Control(I_I_2_999_ind));
C_I_01 = exp(Control(C_I_01_ind));
C_I_1 = exp(Control(C_I_1_ind));
C_I_10 = exp(Control(C_I_10_ind));
C_I_Q1 = exp(Control(C_I_Q1_ind));
C_I_Q2 = exp(Control(C_I_Q2_ind));
C_I_Q3 = exp(Control(C_I_Q3_ind));
C_I_Q4 = exp(Control(C_I_Q4_ind));
C_I_Q5 = exp(Control(C_I_Q5_ind));
C_I_90 = exp(Control(C_I_90_ind));
C_I_99 = exp(Control(C_I_99_ind));
C_I_999 = exp(Control(C_I_999_ind));
I_90to10 = exp(Control(I_90to10_ind));
I_50to10 = exp(Control(I_50to10_ind));
I_90to50 = exp(Control(I_90to50_ind));
W_90to50 = exp(Control(W_90to50_ind));
Inc_1_T10toB10 = exp(Control(Inc_1_T10toB10_ind));
Inc_1_T10toQ3 = exp(Control(Inc_1_T10toQ3_ind));
Inc_1_Q3toB10 = exp(Control(Inc_1_Q3toB10_ind));
Inc_2_T10toB10 = exp(Control(Inc_2_T10toB10_ind));
Inc_2_T10toQ3 = exp(Control(Inc_2_T10toQ3_ind));
Inc_2_Q3toB10 = exp(Control(Inc_2_Q3toB10_ind));
Wlth_T10toQ3 = exp(Control(Wlth_T10toQ3_ind));
C_T10toB10 = exp(Control(C_T10toB10_ind));
C_T10toQ3 = exp(Control(C_T10toQ3_ind));
C_Q3toB10 = exp(Control(C_Q3toB10_ind));
Inc_I_1_T10toB10 = exp(Control(Inc_I_1_T10toB10_ind));
Inc_I_1_T10toQ3 = exp(Control(Inc_I_1_T10toQ3_ind));
Inc_I_1_Q3toB10 = exp(Control(Inc_I_1_Q3toB10_ind));
Inc_I_2_T10toB10 = exp(Control(Inc_I_2_T10toB10_ind));
Inc_I_2_T10toQ3 = exp(Control(Inc_I_2_T10toQ3_ind));
Inc_I_2_Q3toB10 = exp(Control(Inc_I_2_Q3toB10_ind));
Wlth_I_T10toQ3 = exp(Control(Wlth_I_T10toQ3_ind));
C_I_T10toB10 = exp(Control(C_I_T10toB10_ind));
C_I_T10toQ3 = exp(Control(C_I_T10toQ3_ind));
C_I_Q3toB10 = exp(Control(C_I_Q3toB10_ind));


A_hh       = exp(Control(A_hh_ind));
B_hh       = exp(Control(B_hh_ind));
C          = exp(Control(C_ind));
K          = exp(Control(K_ind));
B          = exp(Control(B_ind));
N          = exp(Control(N_ind));
L          = exp(Control(L_ind));
T          = exp(Control(T_ind));
LT         = exp(Control(LT_ind));
G          = exp(Control(G_ind));
MRS        = exp(Control(MRS_ind));
LAMBDA     = exp(Control(Lambda_ind));

PI         = exp(Control(pi_ind));
V          = exp(Control(V_ind));
J          = exp(Control(J_ind));
H          = exp(Control(h_ind));
v          = exp(Control(v_ind));
Y          = exp(Control(Y_ind));
PROFIT     = exp(Control(Profit_ind));
R_K        = exp(Control(r_k_ind));
R_A        = exp(Control(r_a_ind));
MC         = exp(Control(MC_ind));
RRa        = exp(Control(RRa_ind));
RR         = exp(Control(RR_ind));
I          = exp(Control(I_ind));
unemp      = exp(Control(unemp_ind));
x_k        = exp(Control(x_k_ind));
nn         = exp(Control(nn_ind));
M          = exp(Control(M_ind));
f          = exp(Control(f_ind));
B_gov_ncp  = exp(Control(B_gov_ncp_ind));
UB         = exp(Control(UB_ind));

Y_obs      = exp(Control(Y_obs_ind));
C_obs      = exp(Control(C_obs_ind));
I_obs      = exp(Control(I_obs_ind));
w_obs      = exp(Control(w_obs_ind));
PROFIT_obs = exp(Control(PROFIT_obs_ind));
unemp_obs  = exp(Control(unemp_obs_ind));
inf_obs    = exp(Control(inf_obs_ind));
R_obs      = exp(Control(R_obs_ind));
pastw      = exp(Control(pastw_ind));
G_obs      = exp(Control(G_obs_ind));

pastYY      = exp(Control(pastYY_ind));
pastCC      = exp(Control(pastCC_ind));
pastII      = exp(Control(pastII_ind));
pastPPROFIT = exp(Control(pastPPROFIT_ind));
pastuu      = exp(Control(pastuu_ind));
pastGG      = exp(Control(pastGG_ind));
pastA_g     = exp(Control(pastA_g_ind));
pastLT2     = exp(Control(pastLT2_ind));
pastG2      = exp(Control(pastG2_ind));
LT_obs      = exp(Control(LT_obs_ind));
B_gov_ncp2  = exp(Control(B_gov_ncp2_ind));

w2          = exp(Control(w2_ind));
Profit2     = exp(Control(Profit2_ind));
Profit3     = exp(Control(Profit3_ind));
Profit_FI2  = exp(Control(Profit_FI2_ind));
r_a2        = exp(Control(r_a2_ind));

l_lambda    = exp(Control(l_lambda_ind));
pastQ       = exp(Control(pastQ_ind));

x_I       = exp(Control(x_I_ind));
eta2      = exp(Control(eta2_ind));
iota2     = exp(Control(iota2_ind));

% Observables

A_g_obs    = exp(Control(A_g_obs_ind));


tau_w      = param.tau_w;
tau_a      = param.tau_a;

zznext         = exp(Controlnext(zz_ind));
xxnext         = exp(Controlnext(xx_ind));
Profit_FInext  = exp(Controlnext(Profit_FI_ind));
vvnext         = exp(Controlnext(vv_ind));
eenext         = exp(Controlnext(ee_ind));
C_bnext        = exp(Controlnext(C_b_ind));

zz             = exp(Control(zz_ind));
xx             = exp(Control(xx_ind));
Profit_FI      = exp(Control(Profit_FI_ind));
vv             = exp(Control(vv_ind));
ee             = exp(Control(ee_ind));
C_b            = exp(Control(C_b_ind));

%% Write LHS values (t variables)

%  Controls

LHS(nx+Profit_aux_ind) = Profit_aux;
LHS(nx+GiniW_ind) = GiniW;
LHS(nx+GiniI_1_ind) = GiniI_1;
LHS(nx+GiniI_2_ind) = GiniI_2;
LHS(nx+GiniC_ind) = GiniC;
LHS(nx+w_01_ind) = w_01;
LHS(nx+w_1_ind) = w_1;
LHS(nx+w_10_ind) = w_10;
LHS(nx+w_Q1_ind) = w_Q1;
LHS(nx+w_Q2_ind) = w_Q2;
LHS(nx+w_Q3_ind) = w_Q3;
LHS(nx+w_Q4_ind) = w_Q4;
LHS(nx+w_Q5_ind) = w_Q5;
LHS(nx+w_90_ind) = w_90;
LHS(nx+w_99_ind) = w_99;
LHS(nx+w_999_ind) = w_999;
LHS(nx+I_1_01_ind) = I_1_01;
LHS(nx+I_1_1_ind) = I_1_1;
LHS(nx+I_1_10_ind) = I_1_10;
LHS(nx+I_1_Q1_ind) = I_1_Q1;
LHS(nx+I_1_Q2_ind) = I_1_Q2;
LHS(nx+I_1_Q3_ind) = I_1_Q3;
LHS(nx+I_1_Q4_ind) = I_1_Q4;
LHS(nx+I_1_Q5_ind) = I_1_Q5;
LHS(nx+I_1_90_ind) = I_1_90;
LHS(nx+I_1_99_ind) = I_1_99;
LHS(nx+I_1_999_ind) = I_1_999;
LHS(nx+I_2_01_ind) = I_2_01;
LHS(nx+I_2_1_ind) = I_2_1;
LHS(nx+I_2_10_ind) = I_2_10;
LHS(nx+I_2_Q1_ind) = I_2_Q1;
LHS(nx+I_2_Q2_ind) = I_2_Q2;
LHS(nx+I_2_Q3_ind) = I_2_Q3;
LHS(nx+I_2_Q4_ind) = I_2_Q4;
LHS(nx+I_2_Q5_ind) = I_2_Q5;
LHS(nx+I_2_90_ind) = I_2_90;
LHS(nx+I_2_99_ind) = I_2_99;
LHS(nx+I_2_999_ind) = I_2_999;
LHS(nx+C_01_ind) = C_01;
LHS(nx+C_1_ind) = C_1;
LHS(nx+C_10_ind) = C_10;
LHS(nx+C_Q1_ind) = C_Q1;
LHS(nx+C_Q2_ind) = C_Q2;
LHS(nx+C_Q3_ind) = C_Q3;
LHS(nx+C_Q4_ind) = C_Q4;
LHS(nx+C_Q5_ind) = C_Q5;
LHS(nx+C_90_ind) = C_90;
LHS(nx+C_99_ind) = C_99;
LHS(nx+C_999_ind) = C_999;
LHS(nx+w_I_01_ind) = w_I_01;
LHS(nx+w_I_1_ind) = w_I_1;
LHS(nx+w_I_10_ind) = w_I_10;
LHS(nx+w_I_Q1_ind) = w_I_Q1;
LHS(nx+w_I_Q2_ind) = w_I_Q2;
LHS(nx+w_I_Q3_ind) = w_I_Q3;
LHS(nx+w_I_Q4_ind) = w_I_Q4;
LHS(nx+w_I_Q5_ind) = w_I_Q5;
LHS(nx+w_I_90_ind) = w_I_90;
LHS(nx+w_I_99_ind) = w_I_99;
LHS(nx+w_I_999_ind) = w_I_999;
LHS(nx+I_I_1_01_ind) = I_I_1_01;
LHS(nx+I_I_1_1_ind) = I_I_1_1;
LHS(nx+I_I_1_10_ind) = I_I_1_10;
LHS(nx+I_I_1_Q1_ind) = I_I_1_Q1;
LHS(nx+I_I_1_Q2_ind) = I_I_1_Q2;
LHS(nx+I_I_1_Q3_ind) = I_I_1_Q3;
LHS(nx+I_I_1_Q4_ind) = I_I_1_Q4;
LHS(nx+I_I_1_Q5_ind) = I_I_1_Q5;
LHS(nx+I_I_1_90_ind) = I_I_1_90;
LHS(nx+I_I_1_99_ind) = I_I_1_99;
LHS(nx+I_I_1_999_ind) = I_I_1_999;
LHS(nx+I_I_2_01_ind) = I_I_2_01;
LHS(nx+I_I_2_1_ind) = I_I_2_1;
LHS(nx+I_I_2_10_ind) = I_I_2_10;
LHS(nx+I_I_2_Q1_ind) = I_I_2_Q1;
LHS(nx+I_I_2_Q2_ind) = I_I_2_Q2;
LHS(nx+I_I_2_Q3_ind) = I_I_2_Q3;
LHS(nx+I_I_2_Q4_ind) = I_I_2_Q4;
LHS(nx+I_I_2_Q5_ind) = I_I_2_Q5;
LHS(nx+I_I_2_90_ind) = I_I_2_90;
LHS(nx+I_I_2_99_ind) = I_I_2_99;
LHS(nx+I_I_2_999_ind) = I_I_2_999;
LHS(nx+C_I_01_ind) = C_I_01;
LHS(nx+C_I_1_ind) = C_I_1;
LHS(nx+C_I_10_ind) = C_I_10;
LHS(nx+C_I_Q1_ind) = C_I_Q1;
LHS(nx+C_I_Q2_ind) = C_I_Q2;
LHS(nx+C_I_Q3_ind) = C_I_Q3;
LHS(nx+C_I_Q4_ind) = C_I_Q4;
LHS(nx+C_I_Q5_ind) = C_I_Q5;
LHS(nx+C_I_90_ind) = C_I_90;
LHS(nx+C_I_99_ind) = C_I_99;
LHS(nx+C_I_999_ind) = C_I_999;
LHS(nx+I_90to10_ind) = I_90to10;
LHS(nx+I_50to10_ind) = I_50to10;
LHS(nx+I_90to50_ind) = I_90to50;
LHS(nx+W_90to50_ind) = W_90to50;
LHS(nx+Inc_1_T10toB10_ind) = Inc_1_T10toB10;
LHS(nx+Inc_1_T10toQ3_ind) = Inc_1_T10toQ3;
LHS(nx+Inc_1_Q3toB10_ind) = Inc_1_Q3toB10;
LHS(nx+Inc_2_T10toB10_ind) = Inc_2_T10toB10;
LHS(nx+Inc_2_T10toQ3_ind) = Inc_2_T10toQ3;
LHS(nx+Inc_2_Q3toB10_ind) = Inc_2_Q3toB10;
LHS(nx+Wlth_T10toQ3_ind) = Wlth_T10toQ3;
LHS(nx+C_T10toB10_ind) = C_T10toB10;
LHS(nx+C_T10toQ3_ind) = C_T10toQ3;
LHS(nx+C_Q3toB10_ind) = C_Q3toB10;
LHS(nx+Inc_I_1_T10toB10_ind) = Inc_I_1_T10toB10;
LHS(nx+Inc_I_1_T10toQ3_ind) = Inc_I_1_T10toQ3;
LHS(nx+Inc_I_1_Q3toB10_ind) = Inc_I_1_Q3toB10;
LHS(nx+Inc_I_2_T10toB10_ind) = Inc_I_2_T10toB10;
LHS(nx+Inc_I_2_T10toQ3_ind) = Inc_I_2_T10toQ3;
LHS(nx+Inc_I_2_Q3toB10_ind) = Inc_I_2_Q3toB10;
LHS(nx+Wlth_I_T10toQ3_ind) = Wlth_I_T10toQ3;
LHS(nx+C_I_T10toB10_ind) = C_I_T10toB10;
LHS(nx+C_I_T10toQ3_ind) = C_I_T10toQ3;
LHS(nx+C_I_Q3toB10_ind) = C_I_Q3toB10;



LHS(nx+VALUE_ind)     = Control(VALUE_ind);
LHS(nx+mutil_c_ind)   = Control(mutil_c_ind);
LHS(nx+Va_ind)        = Control(Va_ind);

LHS(nx+A_hh_ind)      = A_hh;
LHS(nx+B_hh_ind)      = B_hh;
LHS(nx+C_ind)         = C;
LHS(nx+K_ind)         = K;
LHS(nx+B_ind)         = B;
LHS(nx+N_ind)         = N;
LHS(nx+L_ind)         = L;
LHS(nx+T_ind)         = T;
LHS(nx+LT_ind)        = LT;
LHS(nx+G_ind)         = G;
LHS(nx+MRS_ind)       = MRS; 
LHS(nx+Lambda_ind)    = LAMBDA;

LHS(nx+pi_ind)        = PI;
LHS(nx+V_ind)         = V;       
LHS(nx+J_ind)         = J;      
LHS(nx+h_ind)         = H;      
LHS(nx+v_ind)         = v;      
LHS(nx+Y_ind)         = Y;      
LHS(nx+Profit_ind)    = PROFIT; 
LHS(nx+r_k_ind)       = R_K;    
LHS(nx+r_a_ind)       = R_A;    
LHS(nx+MC_ind)        = MC;
LHS(nx+RRa_ind)       = RRa;
LHS(nx+RR_ind)        = RR;
LHS(nx+I_ind)         = I;
LHS(nx+unemp_ind)     = unemp;
LHS(nx+x_k_ind)       = x_k;

LHS(nx+nn_ind)        = nn;
LHS(nx+M_ind)         = M;
LHS(nx+f_ind)         = f;
LHS(nx+B_gov_ncp_ind) = B_gov_ncp;
LHS(nx+UB_ind)        = UB;

LHS(nx+zz_ind)        = zz;
LHS(nx+xx_ind)        = xx;
LHS(nx+Profit_FI_ind) = Profit_FI;
LHS(nx+vv_ind)        = vv;
LHS(nx+ee_ind)        = ee;
LHS(nx+C_b_ind)       = C_b;
LHS(nx+A_g_obs_ind)   = A_g_obs;

LHS(nx+Y_obs_ind)     = Y_obs;
LHS(nx+C_obs_ind)     = C_obs;
LHS(nx+I_obs_ind)     = I_obs;
LHS(nx+w_obs_ind)     = w_obs;
LHS(nx+PROFIT_obs_ind)= PROFIT_obs;
LHS(nx+unemp_obs_ind) = unemp_obs;
LHS(nx+inf_obs_ind)   = inf_obs;
LHS(nx+R_obs_ind)     = R_obs;
LHS(nx+pastw_ind)     = pastw;
LHS(nx+G_obs_ind)     = G_obs;

LHS(nx+pastYY_ind)      = pastYY;
LHS(nx+pastCC_ind)      = pastCC;
LHS(nx+pastII_ind)      = pastII;
LHS(nx+pastPPROFIT_ind) = pastPPROFIT;
LHS(nx+pastuu_ind)      = pastuu;
LHS(nx+pastGG_ind)      = pastGG;
LHS(nx+pastA_g_ind)     = pastA_g;
LHS(nx+pastLT2_ind)     = pastLT2;
LHS(nx+pastG2_ind)      = pastG2;
LHS(nx+LT_obs_ind)      = LT_obs;
LHS(nx+B_gov_ncp2_ind)  = B_gov_ncp2;

LHS(nx+w2_ind)          = w2; 
LHS(nx+Profit2_ind)     = Profit2; 
LHS(nx+Profit3_ind)     = Profit3; 
LHS(nx+Profit_FI2_ind)  = Profit_FI2; 
LHS(nx+r_a2_ind)        = r_a2; 

LHS(nx+l_lambda_ind)    = l_lambda;
LHS(nx+pastQ_ind)       = pastQ;
LHS(nx+x_I_ind)         = x_I;
LHS(nx+eta2_ind)        = eta2;
LHS(nx+iota2_ind)       = iota2;

%  States
%  Marginal Distributions

LHS(marginal_b_ind)  = Distribution(1:grid.nb-1);
ba = grid.nb;
LHS(marginal_a_ind)  = Distribution(ba+(1:grid.na-1));
ba = grid.nb + grid.na;
LHS(marginal_se_ind) = Distribution(ba+(1:grid.nse-1));

LHS(R_cb_ind)        = (R_cb);
LHS(w_ind)           = (W);
LHS(A_b_ind)         = (A_bnext);
LHS(B_b_ind)         = (B_bnext);
LHS(A_g_ind)         = (A_gauxnext);
LHS(Q_ind)           = (Q);
LHS(lev_ind)         = (lev);
LHS(NW_b_ind)        = (NW_b);
LHS(R_tilde_ind)     = (R_tilde);
LHS(x_cb_ind)        = (x_cb);

LHS(pastpi_ind)      = (pastpi);

LHS(pastY_ind)       = pastY;
LHS(pastC_ind)       = pastC;
LHS(pastI_ind)       = pastI;
LHS(pastPROFIT_ind)  = pastPROFIT;
LHS(pastunemp_ind)   = pastunemp;
LHS(pastG_ind)       = pastG;
LHS(pastLT_ind)      = pastLT;
LHS(R_star_ind)      = R_star;
LHS(B_F_ind)         = B_Fnext;

LHS(Z_ind)           = (Z);
LHS(PSI_RP_ind)      = (PSI_RP);
LHS(eta_ind)         = (eta);
LHS(D_ind)           = (D);
LHS(GG_ind)          = (GG);
LHS(iota_ind)        = (iota);
LHS(BB_ind)          = (BB);
LHS(PSI_W_ind)       = (PSI_W);
LHS(MP_ind)          = (MP);
LHS(p_mark_ind)      = (p_mark);

LHS(pgap_ind)        = (pgap);
LHS(ygap_ind)        = (ygap);


LHS(eps_QE_ind)      = (eps_QEnext);
LHS(eps_RP_ind)      = (eps_RPnext);
LHS(eps_B_F_ind)     = (eps_B_Fnext);
LHS(eps_BB_ind)      = (eps_BBnext);
LHS(eps_Z_ind)       = (eps_Znext);
LHS(eps_G_ind)       = (eps_Gnext);
LHS(eps_D_ind)       = (eps_Dnext);
LHS(eps_R_ind)       = (eps_Rnext);
LHS(eps_iota_ind)    = (eps_iotanext);
LHS(eps_eta_ind)     = (eps_etanext);
LHS(eps_w_ind)       = (eps_wnext);


R_cb        = exp(R_cb)        ; R_cbminus        = exp(R_cbminus);
W           = exp(W)           ; Wminus           = exp(Wminus);
A_bnext     = exp(A_bnext)     ; A_b              = exp(A_b);
B_bnext     = exp(B_bnext)     ; B_b              = exp(B_b);
A_gaux      = exp(A_gaux)      ; A_gauxnext       = exp(A_gauxnext);
Q           = exp(Q)           ; Qminus           = exp(Qminus);
lev         = exp(lev)         ; levminus         = exp(levminus);
NW_b        = exp(NW_b)        ; NW_bminus        = exp(NW_bminus);
R_tilde     = exp(R_tilde)     ; R_tildeminus     = exp(R_tildeminus);
x_cb        = exp(x_cb)        ; x_cbminus        = exp(x_cbminus);
pastpi      = exp(pastpi)      ; pastpiminus      = exp(pastpiminus);
pastY       = exp(pastY)       ; pastYminus       = exp(pastYminus);
pastC       = exp(pastC)       ; pastCminus       = exp(pastCminus);
pastI       = exp(pastI)       ; pastIminus       = exp(pastIminus);
pastPROFIT  = exp(pastPROFIT)  ; pastPROFITminus  = exp(pastPROFITminus);
pastunemp   = exp(pastunemp)   ; pastunempminus   = exp(pastunempminus);
pastG       = exp(pastG)       ; pastGminus       = exp(pastGminus);
pastLT      = exp(pastLT)      ; pastLTminus      = exp(pastLTminus);
B_Fnext     = exp(B_Fnext)     ; B_F              = exp(B_F);
R_star      = exp(R_star)      ; R_starminus      = exp(R_starminus);
Z           = exp(Z)           ; Zminus           = exp(Zminus);
PSI_RP      = exp(PSI_RP)      ; PSI_RPminus      = exp(PSI_RPminus);
eta         = exp(eta)         ; etaminus         = exp(etaminus);
D           = exp(D)           ; Dminus           = exp(Dminus);
GG          = exp(GG)          ; GGminus          = exp(GGminus);
iota        = exp(iota)        ; iotaminus        = exp(iotaminus);
BB          = exp(BB)          ; BBminus          = exp(BBminus);
PSI_W       = exp(PSI_W)       ; PSI_Wminus       = exp(PSI_Wminus);
MP          = exp(MP)          ; MPminus          = exp(MPminus);
p_mark      = exp(p_mark)      ; p_markminus      = exp(p_markminus);
pgap        = exp(pgap)        ; pgapminus        = exp(pgapminus);
ygap        = exp(ygap)        ; ygapminus        = exp(ygapminus);





marginal_b  = Distribution(1:grid.nb)';
marginal_a  = Distribution((1:grid.na)+grid.nb)';
marginal_se = Distribution((1:grid.nse)+grid.nb+grid.na)';

marginal_bminus  = Distributionminus(1:grid.nb)';
marginal_aminus  = Distributionminus((1:grid.na)+grid.nb)';
marginal_seminus = Distributionminus((1:grid.nse)+grid.nb+grid.na)';



A_gnext = A_gauxnext;
A_g     = A_gaux;

cumdist = zeros(grid.nb+1,grid.na+1,grid.nse+1);
cumdist(2:end,2:end,2:end) = Copula({cumsum(marginal_bminus),cumsum(marginal_aminus)',cumsum(marginal_seminus)});
MU = diff(diff(diff(cumdist,1,1),1,2),1,3);

MU_b  = squeeze(sum(sum(MU,2),3));
MU_a  = squeeze(sum(sum(MU,1),3))';
MU_se = squeeze(sum(sum(MU,1),2));


Eshare = param.Eshare;

[meshes.b,meshes.a,meshes.se] = ndgrid(grid.b,grid.a,grid.se);

P_SS3_aux   = kron([1-l_lambda+l_lambda*f,l_lambda*(1-f);f,1-f],eye(grid.ns));  
P_SS3_aux   = [P_SS3_aux, zeros(grid.nse-1,1)];
lastrow     = [zeros(1,2*grid.ns),1];
P_SS3_aux   = [P_SS3_aux; lastrow];

H_tilde  = kron(P_SS3_aux',speye(grid.nb*grid.na));

MU_tilde = H_tilde*MU(:);

MU_tilde = reshape(MU_tilde,[grid.nb grid.na grid.nse]);

MU_tilde_b  = squeeze(sum(sum(MU_tilde,2),3));
MU_tilde_a  = squeeze(sum(sum(MU_tilde,1),3))';
MU_tilde_se = squeeze(sum(sum(MU_tilde,1),2));

marginal_setilde_aux = MU_tilde_se(1:2*grid.ns);

P_SS3next_aux  = kron([1-l_lambdanext+l_lambdanext*fnext,l_lambdanext*(1-fnext);fnext,1-fnext],eye(grid.ns));  
P_SS3next_aux  = [P_SS3next_aux repmat(0, [grid.nse-1 1])];
lastrow        = [zeros(1,2*grid.ns),1];
P_SS3next_aux  = [P_SS3next_aux; lastrow];

P_transition_exp  = param.P_SS2*P_SS3next_aux;
P_transition_dist = param.P_SS2;  % The new transition matrix

%  Income (grids)

auxWW                    = ones(grid.nb,grid.na,grid.nse);
auxWW(:,:,grid.ns+1:end) = zeros(grid.nb,grid.na,grid.ns+1);
auxWW1                   = auxWW;
auxWW                    = auxWW.*(meshes.se.^(param.b_aux-1)); % bouns

auxWW2                    = ones(grid.nb,grid.na,grid.nse);
auxWW2(:,:,1:grid.ns)     = zeros(grid.nb,grid.na,grid.ns);
auxWW2(:,:,end)           = zeros(grid.nb,grid.na,1);

% param.beta_aux = param.beta*gamma_Q;

% param.beta_aux = param.beta;

param.beta_aux = param.beta*D;

% param.beta_aux = param.beta;

R_A_aux = R_A + (param.death_rate/(1-param.death_rate))*(Q);  % taking into aPSI_RPount the annuity payments

% After-tax wage xor unemployment benefit

NW          = param.xi/(1+param.xi)*nn*W;
WW          = (1-tau_w)*(NW*auxWW1 + param.xi/(1+param.xi)*1*param.w_bar.*auxWW2);
WW(:,:,end) = (1-tau_a)*(param.Eratio*PROFIT + Profit_FI - param.fix2)/Eshare*ones(grid.nb,grid.na);
inc.labor    = WW.*(meshes.se);
inc.dividend = meshes.a*R_A_aux;
inc.capital  = meshes.a*Q;
inc.bond     = (R_cbminus/PI).* meshes.b ...
				+ (meshes.b<0).*(param.Rprem/PI).* meshes.b;			
inc.bond     = inc.bond/(1-param.death_rate)*param.b_a_aux;					
inc.transfer = (LT+SS_stats.C_b)/(1+param.frac_b)*ones(grid.nb,grid.na,grid.nse);

%% First Set: Value functions, marginal values
%% Update policies

EVa         = 1*reshape(reshape(Vanext,[grid.nb*grid.na grid.nse])*P_transition_exp',[grid.nb grid.na grid.nse]);
R_tildeaux	= R_tilde/PInext + (meshes.b<0).*(param.Rprem/PInext); 

EVb         = 1*reshape(reshape(R_tildeaux(:)/(1-param.death_rate)*param.b_a_aux.*mutil_cnext,[grid.nb*grid.na grid.nse])*P_transition_exp',[grid.nb grid.na grid.nse]);

% [c_a_star,b_a_star,a_a_star,c_n_star,b_n_star] = policies_update(EVb,EVa,Q,PI,R_tildeminus,1,1,inc,meshes,grid,param);		
[c_a_star,b_a_star,a_a_star,c_n_star,b_n_star] = policies_update(EVb,EVa,Q,PI,R_cbminus,1,1,inc,meshes,grid,param);			

%% Update value functions

meshaux = meshes;


[~,~,meshaux.se_aux] = ndgrid(grid.b,grid.a,1:grid.nse);

EV         = reshape(VALUEnext,[grid.nb*grid.na grid.nse])*P_transition_exp';

VALUE_next = griddedInterpolant(meshaux.b,meshaux.a,meshaux.se_aux,reshape(EV,[grid.nb grid.na grid.nse]));

V_adjust   = util(c_a_star) + param.beta_aux*(1-param.death_rate)*VALUE_next(b_a_star,a_a_star ,meshaux.se_aux);
V_noadjust = util(c_n_star) + param.beta_aux*(1-param.death_rate)*VALUE_next(b_n_star,meshaux.a,meshaux.se_aux);

AProb = max(min(1./(1+exp(-((V_adjust-V_noadjust)-param.mu_chi)./param.sigma_chi)),1-1e-6),1e-6);
AC    = param.sigma_chi.*((1-AProb(:)).*log(1-AProb(:)) + AProb(:).*log(AProb(:)))...
    + param.mu_chi.*AProb(:) -  (param.mu_chi-param.sigma_chi*log(1+exp(param.mu_chi/param.sigma_chi)));
AProb = AProb.*param.nu;
AC    = AC.*param.nu;

VALUEaux  = AProb(:).*V_adjust(:) + (1-AProb(:)).*V_noadjust(:) - AC(:);

%% Update Marginal Value of Bonds
mutil_c_n   = mutil(c_n_star);                         % marginal utility at consumption policy no adjustment
mutil_c_a   = mutil(c_a_star);                         % marginal utility at consumption policy adjustment
mutil_c_aux = AProb.*mutil_c_a + (1-AProb).*mutil_c_n; % Expected marginal utility at consumption policy (w &w/o adjustment)


mutil_cnext_aux = reshape(reshape(mutil_cnext,[grid.nb*grid.na grid.nse])*P_transition_exp',[grid.nb grid.na grid.nse]);
mutil_cnext_aux = griddedInterpolant(meshaux.b,meshaux.a,meshaux.se_aux,reshape(mutil_cnext_aux,[grid.nb grid.na grid.nse]));

MRS_a_aux = param.beta_aux*(1-param.death_rate)*mutil_cnext_aux(b_a_star,a_a_star ,meshaux.se_aux)./mutil_c_a.*  AProb  .*MU_tilde;
MRS_n_aux = param.beta_aux*(1-param.death_rate)*mutil_cnext_aux(b_n_star,meshaux.a,meshaux.se_aux)./mutil_c_n.*(1-AProb).*MU_tilde;

Va_next  = griddedInterpolant(meshaux.b,meshaux.a,meshaux.se_aux,reshape(EVa,[grid.nb grid.na grid.nse]));
Va_aux   = AProb.*((R_A_aux+Q)).*mutil_c_a + (1-AProb).*((R_A_aux)).*mutil_c_n + param.beta_aux*1*(1-param.death_rate).*(1-AProb).*Va_next(b_n_star,meshaux.a,meshaux.se_aux);


[~       ,~       ,meshes.se3_aux] = ndgrid(grid.b,grid.a,grid.se3_aux);


MU_Ent = (meshes.se3_aux==1).*MU;

MU_tilde_Ent   = H_tilde*MU_Ent(:);

MU_tilde_Ent   = reshape(MU_tilde_Ent,   [grid.nb grid.na grid.nse]);

% RHS(nx+C_ind) = q_cons2(c_a_star,c_n_star,Wminus,nnminus,AProb,MU_tilde,meshes,param,grid,tau_wminus);

% Differences for distributions

weight11   = zeros(grid.nb*grid.na, grid.nse, grid.nse);
weight12   = zeros(grid.nb*grid.na, grid.nse, grid.nse);
weight21   = zeros(grid.nb*grid.na, grid.nse, grid.nse);
weight22   = zeros(grid.nb*grid.na, grid.nse, grid.nse);

% find next smallest on-grid value for money and capital choices
[Dist_b_a,idb_a] = genweight(b_a_star,grid.b);
[Dist_b_n,idb_n] = genweight(b_n_star,grid.b);
[Dist_a_a,ida_a] = genweight(a_a_star,grid.a);
[Dist_a_n,ida_n] = genweight(meshaux.a,grid.a);
% [Dist_a_n,ida_n] = genweight(TAU*meshaux.a,grid.a);

[Dist_b_die,idb_die] = genweight(zeros(grid.nb,grid.na,grid.nse),grid.b);
[Dist_a_die,ida_die] = genweight(zeros(grid.nb,grid.na,grid.nse),grid.a);

% Transition matrix for adjustment case
idb_a = repmat(idb_a(:),[1 grid.nse]);
ida_a = repmat(ida_a(:),[1 grid.nse]);
idse  = kron(1:grid.nse,ones(1,grid.nb*grid.na*grid.nse));

index11 = sub2ind([grid.nb grid.na grid.nse],idb_a(:),ida_a(:),idse(:));
index12 = sub2ind([grid.nb grid.na grid.nse],idb_a(:),ida_a(:)+1,idse(:));
index21 = sub2ind([grid.nb grid.na grid.nse],idb_a(:)+1,ida_a(:),idse(:));
index22 = sub2ind([grid.nb grid.na grid.nse],idb_a(:)+1,ida_a(:)+1,idse(:));

for hh = 1:grid.nse

	% Coresponding weights
	weight11_aux = (1-Dist_b_a(:,:,hh)).*(1-Dist_a_a(:,:,hh));
	weight12_aux = (1-Dist_b_a(:,:,hh)).*(Dist_a_a(:,:,hh));
	weight21_aux = (Dist_b_a(:,:,hh)).*(1-Dist_a_a(:,:,hh));
	weight22_aux = (Dist_b_a(:,:,hh)).*(Dist_a_a(:,:,hh));

	% Dimensions (b x a, h', h)
	weight11(:,:,hh) = weight11_aux(:)*P_transition_dist(hh,:);
	weight12(:,:,hh) = weight12_aux(:)*P_transition_dist(hh,:);
	weight21(:,:,hh) = weight21_aux(:)*P_transition_dist(hh,:);
	weight22(:,:,hh) = weight22_aux(:)*P_transition_dist(hh,:);

end

% Dimensions (b x a, h, h')
weight11 = permute(weight11, [1 3 2]);
weight12 = permute(weight12, [1 3 2]);
weight21 = permute(weight21, [1 3 2]);
weight22 = permute(weight22, [1 3 2]);
rowindex = repmat(1:grid.nb*grid.na*grid.nse,[1 4*grid.nse]);

TT_a = sparse(rowindex(:),[index11(:); index21(:); index12(:); index22(:)],...
	  [weight11(:); weight21(:); weight12(:); weight22(:)],grid.nb*grid.na*grid.nse,grid.nb*grid.na*grid.nse);

weight11   = zeros(grid.nb*grid.na, grid.nse, grid.nse);
weight12   = zeros(grid.nb*grid.na, grid.nse, grid.nse);
weight21   = zeros(grid.nb*grid.na, grid.nse, grid.nse);
weight22   = zeros(grid.nb*grid.na, grid.nse, grid.nse);

% Transition matrix for non-adjustment case
idb_n = repmat(idb_n(:),[1 grid.nse]);
ida_n = repmat(ida_n(:),[1 grid.nse]);
idse  = kron(1:grid.nse,ones(1,grid.nb*grid.na*grid.nse));

index11 = sub2ind([grid.nb grid.na grid.nse],idb_n(:),ida_n(:),idse(:));
index12 = sub2ind([grid.nb grid.na grid.nse],idb_n(:),ida_n(:)+1,idse(:));
index21 = sub2ind([grid.nb grid.na grid.nse],idb_n(:)+1,ida_n(:),idse(:));
index22 = sub2ind([grid.nb grid.na grid.nse],idb_n(:)+1,ida_n(:)+1,idse(:));

for hh = 1:grid.nse

	% Coresponding weights
	weight11_aux = (1-Dist_b_n(:,:,hh)).*(1-Dist_a_n(:,:,hh));
	weight12_aux = (1-Dist_b_n(:,:,hh)).*(Dist_a_n(:,:,hh));
	weight21_aux = (Dist_b_n(:,:,hh)).*(1-Dist_a_n(:,:,hh));
	weight22_aux = (Dist_b_n(:,:,hh)).*(Dist_a_n(:,:,hh));

	% Dimensions (b x a, h', h)
	weight11(:,:,hh) = weight11_aux(:)*P_transition_dist(hh,:);
	weight12(:,:,hh) = weight12_aux(:)*P_transition_dist(hh,:);
	weight21(:,:,hh) = weight21_aux(:)*P_transition_dist(hh,:);
	weight22(:,:,hh) = weight22_aux(:)*P_transition_dist(hh,:);

end

% Dimensions (b x a, h, h')
weight11 = permute(weight11, [1 3 2]);
weight12 = permute(weight12, [1 3 2]);
weight21 = permute(weight21, [1 3 2]);
weight22 = permute(weight22, [1 3 2]);
rowindex = repmat(1:grid.nb*grid.na*grid.nse,[1 4*grid.nse]);

TT_n = sparse(rowindex(:),[index11(:); index21(:); index12(:); index22(:)],...
	  [weight11(:); weight21(:); weight12(:); weight22(:)],grid.nb*grid.na*grid.nse,grid.nb*grid.na*grid.nse);

weight11   = zeros(grid.nb*grid.na, grid.nse, grid.nse);
weight12   = zeros(grid.nb*grid.na, grid.nse, grid.nse);
weight21   = zeros(grid.nb*grid.na, grid.nse, grid.nse);
weight22   = zeros(grid.nb*grid.na, grid.nse, grid.nse);

%% Transition matrix die
idb_die =repmat(idb_die(:),[1 grid.nse]);
ida_die =repmat(ida_die(:),[1 grid.nse]);
ids_die =kron(1:grid.nse,ones(1,grid.nb*grid.na*grid.nse));

index11 = sub2ind([grid.nb grid.na grid.nse],idb_die(:),ida_die(:),ids_die(:));
index12 = sub2ind([grid.nb grid.na grid.nse],idb_die(:),ida_die(:)+1,ids_die(:));
index21 = sub2ind([grid.nb grid.na grid.nse],idb_die(:)+1,ida_die(:),ids_die(:));
index22 = sub2ind([grid.nb grid.na grid.nse],idb_die(:)+1,ida_die(:)+1,ids_die(:));

for hh=1:grid.nse

    %Corresponding weights
    weight11_aux = (1-Dist_b_die(:,:,hh)).*(1-Dist_a_die(:,:,hh));
    weight12_aux = (1-Dist_b_die(:,:,hh)).*(Dist_a_die(:,:,hh));
    weight21_aux = Dist_b_die(:,:,hh).*(1-Dist_a_die(:,:,hh));
    weight22_aux = (Dist_b_die(:,:,hh)).*(Dist_a_die(:,:,hh));
    
    % Dimensions (mxk,h',h)
    weight11(:,:,hh) = weight11_aux(:)*P_transition_dist(hh,:);
    weight12(:,:,hh) = weight12_aux(:)*P_transition_dist(hh,:);
    weight21(:,:,hh) = weight21_aux(:)*P_transition_dist(hh,:);
    weight22(:,:,hh) = weight22_aux(:)*P_transition_dist(hh,:);
end
%Dimensions (mxk,h,h')
weight11 = permute(weight11,[1 3 2]);
weight22 = permute(weight22,[1 3 2]);
weight12 = permute(weight12,[1 3 2]);
weight21 = permute(weight21,[1 3 2]);
rowindex = repmat(1:grid.nb*grid.na*grid.nse,[1 4*grid.nse]);

TT_die=sparse(rowindex(:),[index11(:); index21(:); index12(:); index22(:)],...
    [weight11(:); weight21(:); weight12(:); weight22(:)],grid.nb*grid.na*grid.nse,grid.nb*grid.na*grid.nse); % mu'(h',k'), a without interest

% Joint transition matrix and transition

APD_a=sparse(1:length(AProb(:)),1:length(AProb(:)),AProb(:));
APD_n=sparse(1:length(AProb(:)),1:length(AProb(:)),1-AProb(:));

TT = (1-param.death_rate)*(APD_a*TT_a + APD_n*TT_n) + param.death_rate*TT_die;

TT = (H_tilde')*TT;

MUnext = MU(:)'*TT;
MUnext = reshape(MUnext(:),[grid.nb, grid.na, grid.nse]);


aux_b                = squeeze(sum(sum(MUnext,2),3));
aux_a				 = squeeze(sum(sum(MUnext,1),3))';
aux_se               = squeeze(sum(sum(MUnext,1),2));


% State variables

RHS(marginal_b_ind)  = aux_b(1:end-1);
RHS(marginal_a_ind)  = aux_a(1:end-1);
RHS(marginal_se_ind) = aux_se(1:end-1);

RHS(R_cb_ind)        = 0;
RHS(R_star_ind)        = log(param.R_cb) + (1-param.rho_R)*(param.phi_pi*log(PI/param.pi_cb)-param.phi_u*(unemp-SS_stats.u)) + param.rho_R*(log(R_starminus/param.R_cb)) + eps_R;

RHS(w_ind)           = log(param.w_bar) + param.rho_w * log(Wminus/param.w_bar) + param.rho_w*(param.d*log(param.pi_bar/PI)+(1-param.d)*log(pastpiminus/PI)) + (1-param.rho_w)*param.eps_w*log(PSI_W*H/SS_stats.r_l);
RHS(A_b_ind)         = log( lev*NW_b/Q );
RHS(B_b_ind)         = log( T + B_gov_ncpnext - B_gov_ncp*R_cbminus/PI + (Q+R_A-R_cbminus/PI*Qminus)*A_g - UB - LT - param.tau_cp*Q*A_gnext - G + param.b_a_aux2*R_cbminus/PI*B_b - C_b );
        	
RHS(A_g_ind)         = log( ((param.rho_passive_QE)*A_g + (1-param.rho_passive_QE)*SS_stats.A_g) ) ;
RHS(Q_ind)           = log( iota2*(1+param.phi*log(x_k)+param.phi/2*(log(x_k))^2 ) - iota2next*LAMBDA/SS_stats.Lambda*param.phi*log(x_knext)*x_knext );
RHS(lev_ind)         = log( ee/(param.DELTA-vv) );
RHS(NW_b_ind)        = log( (param.theta_b*((RRa-param.b_a_aux2*R_cbminus/PI)*levminus+param.b_a_aux2*R_cbminus/PI)*NW_bminus + param.omega*Qminus*A_b) );
RHS(R_tilde_ind)     = log( R_cb*D );
RHS(x_cb_ind)        = log( 1 +param.rho_X_QE*(x_cbminus-1) + (1-param.rho_X_QE)*(-param.phi_pi_QE*log(PI/param.pi_cb)+param.phi_u_QE*(unemp-SS_stats.u)) ) ;
RHS(pastpi_ind)      = log(PI);

RHS(pastY_ind)       = log(Y);
RHS(pastC_ind)       = log(C);
RHS(pastI_ind)       = log(I);
RHS(pastPROFIT_ind)  = log(PROFIT+param.fix2);
RHS(pastunemp_ind)   = log(unemp);

% switch(param.adjust)
% 	case('G')
% 		RHS(pastG_ind) = log(GG); 
% 	case('LT')
% 		RHS(pastG_ind) = log(G);
% end


RHS(pastG_ind)  = log(w2);
RHS(pastLT_ind) = log(LT);


RHS(B_F_ind)         = param.rho_B_F*log(B_F) + (1-param.rho_B_F)*log(SS_stats.Y/(SS_stats.Y-SS_stats.B_F)) + eps_B_F;
RHS(Z_ind)           = param.rho_Z        * log(Zminus)        + eps_Z;
RHS(PSI_RP_ind)      = param.rho_PSI_RP   * log(PSI_RPminus)   + eps_RP;
RHS(eta_ind)         = log(p_mark/(p_mark-1));
RHS(D_ind )          = param.rho_D        * log(Dminus)        + eps_D;

switch(param.adjust)
	case('G')
		RHS(GG_ind)          = param.rho_G  * log(GGminus) + (1-param.rho_G) * log(SS_stats.Y/(SS_stats.Y-SS_stats.LT-SS_stats.C_b)) + eps_G;
	case('LT')
        RHS(GG_ind)          = param.rho_G  * log(GGminus) + (1-param.rho_G)  * log(SS_stats.Y/(SS_stats.Y-SS_stats.G)) + eps_G;
end



RHS(iota_ind)        = param.rho_iota  * log(iotaminus)  + eps_iota;
RHS(BB_ind)          = param.rho_BB    * log(BBminus)    + eps_BB;
RHS(PSI_W_ind)       = param.rho_PSI_W * log(PSI_Wminus) + eps_w;
RHS(MP_ind)          = param.rho_MP    * log(MPminus)    + eps_R;
RHS(p_mark_ind)      = param.rho_eta   * log(p_markminus) + (1-param.rho_eta)*log(param.eta/(param.eta-1)) + eps_eta;

RHS(pgap_ind)        = log( (PI./param.pi_bar )     *pgapminus^param.rho_pgap  );
RHS(ygap_ind)        = log( exp((SS_stats.u-unemp ))*ygapminus^param.rho_ygap  );


RHS(eps_QE_ind)      = 0;
RHS(eps_B_F_ind)     = 0;
RHS(eps_RP_ind)      = 0;
RHS(eps_eta_ind)     = 0;
RHS(eps_D_ind)       = 0;
RHS(eps_R_ind)       = 0;
RHS(eps_Z_ind)       = 0;
RHS(eps_G_ind)       = 0;
RHS(eps_iota_ind)    = 0;
RHS(eps_w_ind)       = 0;
RHS(eps_BB_ind)      = 0;

% Control variables

RHS(nx+VALUE_ind)   = invutil(VALUEaux(:));
RHS(nx+mutil_c_ind) = invmutil(mutil_c_aux(:));
RHS(nx+Va_ind)      = invmutil(Va_aux(:));

RHS(nx+MRS_ind)     = sum(sum(MRS_a_aux(:,:,end)+MRS_n_aux(:,:,end)))/param.Eshare*param.Lambda_b_aux;

RHS(nx+A_hh_ind)   = sum(grid.a.*marginal_aminus);  
RHS(nx+B_hh_ind)   = sum(grid.b.*marginal_bminus)/(1-param.death_rate)*param.b_a_aux;
RHS(nx+C_ind)      = q_cons2(c_a_star,c_n_star,W,nn,AProb,MU_tilde,meshes,param,grid,tau_w);

RHS(nx+N_ind)      = sum(marginal_seminus(1:grid.ns)); %% beginning of period aggregate employment
RHS(nx+L_ind)      = nn*grid.s*MU_tilde_se(1:grid.ns);
RHS(nx+UB_ind)        = param.w_bar*param.n*param.xi/(1+param.xi)*grid.se(grid.ns+1:2*grid.ns)*MU_tilde_se(grid.ns+1:2*grid.ns);

RHS(nx+K_ind)         = A_hh + A_g  + A_b + SS_stats.A_F;          % beginning of period asset holding
RHS(nx+B_ind)         = B_hh + B_b  + (1-1/1)*SS_stats.Y;   % beginning of period asset holding
RHS(nx+B_gov_ncp_ind) = B_b  + B_hh + (1-1/1)*SS_stats.Y - (Qminus*A_b- NW_bminus) - Qminus *(1+param.tau_cp)* A_g;

RHS(nx+T_ind)      = tau_w.*(W.*L + UB) + tau_a.*(PROFIT);																			

switch(param.adjust)
	case('G')
		RHS(nx+LT_ind) = (1-1/GG)*SS_stats.Y -  SS_stats.C_b;
  	
        RHS(nx+G_ind)  = SS_stats.G + C_b - SS_stats.C_b;										     		           
    case('LT')
		RHS(nx+G_ind)  = (1-1/GG)*SS_stats.Y;
		RHS(nx+LT_ind) = SS_stats.LT;
 		           
end


RHS(nx+Lambda_ind)    = MRS*param.Lambda_aux/param.Lambda_b_aux;																					
RHS(nx+pi_ind)        = param.pi_cb*exp((param.rho_B/(param.rho_B+param.gamma_pi))*log(B_gov_ncp*R_cbminus/(SS_stats.B_gov_ncp*param.R_cb))-(param.gamma_T/(param.rho_B+param.gamma_pi))*log(T/SS_stats.T)-1/(param.rho_B+param.gamma_pi)*log(B_gov_ncpnext/SS_stats.B_gov_ncp));


RHS(nx+V_ind)      = M/param.iota*(J'*grid.s_dist);
RHS(nx+J_ind)      = (H-param.fix_L-W).*(grid.s')*nn + LAMBDA*(1-param.death_rate)*(1-l_lambdanext)*(1-param.in)*param.P_SS*Jnext;

RHS(nx+h_ind)      = Z*MC*(v*K).^param.theta.*(L).^(param.theta_2-1)*(param.theta_2+param.alpha_lk*log(v*K/SS_stats.v/grid.K)+param.alpha_ll*log(L/SS_stats.L)); 
RHS(nx+v_ind)      = min( (R_K/(param.delta_0*param.delta_1))^(1/(param.delta_1-1)),1) ;
RHS(nx+Y_ind)      = Z*(v*K).^(param.theta).*(L).^param.theta_2*((L/SS_stats.L)^(param.alpha_lk*log(v*K/SS_stats.v/grid.K))*(v*K/SS_stats.v/grid.K)^(param.alpha_lk*log(L/SS_stats.L)));
RHS(nx+Profit_ind) = (Y *(1-eta/(2*param.kappa)*(log(PI)-(1-param.GAMMA)*log(param.pi_cb)-param.GAMMA*log(pastpiminus)).^2)+ Y.*(-MC) - param.fix  ...
		                + (H-param.fix_L-W).*L - param.iota.*V)  ...
		                + (R_K*v - param.delta_0.*v.^param.delta_1).*K +  Q*(Knext-K)-(Knext-K) - param.phi/2*(Knext/K-1)^2*K  + Profit_FI - param.fix2 ; 		 		 
	
RHS(nx+r_k_ind)    = Z*MC*(v*K).^(param.theta-1).*(L).^param.theta_2*(param.theta+param.alpha_lk*log(L/SS_stats.L)+param.alpha_kk*log(v*K/SS_stats.v/grid.K));		 
RHS(nx+r_a_ind)    = (1-tau_a)*(1-param.Eratio-param.b_share)*PROFIT/K;

RHS(nx+MC_ind) 	   = 1-1/eta + (log(PI/(pastpiminus^param.GAMMA*param.pi_bar^(1-param.GAMMA)))-LAMBDA*eta2next/eta2*Ynext/Y*log(PInext/(PI^param.GAMMA*param.pi_bar^(1-param.GAMMA)) ))/param.kappa;

RHS(nx+zz_ind)        = (RRa-param.b_a_aux2*R_cbminus/PI)*levminus + param.b_a_aux2*R_cbminus/PI;
RHS(nx+xx_ind)        = (lev/levminus)*zz;
RHS(nx+Profit_FI_ind) = ((1-param.theta_b)*((RRa - param.b_a_aux2*R_cbminus/PI)*levminus+param.b_a_aux2*R_cbminus/PI)*NW_bminus - param.omega*Qminus*A_b);
RHS(nx+vv_ind)        = ((1-param.theta_b)*BB*LAMBDA*(RRanext-R_cb/PInext)+param.theta_b*BB*LAMBDA*xxnext*vvnext);
RHS(nx+ee_ind)        = ((1-param.theta_b)*BB*LAMBDA*R_cb/PInext+param.theta_b*BB*LAMBDA*zznext*eenext); 

RHS(nx+C_b_ind)       = ((1*param.beta_b*R_cb/PInext)/(1+param.phi_B*(B_bnext/B_b-1))*C_bnext^(-param.sigma2))^(-1/param.sigma2);

RHS(nx+RRa_ind)       = (Q+R_A)/Qminus;
RHS(nx+RR_ind)        = R_cbminus/PI;
RHS(nx+I_ind)         = iota*(1+param.phi/2*(log(x_k))^2)*(A_hhnext + A_gnext + A_bnext + SS_stats.A_F) - (1-param.delta_0*v^param.delta_1)*K;
RHS(nx+unemp_ind)     = (1-((1-l_lambda)*N + M)-param.Eshare)/(1-param.Eshare); 


RHS(nx+x_k_ind)       = Knext/K;

RHS(nx+nn_ind)        = ((1-tau_w)*W/(param.psi))^(1/param.xi);
RHS(nx+M_ind)         = ((1-N-Eshare)+l_lambda.*(N)).*V/((((1-N-Eshare)+l_lambda.*(N))^param.alpha+V^param.alpha)^(1/param.alpha));
RHS(nx+f_ind)         = M/((1-N-Eshare)+l_lambda*(N));



RHS(nx+A_g_obs_ind)   = (A_gauxnext);
RHS(nx+Y_obs_ind)     = Y;
RHS(nx+C_obs_ind)     = C;
RHS(nx+I_obs_ind)     = I;
RHS(nx+w_obs_ind)     = W;

RHS(nx+PROFIT_obs_ind)= PROFIT + param.fix2 + log(B_Fnext);
RHS(nx+unemp_obs_ind) = 100*unemp;
RHS(nx+inf_obs_ind)   = PI;
RHS(nx+R_obs_ind)     = R_cb;
RHS(nx+pastw_ind)     = Wminus;

switch(param.adjust)
	case('G')
		RHS(nx+pastGG_ind)      = pastLTminus;
	case('LT')
		RHS(nx+pastGG_ind)      = pastGminus;
end

RHS(nx+pastYY_ind)      = pastYminus;
RHS(nx+pastCC_ind)      = pastCminus;
RHS(nx+pastII_ind)      = pastIminus;
RHS(nx+pastPPROFIT_ind) = pastPROFITminus; 
RHS(nx+pastuu_ind)      = pastunempminus; 
RHS(nx+pastA_g_ind)     = A_g;

RHS(nx+l_lambda_ind)    = param.lambda;
RHS(nx+pastQ_ind)       = Qminus;
RHS(nx+x_I_ind)         = I/pastIminus;
RHS(nx+eta2_ind)        = eta;
RHS(nx+iota2_ind)       = iota;

RHS(nx+G_obs_ind)    = G;
RHS(nx+pastLT2_ind)  = pastLTminus;
RHS(nx+pastG2_ind)   = SS_stats.G;
RHS(nx+LT_obs_ind)   = LT;

RHS(nx+B_gov_ncp2_ind) = B_b  + B_hh - (Qminus*A_b- NW_bminus) - Qminus *(1+param.tau_cp)* A_g;

[mesh.a,mesh.b] = meshgrid(grid.a,grid.b);
mu_redux        = sum(MU,3);

grid.se_aux          = [grid.s,min(param.b_ratio*grid.s,grid.s_bar),1];
[~ ,~,meshes.se_aux] = ndgrid(grid.b,grid.a,grid.se_aux);

total_wealth         = mesh.a(:)*Q+mesh.b(:);
[total_wealth, IX]   = sort(total_wealth);
total_wealth_pdf     = mu_redux(IX);
total_wealth_cdf     = cumsum(total_wealth_pdf);

S_W        = cumsum(total_wealth_pdf.*total_wealth)';
S_W        = [0 S_W];
GiniWealth = 1-(sum(total_wealth_pdf.*(S_W(1:end-1)+S_W(2:end))')/S_W(end));


NW_Gini          = nn*W;
WW_Gini          = (1-param.tau_w)*(NW_Gini*ones(grid.nb,grid.na,grid.nse) + param.b_share*PROFIT/((1-param.u)*(1-param.Eshare)*grid.s2_bar).*auxWW);
WW_Gini(:,:,end) = (1-param.tau_a)*(param.Eratio*PROFIT+Profit_FI-param.fix2)/param.Eshare*ones(grid.nb,grid.na);     

total_wealth         = mesh.a(:)*Q+mesh.b(:);
[total_wealth, IX]   = sort(total_wealth);
total_wealth_pdf     = mu_redux(IX);
total_wealth_cdf     = cumsum(total_wealth_pdf);

income_1_mesh = WW_Gini.*meshes.se_aux                        + R_A_aux.*meshes.a + (R_cbminus/PI-1).*meshes.b.*(meshes.b>=0) + ((R_cbminus+param.Rprem)/PI-1).*meshes.b.*(meshes.b<0) + inc.transfer;
income_2_mesh = WW_Gini.*meshes.se_aux + (Q-Qminus).*meshes.a + R_A_aux.*meshes.a + (R_cbminus/PI-1).*meshes.b.*(meshes.b>=0) + ((R_cbminus+param.Rprem)/PI-1).*meshes.b.*(meshes.b<0) + inc.transfer;

[income_1_sort, IX] = sort(WW_Gini(:).*meshes.se_aux(:)                        + R_A_aux.*meshes.a(:) + (R_cbminus/PI-1).*meshes.b(:).*(meshes.b(:)>=0) + ((R_cbminus+param.Rprem)/PI-1).*meshes.b(:).*(meshes.b(:)<0) + inc.transfer(:));
income_1_pdf        = MU_tilde(IX);
income_1_cdf        = cumsum(income_1_pdf);

S_income_1          = cumsum(income_1_pdf.*income_1_sort)';
S_income_1          = [0 S_income_1];
GiniIncome_1        = 1-(sum(income_1_pdf.*(S_income_1(1:end-1)+S_income_1(2:end))')/S_income_1(end));

[income_2_sort, IX] = sort(WW_Gini(:).*meshes.se_aux(:) + (Q-Qminus).*(meshes.a(:)-(AProb(:).*a_a_star(:)-(1-AProb(:)).*meshes.a(:))) + R_A_aux.*meshes.a(:) + (R_cbminus/PI-1).*meshes.b(:).*(meshes.b(:)>=0) + ((R_cbminus+param.Rprem)/PI-1).*meshes.b(:).*(meshes.b(:)<0) + inc.transfer(:));
income_2_pdf        = MU_tilde(IX);
income_2_cdf        = cumsum(income_2_pdf);

S_income_2          = cumsum(income_2_pdf.*income_2_sort)';
S_income_2          = [0 S_income_2];
GiniIncome_2        = 1-(sum(income_2_pdf.*(S_income_2(1:end-1)+S_income_2(2:end))')/S_income_2(end));

cons_aux          = [c_a_star(:);c_n_star(:)];
mu_dist_tilde_aux = [AProb(:).*MU_tilde(:);(1-AProb(:)).*MU_tilde(:)];

[cons_sort, IX]   = sort(cons_aux);
cons_pdf          = mu_dist_tilde_aux(IX);
cons_cdf          = cumsum(cons_pdf);

S_cons            = cumsum(cons_pdf.*cons_sort)';
S_cons            = [0 S_cons];
GiniCons          = 1-(sum(cons_pdf.*(S_cons(1:end-1)+S_cons(2:end))')/S_cons(end));

meshes.w      = Q*meshes.a + meshes.b;

w0    = total_wealth(1);
w1    = total_wealth(sum(total_wealth_cdf<0.01));
w10   = total_wealth(sum(total_wealth_cdf<0.1));
w20   = total_wealth(sum(total_wealth_cdf<0.2));
w25   = total_wealth(sum(total_wealth_cdf<0.25));
w30   = total_wealth(sum(total_wealth_cdf<0.3));
w40   = total_wealth(sum(total_wealth_cdf<0.4));
w50   = total_wealth(sum(total_wealth_cdf<0.5));
w60   = total_wealth(sum(total_wealth_cdf<0.6));
w70   = total_wealth(sum(total_wealth_cdf<0.7));
w80   = total_wealth(sum(total_wealth_cdf<0.8));
w90   = total_wealth(sum(total_wealth_cdf<0.9));
w99   = total_wealth(sum(total_wealth_cdf<0.99));
w999  = total_wealth(sum(total_wealth_cdf<0.999));
w100  = total_wealth(end);

income0    = income_1_sort(1);
income01   = income_1_sort(sum(income_1_cdf<0.001));
income1    = income_1_sort(sum(income_1_cdf<0.01));
income10   = income_1_sort(sum(income_1_cdf<0.1));
income20   = income_1_sort(sum(income_1_cdf<0.2));
income25   = income_1_sort(sum(income_1_cdf<0.25));
income30   = income_1_sort(sum(income_1_cdf<0.3));
income40   = income_1_sort(sum(income_1_cdf<0.4));
income50   = income_1_sort(sum(income_1_cdf<0.5));
income60   = income_1_sort(sum(income_1_cdf<0.6));
income70   = income_1_sort(sum(income_1_cdf<0.7));
income80   = income_1_sort(sum(income_1_cdf<0.8));
income90   = income_1_sort(sum(income_1_cdf<0.9));
income99   = income_1_sort(sum(income_1_cdf<0.99));
income999  = income_1_sort(sum(income_1_cdf<0.999));
income100  = income_1_sort(end);

mu_dist_B01 = (meshes.w >= w0 ) .*(meshes.w <= w0  ).*MU;
mu_dist_B1  = (meshes.w >= w0 ) .*(meshes.w <= w1  ).*MU;
mu_dist_B10 = (meshes.w >= w0 ) .*(meshes.w <= w10 ).*MU;
mu_dist_Q1  = (meshes.w >= w0 ) .*(meshes.w <= w20 ).*MU;
mu_dist_Q2  = (meshes.w >  w20) .*(meshes.w <= w40 ).*MU;
mu_dist_Q3  = (meshes.w >  w40) .*(meshes.w <= w60 ).*MU;
mu_dist_Q4  = (meshes.w >  w60) .*(meshes.w <= w80 ).*MU;
mu_dist_Q5  = (meshes.w >  w80) .*(meshes.w <= w100).*MU;
mu_dist_P90 = (meshes.w >  w90) .*(meshes.w <= w100).*MU;
mu_dist_P99 = (meshes.w >  w99) .*(meshes.w <= w100).*MU;
mu_dist_P999= (meshes.w >  w999).*(meshes.w <= w100).*MU;


mu_dist_tilde_B01 = reshape(H_tilde*mu_dist_B01(:),[grid.nb grid.na grid.nse]);
mu_dist_tilde_B1  = reshape(H_tilde*mu_dist_B1(:),[grid.nb grid.na grid.nse]);
mu_dist_tilde_B10 = reshape(H_tilde*mu_dist_B10(:),[grid.nb grid.na grid.nse]);
mu_dist_tilde_Q1  = reshape(H_tilde*mu_dist_Q1(:),[grid.nb grid.na grid.nse]);
mu_dist_tilde_Q2  = reshape(H_tilde*mu_dist_Q2(:),[grid.nb grid.na grid.nse]);
mu_dist_tilde_Q3  = reshape(H_tilde*mu_dist_Q3(:),[grid.nb grid.na grid.nse]);
mu_dist_tilde_Q4  = reshape(H_tilde*mu_dist_Q4(:),[grid.nb grid.na grid.nse]);
mu_dist_tilde_Q5  = reshape(H_tilde*mu_dist_Q5(:),[grid.nb grid.na grid.nse]);
mu_dist_tilde_P90 = reshape(H_tilde*mu_dist_P90(:),[grid.nb grid.na grid.nse]);
mu_dist_tilde_P99 = reshape(H_tilde*mu_dist_P99(:),[grid.nb grid.na grid.nse]);
mu_dist_tilde_P999= reshape(H_tilde*mu_dist_P999(:),[grid.nb grid.na grid.nse]);

mu_dist_I_B01 = (income_1_mesh >= income0 ) .*(income_1_mesh <= income01 ).*MU;
mu_dist_I_B1  = (income_1_mesh >= income0 ) .*(income_1_mesh <= income1  ).*MU;
mu_dist_I_B10 = (income_1_mesh >= income0 ) .*(income_1_mesh <= income10 ).*MU;
mu_dist_I_Q1  = (income_1_mesh >= income0 ) .*(income_1_mesh <= income20 ).*MU;
mu_dist_I_Q2  = (income_1_mesh >  income20) .*(income_1_mesh <= income40 ).*MU;
mu_dist_I_Q3  = (income_1_mesh >  income40) .*(income_1_mesh <= income60 ).*MU;
mu_dist_I_Q4  = (income_1_mesh >  income60) .*(income_1_mesh <= income80 ).*MU;
mu_dist_I_Q5  = (income_1_mesh >  income80) .*(income_1_mesh <= income100).*MU;
mu_dist_I_P90 = (income_1_mesh >  income90) .*(income_1_mesh <= income100).*MU;
mu_dist_I_P99 = (income_1_mesh >  income99) .*(income_1_mesh <= income100).*MU;
mu_dist_I_P999= (income_1_mesh >  income999).*(income_1_mesh <= income100).*MU;


mu_dist_tilde_I_B01 = reshape(H_tilde*mu_dist_I_B01(:),[grid.nb grid.na grid.nse]);
mu_dist_tilde_I_B1  = reshape(H_tilde*mu_dist_I_B1(:),[grid.nb grid.na grid.nse]);
mu_dist_tilde_I_B10 = reshape(H_tilde*mu_dist_I_B10(:),[grid.nb grid.na grid.nse]);
mu_dist_tilde_I_Q1  = reshape(H_tilde*mu_dist_I_Q1(:),[grid.nb grid.na grid.nse]);
mu_dist_tilde_I_Q2  = reshape(H_tilde*mu_dist_I_Q2(:),[grid.nb grid.na grid.nse]);
mu_dist_tilde_I_Q3  = reshape(H_tilde*mu_dist_I_Q3(:),[grid.nb grid.na grid.nse]);
mu_dist_tilde_I_Q4  = reshape(H_tilde*mu_dist_I_Q4(:),[grid.nb grid.na grid.nse]);
mu_dist_tilde_I_Q5  = reshape(H_tilde*mu_dist_I_Q5(:),[grid.nb grid.na grid.nse]);
mu_dist_tilde_I_P90 = reshape(H_tilde*mu_dist_I_P90(:),[grid.nb grid.na grid.nse]);
mu_dist_tilde_I_P99 = reshape(H_tilde*mu_dist_I_P99(:),[grid.nb grid.na grid.nse]);
mu_dist_tilde_I_P999= reshape(H_tilde*mu_dist_I_P999(:),[grid.nb grid.na grid.nse]);


w_B1               = sum(sum(sum((AProb.*mu_dist_tilde_B1)   .*meshes.w))) + sum(sum(sum(((1-AProb).*mu_dist_tilde_B1)   .*meshes.w)));
w_B01              = sum(sum(sum((AProb.*mu_dist_tilde_B01)  .*meshes.w))) + sum(sum(sum(((1-AProb).*mu_dist_tilde_B01)  .*meshes.w)));
w_B10              = sum(sum(sum((AProb.*mu_dist_tilde_B10)  .*meshes.w))) + sum(sum(sum(((1-AProb).*mu_dist_tilde_B10)  .*meshes.w)));
w_Q1               = sum(sum(sum((AProb.*mu_dist_tilde_Q1)   .*meshes.w))) + sum(sum(sum(((1-AProb).*mu_dist_tilde_Q1)   .*meshes.w)));
w_Q2               = sum(sum(sum((AProb.*mu_dist_tilde_Q2)   .*meshes.w))) + sum(sum(sum(((1-AProb).*mu_dist_tilde_Q2)   .*meshes.w)));
w_Q3               = sum(sum(sum((AProb.*mu_dist_tilde_Q3)   .*meshes.w))) + sum(sum(sum(((1-AProb).*mu_dist_tilde_Q3)   .*meshes.w)));
w_Q4               = sum(sum(sum((AProb.*mu_dist_tilde_Q4)   .*meshes.w))) + sum(sum(sum(((1-AProb).*mu_dist_tilde_Q4)   .*meshes.w)));
w_Q5               = sum(sum(sum((AProb.*mu_dist_tilde_Q5)   .*meshes.w))) + sum(sum(sum(((1-AProb).*mu_dist_tilde_Q5)   .*meshes.w)));
w_P90              = sum(sum(sum((AProb.*mu_dist_tilde_P90)  .*meshes.w))) + sum(sum(sum(((1-AProb).*mu_dist_tilde_P90)  .*meshes.w)));
w_P99              = sum(sum(sum((AProb.*mu_dist_tilde_P99)  .*meshes.w))) + sum(sum(sum(((1-AProb).*mu_dist_tilde_P99)  .*meshes.w)));
w_P999             = sum(sum(sum((AProb.*mu_dist_tilde_P999) .*meshes.w))) + sum(sum(sum(((1-AProb).*mu_dist_tilde_P999) .*meshes.w)));

income_1_B01 = sum(sum(sum((WW_Gini.*meshes.se_aux + (R_A_aux).*meshes.a                        + ((R_cbminus/(1-param.death_rate)*param.b_a_aux)/PI-1).*meshes.b.*(meshes.b>=0) + (((R_cbminus+param.Rprem)/(1-param.death_rate)*param.b_a_aux)/PI-1).*meshes.b.*(meshes.b<0) + inc.transfer).*mu_dist_tilde_B01))) /sum(sum(sum(mu_dist_tilde_B01)));
income_1_B1  = sum(sum(sum((WW_Gini.*meshes.se_aux + (R_A_aux).*meshes.a                        + ((R_cbminus/(1-param.death_rate)*param.b_a_aux)/PI-1).*meshes.b.*(meshes.b>=0) + (((R_cbminus+param.Rprem)/(1-param.death_rate)*param.b_a_aux)/PI-1).*meshes.b.*(meshes.b<0) + inc.transfer).*mu_dist_tilde_B1)))  /sum(sum(sum(mu_dist_tilde_B1)));
income_1_B10 = sum(sum(sum((WW_Gini.*meshes.se_aux + (R_A_aux).*meshes.a                        + ((R_cbminus/(1-param.death_rate)*param.b_a_aux)/PI-1).*meshes.b.*(meshes.b>=0) + (((R_cbminus+param.Rprem)/(1-param.death_rate)*param.b_a_aux)/PI-1).*meshes.b.*(meshes.b<0) + inc.transfer).*mu_dist_tilde_B10))) /sum(sum(sum(mu_dist_tilde_B10)));
income_1_Q1  = sum(sum(sum((WW_Gini.*meshes.se_aux + (R_A_aux).*meshes.a                        + ((R_cbminus/(1-param.death_rate)*param.b_a_aux)/PI-1).*meshes.b.*(meshes.b>=0) + (((R_cbminus+param.Rprem)/(1-param.death_rate)*param.b_a_aux)/PI-1).*meshes.b.*(meshes.b<0) + inc.transfer).*mu_dist_tilde_Q1)))  /sum(sum(sum(mu_dist_tilde_Q1)));
income_1_Q2  = sum(sum(sum((WW_Gini.*meshes.se_aux + (R_A_aux).*meshes.a                        + ((R_cbminus/(1-param.death_rate)*param.b_a_aux)/PI-1).*meshes.b.*(meshes.b>=0) + (((R_cbminus+param.Rprem)/(1-param.death_rate)*param.b_a_aux)/PI-1).*meshes.b.*(meshes.b<0) + inc.transfer).*mu_dist_tilde_Q2)))  /sum(sum(sum(mu_dist_tilde_Q2)));
income_1_Q3  = sum(sum(sum((WW_Gini.*meshes.se_aux + (R_A_aux).*meshes.a                        + ((R_cbminus/(1-param.death_rate)*param.b_a_aux)/PI-1).*meshes.b.*(meshes.b>=0) + (((R_cbminus+param.Rprem)/(1-param.death_rate)*param.b_a_aux)/PI-1).*meshes.b.*(meshes.b<0) + inc.transfer).*mu_dist_tilde_Q3)))  /sum(sum(sum(mu_dist_tilde_Q3)));
income_1_Q4  = sum(sum(sum((WW_Gini.*meshes.se_aux + (R_A_aux).*meshes.a                        + ((R_cbminus/(1-param.death_rate)*param.b_a_aux)/PI-1).*meshes.b.*(meshes.b>=0) + (((R_cbminus+param.Rprem)/(1-param.death_rate)*param.b_a_aux)/PI-1).*meshes.b.*(meshes.b<0) + inc.transfer).*mu_dist_tilde_Q4)))  /sum(sum(sum(mu_dist_tilde_Q4)));
income_1_Q5  = sum(sum(sum((WW_Gini.*meshes.se_aux + (R_A_aux).*meshes.a                        + ((R_cbminus/(1-param.death_rate)*param.b_a_aux)/PI-1).*meshes.b.*(meshes.b>=0) + (((R_cbminus+param.Rprem)/(1-param.death_rate)*param.b_a_aux)/PI-1).*meshes.b.*(meshes.b<0) + inc.transfer).*mu_dist_tilde_Q5)))  /sum(sum(sum(mu_dist_tilde_Q5)));
income_1_P90 = sum(sum(sum((WW_Gini.*meshes.se_aux + (R_A_aux).*meshes.a                        + ((R_cbminus/(1-param.death_rate)*param.b_a_aux)/PI-1).*meshes.b.*(meshes.b>=0) + (((R_cbminus+param.Rprem)/(1-param.death_rate)*param.b_a_aux)/PI-1).*meshes.b.*(meshes.b<0) + inc.transfer).*mu_dist_tilde_P90))) /sum(sum(sum(mu_dist_tilde_P90)));
income_1_P99 = sum(sum(sum((WW_Gini.*meshes.se_aux + (R_A_aux).*meshes.a                        + ((R_cbminus/(1-param.death_rate)*param.b_a_aux)/PI-1).*meshes.b.*(meshes.b>=0) + (((R_cbminus+param.Rprem)/(1-param.death_rate)*param.b_a_aux)/PI-1).*meshes.b.*(meshes.b<0) + inc.transfer).*mu_dist_tilde_P99))) /sum(sum(sum(mu_dist_tilde_P99)));
income_1_P999= sum(sum(sum((WW_Gini.*meshes.se_aux + (R_A_aux).*meshes.a                        + ((R_cbminus/(1-param.death_rate)*param.b_a_aux)/PI-1).*meshes.b.*(meshes.b>=0) + (((R_cbminus+param.Rprem)/(1-param.death_rate)*param.b_a_aux)/PI-1).*meshes.b.*(meshes.b<0) + inc.transfer).*mu_dist_tilde_P999)))/sum(sum(sum(mu_dist_tilde_P999)));

income_2_B01 = sum(sum(sum((WW_Gini.*meshes.se_aux + (R_A_aux).*meshes.a + (Q-Qminus).*meshes.a + ((R_cbminus/(1-param.death_rate)*param.b_a_aux)/PI-1).*meshes.b.*(meshes.b>=0) + (((R_cbminus+param.Rprem)/(1-param.death_rate)*param.b_a_aux)/PI-1).*meshes.b.*(meshes.b<0) + inc.transfer).*mu_dist_tilde_B01))) /sum(sum(sum(mu_dist_tilde_B01)));
income_2_B1  = sum(sum(sum((WW_Gini.*meshes.se_aux + (R_A_aux).*meshes.a + (Q-Qminus).*meshes.a + ((R_cbminus/(1-param.death_rate)*param.b_a_aux)/PI-1).*meshes.b.*(meshes.b>=0) + (((R_cbminus+param.Rprem)/(1-param.death_rate)*param.b_a_aux)/PI-1).*meshes.b.*(meshes.b<0) + inc.transfer).*mu_dist_tilde_B1)))  /sum(sum(sum(mu_dist_tilde_B1)));
income_2_B10 = sum(sum(sum((WW_Gini.*meshes.se_aux + (R_A_aux).*meshes.a + (Q-Qminus).*meshes.a + ((R_cbminus/(1-param.death_rate)*param.b_a_aux)/PI-1).*meshes.b.*(meshes.b>=0) + (((R_cbminus+param.Rprem)/(1-param.death_rate)*param.b_a_aux)/PI-1).*meshes.b.*(meshes.b<0) + inc.transfer).*mu_dist_tilde_B10))) /sum(sum(sum(mu_dist_tilde_B10)));
income_2_Q1  = sum(sum(sum((WW_Gini.*meshes.se_aux + (R_A_aux).*meshes.a + (Q-Qminus).*meshes.a + ((R_cbminus/(1-param.death_rate)*param.b_a_aux)/PI-1).*meshes.b.*(meshes.b>=0) + (((R_cbminus+param.Rprem)/(1-param.death_rate)*param.b_a_aux)/PI-1).*meshes.b.*(meshes.b<0) + inc.transfer).*mu_dist_tilde_Q1)))  /sum(sum(sum(mu_dist_tilde_Q1)));
income_2_Q2  = sum(sum(sum((WW_Gini.*meshes.se_aux + (R_A_aux).*meshes.a + (Q-Qminus).*meshes.a + ((R_cbminus/(1-param.death_rate)*param.b_a_aux)/PI-1).*meshes.b.*(meshes.b>=0) + (((R_cbminus+param.Rprem)/(1-param.death_rate)*param.b_a_aux)/PI-1).*meshes.b.*(meshes.b<0) + inc.transfer).*mu_dist_tilde_Q2)))  /sum(sum(sum(mu_dist_tilde_Q2)));
income_2_Q3  = sum(sum(sum((WW_Gini.*meshes.se_aux + (R_A_aux).*meshes.a + (Q-Qminus).*meshes.a + ((R_cbminus/(1-param.death_rate)*param.b_a_aux)/PI-1).*meshes.b.*(meshes.b>=0) + (((R_cbminus+param.Rprem)/(1-param.death_rate)*param.b_a_aux)/PI-1).*meshes.b.*(meshes.b<0) + inc.transfer).*mu_dist_tilde_Q3)))  /sum(sum(sum(mu_dist_tilde_Q3)));
income_2_Q4  = sum(sum(sum((WW_Gini.*meshes.se_aux + (R_A_aux).*meshes.a + (Q-Qminus).*meshes.a + ((R_cbminus/(1-param.death_rate)*param.b_a_aux)/PI-1).*meshes.b.*(meshes.b>=0) + (((R_cbminus+param.Rprem)/(1-param.death_rate)*param.b_a_aux)/PI-1).*meshes.b.*(meshes.b<0) + inc.transfer).*mu_dist_tilde_Q4)))  /sum(sum(sum(mu_dist_tilde_Q4)));
income_2_Q5  = sum(sum(sum((WW_Gini.*meshes.se_aux + (R_A_aux).*meshes.a + (Q-Qminus).*meshes.a + ((R_cbminus/(1-param.death_rate)*param.b_a_aux)/PI-1).*meshes.b.*(meshes.b>=0) + (((R_cbminus+param.Rprem)/(1-param.death_rate)*param.b_a_aux)/PI-1).*meshes.b.*(meshes.b<0) + inc.transfer).*mu_dist_tilde_Q5)))  /sum(sum(sum(mu_dist_tilde_Q5)));
income_2_P90 = sum(sum(sum((WW_Gini.*meshes.se_aux + (R_A_aux).*meshes.a + (Q-Qminus).*meshes.a + ((R_cbminus/(1-param.death_rate)*param.b_a_aux)/PI-1).*meshes.b.*(meshes.b>=0) + (((R_cbminus+param.Rprem)/(1-param.death_rate)*param.b_a_aux)/PI-1).*meshes.b.*(meshes.b<0) + inc.transfer).*mu_dist_tilde_P90))) /sum(sum(sum(mu_dist_tilde_P90)));
income_2_P99 = sum(sum(sum((WW_Gini.*meshes.se_aux + (R_A_aux).*meshes.a + (Q-Qminus).*meshes.a + ((R_cbminus/(1-param.death_rate)*param.b_a_aux)/PI-1).*meshes.b.*(meshes.b>=0) + (((R_cbminus+param.Rprem)/(1-param.death_rate)*param.b_a_aux)/PI-1).*meshes.b.*(meshes.b<0) + inc.transfer).*mu_dist_tilde_P99))) /sum(sum(sum(mu_dist_tilde_P99)));
income_2_P999= sum(sum(sum((WW_Gini.*meshes.se_aux + (R_A_aux).*meshes.a + (Q-Qminus).*meshes.a + ((R_cbminus/(1-param.death_rate)*param.b_a_aux)/PI-1).*meshes.b.*(meshes.b>=0) + (((R_cbminus+param.Rprem)/(1-param.death_rate)*param.b_a_aux)/PI-1).*meshes.b.*(meshes.b<0) + inc.transfer).*mu_dist_tilde_P999)))/sum(sum(sum(mu_dist_tilde_P999)));


C_B01    = q_cons2(c_a_star,c_n_star,W,nn,AProb,mu_dist_tilde_B01  ,meshes,param,grid,param.tau_w);
C_B1     = q_cons2(c_a_star,c_n_star,W,nn,AProb,mu_dist_tilde_B1   ,meshes,param,grid,param.tau_w);
C_B10    = q_cons2(c_a_star,c_n_star,W,nn,AProb,mu_dist_tilde_B10  ,meshes,param,grid,param.tau_w);
C_Q1     = q_cons2(c_a_star,c_n_star,W,nn,AProb,mu_dist_tilde_Q1   ,meshes,param,grid,param.tau_w);
C_Q2     = q_cons2(c_a_star,c_n_star,W,nn,AProb,mu_dist_tilde_Q2   ,meshes,param,grid,param.tau_w);
C_Q3     = q_cons2(c_a_star,c_n_star,W,nn,AProb,mu_dist_tilde_Q3   ,meshes,param,grid,param.tau_w);
C_Q4     = q_cons2(c_a_star,c_n_star,W,nn,AProb,mu_dist_tilde_Q4   ,meshes,param,grid,param.tau_w);
C_Q5     = q_cons2(c_a_star,c_n_star,W,nn,AProb,mu_dist_tilde_Q5   ,meshes,param,grid,param.tau_w);
C_P90    = q_cons2(c_a_star,c_n_star,W,nn,AProb,mu_dist_tilde_P90  ,meshes,param,grid,param.tau_w);
C_P99    = q_cons2(c_a_star,c_n_star,W,nn,AProb,mu_dist_tilde_P99  ,meshes,param,grid,param.tau_w);
C_P999   = q_cons2(c_a_star,c_n_star,W,nn,AProb,mu_dist_tilde_P999 ,meshes,param,grid,param.tau_w);

w_I_B1               = sum(sum(sum((AProb.*mu_dist_tilde_I_B1)   .*meshes.w))) + sum(sum(sum(((1-AProb).*mu_dist_tilde_I_B1)   .*meshes.w)));
w_I_B01              = sum(sum(sum((AProb.*mu_dist_tilde_I_B01)  .*meshes.w))) + sum(sum(sum(((1-AProb).*mu_dist_tilde_I_B01)  .*meshes.w)));
w_I_B10              = sum(sum(sum((AProb.*mu_dist_tilde_I_B10)  .*meshes.w))) + sum(sum(sum(((1-AProb).*mu_dist_tilde_I_B10)  .*meshes.w)));
w_I_Q1               = sum(sum(sum((AProb.*mu_dist_tilde_I_Q1)   .*meshes.w))) + sum(sum(sum(((1-AProb).*mu_dist_tilde_I_Q1)   .*meshes.w)));
w_I_Q2               = sum(sum(sum((AProb.*mu_dist_tilde_I_Q2)   .*meshes.w))) + sum(sum(sum(((1-AProb).*mu_dist_tilde_I_Q2)   .*meshes.w)));
w_I_Q3               = sum(sum(sum((AProb.*mu_dist_tilde_I_Q3)   .*meshes.w))) + sum(sum(sum(((1-AProb).*mu_dist_tilde_I_Q3)   .*meshes.w)));
w_I_Q4               = sum(sum(sum((AProb.*mu_dist_tilde_I_Q4)   .*meshes.w))) + sum(sum(sum(((1-AProb).*mu_dist_tilde_I_Q4)   .*meshes.w)));
w_I_Q5               = sum(sum(sum((AProb.*mu_dist_tilde_I_Q5)   .*meshes.w))) + sum(sum(sum(((1-AProb).*mu_dist_tilde_I_Q5)   .*meshes.w)));
w_I_P90              = sum(sum(sum((AProb.*mu_dist_tilde_I_P90)  .*meshes.w))) + sum(sum(sum(((1-AProb).*mu_dist_tilde_I_P90)  .*meshes.w)));
w_I_P99              = sum(sum(sum((AProb.*mu_dist_tilde_I_P99)  .*meshes.w))) + sum(sum(sum(((1-AProb).*mu_dist_tilde_I_P99)  .*meshes.w)));
w_I_P999             = sum(sum(sum((AProb.*mu_dist_tilde_I_P999) .*meshes.w))) + sum(sum(sum(((1-AProb).*mu_dist_tilde_I_P999) .*meshes.w)));

income_1_I_B01 = sum(sum(sum((WW_Gini.*meshes.se_aux + (R_A_aux).*meshes.a                        + ((R_cbminus/(1-param.death_rate)*param.b_a_aux)/PI-1).*meshes.b.*(meshes.b>=0) + (((R_cbminus+param.Rprem)/(1-param.death_rate)*param.b_a_aux)/PI-1).*meshes.b.*(meshes.b<0) + inc.transfer).*mu_dist_tilde_I_B01))) /sum(sum(sum(mu_dist_tilde_I_B01)));
income_1_I_B1  = sum(sum(sum((WW_Gini.*meshes.se_aux + (R_A_aux).*meshes.a                        + ((R_cbminus/(1-param.death_rate)*param.b_a_aux)/PI-1).*meshes.b.*(meshes.b>=0) + (((R_cbminus+param.Rprem)/(1-param.death_rate)*param.b_a_aux)/PI-1).*meshes.b.*(meshes.b<0) + inc.transfer).*mu_dist_tilde_I_B1)))  /sum(sum(sum(mu_dist_tilde_I_B1)));
income_1_I_B10 = sum(sum(sum((WW_Gini.*meshes.se_aux + (R_A_aux).*meshes.a                        + ((R_cbminus/(1-param.death_rate)*param.b_a_aux)/PI-1).*meshes.b.*(meshes.b>=0) + (((R_cbminus+param.Rprem)/(1-param.death_rate)*param.b_a_aux)/PI-1).*meshes.b.*(meshes.b<0) + inc.transfer).*mu_dist_tilde_I_B10))) /sum(sum(sum(mu_dist_tilde_I_B10)));
income_1_I_Q1  = sum(sum(sum((WW_Gini.*meshes.se_aux + (R_A_aux).*meshes.a                        + ((R_cbminus/(1-param.death_rate)*param.b_a_aux)/PI-1).*meshes.b.*(meshes.b>=0) + (((R_cbminus+param.Rprem)/(1-param.death_rate)*param.b_a_aux)/PI-1).*meshes.b.*(meshes.b<0) + inc.transfer).*mu_dist_tilde_I_Q1)))  /sum(sum(sum(mu_dist_tilde_I_Q1)));
income_1_I_Q2  = sum(sum(sum((WW_Gini.*meshes.se_aux + (R_A_aux).*meshes.a                        + ((R_cbminus/(1-param.death_rate)*param.b_a_aux)/PI-1).*meshes.b.*(meshes.b>=0) + (((R_cbminus+param.Rprem)/(1-param.death_rate)*param.b_a_aux)/PI-1).*meshes.b.*(meshes.b<0) + inc.transfer).*mu_dist_tilde_I_Q2)))  /sum(sum(sum(mu_dist_tilde_I_Q2)));
income_1_I_Q3  = sum(sum(sum((WW_Gini.*meshes.se_aux + (R_A_aux).*meshes.a                        + ((R_cbminus/(1-param.death_rate)*param.b_a_aux)/PI-1).*meshes.b.*(meshes.b>=0) + (((R_cbminus+param.Rprem)/(1-param.death_rate)*param.b_a_aux)/PI-1).*meshes.b.*(meshes.b<0) + inc.transfer).*mu_dist_tilde_I_Q3)))  /sum(sum(sum(mu_dist_tilde_I_Q3)));
income_1_I_Q4  = sum(sum(sum((WW_Gini.*meshes.se_aux + (R_A_aux).*meshes.a                        + ((R_cbminus/(1-param.death_rate)*param.b_a_aux)/PI-1).*meshes.b.*(meshes.b>=0) + (((R_cbminus+param.Rprem)/(1-param.death_rate)*param.b_a_aux)/PI-1).*meshes.b.*(meshes.b<0) + inc.transfer).*mu_dist_tilde_I_Q4)))  /sum(sum(sum(mu_dist_tilde_I_Q4)));
income_1_I_Q5  = sum(sum(sum((WW_Gini.*meshes.se_aux + (R_A_aux).*meshes.a                        + ((R_cbminus/(1-param.death_rate)*param.b_a_aux)/PI-1).*meshes.b.*(meshes.b>=0) + (((R_cbminus+param.Rprem)/(1-param.death_rate)*param.b_a_aux)/PI-1).*meshes.b.*(meshes.b<0) + inc.transfer).*mu_dist_tilde_I_Q5)))  /sum(sum(sum(mu_dist_tilde_I_Q5)));
income_1_I_P90 = sum(sum(sum((WW_Gini.*meshes.se_aux + (R_A_aux).*meshes.a                        + ((R_cbminus/(1-param.death_rate)*param.b_a_aux)/PI-1).*meshes.b.*(meshes.b>=0) + (((R_cbminus+param.Rprem)/(1-param.death_rate)*param.b_a_aux)/PI-1).*meshes.b.*(meshes.b<0) + inc.transfer).*mu_dist_tilde_I_P90))) /sum(sum(sum(mu_dist_tilde_I_P90)));
income_1_I_P99 = sum(sum(sum((WW_Gini.*meshes.se_aux + (R_A_aux).*meshes.a                        + ((R_cbminus/(1-param.death_rate)*param.b_a_aux)/PI-1).*meshes.b.*(meshes.b>=0) + (((R_cbminus+param.Rprem)/(1-param.death_rate)*param.b_a_aux)/PI-1).*meshes.b.*(meshes.b<0) + inc.transfer).*mu_dist_tilde_I_P99))) /sum(sum(sum(mu_dist_tilde_I_P99)));
income_1_I_P999= sum(sum(sum((WW_Gini.*meshes.se_aux + (R_A_aux).*meshes.a                        + ((R_cbminus/(1-param.death_rate)*param.b_a_aux)/PI-1).*meshes.b.*(meshes.b>=0) + (((R_cbminus+param.Rprem)/(1-param.death_rate)*param.b_a_aux)/PI-1).*meshes.b.*(meshes.b<0) + inc.transfer).*mu_dist_tilde_I_P999)))/sum(sum(sum(mu_dist_tilde_I_P999)));

income_2_I_B01 = sum(sum(sum((WW_Gini.*meshes.se_aux + (R_A_aux).*meshes.a + (Q-Qminus).*meshes.a + ((R_cbminus/(1-param.death_rate)*param.b_a_aux)/PI-1).*meshes.b.*(meshes.b>=0) + (((R_cbminus+param.Rprem)/(1-param.death_rate)*param.b_a_aux)/PI-1).*meshes.b.*(meshes.b<0) + inc.transfer).*mu_dist_tilde_I_B01))) /sum(sum(sum(mu_dist_tilde_I_B01)));
income_2_I_B1  = sum(sum(sum((WW_Gini.*meshes.se_aux + (R_A_aux).*meshes.a + (Q-Qminus).*meshes.a + ((R_cbminus/(1-param.death_rate)*param.b_a_aux)/PI-1).*meshes.b.*(meshes.b>=0) + (((R_cbminus+param.Rprem)/(1-param.death_rate)*param.b_a_aux)/PI-1).*meshes.b.*(meshes.b<0) + inc.transfer).*mu_dist_tilde_I_B1)))  /sum(sum(sum(mu_dist_tilde_I_B1)));
income_2_I_B10 = sum(sum(sum((WW_Gini.*meshes.se_aux + (R_A_aux).*meshes.a + (Q-Qminus).*meshes.a + ((R_cbminus/(1-param.death_rate)*param.b_a_aux)/PI-1).*meshes.b.*(meshes.b>=0) + (((R_cbminus+param.Rprem)/(1-param.death_rate)*param.b_a_aux)/PI-1).*meshes.b.*(meshes.b<0) + inc.transfer).*mu_dist_tilde_I_B10))) /sum(sum(sum(mu_dist_tilde_I_B10)));
income_2_I_Q1  = sum(sum(sum((WW_Gini.*meshes.se_aux + (R_A_aux).*meshes.a + (Q-Qminus).*meshes.a + ((R_cbminus/(1-param.death_rate)*param.b_a_aux)/PI-1).*meshes.b.*(meshes.b>=0) + (((R_cbminus+param.Rprem)/(1-param.death_rate)*param.b_a_aux)/PI-1).*meshes.b.*(meshes.b<0) + inc.transfer).*mu_dist_tilde_I_Q1)))  /sum(sum(sum(mu_dist_tilde_I_Q1)));
income_2_I_Q2  = sum(sum(sum((WW_Gini.*meshes.se_aux + (R_A_aux).*meshes.a + (Q-Qminus).*meshes.a + ((R_cbminus/(1-param.death_rate)*param.b_a_aux)/PI-1).*meshes.b.*(meshes.b>=0) + (((R_cbminus+param.Rprem)/(1-param.death_rate)*param.b_a_aux)/PI-1).*meshes.b.*(meshes.b<0) + inc.transfer).*mu_dist_tilde_I_Q2)))  /sum(sum(sum(mu_dist_tilde_I_Q2)));
income_2_I_Q3  = sum(sum(sum((WW_Gini.*meshes.se_aux + (R_A_aux).*meshes.a + (Q-Qminus).*meshes.a + ((R_cbminus/(1-param.death_rate)*param.b_a_aux)/PI-1).*meshes.b.*(meshes.b>=0) + (((R_cbminus+param.Rprem)/(1-param.death_rate)*param.b_a_aux)/PI-1).*meshes.b.*(meshes.b<0) + inc.transfer).*mu_dist_tilde_I_Q3)))  /sum(sum(sum(mu_dist_tilde_I_Q3)));
income_2_I_Q4  = sum(sum(sum((WW_Gini.*meshes.se_aux + (R_A_aux).*meshes.a + (Q-Qminus).*meshes.a + ((R_cbminus/(1-param.death_rate)*param.b_a_aux)/PI-1).*meshes.b.*(meshes.b>=0) + (((R_cbminus+param.Rprem)/(1-param.death_rate)*param.b_a_aux)/PI-1).*meshes.b.*(meshes.b<0) + inc.transfer).*mu_dist_tilde_I_Q4)))  /sum(sum(sum(mu_dist_tilde_I_Q4)));
income_2_I_Q5  = sum(sum(sum((WW_Gini.*meshes.se_aux + (R_A_aux).*meshes.a + (Q-Qminus).*meshes.a + ((R_cbminus/(1-param.death_rate)*param.b_a_aux)/PI-1).*meshes.b.*(meshes.b>=0) + (((R_cbminus+param.Rprem)/(1-param.death_rate)*param.b_a_aux)/PI-1).*meshes.b.*(meshes.b<0) + inc.transfer).*mu_dist_tilde_I_Q5)))  /sum(sum(sum(mu_dist_tilde_I_Q5)));
income_2_I_P90 = sum(sum(sum((WW_Gini.*meshes.se_aux + (R_A_aux).*meshes.a + (Q-Qminus).*meshes.a + ((R_cbminus/(1-param.death_rate)*param.b_a_aux)/PI-1).*meshes.b.*(meshes.b>=0) + (((R_cbminus+param.Rprem)/(1-param.death_rate)*param.b_a_aux)/PI-1).*meshes.b.*(meshes.b<0) + inc.transfer).*mu_dist_tilde_I_P90))) /sum(sum(sum(mu_dist_tilde_I_P90)));
income_2_I_P99 = sum(sum(sum((WW_Gini.*meshes.se_aux + (R_A_aux).*meshes.a + (Q-Qminus).*meshes.a + ((R_cbminus/(1-param.death_rate)*param.b_a_aux)/PI-1).*meshes.b.*(meshes.b>=0) + (((R_cbminus+param.Rprem)/(1-param.death_rate)*param.b_a_aux)/PI-1).*meshes.b.*(meshes.b<0) + inc.transfer).*mu_dist_tilde_I_P99))) /sum(sum(sum(mu_dist_tilde_I_P99)));
income_2_I_P999= sum(sum(sum((WW_Gini.*meshes.se_aux + (R_A_aux).*meshes.a + (Q-Qminus).*meshes.a + ((R_cbminus/(1-param.death_rate)*param.b_a_aux)/PI-1).*meshes.b.*(meshes.b>=0) + (((R_cbminus+param.Rprem)/(1-param.death_rate)*param.b_a_aux)/PI-1).*meshes.b.*(meshes.b<0) + inc.transfer).*mu_dist_tilde_I_P999)))/sum(sum(sum(mu_dist_tilde_I_P999)));


C_I_B01    = q_cons2(c_a_star,c_n_star,W,nn,AProb,mu_dist_tilde_I_B01  ,meshes,param,grid,param.tau_w);
C_I_B1     = q_cons2(c_a_star,c_n_star,W,nn,AProb,mu_dist_tilde_I_B1   ,meshes,param,grid,param.tau_w);
C_I_B10    = q_cons2(c_a_star,c_n_star,W,nn,AProb,mu_dist_tilde_I_B10  ,meshes,param,grid,param.tau_w);
C_I_Q1     = q_cons2(c_a_star,c_n_star,W,nn,AProb,mu_dist_tilde_I_Q1   ,meshes,param,grid,param.tau_w);
C_I_Q2     = q_cons2(c_a_star,c_n_star,W,nn,AProb,mu_dist_tilde_I_Q2   ,meshes,param,grid,param.tau_w);
C_I_Q3     = q_cons2(c_a_star,c_n_star,W,nn,AProb,mu_dist_tilde_I_Q3   ,meshes,param,grid,param.tau_w);
C_I_Q4     = q_cons2(c_a_star,c_n_star,W,nn,AProb,mu_dist_tilde_I_Q4   ,meshes,param,grid,param.tau_w);
C_I_Q5     = q_cons2(c_a_star,c_n_star,W,nn,AProb,mu_dist_tilde_I_Q5   ,meshes,param,grid,param.tau_w);
C_I_P90    = q_cons2(c_a_star,c_n_star,W,nn,AProb,mu_dist_tilde_I_P90  ,meshes,param,grid,param.tau_w);
C_I_P99    = q_cons2(c_a_star,c_n_star,W,nn,AProb,mu_dist_tilde_I_P99  ,meshes,param,grid,param.tau_w);
C_I_P999   = q_cons2(c_a_star,c_n_star,W,nn,AProb,mu_dist_tilde_I_P999 ,meshes,param,grid,param.tau_w);  


RHS(nx+GiniW_ind)   = GiniWealth ; 
RHS(nx+GiniI_1_ind) = GiniIncome_1; 
RHS(nx+GiniI_2_ind) = GiniIncome_2; 
RHS(nx+GiniC_ind)   = GiniCons; 

RHS(nx+w_01_ind)  = w_B01+1;
RHS(nx+w_1_ind)   = w_B1+1;
RHS(nx+w_10_ind)  = w_B10+1;
RHS(nx+w_Q1_ind)  = w_Q1+1;
RHS(nx+w_Q2_ind)  = w_Q2;
RHS(nx+w_Q3_ind)  = w_Q3;
RHS(nx+w_Q4_ind)  = w_Q4;
RHS(nx+w_Q5_ind)  = w_Q5;
RHS(nx+w_90_ind)  = w_P90;
RHS(nx+w_99_ind)  = w_P99;
RHS(nx+w_999_ind) = w_P999;


RHS(nx+I_1_01_ind)  = income_1_B01 ; 
RHS(nx+I_1_1_ind)   = income_1_B1 ; 
RHS(nx+I_1_10_ind)  = income_1_B10 ; 
RHS(nx+I_1_Q1_ind)  = income_1_Q1 ; 
RHS(nx+I_1_Q2_ind)  = income_1_Q2 ; 
RHS(nx+I_1_Q3_ind)  = income_1_Q3 ; 
RHS(nx+I_1_Q4_ind)  = income_1_Q4 ; 
RHS(nx+I_1_Q5_ind)  = income_1_Q5 ; 
RHS(nx+I_1_90_ind)  = income_1_P90 ; 
RHS(nx+I_1_99_ind)  = income_1_P99 ; 
RHS(nx+I_1_999_ind) = income_1_P999 ;

RHS(nx+I_2_01_ind)  = income_2_B01 ; 
RHS(nx+I_2_1_ind)   = income_2_B1 ; 
RHS(nx+I_2_10_ind)  = income_2_B10 ; 
RHS(nx+I_2_Q1_ind)  = income_2_Q1 ; 
RHS(nx+I_2_Q2_ind)  = income_2_Q2 ; 
RHS(nx+I_2_Q3_ind)  = income_2_Q3 ; 
RHS(nx+I_2_Q4_ind)  = income_2_Q4 ; 
RHS(nx+I_2_Q5_ind)  = income_2_Q5 ; 
RHS(nx+I_2_90_ind)  = income_2_P90 ; 
RHS(nx+I_2_99_ind)  = income_2_P99 ; 
RHS(nx+I_2_999_ind) = income_2_P999 ; 

RHS(nx+C_01_ind)  = C_B01;
RHS(nx+C_1_ind)   = C_B1;
RHS(nx+C_10_ind)  = C_B10;
RHS(nx+C_Q1_ind)  = C_Q1;
RHS(nx+C_Q2_ind)  = C_Q2;
RHS(nx+C_Q3_ind)  = C_Q3;
RHS(nx+C_Q4_ind)  = C_Q4;
RHS(nx+C_Q5_ind)  = C_Q5;
RHS(nx+C_90_ind)  = C_P90;
RHS(nx+C_99_ind)  = C_P99;
RHS(nx+C_999_ind) = C_P999;

RHS(nx+w_I_01_ind)  = w_I_B01+1;;
RHS(nx+w_I_1_ind)   = w_I_B1+1;;
RHS(nx+w_I_10_ind)  = w_I_B10+1;;
RHS(nx+w_I_Q1_ind)  = w_I_Q1+1;;
RHS(nx+w_I_Q2_ind)  = w_I_Q2;;
RHS(nx+w_I_Q3_ind)  = w_I_Q3;;
RHS(nx+w_I_Q4_ind)  = w_I_Q4;;
RHS(nx+w_I_Q5_ind)  = w_I_Q5;;
RHS(nx+w_I_90_ind)  = w_I_P90;;
RHS(nx+w_I_99_ind)  = w_I_P99;;
RHS(nx+w_I_999_ind) = w_I_P999;;

RHS(nx+I_I_1_01_ind)  = income_1_I_B01; 
RHS(nx+I_I_1_1_ind)   = income_1_I_B1; 
RHS(nx+I_I_1_10_ind)  = income_1_I_B10; 
RHS(nx+I_I_1_Q1_ind)  = income_1_I_Q1; 
RHS(nx+I_I_1_Q2_ind)  = income_1_I_Q2; 
RHS(nx+I_I_1_Q3_ind)  = income_1_I_Q3; 
RHS(nx+I_I_1_Q4_ind)  = income_1_I_Q4; 
RHS(nx+I_I_1_Q5_ind)  = income_1_I_Q5; 
RHS(nx+I_I_1_90_ind)  = income_1_I_P90; 
RHS(nx+I_I_1_99_ind)  = income_1_I_P99; 
RHS(nx+I_I_1_999_ind) = income_1_I_P999; 

RHS(nx+I_I_2_01_ind)  = income_2_I_B01; 
RHS(nx+I_I_2_1_ind)   = income_2_I_B1; 
RHS(nx+I_I_2_10_ind)  = income_2_I_B10; 
RHS(nx+I_I_2_Q1_ind)  = income_2_I_Q1; 
RHS(nx+I_I_2_Q2_ind)  = income_2_I_Q2; 
RHS(nx+I_I_2_Q3_ind)  = income_2_I_Q3; 
RHS(nx+I_I_2_Q4_ind)  = income_2_I_Q4; 
RHS(nx+I_I_2_Q5_ind)  = income_2_I_Q5; 
RHS(nx+I_I_2_90_ind)  = income_2_I_P90; 
RHS(nx+I_I_2_99_ind)  = income_2_I_P99; 
RHS(nx+I_I_2_999_ind) = income_2_I_P999; 

RHS(nx+C_I_01_ind)  = C_I_B01; 
RHS(nx+C_I_1_ind)   = C_I_B1; 
RHS(nx+C_I_10_ind)  = C_I_B10; 
RHS(nx+C_I_Q1_ind)  = C_I_Q1; 
RHS(nx+C_I_Q2_ind)  = C_I_Q2; 
RHS(nx+C_I_Q3_ind)  = C_I_Q3; 
RHS(nx+C_I_Q4_ind)  = C_I_Q4; 
RHS(nx+C_I_Q5_ind)  = C_I_Q5; 
RHS(nx+C_I_90_ind)  = C_I_P90; 
RHS(nx+C_I_99_ind)  = C_I_P99; 
RHS(nx+C_I_999_ind) = C_I_P999; 

RHS(nx+I_90to10_ind) = income90/income10;
RHS(nx+I_50to10_ind) = income50/income10;  
RHS(nx+I_90to50_ind) = income90/income50;
RHS(nx+W_90to50_ind) = w90/w50; 

RHS(nx+Inc_1_T10toB10_ind) = income_1_P90/income_1_B10;
RHS(nx+Inc_1_T10toQ3_ind)  = income_1_P90/income_1_Q3;
RHS(nx+Inc_1_Q3toB10_ind)  = income_1_Q3/income_1_B10;
RHS(nx+Inc_2_T10toB10_ind) = income_2_P90/income_2_B10;
RHS(nx+Inc_2_T10toQ3_ind)  = income_2_P90/income_2_Q3;
RHS(nx+Inc_2_Q3toB10_ind)  = income_2_Q3/income_2_B10;
RHS(nx+Wlth_T10toQ3_ind)   = w_P90/w_Q3;
RHS(nx+C_T10toB10_ind)     = C_P90/C_B10;
RHS(nx+C_T10toQ3_ind)      = C_P90/C_Q3;
RHS(nx+C_Q3toB10_ind)      = C_Q3/C_B10;

RHS(nx+Inc_I_1_T10toB10_ind) = income_1_I_P90/income_1_I_B10;
RHS(nx+Inc_I_1_T10toQ3_ind)  = income_1_I_P90/income_1_I_Q3;
RHS(nx+Inc_I_1_Q3toB10_ind)  = income_1_I_Q3/income_1_I_B10;
RHS(nx+Inc_I_2_T10toB10_ind) = income_2_I_P90/income_2_I_B10;
RHS(nx+Inc_I_2_T10toQ3_ind)  = income_2_I_P90/income_2_I_Q3;
RHS(nx+Inc_I_2_Q3toB10_ind)  = income_2_I_Q3/income_2_I_B10;
RHS(nx+Wlth_I_T10toQ3_ind)   = w_I_P90/w_I_Q3;
RHS(nx+C_I_T10toB10_ind)     = C_I_P90/C_I_B10;
RHS(nx+C_I_T10toQ3_ind)      = C_I_P90/C_I_Q3;
RHS(nx+C_I_Q3toB10_ind)      = C_I_Q3/C_I_B10;


RHS(nx+Profit_aux_ind) = (Y *(1-eta/(2*param.kappa)*(log(PI)-(1-param.GAMMA)*log(param.pi_cb)-param.GAMMA*log(pastpiminus)).^2)+ Y.*(-MC) - (param.fix_ratio+log(B_Fnext))*SS_stats.Y  ...
		                + (H-param.fix_L-W).*L - param.iota.*V)  ...
		                + (R_K*v - param.delta_0.*v.^param.delta_1).*K +  Q*(Knext-K)-(Knext-K) - param.phi/2*(Knext/K-1)^2*K; 		 		 


RHS(nx+w2_ind)          = exp(log(param.w_bar) + param.rho_w * log(Wminus/param.w_bar) + param.rho_w*(param.d*log(param.pi_bar/PI)+(1-param.d)*log(pastpiminus/PI)) + (1-param.rho_w)*param.eps_w*log(PSI_W*H/SS_stats.r_l));         
RHS(nx+Profit2_ind)     = (Y *(1-eta/(2*param.kappa)*(log(PI)-(1-param.GAMMA)*log(param.pi_cb)-param.GAMMA*log(pastpiminus)).^2)+ Y.*(-MC) - param.fix  ...
		                   + (H-param.fix_L-W).*L - param.iota.*V)  ...
		                   + (R_K*v - param.delta_0.*v.^param.delta_1).*K +  Q*(Knext-K)-(Knext-K) - param.phi/2*(Knext/K-1)^2*K; 		 		          
RHS(nx+Profit3_ind)     = (Y *(1-eta/(2*param.kappa)*(log(PI)-(1-param.GAMMA)*log(param.pi_cb)-param.GAMMA*log(pastpiminus)).^2)+ Y.*(-MC) - param.fix  ...
		                   + (H-param.fix_L-W).*L - param.iota.*V)  ...
		                   + (R_K*v - param.delta_0.*v.^param.delta_1).*K +  Q*(Knext-K)-(Knext-K) - param.phi/2*(Knext/K-1)^2*K  + Profit_FI - param.fix2 ; 		 		          
RHS(nx+Profit_FI2_ind)  = ((1-param.theta_b)*((RRa - param.b_a_aux2*R_cbminus/PI)*levminus+param.b_a_aux2*R_cbminus/PI)*NW_bminus - param.omega*Qminus*A_b);         
RHS(nx+r_a2_ind)        = (1-tau_a)*(1-param.Eratio-param.b_share)*PROFIT/K;


%% Difference
Difference=InvGamma(:,:)*((LHS-RHS)./[ones(nx,1);(ControlSS(1:end-grid.oc));ones(grid.oc,1)]);


%% Calculate check on quality of approximation of Distribution
cumdist = zeros(grid.nb+1,grid.na+1,grid.nse+1);
cumdist(2:end,2:end,2:end) = Copula({cumsum(marginal_b),cumsum(marginal_a)',cumsum(marginal_se)});
MUMU    = diff(diff(diff(cumdist,1,1),1,2),1,3);
IMSE    = sum((MUnext(:)-MUMU(:)).^2);
R2like  = 1-IMSE/sum((MUnext(:)-1./numel(MUnext)).^2);


end

function [ weight,index ] = genweight( x,xgrid )
% function: GENWEIGHT generates weights and indexes used for linear interpolation
%
% X: Points at which function is to be interpolated.
% xgrid: grid points at which function is measured
% no extrapolation allowed
[~,index] = histc(x,xgrid);
index(x<=xgrid(1))=1;
index(x>=xgrid(end))=length(xgrid)-1;

weight = (x-xgrid(index))./(xgrid(index+1)-xgrid(index)); % weight xm of higher gridpoint
weight(weight<=0) = 1.e-16; % no extrapolation
weight(weight>=1) = 1-1.e-16; % no extrapolation

end  % function


