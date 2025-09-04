% TFP

param.kappa   = 0.0275;

param.GAMMA   = 0.1454;                               % degree of indexation to previous period's inflation


param.phi     = 47.847;
param.phi_x   = 0.01;       % sales-to-stock ratio adjustment costs


param.rho_w   = 0.8116;
param.eps_w   = 1.0045;      % real wage rigidity
param.d       = 1-0.0806;  % degree of partial indexation of the nominal wage


param.rho_B    = 0.6323;
param.gamma_pi = 0;
param.gamma_T  = 0;

param.rho_R   = 0.8307; 
param.sigma_R = 0.000625;                         % standatd deviation of MP shock
param.phi_pi  = 1.9853;                              % response to inflation
param.phi_y   = 0.4511;                              % response to output gap
param.phi_u   = 0.4511;                              % response to output gap
param.rho_m   = 0;                                % persistence of MP shock

param.R_ZLB   = 1;

param.alpha_lk = 0;


%  9) Capital quality shock


param.nu_cp    = 1.5;    % GK (2011): moderate intervention


param.phi_A_1 = 0;
param.phi_A_2 = 0;
param.phi_B   = 0;

param.tau_FG = 0;

param.rho_MRS_ent = 0;

%% QE


param.rho_passive_QE = 0.7;
param.rho_X_QE       = 0.7;

load('QE_param.mat');

param.phi_pi_QE = 10*param.phi_pi;
param.phi_u_QE  = 10*param.phi_y;

param.rho_A_g 	= 0.7;


param.rho_psi_cp = 0;

param.sigma2 = 1.5;

param.rho_Z       = 0.9;
param.rho_gamma_Q = 0.9;
param.rho_PSI_RP  = 0.9;
param.rho_eta     = 0.9;
param.rho_D       = 0.9;
param.rho_G       = 0.9;
param.rho_LT      = 0.9;
param.rho_TAU     = 0.9;
param.rho_varphi  = 0.9;
param.rho_iota    = 0.9;
param.rho_MP      = 0.8;
param.rho_BB      = 0.99995;
param.rho_PSI_W   = 0.8;
param.rho_MP      = 0;

param.rho_B_F     = 0.9;



param.delta_1 = SS_stats.r_k/param.delta_ss*param.v;
param.delta_0 = param.delta_ss/param.v^param.delta_1;

param.alpha_ll = 0;
param.alpha_kk = 0;
