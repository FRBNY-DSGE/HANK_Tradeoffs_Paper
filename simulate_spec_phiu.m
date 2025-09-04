%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% simulate_spec.m
%
% Description: This file is used to run simulations with a specified set of
% shocks for a user-defined set of phi_pi values
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all;
close all;
clc;

addpath('auxfiles/')
addpath('helper_functions/')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% (1) Simulation Settings:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


spec_settings.policy_rule = 'Taylor';
spec_settings.simul_ELB = true; %Use the ELB?
spec_settings.simul_disagg = false; %Don't change unless needed
spec_settings.use_covidGR = true;

% Kappa Settings
spec_settings.high_kappa = false;
spec_settings.high_kappa_val = 0.1; %Don't change unless needed


% Load in file with saved shock series

if spec_settings.use_covidGR
    load('shock_series_combined_hk_5.mat')
    % Shock label: Shock identifier in simulation results file
    spec_settings.shock_label = "shock_series_combined_hk_5";
    spec_settings.shock_series = shock_series_combined_hk_5(2:end,:)';
else
    load('shock_series_nocovidGR.mat')
    % Shock label: Shock identifier in simulation results file
    spec_settings.shock_label = "shock_series_nocovidGR";
    spec_settings.shock_series = shock_series_rand_new(2:end,:)';
end
spec_settings.num_simul_period = size(spec_settings.shock_series,1);




% Parameter estimation string
spec_settings.param_input = 'baseline';

% Identifier for the simulation
spec_settings.simul_string = strcat(spec_settings.param_input,'_pol=', spec_settings.policy_rule, '_shocks=', spec_settings.shock_label, '_ELB=', string(spec_settings.simul_ELB));
spec_settings.sim_save_folder = 'paper_simulations_replication/';

% Parallelization
pool_size = 100; % Number of parallel workers (cores)
%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 2. Parameter Loading under Specified Policy Rule
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Determine which values to loop over
phi_pi_vals = [1.5, 1.8533, 2, 2.5, 3, 3.5, 4];



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% (2.A) Parallelization:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Begin outer loop over phi_pi values:
parpool(pool_size);
for phi_u_tmp = [0.2, 0.36, 0.5]
    spec_settings.simul_string = strcat(spec_settings.param_input,'_pol=', spec_settings.policy_rule, '_shocks=', spec_settings.shock_label, '_ELB=', string(spec_settings.simul_ELB), '_phiu=',strrep(phi_u_tmp, '.', '_'));
    parfor phi_pi_ind = 1:length(phi_pi_vals)
        run_simul(spec_settings, phi_pi_vals, phi_pi_ind, phi_u_tmp);
    end
end



%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% (2.A.1) Simulation Loading
%%%%%%%%%%

function run_simul(spec_settings, phi_pi_vals, phi_pi_ind, phi_u_val)

    grid = 0;

    if strcmp(spec_settings.param_input, 'baseline')
        load('MCMC_posterior_mode.mat')
    else
        warn('Invalid parameters')
        assert(false);
    end

    INIT = xhat;

    % Set kappa based on spec settings:
    if spec_settings.high_kappa == true
        param_update.kappa = spec_settings.high_kappa_val;
    else
        param_update.kappa     = INIT(1 );
    end

    param_update.rho_w     = INIT(2 );
    param_update.d         = 1-INIT(3);
    param_update.phi       = INIT(4 )*10;

    % Assign phi_pi according to loop
    param_update.phi_pi    = phi_pi_vals(phi_pi_ind);

    param_update.phi_u     = phi_u_val;
    param_update.rho_BB    = INIT(7 );
    param_update.rho_B     = INIT(8 );
    param_update.rho_D     = INIT(9 ) ;
    param_update.rho_G     = INIT(10) ;
    param_update.rho_R     = INIT(11) ;
    param_update.rho_Z     = INIT(12) ;
    param_update.rho_eta   = INIT(13) ;
    param_update.rho_iota  = INIT(14) ;
    param_update.sig_D     = INIT(15)/100 ;
    param_update.sig_G     = INIT(16)/100 ;
    param_update.sig_R     = INIT(17)/100 ;
    param_update.sig_Z     = INIT(18)/100 ;
    param_update.sig_eta   = INIT(19)/100 ;
    param_update.sig_iota  = INIT(20)/100 ;
    param_update.sig_w     = INIT(21)/100 ;
    param_update.sig_BB    = INIT(22)/100 ;
    param_update.GAMMA     = INIT(23);
    param_update.sig_B_F   = INIT(24)/100;
    param_update.rho_B_F   = INIT(25);
    param_update.iota      = INIT(26);

    param_update.eps_w     = 1;
    param_update.phi_x     = 0;
    param_update.rho_PSI_W = 0;
    param_update.varphi    = 0;
    param_update.alpha_lk  = 0;
    param_update.theta_b   = 0.97;

    param_update.rho_X_QE   = 0.95;
    param_update.phi_pi_QE  = 0.5;
    param_update.phi_u_QE   = 10;

    param_update.gamma_B   = 0;
    param_update.gamma_pi  = 0;
    param_update.gamma_Y   = 0;

    param_update.rho_MP    = 0;
    param_update.gamma_pi  = 0;
    param_update.gamma_T   = 0;

    param_update.rho_PSI_RP = 0;
    param_update.sig_RP     = 0.1/100;
    param_update.Calvo      = 0.7;

    param.QE     = 'QE';
    param.regime = 'Taylor';

    update_ss_v5;

    parameters_agg_EST_v3;

    state_reduc_analysis;

    indexes_analysis;

    param.rho_pgap = 0.93;
    param.rho_ygap = 0.0;

    param.phi_pgap = param.phi_pi;


    param.phi_ygap = param.phi_u;
    param.rho_R_AIT = param.rho_R;

    % Load jacobian base depending on whether simulation is run with ELB

    load('Jacob_base_Taylor_G');
    Jacob_base = Jacob_base_Taylor_G;
    out_Jacob   = compute_Jacob_base_Taylor_v1(param,grid,SS_stats,idx,Xss,Yss);

    load('Jacob_base_ELB_G');
    Jacob_ZLB = Jacob_base_ELB_G;
    out_Jacob_ELB      = compute_Jacob_GR_Taylor_v1(param,grid,SS_stats,idx,Xss,Yss);


    [hx, gx, F1_aux, F2_aux, F3_aux, F4_aux, param, indicator_1] = SGU_EST_ref_v2(param,grid,Jacob_base,idx,out_Jacob);


    [P_ref,Q_ref,cof_ref,J_ref] = compute_lin_sol(hx, gx, F1_aux, F2_aux, F3_aux, F4_aux, grid);

    F_ZLB = @(a,b,c,d)F_sys_ref_anal_ELB(a,b,c,d,Xss,Yss,Gamma_state,Gamma_control,InvGamma,param,grid,SS_stats,Copula,P_SE);

    Fss_ZLB_old = Fss_ZLB(:);
    [Fss_ZLB,LHS_ZLB,RHS_ZLB,Distr_ZLB] = F_ZLB(State,State_m,Contr,Contr_m);

    [F1_ZLB_aux,F2_ZLB_aux,F3_ZLB_aux,F4_ZLB_aux,param]     = SGU_EST_ZLB_v2(param,grid,Jacob_ZLB,idx,out_Jacob_ELB);

    F1_ZLB = [F1_ZLB_aux(1:grid.numstates_endo,:);F1_ZLB_aux(grid.numstates+1:end,:)];
    F2_ZLB = [F2_ZLB_aux(1:grid.numstates_endo,:);F2_ZLB_aux(grid.numstates+1:end,:)];
    F3_ZLB = [F3_ZLB_aux(1:grid.numstates_endo,:);F3_ZLB_aux(grid.numstates+1:end,:)];
    F4_ZLB = [F4_ZLB_aux(1:grid.numstates_endo,:);F4_ZLB_aux(grid.numstates+1:end,:)];

    F11_ZLB = F1_ZLB(:,1:grid.numstates_endo);
    F12_ZLB = F1_ZLB(:,grid.numstates_endo+1:end);

    F31_ZLB = F3_ZLB(:,1:grid.numstates_endo);
    F32_ZLB = F3_ZLB(:,grid.numstates_endo+1:end);

    % Construct a coefficient matrix

    A_ZLB = [zeros(grid.numstates_endo+grid.numcontrols,grid.numstates_endo), F2_ZLB];
    B_ZLB = [F11_ZLB, F4_ZLB];
    C_ZLB = [F31_ZLB, zeros(grid.numstates_endo+grid.numcontrols,grid.numcontrols)];
    E_ZLB = F32_ZLB;
    D_ZLB = [Fss_ZLB(1:grid.numstates_endo,:);Fss_ZLB(grid.numstates+1:end,:)];

    cof_ZLB = [C_ZLB,B_ZLB,A_ZLB];
    J_ZLB   = E_ZLB;

    load('Jacob_base_alt_G_EST');

    Jacob_base_alt = Jacob_base_alt_G_EST;

    load('data_est_covid');

    PD_filter = @(x)post_density_covid_filter_v10(x,prior_info,data_est_covid,ZLB_duration_1,ZLB_duration_2,Jacob_base_alt,Jacob_base_alt,Jacob_base_alt,SS_stats,param,grid,mu_dist,Value,mutil_c,Va,Copula,Fss_ZLB_old,H_aux);

    est_shock_series = PD_filter(xhat);

    % Run simulation file
    simul_main_ELB_v1;

end
