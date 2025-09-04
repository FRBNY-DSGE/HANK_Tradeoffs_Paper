
	CEs_shocks        = zeros(15,9);
	CEs_RR_shocks     = zeros(15,9);
	CEs_PROFIT_shocks = zeros(15,9);
	CEs_W_shocks      = zeros(15,9);
	CEs_f_shocks      = zeros(15,9);
	CEs_GOV_shocks    = zeros(15,9);
	CEs_Q_shocks      = zeros(15,9);

	%fig_names = {'FC','RP','Z','LT','LP','MP','IT','PM','WM'};
    
    %for lllllll = 6:6
        
        close all;
        
 load('MCMC_posterior_precovid_ugap_v4.mat');
 xhat = out.post_mode;
clear title

        
        
        maxiter = 10;
        tol     = 1e-8;
        
        param.adjust     = 'G';
        




        
        parameters_agg;
        
 
        param_update.kappa     = xhat(1 );
        param_update.rho_w     = xhat(2 );

        param_update.d         = 1-xhat(3);

        param_update.eps_w     = 1;
        param_update.phi       = xhat(4 )*10;
        param_update.phi_x     = 0;
        param_update.phi_pi    = xhat(5 );

        param_update.phi_u     = xhat(6 );
        param_update.rho_BB    = xhat(7 );
        param_update.rho_B     = xhat(8 );
        param_update.rho_D     = xhat(9 ) ;
        param_update.rho_G     = xhat(10) ;
        param_update.rho_R     = xhat(11) ;
        param_update.rho_Z     = xhat(12) ;
        param_update.rho_eta   = xhat(13) ;
        param_update.rho_iota  = xhat(14) ;

        param_update.rho_PSI_W = 0 ;
        param_update.sig_D     = xhat(15)/100 ;
        param_update.sig_G     = xhat(16)/100 ;
        param_update.sig_R     = xhat(17)/100 ;
        param_update.sig_Z     = xhat(18)/100 ;
        param_update.sig_eta   = xhat(19)/100 ;
        param_update.sig_iota  = xhat(20)/100 ;
        param_update.sig_w     = xhat(21)/100 ;
        param_update.sig_BB    = xhat(22)/100 ;

        param_update.iota      = xhat(26);

        param_update.varphi    = 0;

        param_update.alpha_lk  = 0;
        param_update.theta_b   = 0.97;

        param_update.GAMMA     = xhat(23);
        param_update.rho_B_F   = xhat(25);

        param_update.sig_B_F   = xhat(24)/100;

        
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
  

        fprintf('Phi pi = %.2f ... \n', param_update.phi_pi);
        
        
        update_ss_v5;
        
        parameters_agg_EST_v3;



        param.rho_pgap = 0.93; 
        param.rho_ygap = 0.93; 

        param.phi_pgap = param.phi_pi;
        param.phi_ygap = param.phi_u;

        param.rho_R_AIT = param.rho_R;
        

        shock_sets = [param.sig_B_F;-param.sig_BB;-param.sig_Z;param.sig_G;param.sig_D;param.sig_R;-param.sig_iota;param.sig_eta;param.sig_w];

        %state_reduc_anal_AIT;
        state_reduc_analysis;
        
        param.alpha_lk = 0;
        
        param.FEAR_shock = 'both';  % 'both','unemp','fin'
        
        param.QE = 'No QE';
        
        param.regime = 'Taylor';
        
        F_ref = @(a,b,c,d)F_sys_ref_anal_Taylor(a,b,c,d,Xss,Yss,Gamma_state,Gamma_control,InvGamma,param,grid,SS_stats,Copula,P_SE);

        [Fss_ref,LHS_ref,RHS_ref,Distr_ref] = F_ref(State,State_m,Contr,Contr_m);
        
        disp(max(abs(Fss_ref)));
        
        %% Solve RE via Schmitt-Grohe-Uribe Form
        %  Initialize parallel pool
        
        p = gcp('nocreate');
        pool='local'; % select machine
        
        if isempty(p)
            c = parcluster;
            parpool(pool,c.NumWorkers-1)
            p = gcp('nocreate');
        end
        

        %indexes_anal_AIT;
        indexes_analysis;

        out_Jacob = compute_Jacob_base_Taylor_v1(param,grid,SS_stats,idx,Xss,Yss);
       

        if spec_settings.fix_hh_bins
            load('Jacob_base_Taylor_G_v2')
            Jacob_base = Jacob_base_Taylor_G_v2;
        else
            load('Jacob_base_Taylor_G')
            Jacob_base = Jacob_base_Taylor_G;
        end


        param.overrideEigen = true;  % Warning appears, but critical Eigenvalue shifted
        tic;


        [hx,gx,F1_aux,F2_aux,F3_aux,F4_aux,param] = SGU_EST2(F_ref,param,grid,Jacob_base,idx,out_Jacob,p);
        toc;
        % Compute matrices for the linearized equilbrium system
        
        F1 = [F1_aux(1:grid.numstates_endo,:);F1_aux(grid.numstates+1:end,:)];
        F2 = [F2_aux(1:grid.numstates_endo,:);F2_aux(grid.numstates+1:end,:)];
        F3 = [F3_aux(1:grid.numstates_endo,:);F3_aux(grid.numstates+1:end,:)];
        F4 = [F4_aux(1:grid.numstates_endo,:);F4_aux(grid.numstates+1:end,:)];
        
        F11 = F1(:,1:grid.numstates_endo);
        F12 = F1(:,grid.numstates_endo+1:end);
        
        F31 = F3(:,1:grid.numstates_endo);
        F32 = F3(:,grid.numstates_endo+1:end);
        
        % Construct a coefficient matrix
        
        A_ref = [zeros(grid.numstates_endo+grid.numcontrols,grid.numstates_endo), F2];
        B_ref = [F11, F4];
        C_ref = [F31, zeros(grid.numstates_endo+grid.numcontrols,grid.numcontrols)];
        E_ref = F32;
        
        % Constuct P, Q
        

        hx_aux = hx(1:grid.numstates_endo,:);
        H_aux  = [hx_aux;gx];
        H1_aux = H_aux(:,1:grid.numstates_endo);
        H2_aux = H_aux(:,end-grid.numstates_shocks+1:end);
       

        P_ref  = [H1_aux,zeros(grid.numstates_endo+grid.numcontrols,grid.numcontrols)];
        Q_ref  = H2_aux;
        
        cof_ref = [C_ref,B_ref,A_ref];
        J_ref   = E_ref;

        %IRFs_shocks_anal_AIT;
        IRFs_shocks_analysis;

        IRFs_taylor = IRFs;

        if spec_settings.calc_consumption_equivalence

            CEs_shocks_anal;
    
            CEs_shocks(:,lllllll)        = CE_total';
            CEs_PROFIT_shocks(:,lllllll) = CE_PROFIT_total'; %%Profit
            CEs_Q_shocks(:,lllllll)      = CE_Q_total'; %% Equity price
            CEs_W_shocks(:,lllllll)      = CE_W_total'; %%Wage
            CEs_f_shocks(:,lllllll)      = CE_f_total'; %%job finding rate
            CEs_RR_shocks(:,lllllll)     = CE_RR_total'; %%Real rate
            CEs_GOV_shocks(:,lllllll)    = CE_GOV_total'; %%Gov transfer
    
            nperiods_plot = 20;
    
            %Save CE results
            CEs_shocks_anal_collect.CEs_collect_shocks = CEs_shocks;
            CEs_shocks_anal_collect.CEs_collect_PROFIT_shocks = CEs_PROFIT_shocks;
            CEs_shocks_anal_collect.CEs_collect_Q_shocks = CEs_Q_shocks;
            CEs_shocks_anal_collect.CEs_collect_W_shocks = CEs_W_shocks;
            CEs_shocks_anal_collect.CEs_collect_f_shocks = CEs_f_shocks;
            CEs_shocks_anal_collect.CEs_collect_RR_shocks = CEs_RR_shocks;
            CEs_shocks_anal_collect.CEs_collect_GOV_shocks = CEs_GOV_shocks;
    
            CEs_baseline = CEs_shocks_anal_collect;
            com_CE;
            CE_PI_bar_base = CE_PI_bar;
            CE_bar_base = CE_bar;
    
            CEs_shocks;

        end





        %%%%% Repeat for high Phi %%%%%%%
if spec_settings.calculate_high_phi_pi

        % Use phi_pi = 3 
         param_update.phi_pi = spec_settings.high_phi_pi; 

        % If high kappa activated
        if spec_settings.calculate_high_kappa
            param_update.kappa = spec_settings.high_kappa;
            param.kappa = 0.1;
            fprintf('kappa = %.2f ... \n', param_update.kappa);
        else
            param_update.kappa = xhat(1);
            param.kappa = xhat(1);
            fprintf('kappa = %.2f ... \n', param.kappa);
        end
         
        update_ss_v5;
        
        parameters_agg_EST_v3;
        param.FEAR_shock = 'both';  % 'both','unemp','fin'
        param.QE = 'No QE';
        param.regime = 'Taylor';
        

        F_ref = @(a,b,c,d)F_sys_ref_anal_Taylor(a,b,c,d,Xss,Yss,Gamma_state,Gamma_control,InvGamma,param,grid,SS_stats,Copula,P_SE);
        
        %F_ref = @(a,b,c,d)F_sys_ref_anal_Taylor(a,b,c,d,Xss,Yss,Gamma_state,Gamma_control,InvGamma,param,grid,SS_stats,Copula,P_SE);
        [Fss_ref,LHS_ref,RHS_ref,Distr_ref] = F_ref(State,State_m,Contr,Contr_m);
        disp(max(abs(Fss_ref)));
        
        %% Solve RE via Schmitt-Grohe-Uribe Form
        %  Initialize parallel pool
        
        p = gcp('nocreate');
        pool='local'; % select machine
        
        if isempty(p)
            c = parcluster;
            parpool(pool,c.NumWorkers-1)
            p = gcp('nocreate');
        end
        

        %indexes_anal_AIT;
        indexes_analysis;
        
        if spec_settings.fix_hh_bins
            disp('Using Jacob_base_Taylor_v2')
            out_Jacob = compute_Jacob_base_Taylor_v2(param,grid,SS_stats,idx,Xss,Yss);
        else
            out_Jacob = compute_Jacob_base_Taylor_v1(param,grid,SS_stats,idx,Xss,Yss);
        end
        
        if spec_settings.fix_hh_bins
            load('Jacob_base_Taylor_G_v2')
            Jacob_base = Jacob_base_Taylor_G_v2;
        else
            load('Jacob_base_Taylor_G')
            Jacob_base = Jacob_base_Taylor_G;
        end

        param.overrideEigen = true;  % Warning appears, but critical Eigenvalue shifted
        tic;
        [hx,gx,F1_aux,F2_aux,F3_aux,F4_aux,param] = SGU_EST2(F_ref,param,grid,Jacob_base,idx,out_Jacob,p);
        toc;
        % Compute matrices for the linearized equilbrium system
        
        F1 = [F1_aux(1:grid.numstates_endo,:);F1_aux(grid.numstates+1:end,:)];
        F2 = [F2_aux(1:grid.numstates_endo,:);F2_aux(grid.numstates+1:end,:)];
        F3 = [F3_aux(1:grid.numstates_endo,:);F3_aux(grid.numstates+1:end,:)];
        F4 = [F4_aux(1:grid.numstates_endo,:);F4_aux(grid.numstates+1:end,:)];
        
        F11 = F1(:,1:grid.numstates_endo);
        F12 = F1(:,grid.numstates_endo+1:end);
        
        F31 = F3(:,1:grid.numstates_endo);
        F32 = F3(:,grid.numstates_endo+1:end);
        
        % Construct a coefficient matrix
        
        A_ref = [zeros(grid.numstates_endo+grid.numcontrols,grid.numstates_endo), F2];
        B_ref = [F11, F4];
        C_ref = [F31, zeros(grid.numstates_endo+grid.numcontrols,grid.numcontrols)];
        E_ref = F32;
        
        % Constuct P, Q
        
        hx_aux = hx(1:grid.numstates_endo,:);
        H_aux  = [hx_aux;gx];
        H1_aux = H_aux(:,1:grid.numstates_endo);
        H2_aux = H_aux(:,end-grid.numstates_shocks+1:end);
        
        P_ref  = [H1_aux,zeros(grid.numstates_endo+grid.numcontrols,grid.numcontrols)];
        Q_ref  = H2_aux;
        
        cof_ref = [C_ref,B_ref,A_ref];
        J_ref   = E_ref;

        IRFs_shocks_analysis;

        IRFs_high_phi = IRFs;


        if spec_settings.calc_consumption_equivalence
            CEs_shocks_anal;
    
            CEs_collect_shocks(:,lllllll)        = CE_total';
            CEs_collect_PROFIT_shocks(:,lllllll) = CE_PROFIT_total';
            CEs_collect_Q_shocks(:,lllllll)      = CE_Q_total';
            CEs_collect_W_shocks(:,lllllll)      = CE_W_total';
            CEs_collect_f_shocks(:,lllllll)      = CE_f_total';
            CEs_collect_RR_shocks(:,lllllll)     = CE_RR_total';
            CEs_collect_GOV_shocks(:,lllllll)    = CE_GOV_total';
    
    
            CEs_shocks_anal_collect.CEs_collect_shocks = CEs_collect_shocks;
            CEs_shocks_anal_collect.CEs_collect_PROFIT_shocks = CEs_collect_PROFIT_shocks;
            CEs_shocks_anal_collect.CEs_collect_Q_shocks = CEs_collect_Q_shocks;
            CEs_shocks_anal_collect.CEs_collect_W_shocks = CEs_collect_W_shocks;
            CEs_shocks_anal_collect.CEs_collect_f_shocks = CEs_collect_f_shocks;
            CEs_shocks_anal_collect.CEs_collect_RR_shocks = CEs_collect_RR_shocks;
            CEs_shocks_anal_collect.CEs_collect_GOV_shocks = CEs_collect_GOV_shocks;

            com_CE;
            CE_PI_bar_compare = CE_PI_bar;
            CE_bar_compare = CE_bar;

        end

elseif spec_settings.calculate_high_kappa


    param_update.kappa = spec_settings.high_kappa_val;
    fprintf('kappa = %.2f ... \n', param_update.kappa);

    update_ss_v5;

    parameters_agg_EST_v3;
    param.FEAR_shock = 'both';  % 'both','unemp','fin'
    param.QE = 'No QE';
    param.regime = 'Taylor';


    F_ref = @(a,b,c,d)F_sys_ref_anal_Taylor(a,b,c,d,Xss,Yss,Gamma_state,Gamma_control,InvGamma,param,grid,SS_stats,Copula,P_SE);

    [Fss_ref,LHS_ref,RHS_ref,Distr_ref] = F_ref(State,State_m,Contr,Contr_m);
    disp(max(abs(Fss_ref)));

    %%% Solve RE via Schmitt-Grohe-Uribe Form
    %  Initialize parallel pool

    p = gcp('nocreate');
    pool='local'; % select machine

    if isempty(p)
        c = parcluster;
        parpool(pool,c.NumWorkers-1)
        p = gcp('nocreate');
    end


    indexes_analysis;
    out_Jacob = compute_Jacob_base_Taylor_v1(param,grid,SS_stats,idx,Xss,Yss);


    if spec_settings.fix_hh_bins
        load('Jacob_base_Taylor_G_v2')
        Jacob_base = Jacob_base_Taylor_G_v2;
    else
        load('Jacob_base_Taylor_G')
        Jacob_base = Jacob_base_Taylor_G;
    end


    %
    param.overrideEigen = true;  % Warning appears, but critical Eigenvalue shifted
    tic;
    [hx,gx,F1_aux,F2_aux,F3_aux,F4_aux,param] = SGU_EST2(F_ref,param,grid,Jacob_base,idx,out_Jacob,p);
    toc;
    % Compute matrices for the linearized equilbrium system

    F1 = [F1_aux(1:grid.numstates_endo,:);F1_aux(grid.numstates+1:end,:)];
    F2 = [F2_aux(1:grid.numstates_endo,:);F2_aux(grid.numstates+1:end,:)];
    F3 = [F3_aux(1:grid.numstates_endo,:);F3_aux(grid.numstates+1:end,:)];
    F4 = [F4_aux(1:grid.numstates_endo,:);F4_aux(grid.numstates+1:end,:)];

    F11 = F1(:,1:grid.numstates_endo);
    F12 = F1(:,grid.numstates_endo+1:end);

    F31 = F3(:,1:grid.numstates_endo);
    F32 = F3(:,grid.numstates_endo+1:end);

    % Construct a coefficient matrix

    A_ref = [zeros(grid.numstates_endo+grid.numcontrols,grid.numstates_endo), F2];
    B_ref = [F11, F4];
    C_ref = [F31, zeros(grid.numstates_endo+grid.numcontrols,grid.numcontrols)];
    E_ref = F32;

    % Constuct P, Q

    hx_aux = hx(1:grid.numstates_endo,:);
    H_aux  = [hx_aux;gx];
    H1_aux = H_aux(:,1:grid.numstates_endo);
    H2_aux = H_aux(:,end-grid.numstates_shocks+1:end);

    P_ref  = [H1_aux,zeros(grid.numstates_endo+grid.numcontrols,grid.numcontrols)];
    Q_ref  = H2_aux;

    cof_ref = [C_ref,B_ref,A_ref];
    J_ref   = E_ref;

    IRFs_shocks_analysis;
    IRFs_high_phi = IRFs;


    if spec_settings.calc_consumption_equivalence
        CEs_shocks_anal;

        CEs_collect_shocks(:,lllllll)        = CE_total';
        CEs_collect_PROFIT_shocks(:,lllllll) = CE_PROFIT_total';
        CEs_collect_Q_shocks(:,lllllll)      = CE_Q_total';
        CEs_collect_W_shocks(:,lllllll)      = CE_W_total';
        CEs_collect_f_shocks(:,lllllll)      = CE_f_total';
        CEs_collect_RR_shocks(:,lllllll)     = CE_RR_total';
        CEs_collect_GOV_shocks(:,lllllll)    = CE_GOV_total';

        %plot_CEs_GR_mainQE_v4;
        CEs_shocks_anal_collect.CEs_collect_shocks = CEs_collect_shocks;
        CEs_shocks_anal_collect.CEs_collect_PROFIT_shocks = CEs_collect_PROFIT_shocks;
        CEs_shocks_anal_collect.CEs_collect_Q_shocks = CEs_collect_Q_shocks;
        CEs_shocks_anal_collect.CEs_collect_W_shocks = CEs_collect_W_shocks;
        CEs_shocks_anal_collect.CEs_collect_f_shocks = CEs_collect_f_shocks;
        CEs_shocks_anal_collect.CEs_collect_RR_shocks = CEs_collect_RR_shocks;
        CEs_shocks_anal_collect.CEs_collect_GOV_shocks = CEs_collect_GOV_shocks;

        com_CE;
        CE_PI_bar_compare = CE_PI_bar;
        CE_bar_compare = CE_bar;

    end


end


