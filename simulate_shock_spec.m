%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% simulate_shock_spec.m
%
% Description: Code for generating shock series similar to
% shock_series_combined.mat and shock_series_nocovidGR.mat
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all;
close all;
clc;

addpath('auxfiles/')
addpath('helper_functions/')
addpath("shockfiles") % Load in shock trimming helper functions
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% (1) General Settings:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Toggle for generating shock_series_combined_hk_5 (true) or shock_series_noCovidGR
% (false)
generate_CovidGR_shocks = true;

% Filepath to store shock series
shock_name = "shock_series_combined_custom"; % Change name here
shock_folder = "shocks_replication/" + shock_name + "/";

% Trim settings
kappa_val = 0.1;
phi_pi_seq = [4, 1.5, 4];
ELB_threshold = 101;
num_simul_thresholds = 60000; % Number of simulation periods with at least "ELB_threshold" number of ELB periods
sh_series_periods = 600000; % Initial number of shock series periods to trim from 

% Parallelization (# of cores)
pool_size = 100;

%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% (2) Run simulation:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%
% (3.1) Run in-simulation trimming where we augment shock series with
% opposite-signed shocks
%%%%%%%%%%%%%%%%%%%%

parpool(pool_size);

for trim_num = 1:length(phi_pi_seq)
    % Loads parameter values into workspace

    load('mode_find_results_precovid_ugap_v4.mat')
    load('precovid_ugap_v4_post_mode.mat')

    % Adds filtered shocks, jacobian bases, transition matricies, steady state
    % variables

    INIT = precovid_ugap_v4_post_mode;

    param_update.kappa     = INIT(1 );
    param_update.rho_w     = INIT(2 );
    param_update.d         = 1-INIT(3);
    param_update.phi       = INIT(4 )*10;
    param_update.phi_pi    = INIT(5 );
    param_update.phi_u     = INIT(6 );
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

    param_update.phi_pi = 1.5;

    param.QE     = 'QE';
    param.regime = 'Taylor';

    update_ss_v5;

    parameters_agg_EST_v3;

    state_reduc_analysis;

    indexes_analysis;

    param.rho_pgap = 0.93;
    param.rho_ygap = 0.0;

    param.phi_pi = phi_pi_seq(trim_num);

    param.kappa = kappa_val;

    param.phi_pgap = 1 * (1/(1-param.rho_pgap));
 
    param.phi_ygap = param.phi_u;

    param.rho_R_AIT = param.rho_R;

    % Load jacob bases into workspace
    load('Jacob_base_Taylor_G');
    Jacob_base_Taylor = Jacob_base_Taylor_G;

    load('Jacob_base_ELB_G');
    Jacob_ZLB = Jacob_base_ELB_G;

    load('Jacob_base_alt_G_EST');
    Jacob_base = Jacob_base_alt_G_EST;

    out_Jacob_Taylor   = compute_Jacob_base_Taylor_v1(param,grid,SS_stats,idx,Xss,Yss);

    [hx    ,gx    ,F1_aux    ,F2_aux    ,F3_aux    ,F4_aux    ,param,indicator_1] = SGU_EST_ref_v2(param,grid,Jacob_base_Taylor,idx,out_Jacob_Taylor);
    [P_ref,Q_ref,cof_ref,J_ref] = compute_lin_sol(hx    ,gx    ,F1_aux    ,F2_aux    ,F3_aux    ,F4_aux    ,grid);

    F_ZLB = @(a,b,c,d)F_sys_ref_anal_ELB(a,b,c,d,Xss,Yss,Gamma_state,Gamma_control,InvGamma,param,grid,SS_stats,Copula,P_SE);
    Fss_ZLB_old = Fss_ZLB(:);
    [Fss_ZLB,LHS_ZLB,RHS_ZLB,Distr_ZLB] = F_ZLB(State,State_m,Contr,Contr_m);

    % Load in ZLB jacobian 
    load("out_Jacob_ZLB_paper.mat")

    [F1_ZLB_aux,F2_ZLB_aux,F3_ZLB_aux,F4_ZLB_aux,param] = SGU_EST_ZLB_v2(param,grid,Jacob_ZLB,idx,out_Jacob_ELB);

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

    load('data_est_covid');

    % Read in filtered covid and great recession shocks
    load("covid_shocks.mat")
    load("gr_shocks.mat")

    SIGMA = eye(grid.numstates_shocks-1);
    SIGMA(1 ,1 ) = param.sig_B_F.^2;
    SIGMA(2 ,2 ) = param.sig_BB.^2;
    SIGMA(3 ,3 ) = param.sig_Z.^2;
    SIGMA(4 ,4 ) = param.sig_G.^2;
    SIGMA(5 ,5 ) = param.sig_D.^2;
    SIGMA(6 ,6 ) = param.sig_R.^2;
    SIGMA(7 ,7 ) = param.sig_iota.^2;
    SIGMA(8 ,8 ) = param.sig_eta.^2;
    SIGMA(9 ,9 ) = param.sig_w.^2;

    pull_rec_shocks = -1;

    covid_shocks(1,:) = 0;
    gr_shocks(1,:) = 0;

    covid_shocks(7,:) = 0;
    gr_shocks(7,:) = 0;
%%

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % (3.2) Create or load in previous shock series depending on
    % trimming step

    if generate_CovidGR_shocks
        % Trimming Step 1: Generate main (original) shock series (under
        % high phi_pi)
        if trim_num == 1
            num_simul_thresholds = 60000;
            shock_series = mvnrnd(zeros(1, grid.numstates_shocks-1),SIGMA,sh_series_periods);

            % Trimming Step 2: Trim main shock series under low phi_pi
        elseif trim_num == 2
            url = shock_folder + "trim1=high_phipi.mat";
            load(url)
            shock_series = shock_series_trim_1(2:end, :)';

            % Trimming Step 3: Generate sign-reversed shock series trimmed
            % under high phi_pi
        else
            num_simul_thresholds = 30000;
            url = shock_folder + "trim2=low_phipi.mat";
            load(url)
            shock_series = -shock_series_trim_2(2:end, :)';
        end

    else
        % Trimming Step 1: Generate main shock series (under high phi_pi)
        if trim_num == 1
            num_simul_thresholds = 60000;
            shock_series = mvnrnd(zeros(1, grid.numstates_shocks-1),SIGMA,sh_series_periods);

        % Trimming Step 2: Trim main shock series under low phi_pi
        elseif trim_num == 2
     
            url = shock_folder + "trim1=noCovidGR_high_phipi.mat"
            load(url)
            shock_series = shock_series_trim_1(2:end, :)';

        % Trimming Step 3: Generate sign-reversed shock series trimmed
        % under high phi_pi
        else
            num_simul_thresholds = 30000;
            url = shock_folder + "trim2=noCovidGR_low_phipi.mat"
            load(url)
            shock_series = -shock_series_trim_2(2:end, :)';
        end
    end

    sh_series_size = size(shock_series, 1);

    refresh_size = 300;

    num_set_periods = floor(sh_series_size/refresh_size);

    % Pre-defined struct to store blocks that don't result in explosive
    % behavior
    X_vec_aux = struct;

    X_vec_aux(num_set_periods).X_Taylor_series = [];
    X_vec_aux(num_set_periods).ELB_Taylor_hits = [];
    X_vec_aux(num_set_periods).ELB_Taylor_expected = [];

    Explosive_check = zeros(num_set_periods,1);

    simul_shock_series_aux = struct;
    simul_shock_series_raw_aux = struct;

    simul_shock_series_aux(num_set_periods).shock_series = [];
    simul_shock_series_raw_aux(num_set_periods).shock_series = [];

    num_sets_noexplosive = 1;

    GR_shock_frequency = 1/20;

    num_cooldown_periods = 6; % 3 years of cool down period

    scale       = 1;
    scale_first = 1;

    gr_shocks = scale*gr_shocks;

    gr_shocks(:,1) = scale_first*gr_shocks(:,1);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % (3.3) Run trimming scheme

    %%
    parfor kk = 1:num_set_periods
        output = run_simul(kk, ELB_threshold, num_simul_thresholds, shock_series, sh_series_size, ...
            refresh_size, num_set_periods, ...
            num_sets_noexplosive, GR_shock_frequency, ...
            num_cooldown_periods, scale, gr_shocks, covid_shocks, ...
            P_ref, Q_ref, cof_ref, J_ref, cof_ZLB, J_ZLB, D_ZLB, grid, ...
            param, irfshock, trim_num, idx, generate_CovidGR_shocks)

        X_vec_aux(kk).X_Taylor_series = output.X_Taylor_series
        X_vec_aux(kk).ELB_Taylor_hits = output.ELB_Taylor_hits
        X_vec_aux(kk).ELB_Taylor_expected = output.ELB_Taylor_expected
        simul_shock_series_aux(kk).shock_series = output.shock_series
        simul_shock_series_raw_aux(kk).shock_series = output.shock_series_raw
        Explosive_check(kk, 1) = output.Explosive_check
    end

    num_simul_noexplosive = num_set_periods-sum(Explosive_check);
    %%
    % For first two trims
    if trim_num ~= 3

        % Retain order of X_vec and shock series based non-explosive blocks
        num_set_periods_vec = 1:num_set_periods;
        non_exp_block_inds = num_set_periods_vec(logical(~Explosive_check)); % Find non-explosive blocks indices (in order)
        non_exp_order_iter = 1;

        for non_exp_ind = non_exp_block_inds % Re-index structs from lowest index (kk) to highest index for non-explosive blocks
            X_vec(non_exp_order_iter).X_Taylor_series = X_vec_aux(non_exp_ind).X_Taylor_series;
            X_vec(non_exp_order_iter).ELB_Taylor_hits = X_vec_aux(non_exp_ind).ELB_Taylor_hits;
            X_vec(non_exp_order_iter).ELB_Taylor_expected = X_vec_aux(non_exp_ind).ELB_Taylor_expected;
            simul_shock_series(non_exp_order_iter).shock_series =  simul_shock_series_aux(non_exp_ind).shock_series;
            simul_shock_series_raw(non_exp_order_iter).shock_series = simul_shock_series_raw_aux(non_exp_ind).shock_series;
            non_exp_order_iter = non_exp_order_iter + 1;
        end

        % For first trim, we only save the shock series
        if trim_num == 1

            find_ELB_threshold;
            X_vec_trim;

            shock_series_trim_1 = shock_series_trim;

            if generate_CovidGR_shocks
                url = shock_folder + "trim1=high_phipi";
                save(url, "shock_series_trim_1")
            else
                url = shock_folder + "trim1=noCovidGR_high_phipi";
                save(url, "shock_series_trim_1")
            end

        else
            %Turn 300-period blocks into continuous sequence
            for jj = 1:num_simul_noexplosive
               shock_series_trim(:,(jj-1)*refresh_size+1:jj*refresh_size) = simul_shock_series(jj).shock_series;
               shock_series_raw_trim(:,(jj-1)*refresh_size+1:jj*refresh_size) = simul_shock_series_raw(jj).shock_series;
            end

            % For second trim, we save shock series and results
            results_trim.X_vec_trim               = X_vec;
            results_trim.simul_shock_series_trim  = simul_shock_series_trim;
            results_trim.simul_shock_series_raw_trim = simul_shock_series_raw_trim;
            results_trim.shock_series_trim  = shock_series_trim;
            results_trim.shock_series_raw_trim = shock_series_raw_trim;
            results_trim.Explosive_check          = Explosive_check;

            shock_series_trim_2 = results_trim.shock_series_trim;

            if generate_CovidGR_shocks
                url1 = shock_folder + 'results_trim';
                url2 = shock_folder + 'trim2=low_phipi';
                save(url1,'results_trim');
                save(url2, 'shock_series_trim_2')
            else
                url1 = shock_folder + 'results_trim_noCovidGR';
                url2 = shock_folder + 'trim2=noCovidGR_low_phipi';
                save(url1, 'results_trim')
                save(url2, 'shock_series_trim_2')
            end

        end

        % Sign-reversed shock series (third trim)
    else
        % Save results
        results_trim_reverse.X_vec_trim               = X_vec_aux;
        results_trim_reverse.simul_shock_series_trim  = simul_shock_series_aux;
        results_trim_reverse.simul_shock_series_raw_trim = simul_shock_series_raw_aux;
        results_trim_reverse.Explosive_check          = Explosive_check;

        if generate_CovidGR_shocks
            url = shock_folder + 'results_trim_reverse';
            save(url, 'results_trim_reverse')
        else
            url = shock_folder + 'results_trim_reverse_noCovidGR'
            save(url, 'results_trim_reverse')
        end
    end
end

%%

% Augment main shock series with its sign-reversed counterpart
X_vec_finalize;

% Save shock series
shock_series = shock_series_final;
url = shock_folder + shock_name;

save(url, "shock_series")

%%

function out = run_simul(kk, ELB_threshold, num_simul_thresholds, shock_series, sh_series_size, ...
                refresh_size, num_set_periods,...
                num_sets_noexplosive, GR_shock_frequency, ...
                num_cooldown_periods, scale, gr_shocks, covid_shocks, ...
                P_ref, Q_ref, cof_ref, J_ref, cof_ZLB, J_ZLB, D_ZLB, grid, ...
                param, irfshock, trim_num, idx, generate_CovidGR_shocks)


            shock_block = shock_series((kk-1)*refresh_size+1:kk*refresh_size,:);

            block_size = refresh_size;

            X_Taylor = zeros(grid.num_endo, 1);

            ELB_Taylor_hits     = zeros(block_size,1);
            ELB_Taylor_expected = zeros(block_size,1);
            Taylor_ELB_converge = zeros(block_size,1);

            X_Taylor_series = zeros(grid.num_endo,block_size);

            EXPLOSIVE_cons  = false;
            EXPLOSIVE_unemp = false;
            EXPLOSIVE_phipi = 0;

            bb = 1;

            idx_GR = 0;

            shock_series_aux = zeros(grid.numstates_shocks,refresh_size);

            % Function output
            out = struct;
            out.X_Taylor_series = zeros(grid.num_endo,block_size);
            out.ELB_Taylor_hits = zeros(block_size,1);
            out.ELB_Taylor_expected = zeros(block_size,1);
            out.shock_series = zeros(grid.numstates_shocks,refresh_size);
            out.shock_series_raw = zeros(grid.numstates_shocks,refresh_size);
            out.Explosive_check = 0;

            while bb <= block_size && EXPLOSIVE_phipi == 0

                if generate_CovidGR_shocks && trim_num == 1
                    % First iteration under high phi_pi generates original shock series (main
                    % series)
                    if bb <= refresh_size - (6+num_cooldown_periods)+1

                        if idx_GR == 0

                            eps_t = shock_block(bb, :);
                            eps_t = [0; eps_t(:)];
                            eps_t(7) = 0;

                            xx_aux = unifrnd(0,1);

                            if xx_aux <= GR_shock_frequency

                                if num_cooldown_periods == 0

                                    xx_aux_2 = unifrnd(0,1);

                                    if xx_aux_2 <= 0.5

                                        new_gr_shocks = gr_shocks;

                                    else

                                        new_gr_shocks = covid_shocks;

                                    end


                                else

                                    xx_aux_2 = unifrnd(0,1);

                                    if xx_aux_2 <= 0.5


                                        gen_new_GR_shocks;

                                    else

                                        gen_new_covid_shocks;

                                        new_gr_shocks = new_covid_shocks;

                                    end

                                end

                                idx_GR = 1;

                                idx_GR_shocks = 1;

                                eps_t  = new_gr_shocks(:,idx_GR_shocks);

                                idx_GR_shocks = idx_GR_shocks + 1;

                            end

                        elseif idx_GR == 1

                            if idx_GR_shocks < 6+num_cooldown_periods

                                eps_t  = new_gr_shocks(:,idx_GR_shocks);

                                idx_GR_shocks = idx_GR_shocks + 1;

                            elseif idx_GR_shocks == 6+num_cooldown_periods

                                eps_t  = new_gr_shocks(:,idx_GR_shocks);
                                clear idx_GR_shocks
                                idx_GR = 0; % reset

                            end

                        end

                    elseif bb > refresh_size - (6+num_cooldown_periods)+1

                        eps_t = shock_block(bb, :);
                        eps_t = [0; eps_t(:)];
                        eps_t(7) = 0;

                    end

                else
                    % For all trimming schemes but the first, load in
                    % shock_series
                    eps_t = shock_block(bb, :);
                    eps_t = [0; eps_t(:)];
                    eps_t(7) = 0;
                end

                changes_Taylor = 0;

                shock_series_aux(:,bb) = eps_t;

                X_Taylor_aug = P_ref*X_Taylor + Q_ref*eps_t;

                ELB_indicator_Taylor = exp(X_Taylor_aug(grid.numstates-grid.os+idx.R_cb_ind))*param.R_cb < 1+1e-8;

                if ELB_indicator_Taylor
                    OccBin_one_con_endo_ZLB_shocks;
                    X_Taylor_aug = zdatapiecewise(1,:)';
                    ELB_Taylor_hits(bb)     = 1;
                    ELB_Taylor_expected(bb) = track_Tmax;
                    Taylor_ELB_converge(bb) = changes;
                    changes_Taylor = changes;
                end

                % Check for explosiveness live within block rather than ex-post

                EXPLOSIVE_cons_Taylor  = abs(X_Taylor_aug(end-grid.oc+idx.C_10_ind))  > 3*5;

                EXPLOSIVE_unemp_Taylor = abs(X_Taylor_aug(end-grid.oc+idx.unemp_ind)) > 3*10;

                if EXPLOSIVE_cons_Taylor

                    EXPLOSIVE_cons = 1;

                else

                    EXPLOSIVE_cons = 0;

                end

                if EXPLOSIVE_unemp_Taylor

                    EXPLOSIVE_unemp = 1;

                else

                    EXPLOSIVE_unemp = 0;

                end

                changes_sum = changes_Taylor;

                if changes_sum > 0

                    EXPLOSIVE_cons = 1;

                end

                if EXPLOSIVE_unemp || EXPLOSIVE_cons
                    EXPLOSIVE_phipi = 1;
                    out.Explosive_check = 1;
                    bb = block_size;

                else

                    X_Taylor = X_Taylor_aug;
                    X_Taylor_series(:,bb) = X_Taylor;

                    bb = bb+1;

                end

            end

            if EXPLOSIVE_phipi < 1

                % We cannot keep track of and index the non-explosive
                % blocks within each process (workers would need to share
                % information) nor can we just append to end of struct
                % array. Instead, place non-explosive info in indices of
                % struct array where no explosive behavior occurs.
                kk;

                [sum(ELB_Taylor_hits)/refresh_size];

                [max(abs(X_Taylor_series(end-grid.oc+idx.C_10_ind,:)))];
                [max(abs(X_Taylor_series(end-grid.oc+idx.unemp_ind,:)))];

                out.X_Taylor_series     = X_Taylor_series;
                out.ELB_Taylor_hits     = ELB_Taylor_hits;
                out.ELB_Taylor_expected = ELB_Taylor_expected;

                out.shock_series = shock_series_aux;
                out.shock_series_raw = [zeros(1,refresh_size); shock_block'];

            else
                out.Explosive_check = 1;
            end

end
