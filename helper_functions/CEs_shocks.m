tic;

Tmax = 1001;

nperiods_decomp_IRF = Tmax;
num_decomp_case     = 8;
compute_start_period = 401;

H_tilde_ss = SS_stats.H_tilde;

clear input

x0 = zeros(grid.numstates,1);

x0(end-9+lllllll) = shock_sets(lllllll);

clear MX IRF_state_sparse x

MX               = [eye(length(x0));gx];
IRF_state_sparse = [];
x                = x0;

for t=1:Tmax
    IRF_state_sparse(:,t)=(MX*x)';
    x=hx*x;
    
end

[meshes.b,meshes.a,meshes.se] = ndgrid(grid.b,grid.a,grid.se);

util     = @(c)  (c.^(1-param.sigma))./(1-param.sigma);
mutil    = @(c)  (1./(c.^param.sigma));
invutil  = @(u)  (((1-param.sigma).*u).^(1/(1-param.sigma)));
invmutil = @(mu) ((1./mu).^(1/param.sigma));

NN   = grid.nb*grid.na*grid.nse; % Number of points in the full gri

VALUE_ind   = 1:NN;
mutil_c_ind = NN + (1:NN);
Va_ind      = 2*NN + (1:NN);

IRF_distr=Gamma_state*IRF_state_sparse(1:grid.numstates-grid.os,1:Tmax);

Statenext_IRFs   = repmat(Xss(end-grid.os+1:end),[1 Tmax]) + IRF_state_sparse(grid.numstates-grid.os+1:grid.numstates,:);
Control_IRFs     = repmat(Yss(end-grid.oc+1:end),[1 Tmax]) + Gamma_control(end-grid.oc+1:end,:)*IRF_state_sparse(grid.numstates+1:end,:);

Statenext_IRFs   = exp(Statenext_IRFs);
Control_IRFs     = exp(Control_IRFs);

% preparation

dist_a  = diff([0;marginal_a']);
dist_b  = diff([0;marginal_b]);
dist_se = diff([0;marginal_se]);

IRF_distr_full = repmat([dist_b;dist_a;dist_se],[1 Tmax]) + IRF_distr;

IRF_distr_full_b  = IRF_distr_full(1:grid.nb,:);
IRF_distr_full_a  = IRF_distr_full(grid.nb+(1:grid.na),:);
IRF_distr_full_se = IRF_distr_full(grid.nb+grid.na+(1:grid.nse),:);

indexes_v2_anal;

compute_history_Taylor_QE_compare;

param.beta_aux = param.beta;

% for kk = 1:2
    
    clear input out
        
    Controlnext_sparse_aux = IRF_state_sparse(grid.numstates+1:end,nperiods_decomp_IRF);
    
    Controlnext_aux    = Yss .* (1+Gamma_control*(Controlnext_sparse_aux));
    
    VALUE_next_aux2    = util(Controlnext_aux(VALUE_ind));
    Va_next_aux2       = mutil(Controlnext_aux(Va_ind));
    mutil_c_next_aux2  = mutil(Controlnext_aux(mutil_c_ind));
    
    for jjjj = 1:num_decomp_case
        
        input(jjjj).VALUE_next_input    = VALUE_next_aux2(:);
        input(jjjj).Va_next_input       = Va_next_aux2(:);
        input(jjjj).mutil_c_next_input  = mutil_c_next_aux2(:);
        input(jjjj).VALUE_next_input2   = VALUE_next_aux2(:);

        EVs(jjjj).VALUEs_next   = zeros(grid.nb*grid.na*grid.nse,40);
        EVs(jjjj).Vas_next      = zeros(grid.nb*grid.na*grid.nse,40);
        EVs(jjjj).mutil_cs_next = zeros(grid.nb*grid.na*grid.nse,40);
        
    end

    input(1).VALUE_next_input    = Value(:);
    input(1).Va_next_input       = Va(:);
    input(1).mutil_c_next_input  = mutil_c(:);
    input(1).VALUE_next_input2   = Value(:);

    decomp_agg(lllllll).a_decomp_IRF = zeros(num_decomp_case,nperiods_decomp_IRF-1);
    decomp_agg(lllllll).b_decomp_IRF = zeros(num_decomp_case,nperiods_decomp_IRF-1);
    decomp_agg(lllllll).c_decomp_IRF = zeros(num_decomp_case,nperiods_decomp_IRF-1);

    for jj =1:compute_start_period-1
        
        tt = compute_start_period-jj;

        marginal_bminus  = IRF_distr_full_b(:,tt);
        marginal_aminus  = IRF_distr_full_a(:,tt);
        marginal_seminus = IRF_distr_full_se(:,tt);
        
        cumdist = zeros(grid.nb+1,grid.na+1,grid.nse+1);
        cumdist(2:end,2:end,2:end) = Copula({cumsum(marginal_bminus),cumsum(marginal_aminus)',cumsum(marginal_seminus)});
        MU = diff(diff(diff(cumdist,1,1),1,2),1,3);
    
        N        = N_IRFs(1,tt);
        lambda   = param.lambda;
        V        = V_IRFs(1,tt);
        Vnext    = V_IRFs(1,tt+1);
        W        = W_IRFs(1,tt+1);
        nn       = ((1-param.tau_w)*W/param.psi)^(1/param.xi);
        lambda_e = param.lambda;
        vartheta = 1;
        PROFIT   = Profit_IRFs(1,tt);
        Profit_FI= Profit_FI_IRFs_p(1,tt);
        C_b      = C_b_IRFs(1,tt);
        R_A      = r_a_IRFs(1,tt);
        Q        = Q_IRFs(1,tt+1);
        R_cb     = Rcb_IRFs(1,tt+1);
        LT       = LT_IRFs(1,tt);
        
        R_cbminus = Rcb_IRFs(1,tt);

        PI       = pi_IRFs(1,tt);
        PInext   = pi_IRFs(1,tt+1);
        UU       = 1;

        AA_aux   = 1;

        U      = 1-N-param.Eshare;
        MATCH  = (U+lambda.*(N)).*V/(((U+lambda.*(N))^param.alpha+V^param.alpha)^(1/param.alpha));
        
        Ntilde = (1-lambda)*N + MATCH;
        Utilde = (1-Ntilde-param.Eshare)/(1-param.Eshare);
        
        f_rate = MATCH/(U+lambda*(N));
        
        P_SS3_aux   = kron([1-lambda+lambda*f_rate,lambda*(1-f_rate);f_rate,1-f_rate],eye(grid.ns));
        P_SS3_aux   = [P_SS3_aux, zeros(grid.nse-1,1)];
        lastrow     = [zeros(1,2*grid.ns),1];
        P_SS3_aux   = [P_SS3_aux; lastrow];
        
        H_tilde_ss  = kron(param.P_SS3',speye(grid.nb*grid.na));

        H_tilde     = kron(P_SS3_aux',speye(grid.nb*grid.na));
        
        MU_tilde    = H_tilde*MU(:);
        
        MU_tilde = reshape(MU_tilde,[grid.nb grid.na grid.nse]);
        
        MU_tilde_ss = H_tilde_ss*MU(:);
        
        MU_tilde_ss = reshape(MU_tilde_ss,[grid.nb grid.na grid.nse]);
        
        N_aux  = (1-param.in)*Ntilde;
        U_aux  = (1-param.in)*(1-Ntilde-param.Eshare) + param.out*param.Eshare;
        
        MATCHnext = (U_aux+lambda.*(N_aux)).*Vnext/(((U_aux+lambda.*(N_aux))^param.alpha+Vnext^param.alpha)^(1/param.alpha));
        
        f_ratenext = MATCHnext/(U_aux+lambda*(N_aux));
        
        P_SS3next_aux  = kron([1-lambda+lambda*f_ratenext,lambda*(1-f_ratenext);f_ratenext,1-f_ratenext],eye(grid.ns));
        P_SS3next_aux  = [P_SS3next_aux repmat(0, [grid.nse-1 1])];
        lastrow        = [zeros(1,2*grid.ns),1];
        P_SS3next_aux  = [P_SS3next_aux; lastrow];
        
        P_transition_exp = param.P_SS2*P_SS3next_aux;

        P_transition_exp_ss = param.P_SS2*param.P_SS3;
        
        switch(param.FEAR_shock)
            case('unemp')
                lambda2_aux = lambda_e;
            case('fin')
                lambda2_aux = param.lambda;
            case('both')
                lambda2_aux = lambda_e;
        end
        
        alpha_aux   = param.alpha;
        MATCHnext2  = (U_aux+lambda2_aux.*(N_aux)).*Vnext/(((U_aux+lambda2_aux.*(N_aux))^alpha_aux+Vnext^alpha_aux)^(1/alpha_aux));
        f_ratenext2 = MATCHnext2/(U_aux+lambda2_aux*(N_aux));
        
        P_SS3next_aux2  = kron([1-lambda2_aux+lambda2_aux*f_ratenext2,lambda2_aux*(1-f_ratenext2);f_ratenext2,1-f_ratenext2],eye(grid.ns));
        P_SS3next_aux2  = [P_SS3next_aux2 repmat(0, [grid.nse-1 1])];
        lastrow        = [zeros(1,2*grid.ns),1];
        P_SS3next_aux2  = [P_SS3next_aux2; lastrow];
        
        P_transition_exp2  = param.P_SS2*P_SS3next_aux2;

        % 1) base : every price variables take on their ss value
        
        input(1).nn                 = SS_stats.nn;
        input(1).N                  = SS_stats.N;
        input(1).W                  = SS_stats.W;
        input(1).tau_w              = param.tau_w;
        input(1).tau_a              = param.tau_a;
        input(1).PROFIT             = SS_stats.PROFIT;
        input(1).C_b                = SS_stats.C_b;
        input(1).Ntilde             = SS_stats.Ntilde;
        input(1).R_A                = SS_stats.R_A;
        input(1).Q                  = SS_stats.Q;
        input(1).R_cbminus          = SS_stats.R_cbminus;
        input(1).R_cb               = SS_stats.R_cb;
        input(1).PI                 = SS_stats.PI;
        input(1).PInext             = SS_stats.PInext;
        input(1).P_transition_exp2  = P_transition_exp_ss;
        input(1).MU_tilde           = MU_tilde_ss;
        input(1).UU                 = 1;
        input(1).LT                 = SS_stats.LT;
        input(1).vartheta           = 1;
        input(1).P_transition_exp   = P_transition_exp_ss;
        input(1).AA_aux             = 1;
        input(1).Profit_FI          = SS_stats.Profit_FI;
        
        % 2) full : every price variables respond
        
        input(2).nn                 = nn;
        input(2).N                  = N;
        input(2).W                  = W;
        input(2).tau_w              = param.tau_w;
        input(2).tau_a              = param.tau_a;
        input(2).PROFIT             = PROFIT;
        input(2).C_b                = C_b;
        input(2).Ntilde             = Ntilde;
        input(2).R_A                = R_A;
        input(2).Q                  = Q;
        input(2).R_cbminus          = R_cbminus;
        input(2).R_cb               = R_cb;
        input(2).PI                 = PI;
        input(2).PInext             = PInext;
        input(2).P_transition_exp2  = P_transition_exp2;
        input(2).MU_tilde           = MU_tilde;
        input(2).UU                 = 1;
        input(2).LT                 = LT;
        input(2).vartheta           = 1;
        input(2).P_transition_exp   = P_transition_exp;
        input(2).AA_aux             = 1;
        input(2).Profit_FI          = Profit_FI;

        % 3) Profit + R_A
        
        input(3).nn                 = SS_stats.nn;
        input(3).N                  = SS_stats.N;
        input(3).W                  = SS_stats.W;
        input(3).tau_w              = param.tau_w;
        input(3).tau_a              = param.tau_a;
        input(3).PROFIT             = PROFIT;
        input(3).C_b                = SS_stats.C_b;
        input(3).Ntilde             = Ntilde;
        input(3).R_A                = R_A;
        input(3).Q                  = 1;
        input(3).R_cbminus          = SS_stats.R_cbminus;
        input(3).R_cb               = SS_stats.R_cb;
        input(3).PI                 = SS_stats.PI;
        input(3).PInext             = SS_stats.PInext;
        input(3).P_transition_exp2  = P_transition_exp_ss;
        input(3).MU_tilde           = MU_tilde_ss;
        input(3).UU                 = 1;
        input(3).LT                 = SS_stats.LT;
        input(3).vartheta           = 1;
        input(3).P_transition_exp   = P_transition_exp_ss;
        input(3).AA_aux             = 1;
        input(3).Profit_FI          = Profit_FI;

        % 4) R_cb + PI
        
        input(4).nn                 = SS_stats.nn;
        input(4).N                  = SS_stats.N;
        input(4).W                  = SS_stats.W;
        input(4).tau_w              = param.tau_w;
        input(4).tau_a              = param.tau_a;
        input(4).PROFIT             = SS_stats.PROFIT;
        input(4).C_b                = SS_stats.C_b;
        input(4).Ntilde             = SS_stats.Ntilde;
        input(4).R_A                = SS_stats.R_A;
        input(4).Q                  = SS_stats.Q;
        input(4).R_cbminus          = R_cbminus;
        input(4).R_cb               = R_cb;
        input(4).PI                 = PI;
        input(4).PInext             = PInext;
        input(4).P_transition_exp2  = P_transition_exp_ss;
        input(4).MU_tilde           = MU_tilde_ss;
        input(4).UU                 = 1;
        input(4).LT                 = SS_stats.LT;
        input(4).vartheta           = 1;
        input(4).P_transition_exp   = P_transition_exp_ss;
        input(4).AA_aux             = 1;
        input(4).Profit_FI          = SS_stats.Profit_FI;

        % 5) W
        
        input(5).nn                 = nn;
        input(5).N                  = N;
        input(5).W                  = W;
        input(5).tau_w              = param.tau_w;
        input(5).tau_a              = param.tau_a;
        input(5).PROFIT             = SS_stats.PROFIT;
        input(5).C_b                = SS_stats.C_b;
        input(5).Ntilde             = SS_stats.Ntilde;
        input(5).R_A                = SS_stats.R_A;
        input(5).Q                  = SS_stats.Q;
        input(5).R_cbminus          = SS_stats.R_cbminus;
        input(5).R_cb               = SS_stats.R_cb;
        input(5).PI                 = SS_stats.PI;
        input(5).PInext             = SS_stats.PInext;
        input(5).P_transition_exp2  = P_transition_exp_ss;
        input(5).MU_tilde           = MU_tilde_ss;
        input(5).UU                 = 1;
        input(5).LT                 = SS_stats.LT;
        input(5).vartheta           = 1;
        input(5).P_transition_exp   = P_transition_exp_ss;
        input(5).AA_aux             = 1;
        input(5).Profit_FI          = SS_stats.Profit_FI;

        % 6) f
        
        input(6).nn                 = SS_stats.nn;
        input(6).N                  = SS_stats.N;
        input(6).W                  = SS_stats.W;
        input(6).tau_w              = param.tau_w;
        input(6).tau_a              = param.tau_a;
        input(6).PROFIT             = SS_stats.PROFIT;
        input(6).C_b                = SS_stats.C_b;
        input(6).Ntilde             = SS_stats.Ntilde;
        input(6).R_A                = SS_stats.R_A;
        input(6).Q                  = SS_stats.Q;
        input(6).R_cbminus          = SS_stats.R_cbminus;
        input(6).R_cb               = SS_stats.R_cb;
        input(6).PI                 = SS_stats.PI;
        input(6).PInext             = SS_stats.PInext;
        input(6).P_transition_exp2  = P_transition_exp2;
        input(6).MU_tilde           = MU_tilde;
        input(6).UU                 = 1;
        input(6).LT                 = SS_stats.LT;
        input(6).vartheta           = 1;
        input(6).P_transition_exp   = P_transition_exp;
        input(6).AA_aux             = 1;
        input(6).Profit_FI          = SS_stats.Profit_FI;

        % 7) Gov
        
        input(7).nn                 = SS_stats.nn;
        input(7).N                  = SS_stats.N;
        input(7).W                  = SS_stats.W;
        input(7).tau_w              = param.tau_w;
        input(7).tau_a              = param.tau_a;
        input(7).PROFIT             = SS_stats.PROFIT;
        input(7).C_b                = C_b;
        input(7).Ntilde             = SS_stats.Ntilde;
        input(7).R_A                = SS_stats.R_A;
        input(7).Q                  = SS_stats.Q;
        input(7).R_cbminus          = SS_stats.R_cbminus;
        input(7).R_cb               = SS_stats.R_cb;
        input(7).PI                 = SS_stats.PI;
        input(7).PInext             = SS_stats.PInext;
        input(7).P_transition_exp2  = P_transition_exp_ss;
        input(7).MU_tilde           = MU_tilde_ss;
        input(7).UU                 = 1;
        input(7).LT                 = LT;
        input(7).vartheta           = 1;
        input(7).P_transition_exp   = P_transition_exp_ss;
        input(7).AA_aux             = 1;
        input(7).Profit_FI          = SS_stats.Profit_FI;

        % 8) Q
        
        input(8).nn                 = SS_stats.nn;
        input(8).N                  = SS_stats.N;
        input(8).W                  = SS_stats.W;
        input(8).tau_w              = param.tau_w;
        input(8).tau_a              = param.tau_a;
        input(8).PROFIT             = SS_stats.PROFIT;
        input(8).C_b                = SS_stats.C_b;
        input(8).Ntilde             = SS_stats.Ntilde;
        input(8).R_A                = SS_stats.R_A;
        input(8).Q                  = Q;
        input(8).R_cbminus          = SS_stats.R_cbminus;
        input(8).R_cb               = SS_stats.R_cb;
        input(8).PI                 = SS_stats.PI;
        input(8).PInext             = SS_stats.PInext;
        input(8).P_transition_exp2  = P_transition_exp_ss;
        input(8).MU_tilde           = MU_tilde_ss;
        input(8).UU                 = 1;
        input(8).LT                 = SS_stats.LT;
        input(8).vartheta           = 1;
        input(8).P_transition_exp   = P_transition_exp_ss;
        input(8).AA_aux             = 1;
        input(8).Profit_FI          = SS_stats.Profit_FI;

        parfor kkkk = 1:num_decomp_case
            
            out(kkkk)     = compute_CE(input(kkkk),param,grid);
            
        end
        
        for jjjj = 1:num_decomp_case

            input(jjjj).VALUE_next_input   = out(jjjj).VALUE(:);
            input(jjjj).Va_next_input      = out(jjjj).Va(:);
            input(jjjj).mutil_c_next_input = out(jjjj).mutil_c(:);
            input(jjjj).VALUE_next_input2  = out(jjjj).VALUE_aux3(:);

            decomp_agg(lllllll).a_decomp_IRF(jjjj,tt) = out(jjjj).A;
            decomp_agg(lllllll).b_decomp_IRF(jjjj,tt) = out(jjjj).B;
            decomp_agg(lllllll).c_decomp_IRF(jjjj,tt) = out(jjjj).C;
    
            if tt <= 41

                EVs(jjjj).VALUEs_next(:,tt)   = out(jjjj).VALUE(:);
                EVs(jjjj).Vas_next(:,tt)      = out(jjjj).Va(:);
                EVs(jjjj).mutil_cs_next(:,tt) = out(jjjj).mutil_c(:);
            
            end

        end
        
        if tt == 1
            
            VALUE_base = H_tilde_ss'*out(1).VALUE_aux2(:);
            AC_base    = H_tilde_ss'*out(1).AC(:);
            VALUE_base = reshape(VALUE_base,[grid.nb grid.na grid.nse]);
            AC_base    = reshape(AC_base,[grid.nb grid.na grid.nse]);
            
            Vshocks.VALUE_shock = H_tilde'*out(2).VALUE_aux2(:);
            ACshocks.AC_shock   = H_tilde'*out(2).AC(:);
            Vshocks.VALUE_shock = reshape(Vshocks.VALUE_shock,[grid.nb grid.na grid.nse]);
            ACshocks.AC_shock   = reshape(ACshocks.AC_shock,[grid.nb grid.na grid.nse]);

            Vshocks.VALUE_PROFIT = H_tilde_ss'*out(3).VALUE_aux2(:);
            ACshocks.AC_PROFIT   = H_tilde_ss'*out(3).AC(:);
            Vshocks.VALUE_PROFIT = reshape(Vshocks.VALUE_PROFIT,[grid.nb grid.na grid.nse]);
            ACshocks.AC_PROFIT   = reshape(ACshocks.AC_PROFIT,[grid.nb grid.na grid.nse]);

            Vshocks.VALUE_RR = H_tilde_ss'*out(4).VALUE_aux2(:);
            ACshocks.AC_RR   = H_tilde_ss'*out(4).AC(:);
            Vshocks.VALUE_RR = reshape(Vshocks.VALUE_RR,[grid.nb grid.na grid.nse]);
            ACshocks.AC_RR   = reshape(ACshocks.AC_RR,[grid.nb grid.na grid.nse]);

            Vshocks.VALUE_W = H_tilde_ss'*out(5).VALUE_aux2(:);
            ACshocks.AC_W   = H_tilde_ss'*out(5).AC(:);
            Vshocks.VALUE_W = reshape(Vshocks.VALUE_W,[grid.nb grid.na grid.nse]);
            ACshocks.AC_W   = reshape(ACshocks.AC_W,[grid.nb grid.na grid.nse]);

            Vshocks.VALUE_f = H_tilde'*out(6).VALUE_aux2(:);
            ACshocks.AC_f   = H_tilde'*out(6).AC(:);
            Vshocks.VALUE_f = reshape(Vshocks.VALUE_f,[grid.nb grid.na grid.nse]);
            ACshocks.AC_f   = reshape(ACshocks.AC_f,[grid.nb grid.na grid.nse]);

            Vshocks.VALUE_GOV = H_tilde_ss'*out(7).VALUE_aux2(:);
            ACshocks.AC_GOV   = H_tilde_ss'*out(7).AC(:);
            Vshocks.VALUE_GOV = reshape(Vshocks.VALUE_GOV,[grid.nb grid.na grid.nse]);
            ACshocks.AC_GOV   = reshape(ACshocks.AC_GOV,[grid.nb grid.na grid.nse]);

            Vshocks.VALUE_Q   = H_tilde_ss'*out(8).VALUE_aux2(:);
            ACshocks.AC_Q     = H_tilde_ss'*out(8).AC(:);
            Vshocks.VALUE_Q   = reshape(Vshocks.VALUE_Q,[grid.nb grid.na grid.nse]);
            ACshocks.AC_Q     = reshape(ACshocks.AC_Q,[grid.nb grid.na grid.nse]);
            
            CEs         = ((Vshocks.VALUE_shock  - ACshocks.AC_shock  + AC_base)./VALUE_base).^(1/(1-param.sigma))-1;
            CEs_PROFIT  = ((Vshocks.VALUE_PROFIT - ACshocks.AC_PROFIT + AC_base)./VALUE_base).^(1/(1-param.sigma))-1;
            CEs_RR      = ((Vshocks.VALUE_RR     - ACshocks.AC_RR     + AC_base)./VALUE_base).^(1/(1-param.sigma))-1;
            CEs_W       = ((Vshocks.VALUE_W      - ACshocks.AC_W      + AC_base)./VALUE_base).^(1/(1-param.sigma))-1;
            CEs_f       = ((Vshocks.VALUE_f      - ACshocks.AC_f      + AC_base)./VALUE_base).^(1/(1-param.sigma))-1;
            CEs_GOV     = ((Vshocks.VALUE_GOV    - ACshocks.AC_GOV    + AC_base)./VALUE_base).^(1/(1-param.sigma))-1;
            CEs_Q       = ((Vshocks.VALUE_Q      - ACshocks.AC_Q      + AC_base)./VALUE_base).^(1/(1-param.sigma))-1;

        end
        
        clear out
        
    end
    
    
% end

com_CE;

decomp_CE      = abs(CE      - CE_PROFIT      - CE_RR      - CE_Q      - CE_W      - CE_f      - CE_GOV)      ./abs(CE);
decomp_CE_Q1   = abs(CE_Q1   - CE_PROFIT_Q1   - CE_RR_Q1   - CE_Q_Q1   - CE_W_Q1   - CE_f_Q1   - CE_GOV_Q1)   ./abs(CE_Q1);
decomp_CE_Q2   = abs(CE_Q2   - CE_PROFIT_Q2   - CE_RR_Q2   - CE_Q_Q2   - CE_W_Q2   - CE_f_Q2   - CE_GOV_Q2)   ./abs(CE_Q2);
decomp_CE_Q3   = abs(CE_Q3   - CE_PROFIT_Q3   - CE_RR_Q3   - CE_Q_Q3   - CE_W_Q3   - CE_f_Q3   - CE_GOV_Q3)   ./abs(CE_Q3);
decomp_CE_Q4   = abs(CE_Q4   - CE_PROFIT_Q4   - CE_RR_Q4   - CE_Q_Q4   - CE_W_Q4   - CE_f_Q4   - CE_GOV_Q4)   ./abs(CE_Q4);
decomp_CE_Q5   = abs(CE_Q5   - CE_PROFIT_Q5   - CE_RR_Q5   - CE_Q_Q5   - CE_W_Q5   - CE_f_Q5   - CE_GOV_Q5)   ./abs(CE_Q5);
decomp_CE_P90  = abs(CE_P90  - CE_PROFIT_P90  - CE_RR_P90  - CE_Q_P90  - CE_W_P90  - CE_f_P90  - CE_GOV_P90)  ./abs(CE_P90);
decomp_CE_P99  = abs(CE_P99  - CE_PROFIT_P99  - CE_RR_P99  - CE_Q_P99  - CE_W_P99  - CE_f_P99  - CE_GOV_P99)  ./abs(CE_P99);
decomp_CE_P999 = abs(CE_P999 - CE_PROFIT_P999 - CE_RR_P999 - CE_Q_P999 - CE_W_P999 - CE_f_P999 - CE_GOV_P999) ./abs(CE_P999);


disp('-------------------------');
disp(' CE (Total effect) ');
disp('-------------------------');
disp('')
fprintf(' CE for B01  : %1.2f \n', 10000*CE_B01);
fprintf(' CE for B1   : %1.2f \n', 10000*CE_B1);
fprintf(' CE for B10  : %1.2f \n', 10000*CE_B10);
fprintf(' CE for Q1   : %1.2f \n', 10000*CE_Q1);
fprintf(' CE for Q2   : %1.2f \n', 10000*CE_Q2);
fprintf(' CE for Q3   : %1.2f \n', 10000*CE_Q3);
fprintf(' CE for Q4   : %1.2f \n', 10000*CE_Q4);
fprintf(' CE for Q5   : %1.2f \n', 10000*CE_Q5);
fprintf(' CE for P90  : %1.2f \n', 10000*CE_P90);
fprintf(' CE for P99  : %1.2f \n', 10000*CE_P99);
fprintf(' CE for P999 : %1.2f \n', 10000*CE_P999);
fprintf(' CE Ent      : %1.2f \n', 10000*CE_Ent);
fprintf(' CE Emp      : %1.2f \n', 10000*CE_Emp);
fprintf(' CE Unemp    : %1.2f \n', 10000*CE_Unemp);
fprintf(' CE Avg      : %1.2f \n', 10000*CE);


disp('-------------------------');
disp(' CE (RR effect) ');
disp('-------------------------');
disp('')
fprintf(' CE for B01  : %1.2f \n', 10000*CE_RR_B01);
fprintf(' CE for B1   : %1.2f \n', 10000*CE_RR_B1);
fprintf(' CE for B10  : %1.2f \n', 10000*CE_RR_B10);
fprintf(' CE for Q1   : %1.2f \n', 10000*CE_RR_Q1);
fprintf(' CE for Q2   : %1.2f \n', 10000*CE_RR_Q2);
fprintf(' CE for Q3   : %1.2f \n', 10000*CE_RR_Q3);
fprintf(' CE for Q4   : %1.2f \n', 10000*CE_RR_Q4);
fprintf(' CE for Q5   : %1.2f \n', 10000*CE_RR_Q5);
fprintf(' CE for P90  : %1.2f \n', 10000*CE_RR_P90);
fprintf(' CE for P99  : %1.2f \n', 10000*CE_RR_P99);
fprintf(' CE for P999 : %1.2f \n', 10000*CE_RR_P999);
fprintf(' CE Ent      : %1.2f \n', 10000*CE_RR_Ent);
fprintf(' CE Emp      : %1.2f \n', 10000*CE_RR_Emp);
fprintf(' CE Unemp    : %1.2f \n', 10000*CE_RR_Unemp);
fprintf(' CE Avg      : %1.2f \n', 10000*CE_RR);

disp('-------------------------');
disp(' CE (PROFIT effect) ');
disp('-------------------------');
disp('')
fprintf(' CE for B01  : %1.2f \n', 10000*CE_PROFIT_B01);
fprintf(' CE for B1   : %1.2f \n', 10000*CE_PROFIT_B1);
fprintf(' CE for B10  : %1.2f \n', 10000*CE_PROFIT_B10);
fprintf(' CE for Q1   : %1.2f \n', 10000*CE_PROFIT_Q1);
fprintf(' CE for Q2   : %1.2f \n', 10000*CE_PROFIT_Q2);
fprintf(' CE for Q3   : %1.2f \n', 10000*CE_PROFIT_Q3);
fprintf(' CE for Q4   : %1.2f \n', 10000*CE_PROFIT_Q4);
fprintf(' CE for Q5   : %1.2f \n', 10000*CE_PROFIT_Q5);
fprintf(' CE for P90  : %1.2f \n', 10000*CE_PROFIT_P90);
fprintf(' CE for P99  : %1.2f \n', 10000*CE_PROFIT_P99);
fprintf(' CE for P999 : %1.2f \n', 10000*CE_PROFIT_P999);
fprintf(' CE Ent      : %1.2f \n', 10000*CE_PROFIT_Ent);
fprintf(' CE Emp      : %1.2f \n', 10000*CE_PROFIT_Emp);
fprintf(' CE Unemp    : %1.2f \n', 10000*CE_PROFIT_Unemp);
fprintf(' CE Avg      : %1.2f \n', 10000*CE_PROFIT);

disp('-------------------------');
disp(' CE (Q effect) ');
disp('-------------------------');
disp('')
fprintf(' CE for B01  : %1.2f \n', 10000*CE_Q_B01);
fprintf(' CE for B1   : %1.2f \n', 10000*CE_Q_B1);
fprintf(' CE for B10  : %1.2f \n', 10000*CE_Q_B10);
fprintf(' CE for Q1   : %1.2f \n', 10000*CE_Q_Q1);
fprintf(' CE for Q2   : %1.2f \n', 10000*CE_Q_Q2);
fprintf(' CE for Q3   : %1.2f \n', 10000*CE_Q_Q3);
fprintf(' CE for Q4   : %1.2f \n', 10000*CE_Q_Q4);
fprintf(' CE for Q5   : %1.2f \n', 10000*CE_Q_Q5);
fprintf(' CE for P90  : %1.2f \n', 10000*CE_Q_P90);
fprintf(' CE for P99  : %1.2f \n', 10000*CE_Q_P99);
fprintf(' CE for P999 : %1.2f \n', 10000*CE_Q_P999);
fprintf(' CE Ent      : %1.2f \n', 10000*CE_Q_Ent);
fprintf(' CE Emp      : %1.2f \n', 10000*CE_Q_Emp);
fprintf(' CE Unemp    : %1.2f \n', 10000*CE_Q_Unemp);
fprintf(' CE Avg      : %1.2f \n', 10000*CE_Q);

disp('-------------------------');
disp(' CE (W effect) ');
disp('-------------------------');
disp('')
fprintf(' CE for B01  : %1.2f \n', 10000*CE_W_B01);
fprintf(' CE for B1   : %1.2f \n', 10000*CE_W_B1);
fprintf(' CE for B10  : %1.2f \n', 10000*CE_W_B10);
fprintf(' CE for Q1   : %1.2f \n', 10000*CE_W_Q1);
fprintf(' CE for Q2   : %1.2f \n', 10000*CE_W_Q2);
fprintf(' CE for Q3   : %1.2f \n', 10000*CE_W_Q3);
fprintf(' CE for Q4   : %1.2f \n', 10000*CE_W_Q4);
fprintf(' CE for Q5   : %1.2f \n', 10000*CE_W_Q5);
fprintf(' CE for P90  : %1.2f \n', 10000*CE_W_P90);
fprintf(' CE for P99  : %1.2f \n', 10000*CE_W_P99);
fprintf(' CE for P999 : %1.2f \n', 10000*CE_W_P999);
fprintf(' CE Ent      : %1.2f \n', 10000*CE_W_Ent);
fprintf(' CE Emp      : %1.2f \n', 10000*CE_W_Emp);
fprintf(' CE Unemp    : %1.2f \n', 10000*CE_W_Unemp);
fprintf(' CE Avg      : %1.2f \n', 10000*CE_W);

disp('-------------------------');
disp(' CE (f effect) ');
disp('-------------------------');
disp('')
fprintf(' CE for B01  : %1.2f \n', 10000*CE_f_B01);
fprintf(' CE for B1   : %1.2f \n', 10000*CE_f_B1);
fprintf(' CE for B10  : %1.2f \n', 10000*CE_f_B10);
fprintf(' CE for Q1   : %1.2f \n', 10000*CE_f_Q1);
fprintf(' CE for Q2   : %1.2f \n', 10000*CE_f_Q2);
fprintf(' CE for Q3   : %1.2f \n', 10000*CE_f_Q3);
fprintf(' CE for Q4   : %1.2f \n', 10000*CE_f_Q4);
fprintf(' CE for Q5   : %1.2f \n', 10000*CE_f_Q5);
fprintf(' CE for P90  : %1.2f \n', 10000*CE_f_P90);
fprintf(' CE for P99  : %1.2f \n', 10000*CE_f_P99);
fprintf(' CE for P999 : %1.2f \n', 10000*CE_f_P999);
fprintf(' CE Ent      : %1.2f \n', 10000*CE_f_Ent);
fprintf(' CE Emp      : %1.2f \n', 10000*CE_f_Emp);
fprintf(' CE Unemp    : %1.2f \n', 10000*CE_f_Unemp);
fprintf(' CE Avg      : %1.2f \n', 10000*CE_f);

disp('-------------------------');
disp(' CE (GOV effect) ');
disp('-------------------------');
disp('')
fprintf(' CE for B01  : %1.2f \n', 10000*CE_GOV_B01);
fprintf(' CE for B1   : %1.2f \n', 10000*CE_GOV_B1);
fprintf(' CE for B10  : %1.2f \n', 10000*CE_GOV_B10);
fprintf(' CE for Q1   : %1.2f \n', 10000*CE_GOV_Q1);
fprintf(' CE for Q2   : %1.2f \n', 10000*CE_GOV_Q2);
fprintf(' CE for Q3   : %1.2f \n', 10000*CE_GOV_Q3);
fprintf(' CE for Q4   : %1.2f \n', 10000*CE_GOV_Q4);
fprintf(' CE for Q5   : %1.2f \n', 10000*CE_GOV_Q5);
fprintf(' CE for P90  : %1.2f \n', 10000*CE_GOV_P90);
fprintf(' CE for P99  : %1.2f \n', 10000*CE_GOV_P99);
fprintf(' CE for P999 : %1.2f \n', 10000*CE_GOV_P999);
fprintf(' CE Ent      : %1.2f \n', 10000*CE_GOV_Ent);
fprintf(' CE Emp      : %1.2f \n', 10000*CE_GOV_Emp);
fprintf(' CE Unemp    : %1.2f \n', 10000*CE_GOV_Unemp);
fprintf(' CE Avg      : %1.2f \n', 10000*CE_GOV);


CE_Taylor.Vshocks  = Vshocks;
CE_Taylor.ACshocks = ACshocks;

CE_total        = [CE_B01,CE_B1,CE_B10,CE_Q1,CE_Q2,CE_Q3,CE_Q4,CE_Q5,CE_P90,CE_P99,CE_P999,CE_Ent,CE_Emp,CE_Unemp,CE];
CE_PROFIT_total = [CE_PROFIT_B01,CE_PROFIT_B1,CE_PROFIT_B10,CE_PROFIT_Q1,CE_PROFIT_Q2,CE_PROFIT_Q3,CE_PROFIT_Q4,CE_PROFIT_Q5,CE_PROFIT_P90,CE_PROFIT_P99,CE_PROFIT_P999,CE_PROFIT_Ent,CE_PROFIT_Emp,CE_PROFIT_Unemp,CE_PROFIT];
CE_Q_total      = [CE_Q_B01,CE_Q_B1,CE_Q_B10,CE_Q_Q1,CE_Q_Q2,CE_Q_Q3,CE_Q_Q4,CE_Q_Q5,CE_Q_P90,CE_Q_P99,CE_Q_P999,CE_Q_Ent,CE_Q_Emp,CE_Q_Unemp,CE_Q];
CE_W_total      = [CE_W_B01,CE_W_B1,CE_W_B10,CE_W_Q1,CE_W_Q2,CE_W_Q3,CE_W_Q4,CE_W_Q5,CE_W_P90,CE_W_P99,CE_W_P999,CE_W_Ent,CE_W_Emp,CE_W_Unemp,CE_W];
CE_f_total      = [CE_f_B01,CE_f_B1,CE_f_B10,CE_f_Q1,CE_f_Q2,CE_f_Q3,CE_f_Q4,CE_f_Q5,CE_f_P90,CE_f_P99,CE_f_P999,CE_f_Ent,CE_f_Emp,CE_f_Unemp,CE_f];
CE_RR_total     = [CE_RR_B01,CE_RR_B1,CE_RR_B10,CE_RR_Q1,CE_RR_Q2,CE_RR_Q3,CE_RR_Q4,CE_RR_Q5,CE_RR_P90,CE_RR_P99,CE_RR_P999,CE_RR_Ent,CE_RR_Emp,CE_RR_Unemp,CE_RR];
CE_GOV_total    = [CE_GOV_B01,CE_GOV_B1,CE_GOV_B10,CE_GOV_Q1,CE_GOV_Q2,CE_GOV_Q3,CE_GOV_Q4,CE_GOV_Q5,CE_GOV_P90,CE_GOV_P99,CE_GOV_P999,CE_GOV_Ent,CE_GOV_Emp,CE_GOV_Unemp,CE_GOV];

toc;
