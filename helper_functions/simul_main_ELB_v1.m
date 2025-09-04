num_simul = 1;

if spec_settings.simul_disagg
    sh_type_iter = ["demand", "supply", "all"];
else
    sh_type_iter = ["all"];
end

for sh_type = sh_type_iter

    ELB_hits = [];
    ELB_expected = [];


    ELB_converge = zeros(spec_settings.num_simul_period, 1);

    X_t = zeros(grid.num_endo,1);

    X_series = zeros(grid.num_endo,spec_settings.num_simul_period);

    tic;

    for tt = 1:spec_settings.num_simul_period


        tt

        if tt ~= 1 && mod(tt, 300) == 1
            X_t = zeros(grid.num_endo, 1);

        end

        eps_t = spec_settings.shock_series(tt, :);

        eps_t = [0; eps_t(:)];
        eps_t(7) = 0;


        if sh_type == "demand"
            %If we are looking at demand shocks, we zero out all supply
            %shocks -- these are Productivity (z) and Price Markup
            eps_t(4) = 0; %Producitivity
            eps_t(9) = 0; %Price Markup
        elseif sh_type == "supply"
            % Else, we zero out all demand shocks (the rest)
            eps_t(2) = 0; %Fixed Cost
            eps_t(3) = 0; %Risk Premia
            eps_t(5) = 0; %Lump-Sum Transfer
            eps_t(6) = 0; %Discount Factor
            eps_t(8) = 0; %Investment Technology
            eps_t(10) = 0; %Wage Markup
        end

        % Iterate transition Eq.
        X_t_aug = P_ref*X_t + Q_ref*eps_t;

    if spec_settings.simul_ELB
        ELB_indicator = exp(X_t_aug(grid.numstates-grid.os+idx.R_cb_ind))*param.R_cb < 1+1e-8;



        if ELB_indicator
            OccBin_one_con_endo_ZLB_v4_pg;
            X_t_aug = zdatapiecewise(1,:)';
            ELB_hits(end+1) = tt;
            ELB_expected(end+1) = track_Tmax;

            ELB_converge(tt) = changes;
        end
    end

        % New X_t, next period.
        X_t = X_t_aug;
        X_series(:,tt)    = X_t;
    end

    toc;
%%

% Compute deviation of variables from steady state and store simulation
% results:


    IRFs(1) = compute_IRFs_simul(X_series,Xss,Yss,grid,param,SS_stats,idx);

    [mean_simul,std_simul,std_simul_annual, std_simul_quarterly] = simul_vars_compute(IRFs,num_simul,SS_stats,param);

    simul_results.mean_simul = mean_simul;

    simul_results.std_simul  = std_simul;

    simul_results.std_simul_annual = std_simul_annual;

    simul_results.std_simul_quarterly = std_simul_quarterly;

    simul_results.shocks_simul   = spec_settings.shock_series;

    simul_results.IRFs = IRFs;

    results = simul_results;

    Policy_Rule = [spec_settings.policy_rule; spec_settings.policy_rule];
    Aggregation_Time = ["Annual"; "Quarterly"];
    Output = [std_simul_annual.Y; std_simul_quarterly.Y];
    AggConsumption = [std_simul_annual.C; std_simul_quarterly.C];
    Inflation = [std_simul_annual.pi; std_simul_quarterly.pi];
    AggIncome = [std_simul_annual.I; std_simul_quarterly.I];
    NominalRate = [std_simul_annual.R_cb; std_simul_quarterly.R_cb];
    RealRate = [std_simul_annual.R; std_simul_quarterly.R];
    Wage = [std_simul_annual.W; std_simul_quarterly.W];
    Unemployment = [std_simul_annual.unemp * 100; std_simul_quarterly.unemp * 100];
    Profit = [std_simul_annual.Profit; std_simul_quarterly.Profit];

    Output(isnan(Output)) = -1;
    AggConsumption(isnan(AggConsumption)) = -1;
    Inflation(isnan(Inflation)) = -1;
    AggIncome(isnan(AggIncome)) = -1;
    NominalRate(isnan(NominalRate)) = -1;
    RealRate(isnan(RealRate)) = -1;
    Wage(isnan(Wage)) = -1;
    Unemployment(isnan(Unemployment)) = -1;
    Profit(isnan(Profit)) = -1;

    ConsumptionB1_W = [std_simul_annual.C_1; std_simul_quarterly.C_1];
    ConsumptionB10_W = [std_simul_annual.C_10; std_simul_quarterly.C_10];
    ConsumptionQ1_W = [std_simul_annual.C_Q1; std_simul_quarterly.C_Q1];
    ConsumptionQ2_W = [std_simul_annual.C_Q2; std_simul_quarterly.C_Q2];
    ConsumptionQ3_W = [std_simul_annual.C_Q3; std_simul_quarterly.C_Q3];
    ConsumptionQ4_W = [std_simul_annual.C_Q4; std_simul_quarterly.C_Q4];
    ConsumptionQ5_W = [std_simul_annual.C_Q5; std_simul_quarterly.C_Q5];
    ConsumptionT10_W = [std_simul_annual.C_90; std_simul_quarterly.C_90];
    ConsumptionT1_W = [std_simul_annual.C_99; std_simul_quarterly.C_99];

    ConsumptionB1_W(isnan(ConsumptionB1_W)) = -1;
    ConsumptionB10_W(isnan(ConsumptionB10_W)) = -1;
    ConsumptionQ1_W(isnan(ConsumptionQ1_W)) = -1;
    ConsumptionQ2_W(isnan(ConsumptionQ2_W)) = -1;
    ConsumptionQ3_W(isnan(ConsumptionQ3_W)) = -1;
    ConsumptionQ4_W(isnan(ConsumptionQ4_W)) = -1;
    ConsumptionQ5_W(isnan(ConsumptionQ5_W)) = -1;
    ConsumptionT10_W(isnan(ConsumptionT10_W)) = -1;
    ConsumptionT1_W(isnan(ConsumptionT1_W)) = -1;

    ConsumptionB1_I = [std_simul_annual.C_I_1; std_simul_quarterly.C_I_1];
    ConsumptionB10_I = [std_simul_annual.C_I_10; std_simul_quarterly.C_I_10];
    ConsumptionQ1_I = [std_simul_annual.C_I_Q1; std_simul_quarterly.C_I_Q1];
    ConsumptionQ2_I = [std_simul_annual.C_I_Q2; std_simul_quarterly.C_I_Q2];
    ConsumptionQ3_I = [std_simul_annual.C_I_Q3; std_simul_quarterly.C_I_Q3];
    ConsumptionQ4_I = [std_simul_annual.C_I_Q4; std_simul_quarterly.C_I_Q4];
    ConsumptionQ5_I = [std_simul_annual.C_I_Q5; std_simul_quarterly.C_I_Q5];
    ConsumptionT10_I = [std_simul_annual.C_I_90; std_simul_quarterly.C_I_90];
    ConsumptionT1_I = [std_simul_annual.C_I_99; std_simul_quarterly.C_I_99];

    ConsumptionB1_I(isnan(ConsumptionB1_I)) = -1;
    ConsumptionB10_I(isnan(ConsumptionB10_I)) = -1;
    ConsumptionQ1_I(isnan(ConsumptionQ1_I)) = -1;
    ConsumptionQ2_I(isnan(ConsumptionQ2_I)) = -1;
    ConsumptionQ3_I(isnan(ConsumptionQ3_I)) = -1;
    ConsumptionQ4_I(isnan(ConsumptionQ4_I)) = -1;
    ConsumptionQ5_I(isnan(ConsumptionQ5_I)) = -1;
    ConsumptionT10_I(isnan(ConsumptionT10_I)) = -1;
    ConsumptionT1_I(isnan(ConsumptionT1_I)) = -1;

    IncomeB1_I = [std_simul_annual.I_I_1_1; std_simul_quarterly.I_I_1_1];
    IncomeB10_I = [std_simul_annual.I_I_1_10; std_simul_quarterly.I_I_1_10];
    IncomeQ1_I = [std_simul_annual.I_I_1_Q1; std_simul_quarterly.I_I_1_Q1];
    IncomeQ2_I = [std_simul_annual.I_I_1_Q2; std_simul_quarterly.I_I_1_Q2];
    IncomeQ3_I = [std_simul_annual.I_I_1_Q3; std_simul_quarterly.I_I_1_Q3];
    IncomeQ4_I = [std_simul_annual.I_I_1_Q4; std_simul_quarterly.I_I_1_Q4];
    IncomeQ5_I = [std_simul_annual.I_I_1_Q5; std_simul_quarterly.I_I_1_Q5];
    IncomeT10_I = [std_simul_annual.I_I_1_90; std_simul_quarterly.I_I_1_90];
    IncomeT1_I = [std_simul_annual.I_I_1_99; std_simul_quarterly.I_I_1_99];

    IncomeB1_I(isnan(IncomeB1_I)) = -1;
    IncomeB10_I(isnan(IncomeB10_I)) = -1;
    IncomeQ1_I(isnan(IncomeQ1_I)) = -1;
    IncomeQ2_I(isnan(IncomeQ2_I)) = -1;
    IncomeQ3_I(isnan(IncomeQ3_I)) = -1;
    IncomeQ4_I(isnan(IncomeQ4_I)) = -1;
    IncomeQ5_I(isnan(IncomeQ5_I)) = -1;
    IncomeT10_I(isnan(IncomeT10_I)) = -1;
    IncomeT1_I(isnan(IncomeT1_I)) = -1;

    IncomeB1_W = [std_simul_annual.I_1_1; std_simul_quarterly.I_1_1];
    IncomeB10_W = [std_simul_annual.I_1_10; std_simul_quarterly.I_1_10];
    IncomeQ1_W = [std_simul_annual.I_1_Q1; std_simul_quarterly.I_1_Q1];
    IncomeQ2_W = [std_simul_annual.I_1_Q2; std_simul_quarterly.I_1_Q2];
    IncomeQ3_W = [std_simul_annual.I_1_Q3; std_simul_quarterly.I_1_Q3];
    IncomeQ4_W = [std_simul_annual.I_1_Q4; std_simul_quarterly.I_1_Q4];
    IncomeQ5_W = [std_simul_annual.I_1_Q5; std_simul_quarterly.I_1_Q5];
    IncomeT10_W = [std_simul_annual.I_1_90; std_simul_quarterly.I_1_90];
    IncomeT1_W = [std_simul_annual.I_1_99; std_simul_quarterly.I_1_99];

    IncomeB1_W(isnan(IncomeB1_W)) = -1;
    IncomeB10_W(isnan(IncomeB10_W)) = -1;
    IncomeQ1_W(isnan(IncomeQ1_W)) = -1;
    IncomeQ2_W(isnan(IncomeQ2_W)) = -1;
    IncomeQ3_W(isnan(IncomeQ3_W)) = -1;
    IncomeQ4_W(isnan(IncomeQ4_W)) = -1;
    IncomeQ5_W(isnan(IncomeQ5_W)) = -1;
    IncomeT10_W(isnan(IncomeT10_W)) = -1;
    IncomeT1_W(isnan(IncomeT1_W)) = -1;


    rates_of_change_agg = table(Policy_Rule, Aggregation_Time, Output, AggConsumption, Inflation, AggIncome, NominalRate, RealRate, Wage, Unemployment, Profit);

    rates_of_change_C_W = table(Policy_Rule, Aggregation_Time, AggConsumption, ConsumptionT1_W, ConsumptionT10_W, ConsumptionQ1_W, ...
        ConsumptionQ2_W, ConsumptionQ3_W, ConsumptionQ4_W, ConsumptionQ5_W, ConsumptionB10_W, ConsumptionB1_W);

    rates_of_change_C_I = table(Policy_Rule, Aggregation_Time, AggConsumption, ConsumptionT1_I, ConsumptionT10_I, ConsumptionQ1_I, ...
        ConsumptionQ2_I, ConsumptionQ3_I, ConsumptionQ4_I, ConsumptionQ5_I, ConsumptionB10_I, ConsumptionB1_I);

    rates_of_change_I_W = table(Policy_Rule, AggIncome, IncomeT1_W, IncomeT10_W, IncomeQ1_W, ...
        IncomeQ2_W, IncomeQ3_W, IncomeQ4_W, IncomeQ5_W, IncomeB10_W, IncomeB1_W);

    rates_of_change_I_I = table(Policy_Rule, AggIncome, IncomeT1_I, IncomeT10_I, IncomeQ1_I, ...
        IncomeQ2_I, IncomeQ3_I, IncomeQ4_I, IncomeQ5_I, IncomeB10_I, IncomeB1_I);


    IRFs_R_cb_all = [IRFs(1).R_cb_IRFs_p];
    ZLB_freq = sum(IRFs_R_cb_all<=1+1e-8)./length(IRFs_R_cb_all);


    IRFs_pi_all = [IRFs(1).pi_IRFs_p];
    hipi5_freq = sum(IRFs_pi_all>=1.0125)./length(IRFs_pi_all);
    hipi7_freq = sum(IRFs_pi_all>=1.0175)./length(IRFs_pi_all);
    hipi10_freq = sum(IRFs_pi_all>=1.025)./length(IRFs_pi_all);


    Policy_Rule_freq = [spec_settings.policy_rule; spec_settings.policy_rule; spec_settings.policy_rule; spec_settings.policy_rule];
    Frequency_Category = ["High Inflation 5"; "High Inflation 7"; "High Inflation 10"; "ZLB"];
    Frequency = [hipi5_freq; hipi7_freq; hipi10_freq; ZLB_freq];

    freq_table = table(Policy_Rule_freq, Frequency_Category, Frequency);

    results.freq_table = freq_table;
    results.rates_table_agg = rates_of_change_agg;
    results.rates_table_C_W = rates_of_change_C_W;
    results.rates_table_C_I = rates_of_change_C_I;
    results.rates_table_I_W = rates_of_change_I_W;
    results.rates_table_I_I = rates_of_change_I_I;

    merged_table = table(Policy_Rule, Aggregation_Time, Output, AggConsumption, Inflation, AggIncome, NominalRate, RealRate, Wage, Unemployment, Profit, ...
        ConsumptionT1_W, ConsumptionT10_W, ConsumptionQ1_W, ...
        ConsumptionQ2_W, ConsumptionQ3_W, ConsumptionQ4_W, ConsumptionQ5_W, ConsumptionB10_W, ConsumptionB1_W, ...
        ConsumptionT1_I, ConsumptionT10_I, ConsumptionQ1_I, ...
        ConsumptionQ2_I, ConsumptionQ3_I, ConsumptionQ4_I, ConsumptionQ5_I, ConsumptionB10_I, ConsumptionB1_I, ...
        IncomeT1_W, IncomeT10_W, IncomeQ1_W, ...
        IncomeQ2_W, IncomeQ3_W, IncomeQ4_W, IncomeQ5_W, IncomeB10_W, IncomeB1_W, ...
        IncomeT1_I, IncomeT10_I, IncomeQ1_I, ...
        IncomeQ2_I, IncomeQ3_I, IncomeQ4_I, IncomeQ5_I, IncomeB10_I, IncomeB1_I);

    results.merged_table = merged_table;

    results.param = param;
    results.X_series = X_series;

    tag = strcat(spec_settings.simul_string);
    file_label = strcat(spec_settings.sim_save_folder,tag,'_',sh_type, '_phipi=', strrep(num2str(param.phi_pi), ".", "_"), "_kappa=", strrep(num2str(round(param.kappa,3)), ".", "_"));

    save(file_label, 'results');

    clear results simul_results merged_table freq_table


end