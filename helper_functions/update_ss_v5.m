
SS_stats_base = SS_stats;
param_base    = param;

param.iota    = param_update.iota;
param.theta_b = param_update.theta_b;

% 1) iota

param.fix_L   = (SS_stats.J_bar1-param.iota*SS_stats.V/SS_stats.M)/SS_stats.J_bar2;

SS_stats.J     = inv(eye(grid.ns)-param.beta_L*(1-param.death_rate)*(1-param.lambda)*(1-param.in)*param.P_SS)*((SS_stats.r_l-param.fix_L-param.w_bar)*param.n*(grid.s'));  % value of matching
SS_stats.J_bar = (SS_stats.J')*grid.s_dist;     

SS_stats.Profit_L = (SS_stats.r_l-param.fix_L-param.w_bar)*SS_stats.L - param.iota*SS_stats.V;

% % 2) varphi

param.fix   = SS_stats.Y - SS_stats_base.mc*SS_stats_base.Y + SS_stats.Profit_L + SS_stats.Profit_K - SS_stats_base.Profit;

SS_stats.Profit = SS_stats.Y - SS_stats_base.mc*SS_stats_base.Y + SS_stats.Profit_L + SS_stats.Profit_K - param.fix;

SS_stats.Profit_int = (1-SS_stats.mc)*SS_stats.Y - param.fix;

% 3) theta_b

param.omega = (1-param.theta_b*((SS_stats_base.R_a-SS_stats_base.R)*SS_stats_base.leverage+SS_stats_base.R))/SS_stats_base.leverage;

SS_stats.vv   	  = ((1-param.theta_b)*SS_stats.Lambda_b*(SS_stats.R_a-param.b_a_aux2*SS_stats.R))/(1-param.theta_b*SS_stats.Lambda_b*SS_stats.xx);
SS_stats.ee   	  = (1-param.theta_b)*SS_stats.Lambda_b*param.b_a_aux2*SS_stats.R/(1-param.theta_b*SS_stats.Lambda_b*SS_stats.zz);

param.DELTA = SS_stats.vv + SS_stats.ee/SS_stats_base.leverage;

SS_stats.Profit_FI = (1-param.theta_b)*((SS_stats_base.R_a-SS_stats_base.R)*SS_stats_base.leverage+SS_stats_base.R)*SS_stats_base.NW_b-param.omega*param.q*SS_stats_base.A_b;  % Profits of FI distributed to the bankers family from existing bankers

param.fix2 = SS_stats.Profit_FI;

SS_stats.AvgC = (param.delta_ss*grid.K + param.w_bar*grid.L + param.iota*SS_stats.V + param.fix)/SS_stats.Y;

param.fix_ratio = param.fix/SS_stats.Y;