function out = compute_CE(input,param,grid)

util     = @(c)  (c.^(1-param.sigma))./(1-param.sigma);
mutil    = @(c)  (1./(c.^param.sigma));
invutil  = @(u)  (((1-param.sigma).*u).^(1/(1-param.sigma)));
invmutil = @(mu) ((1./mu).^(1/param.sigma));

nn                 = input.nn;
N                  = input.N;
W 	               = input.W;
tau_w              = input.tau_w;
tau_a              = input.tau_a;
PROFIT             = input.PROFIT;
C_b                = input.C_b;
Ntilde             = input.Ntilde;
R_A                = input.R_A;
Q                  = input.Q;
R_cbminus          = input.R_cbminus;
R_cb               = input.R_cb;
PI                 = input.PI;
PInext             = input.PInext;
Va_next_input      = input.Va_next_input;
P_transition_exp2  = input.P_transition_exp2;
P_transition_exp   = input.P_transition_exp;
mutil_c_next_input = input.mutil_c_next_input;
MU_tilde           = input.MU_tilde;
UU                 = input.UU;
VALUEnext          = input.VALUE_next_input;
VALUEnext2         = input.VALUE_next_input2;
LT                 = input.LT;
vartheta           = input.vartheta;
AA_aux             = input.AA_aux;
Profit_FI          = input.Profit_FI;
D                  = input.D; 

Eshare = param.Eshare;

[meshes.b,meshes.a,meshes.se] = ndgrid(grid.b,grid.a,grid.se);

meshaux = meshes;

[~,~,meshaux.se_aux] = ndgrid(grid.b,grid.a,1:grid.nse);

auxWW                    = ones(grid.nb,grid.na,grid.nse);
auxWW(:,:,grid.ns+1:end) = zeros(grid.nb,grid.na,grid.ns+1);
auxWW1                   = auxWW;
auxWW                    = auxWW.*(meshes.se.^(param.b_aux-1));

auxWW2                    = ones(grid.nb,grid.na,grid.nse);
auxWW2(:,:,1:grid.ns)     = zeros(grid.nb,grid.na,grid.ns);
auxWW2(:,:,end)           = zeros(grid.nb,grid.na,1);

R_A_aux     = R_A + (param.death_rate/(1-param.death_rate))*(Q);

NW          = param.xi/(1+param.xi)*nn*W;
WW          = (1-tau_w)*(NW*auxWW1 + param.xi/(1+param.xi)*1*param.w_bar.*auxWW2  + param.b_share*PROFIT/(Ntilde*grid.s2_bar).*auxWW);
WW(:,:,end) = (1-tau_a)*(param.Eratio*PROFIT+Profit_FI-param.fix2)/Eshare*ones(grid.nb,grid.na);
inc.labor    = WW.*(meshes.se);
inc.dividend = meshes.a*R_A_aux;
inc.capital  = meshes.a*Q;
inc.bond     = (R_cbminus/PI).* meshes.b ...
				+ (meshes.b<0).*(param.Rprem/PI).* meshes.b;
inc.bond     = inc.bond/(1-param.death_rate)*param.b_a_aux;					
inc.transfer = (LT+C_b)/(1+param.frac_b)*ones(grid.nb,grid.na,grid.nse);	

param.beta_aux = D*param.beta;					

Vanext_aux1 = permute(interp1(grid.a,permute(reshape(Va_next_input, [grid.nb grid.na grid.nse]),[2 1 3]),min(vartheta,1)*grid.a),[2 1 3]);
Vanext_aux1 = Vanext_aux1(:);

Vanext_aux2 = interp1(grid.b,reshape(Vanext_aux1, [grid.nb grid.na grid.nse]),AA_aux*grid.b,'linear','extrap');
Vanext_aux2 = Vanext_aux2(:);

Vanext_aux  = Vanext_aux2;

EVa     = vartheta*reshape(reshape(Vanext_aux,[grid.nb*grid.na grid.nse])*P_transition_exp2',[grid.nb grid.na grid.nse]);
R_cbaux	= R_cb/PInext + (meshes.b<0).*(param.Rprem/PInext); 


mutil_cnext_aux1 = permute(interp1(grid.a,permute(reshape(mutil_c_next_input, [grid.nb grid.na grid.nse]),[2 1 3]),min(vartheta,1)*grid.a),[2 1 3]);
mutil_cnext_aux1 = mutil_cnext_aux1(:);

mutil_cnext_aux2 = interp1(grid.b,reshape(mutil_cnext_aux1, [grid.nb grid.na grid.nse]),AA_aux*grid.b,'linear','extrap');
mutil_cnext_aux2 = mutil_cnext_aux2(:);

mutil_cnext_aux  = mutil_cnext_aux2;

EVb     = AA_aux*reshape(reshape(R_cbaux(:)/(1-param.death_rate)*param.b_a_aux.*mutil_cnext_aux,[grid.nb*grid.na grid.nse])*P_transition_exp2',[grid.nb grid.na grid.nse]);

[c_a_star,b_a_star,a_a_star,c_n_star,b_n_star] = policies_update(EVb,EVa,Q,PI,R_cbminus,UU,1,inc,meshes,grid,param);			

EV         = reshape(VALUEnext ,[grid.nb*grid.na grid.nse])*P_transition_exp2';
EV_actual  = reshape(VALUEnext2,[grid.nb*grid.na grid.nse])*P_transition_exp';

VALUE_next        = griddedInterpolant(meshaux.b,meshaux.a,meshaux.se_aux,reshape(EV,[grid.nb grid.na grid.nse]));
VALUE_actual_next = griddedInterpolant(meshaux.b,meshaux.a,meshaux.se_aux,reshape(EV_actual,[grid.nb grid.na grid.nse]));

V_adjust   = util(c_a_star) + param.beta_aux*UU*(1-param.death_rate)*VALUE_next(AA_aux*b_a_star,vartheta*a_a_star ,meshaux.se_aux);
V_noadjust = util(c_n_star) + param.beta_aux*UU*(1-param.death_rate)*VALUE_next(AA_aux*b_n_star,vartheta*meshaux.a,meshaux.se_aux);

V_a_actual  = util(c_a_star) + param.beta_aux*UU*(1-param.death_rate)*VALUE_actual_next(b_a_star,a_a_star ,meshaux.se_aux);
V_na_actual = util(c_n_star) + param.beta_aux*UU*(1-param.death_rate)*VALUE_actual_next(b_n_star,meshaux.a,meshaux.se_aux);

mu_chi = param.mu_chi;

AProb = max(min(1./(1+exp(-((V_adjust-V_noadjust)-mu_chi)./param.sigma_chi)),1-1e-6),1e-6);
AC    = param.sigma_chi.*((1-AProb(:)).*log(1-AProb(:)) + AProb(:).*log(AProb(:)))...
    + mu_chi.*AProb(:) -  (mu_chi-param.sigma_chi*log(1+exp(mu_chi/param.sigma_chi)));
AProb = AProb.*param.nu;
AC    = AC.*param.nu;

VALUE_aux  = AProb(:).*V_adjust(:) + (1-AProb(:)).*V_noadjust(:) - AC(:);

out.VALUE_aux2  = AProb(:).*V_a_actual(:) + (1-AProb(:)).*V_na_actual(:);
out.AC          = AC;
out.VALUE_aux3  = AProb(:).*V_a_actual(:) + (1-AProb(:)).*V_na_actual(:) - AC(:);

mutil_c_n   = mutil(c_n_star);                         % marginal utility at consumption policy no adjustment
mutil_c_a   = mutil(c_a_star);                         % marginal utility at consumption policy adjustment
mutil_c_aux = AProb.*mutil_c_a + (1-AProb).*mutil_c_n; % Expected marginal utility at consumption policy (w &w/o adjustment)

Va_next  = griddedInterpolant(meshaux.b,meshaux.a,meshaux.se_aux,reshape(EVa,[grid.nb grid.na grid.nse]));
Va_aux   = AProb.*((R_A_aux+Q)).*mutil_c_a + (1-AProb).*((R_A_aux)).*mutil_c_n + param.beta_aux*UU*(1-param.death_rate).*(1-AProb).*Va_next(b_n_star,meshaux.a,meshaux.se_aux);

out.VALUE   = VALUE_aux;
out.mutil_c = mutil_c_aux;
out.Va      = Va_aux;

C_hh = q_cons2(c_a_star,c_n_star,W,nn,AProb,MU_tilde,meshes,param,grid,tau_w);

A_hh = sum(sum(sum((AProb.*MU_tilde).*a_a_star))) + sum(sum(sum(((1-AProb).*MU_tilde).*meshes.a)));

B_hh = sum(sum(sum((AProb.*MU_tilde).*b_a_star))) + sum(sum(sum(((1-AProb).*MU_tilde).*b_n_star)));

out.A = (1-param.death_rate)*A_hh;
out.B = (1-param.death_rate)*B_hh;
out.C = C_hh;


end

