eps_t_aux = eps_t;
SHOCKS  = [ [eps_t_aux'] ; zeros(49,length(eps_t))];

% set some initial conditions and loop through the shocks period by period

shockssequence = SHOCKS;
nperiods = length(SHOCKS)+1;

nshocks = size(shockssequence,1);
init    = X_Taylor;


if ~exist('maxiter')
    maxiter = 30;
end

if ~exist('nperiods')
    nperiods = 50;
end

init_orig      = init;
zdatapiecewise = zeros(nperiods,grid.num_endo);
violvecbool    = zeros(nperiods+1,1);
track_Tmax = 0;

for ishock = 1: nshocks

    %ishock

	changes = 1;
	iter    = 0;

	while (changes && iter < maxiter)

		iter = iter + 1;

		% analyze when each regime starts based on current guess

		[regime regimestart] = map_regime(violvecbool);

		% get the hypothesized piece wise linear solution

		[zdatalinear, Tmax] = mkdatap_anticipated(nperiods,P_ref, Q_ref,...
            cof_ref,J_ref,cof_ZLB,J_ZLB,D_ZLB,...
            regime,regimestart,violvecbool,...
            grid,irfshock,shockssequence(ishock,:),init);

        if iter == 2
            track_Tmax = Tmax;
        end
  		r_difference = zdatalinear(:,grid.nb+grid.na+grid.nse-3+1);
        newviolvecbool  = exp(r_difference)*param.R_cb <= param.R_ZLB + 1e-8;
        relaxconstraint = exp(r_difference)*param.R_cb  > param.R_ZLB + 1e-8;

  		% check if changes to the hypothesis of the duration for each regime
  		% If a guess for regimes in each period has been changed or the consraint is slack in a period in which it is guessed not previously

  		if (max(newviolvecbool - violvecbool>0)) || sum(relaxconstraint(find(violvecbool==1))>0)
  			changes = 1;
  		else
  			changes = 0;
  		end

  		% Update the guess
  		
  		violvecbool = (violvecbool|newviolvecbool) - (relaxconstraint & violvecbool);

	end

	init                     = zdatalinear(1,:);
	zdatapiecewise(ishock,:) = init;
	init   				           = init';

    

	% reset violvecbool_ for next period's shock -- this resetting is 
    % consistent with expecting no additional shocks
    violvecbool = [violvecbool(2:end);0];

end

zdatapiecewise(ishock+1:end,:)=zdatalinear(2:nperiods-ishock+1,:);

% get the linear responses
%zdatalinear = mkdata(max(nperiods,size(shockssequence,1)),...
%                  P_ref,Q_ref,grid,...
%                  irfshock,shockssequence,init_orig);


if changes ==1
    display('Did not converge -- increase maxiter_')
end 



