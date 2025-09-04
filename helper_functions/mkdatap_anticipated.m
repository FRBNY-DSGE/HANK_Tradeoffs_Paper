function [zdata, Tmax]=mkdatap_anticipated(nperiods,decrulea,decruleb,cof,Jbarmat,cofstar,Jstarbarmat,Dstarbarmat,...
    regime,regimestart,violvecbool,grid,irfshock,scalefactormod,init)

nvars = grid.num_endo;

if nargin<15
	init=zeros(nvars,1);
end

if nargin<14
	scalefactormod=1;
end

nshocks  = grid.numstates_shocks;

for i = 1:nshocks
    shockpos = strmatch(irfshock(i,:),grid.exog_names,'exact');
    if ~isempty(shockpos)
        irfshockpos(i) = shockpos;
    else
        error(['Shock ',irfshock(i,:),' is not in the model']);
    end
end

nregimes = length(regime); % number of regimes in the periods considered

Cbarmat = cof(:,1:nvars);
Bbarmat = cof(:,nvars+1:2*nvars);
Abarmat = cof(:,2*nvars+1:3*nvars);

% cofstar contains the system for the model when the constraint binds

Cstarbarmat = cofstar(:,1:nvars);
Bstarbarmat = cofstar(:,nvars+1:2*nvars);
Astarbarmat = cofstar(:,2*nvars+1:3*nvars);

% get the time-dependent decision rules

Tmax = regimestart(nregimes)-1;  % Tmax is the position of the last period when the constraint binds

if Tmax > 0

	% At Tmax+1, the model is at the reference regime, so decrulea and decruleb are coeffecients of the linearized solution
	% At Tmax,   the model is at the alternative regime

	P = zeros(nvars,nvars,Tmax);
	D = zeros(nvars,Tmax);

	invmat      = inv((Astarbarmat*decrulea+Bstarbarmat));
	P(:,:,Tmax) = -invmat*Cstarbarmat;
	D(:,Tmax)   = -invmat*Dstarbarmat;

	% equivalent to pre-multiplying by the inverse above if the target
    % matrix is invertible. Otherwise it yields the minimum state solution
    % P(:,:,Tmax) = -(Astarbarmat*decrulea+Bstarbarmat)\Cstarbarmat;
    % D(:,Tmax)   = -(Astarbarmat*decrulea+Bstarbarmat)\Dstarbarmat;

    % Iterate backwardly to the first period

    for i = Tmax-1:-1:1

    	if violvecbool(i) % if violvecbool(i) == 1, then the model is still at the alternative regime
    		invmat   = inv(Astarbarmat*P(:,:,i+1)+Bstarbarmat);
    		P(:,:,i) = -invmat*Cstarbarmat;
    		D(:,i) = -invmat*(Astarbarmat*D(:,i+1)+Dstarbarmat);
    	else % the model is at the reference regime
    		invmat   = inv(Abarmat*P(:,:,i+1)+Bbarmat);
    		P(:,:,i) = -invmat*Cbarmat;
    		D(:,i) = -invmat*(Abarmat*D(:,i+1));
    	end

    end

    % Check the first perid. At the first period, there is a shock

    if Tmax > 1
    	if violvecbool(1) % at the first period, the alternative regime applies
    		E = -invmat*Jstarbarmat;
    	else
    		E = -invmat*Jbarmat;
    	end
    else % Tmax == 1
    	invmat = inv(Astarbarmat*decrulea+Bstarbarmat);
    	E      = -invmat*Jstarbarmat;
    end

end

% generate data
% history will contain data, the state vector at each period in time will be stored columnwise

history      = zeros(nvars,nperiods+1);
history(:,1) = init;
errvec       = zeros(grid.numstates_shocks,1);

% deal with predetermined conditions
for i = 1:nshocks
	errvec(irfshockpos(i)) = scalefactormod(i);
end

% deal with shocks
irfpos = 1;
if irfpos <= Tmax
	history(:,irfpos+1) = P(:,:,irfpos)*history(:,irfpos)+...
		D(:,irfpos) + E*errvec;
else
	history(:,irfpos+1) = decrulea*history(:,irfpos) + decruleb*errvec;
end

% all other periods
for irfpos=2:nperiods+1
	if irfpos <= Tmax
			history(:,irfpos+1) = P(:,:,irfpos)*history(:,irfpos)+...
		D(:,irfpos);
	else
		history(:,irfpos+1) = decrulea*history(:,irfpos);
	end
end

history = history';
zdata   = history(2:end,:); %zdata   = history(2:end,:);

end