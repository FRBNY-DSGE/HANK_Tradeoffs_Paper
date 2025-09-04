function [c_a_star,b_a_star,a_a_star,c_n_star,b_n_star] = policies_update(EVb,EVa,Qminus,PIminus,R_cbminus,UU,KK,inc,meshes,grid,param)

%% EGM Step 1: Non-adjustment case
EMU     = param.beta_aux*UU*(1-param.death_rate)*reshape(EVb,[grid.nb grid.na grid.nse]);
%  Linear interpolation (taking into account the financial shock)
% EMU_aux = permute(interp1(grid.a,permute(EMU,[2 1 3]),min(KK,1)*grid.a),[2 1 3]);
EMU_aux = EMU;

c_new = 1./(EMU_aux.^(1/param.sigma));

b_star_n = (c_new + meshes.b - inc.labor - inc.dividend - inc.transfer);
b_star_n = (b_star_n<0).*b_star_n./((R_cbminus+param.Rprem)/PIminus) ...
			+ (b_star_n>=0).*b_star_n./(R_cbminus/PIminus);
b_star_n = (1-param.death_rate)*b_star_n/param.b_a_aux;						

binding_constraints = meshes.b < repmat(b_star_n(1,:,:),[grid.nb 1 1]);

Resource 	 = inc.labor + inc.dividend + inc.bond + inc.transfer;
Resource_aux = reshape(Resource,[grid.nb, grid.na*grid.nse]); 

b_star_n = reshape(b_star_n,[grid.nb grid.na*grid.nse]);
c_n_aux  = reshape(c_new,[grid.nb grid.na*grid.nse]);

% Check monotonicity of b_star_n
if max(sum(abs(diff(sign(diff(b_star_n))))))~=0
    warning('non monotone future liquid asset choice encountered')
end			

c_update = zeros(grid.nb, grid.na*grid.nse);
b_update = zeros(grid.nb, grid.na*grid.nse);

for hh = 1:grid.na*grid.nse
	% Savings        = griddedInterpolant(b_star_n(:,hh),grid.b);
	% b_update(:,hh) = Savings(grid.b);
	Consumption    = griddedInterpolant(b_star_n(:,hh),c_n_aux(:,hh));
	c_update(:,hh) = Consumption(grid.b);
	b_update(:,hh) = Resource_aux(:,hh) - c_update(:,hh);
end

c_n_star = reshape(c_update,[grid.nb, grid.na, grid.nse]);
b_n_star = reshape(b_update,[grid.nb, grid.na, grid.nse]);

c_n_star(binding_constraints) = Resource(binding_constraints) - grid.b(1);
b_n_star(binding_constraints) = min(grid.b);

b_n_star(b_n_star>grid.b(end)) = grid.b(end);
% c_n_star(b_n_star>grid.b(end)) = Resource(b_n_star>grid.b(end)) - grid.b(end);

%% EGM Step 2: Find optimal portfolio compositions

term1 = param.beta_aux*UU*(1-param.death_rate)*reshape(EVa,[grid.nb,grid.na,grid.nse]);

E_return_diff = term1./Qminus - EMU;

% Check quasi-monotonicity of E_return_diff
if max(sum(abs(diff(sign(E_return_diff)))))>2
    warning('multiple roots of portfolio choice encountered')
end

% Find b_a* for given k' that solves the difference equation

b_a_aux = Fastroot(grid.b,E_return_diff);
b_a_aux = max(b_a_aux,grid.b(1));
b_a_aux = min(b_a_aux,grid.b(end));
b_a_aux = reshape(b_a_aux,[grid.na grid.nse]);   % b(a')

%% EGM Step 3: Construct cons_list, res_list, bond_list
EMU = reshape(EMU,[grid.nb, grid.na*grid.nse]);

[~,idx] = histc(b_a_aux,grid.b); 
idx(b_a_aux<=grid.b(1))   = 1;             % if below minimum
idx(b_a_aux>=grid.b(end)) = grid.nb-1;     % if above minimum
step = diff(grid.b);                       % stepzie on grid
s    = (b_a_aux - grid.b(idx))./step(idx); % compute the slope

aux_index = (0:(grid.na*grid.nse)-1)*grid.nb;  % aux for linear indices
aux3      = EMU(idx(:)+aux_index(:));          % vectorized EMU

% interpolate EMU
EMU_star        = aux3 + s(:).*(EMU(idx(:)+ aux_index(:)+1)-aux3);
c_a_aux         = 1./(EMU_star.^(1/param.sigma));   % c(k')
cap_expenditure = shiftdim(inc.capital(1,:,:))/KK;     % size na * nse
auxL            = shiftdim(inc.labor(1,:,:));          % size na * nse
auxLT           = shiftdim(inc.transfer(1,:,:));       % size na * nse

% Resources that lead to capital choice k'
% = c + b*(a') + q*a' - w*s
Resource = c_a_aux + b_a_aux(:) + cap_expenditure(:) - auxL(:) - auxLT(:);

c_a_aux  = reshape(c_a_aux, [grid.na grid.nse]);
Resource = reshape(Resource, [grid.na grid.nse]);

% 3.1) Money constraint is not binding, but capital constraint is binding
b_star_zero = shiftdim(b_a_aux(1,:));   % b*(a'=0)

% Use c_n(a'=0) from constrained problem, when b' is on grid

aux_c     = reshape(c_new(:,1,:),[grid.nb grid.nse]);
aux_inc   = reshape(inc.labor(1,1,:)+inc.transfer(1,1,:),[1 grid.nse]);
cons_list = cell(grid.nse,1);
res_list  = cell(grid.nse,1);
b_list    = cell(grid.nse,1);
a_list    = cell(grid.nse,1);

for j=1:grid.nse
	% When choosing a' = 0, HHs might still want to choose b' smaller thatn b*(a'=0)
	if b_star_zero(j)>grid.b(1)
		log_index    = grid.b<b_star_zero(j);
		c_a_cons     = aux_c(log_index,j);  
		cons_list{j} = c_a_cons;
		res_list{j}  = c_a_cons + grid.b(log_index)' - aux_inc(j);
		b_list{j}    = grid.b(log_index)';
		a_list{j}    = zeros(sum(log_index),1);
	end
end

% Merge lists
c_a_aux  = reshape(c_a_aux,[grid.na grid.nse]);
b_a_aux  = reshape(b_a_aux,[grid.na grid.nse]);
Resource = reshape(Resource,[grid.na grid.nse]);

for j=1:grid.nse
	cons_list{j} = [cons_list{j}; c_a_aux(:,j)];
	res_list{j}  = [res_list{j};  Resource(:,j)];
	b_list{j}    = [b_list{j};    b_a_aux(:,j)];
	a_list{j}    = [a_list{j};    grid.a'];
end

%% EGM Step 4: Inteolate to fixed grid
c_a_star = NaN([grid.nb*grid.na grid.nse]);
b_a_star = NaN([grid.nb*grid.na grid.nse]);
a_a_star = NaN([grid.nb*grid.na grid.nse]);

Resource_grid  = reshape(inc.capital+inc.bond+inc.dividend,[grid.nb*grid.na, grid.nse]);
labor_inc_grid = reshape(inc.labor + inc.transfer,[grid.nb*grid.na grid.nse]);

for j = 1:grid.nse
	log_index = Resource_grid(:,j)<res_list{j}(1);

	% Check monotonicity of resources
	if max(sum(abs(diff(sign(diff(res_list{j}))))))~=0
        warning('non monotone resource list encountered')
    end

    cons          = griddedInterpolant(res_list{j},cons_list{j});
    c_a_star(:,j) = cons(Resource_grid(:,j));
    bond          = griddedInterpolant(res_list{j},b_list{j});
    b_a_star(:,j) = bond(Resource_grid(:,j));
    equity        = griddedInterpolant(res_list{j},a_list{j});
    a_a_star(:,j) = equity(Resource_grid(:,j));

    % Lowest value of res_list corresponds to b_a'=0 and a_a'=0
    % If resources are smaller than the lowest value of res_list, then HHs consume everything that they have

    c_a_star(log_index,j) = Resource_grid(log_index,j) + labor_inc_grid(log_index,j)-grid.b(1);
    b_a_star(log_index,j) = grid.b(1);
    a_a_star(log_index,j) = 0;
end

c_a_star = reshape(c_a_star,[grid.nb, grid.na, grid.nse]);
b_a_star = reshape(b_a_star,[grid.nb, grid.na, grid.nse]);
a_a_star = reshape(a_a_star,[grid.nb, grid.na, grid.nse]);

b_a_star(b_a_star>grid.b(end)) = grid.b(end);
a_a_star(a_a_star>grid.a(end)) = grid.a(end);

% c_a_star = inc.labor + inc.capital + inc.dividend + inc.bond - Qminus*a_a_star - b_a_star;

end

function roots = Fastroot(xgrid,fx)
%fast linear interpolation root finding
%(=one Newton step at largest negative function value)
%   stripped down version of interp1 that accepts multiple inputs (max 3)
%   that are interpolated over the same grids x & xi
xgrid=xgrid(:);
fx=reshape(fx,[numel(xgrid),numel(fx)/numel(xgrid)]);

dxgrid=diff(xgrid);
dfx=diff(fx);
idx=ones(1,numel(fx)/numel(xgrid));

% Make use of the fact that the difference equation is monotonically
% increasing in m
idx_min=(fx(1,:)>0); %Corner solutions left (if no solution x* to f(x)=0 exists)
idx_max=(fx(end,:)<0); %Corner solutions right (if no solution x* to f(x)=0 exists)
index=find(and(not(idx_min),not(idx_max))); %interior solutions (if solution x* to f(x)=0 exists)

% Find index of two gridpoints where sign of fx changes from positive to negative,
[~,idx(index)]=max(diff(sign(fx(:,index))));

aux_index  = (0:numel(fx)/numel(xgrid)-1)*numel(xgrid); %aux for linear indexes
aux_index2 = (0:numel(fx)/numel(xgrid)-1)*(numel(xgrid)-1);
fxx  = fx(idx+aux_index);
xl   = xgrid(idx)';
dx   = dxgrid(idx)';
dfxx = dfx(idx+aux_index2);
% Because function is piecewise linear in gridpoints, one newton step is
% enough to find the solution
roots = xl-fxx.*dx./dfxx;

roots(idx_min)=xgrid(1);   % constrained choice
roots(idx_max)=xgrid(end); % no-extrapolation
end
