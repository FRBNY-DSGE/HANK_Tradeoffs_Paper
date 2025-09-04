% Tmax = 5;


% At Tmax+1, the model is at the reference regime, so decrulea and decruleb are coeffecients of the linearized solution
% At Tmax,   the model is at the alternative regime

Ps = zeros(nvars,nvars,Tmax);
Ds = zeros(nvars,Tmax);


invmat       = inv((Astarbarmat*decrulea+Bstarbarmat));
Ps(:,:,Tmax) = -invmat*Cstarbarmat;
Ds(:,Tmax)   = -invmat*D_ZLB;
    
for ii = Tmax-1:-1:1
        

    invmat      = inv(Astarbarmat*Ps(:,:,ii+1)+Bstarbarmat);
    Ps(:,:,ii)  = -invmat*Cstarbarmat;
    Ds(:,ii)    = -invmat*(Astarbarmat*Ds(:,ii+1)+D_ZLB);
        
end
 
E = -invmat*J_ZLB;
       
P_1 = Ps(:,:,1);
D_1 = Ds(:,1);
E_1 = E;
    