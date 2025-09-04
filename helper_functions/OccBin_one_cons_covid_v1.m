
Ps_afterGR = zeros(nvars,nvars,Tmax);
Ds_afterGR = zeros(nvars,Tmax);

invmat               = inv((Astarbarmat*decrulea_afterGR+Bstarbarmat));
Ps_afterGR(:,:,Tmax) = -invmat*Cstarbarmat;
Ds_afterGR(:,Tmax)   = -invmat*D_ZLB;
    
for ii = Tmax-1:-1:1
        
    invmat              = inv(Astarbarmat*Ps_afterGR(:,:,ii+1)+Bstarbarmat);
    Ps_afterGR(:,:,ii)  = -invmat*Cstarbarmat;
    Ds_afterGR(:,ii)    = -invmat*(Astarbarmat*Ds_afterGR(:,ii+1)+D_ZLB);
        
end

E_afterGR = -invmat*J_ZLB;
       
P_1 = Ps_afterGR(:,:,1);
D_1 = Ds_afterGR(:,1);
E_1 = E_afterGR;
    
