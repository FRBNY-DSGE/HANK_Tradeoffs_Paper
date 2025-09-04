function [P_ref,Q_ref,cof_ref,J_ref] = compute_lin_sol(hx,gx,F1_aux,F2_aux,F3_aux,F4_aux,grid)


F1 = [F1_aux(1:grid.numstates_endo,:);F1_aux(grid.numstates+1:end,:)];
F2 = [F2_aux(1:grid.numstates_endo,:);F2_aux(grid.numstates+1:end,:)];
F3 = [F3_aux(1:grid.numstates_endo,:);F3_aux(grid.numstates+1:end,:)];
F4 = [F4_aux(1:grid.numstates_endo,:);F4_aux(grid.numstates+1:end,:)];

F11 = F1(:,1:grid.numstates_endo);
F12 = F1(:,grid.numstates_endo+1:end);

F31 = F3(:,1:grid.numstates_endo);
F32 = F3(:,grid.numstates_endo+1:end);

% Construct a coefficient matrix

A_ref = [zeros(grid.numstates_endo+grid.numcontrols,grid.numstates_endo), F2];
B_ref = [F11, F4];
C_ref = [F31, zeros(grid.numstates_endo+grid.numcontrols,grid.numcontrols)];
E_ref = F32;


% Constuct P, Q

hx_aux = hx(1:grid.numstates_endo,:);
HH_aux = [hx_aux;gx];
H1_aux = HH_aux(:,1:grid.numstates_endo);
H2_aux = HH_aux(:,end-grid.numstates_shocks+1:end);

P_ref  = [H1_aux,zeros(grid.numstates_endo+grid.numcontrols,grid.numcontrols)];
Q_ref  = H2_aux;

cof_ref = [C_ref,B_ref,A_ref];
J_ref   = E_ref;	