function out = find_f(in,uminus,u,param,grid)

f = exp(in);

P_SS2_reduc_aux = [1-param.in,0,param.in;0,1-param.in,param.in;0,param.out,1-param.out];

dist_minus = [(1-uminus)*(1-param.Eshare);uminus*(1-param.Eshare);param.Eshare];	

dist = dist_minus(:)'*P_SS2_reduc_aux;

dist = dist(:);

P_SS3_reduc_aux   = kron([1-param.lambda+param.lambda*f,param.lambda*(1-f);f,1-f],eye(1));
P_SS3_reduc_aux   = [P_SS3_reduc_aux, zeros(3-1,1)];
lastrow     = [zeros(1,2),1];
P_SS3_reduc_aux   = [P_SS3_reduc_aux; lastrow];

dist_tilde_aux = dist(:)'*P_SS3_reduc_aux;

u_aux = (dist_tilde_aux(2))./(dist_tilde_aux(1)+dist_tilde_aux(2));

out = u - u_aux;