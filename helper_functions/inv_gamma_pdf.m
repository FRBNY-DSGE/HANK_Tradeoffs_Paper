function out = inv_gamma_pdf(x,m,sig)

% m    : mean
% sig  : standard deviation

alpha_aux = m.^2./sig.^2 + 2;    % location parameter
% alpha_aux = m.^2 + 2;
beta_aux  = m*(alpha_aux-1);     % shape parameter
% beta_aux  = beta_aux/100;

out = beta_aux.^alpha_aux./gamma(alpha_aux).*x.^(-alpha_aux-1).*exp(-beta_aux./x);