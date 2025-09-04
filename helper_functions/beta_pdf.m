function out = beta_pdf(x, mean, std)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  Description : This code gives probability density of gamma at x
%  with given mean and variance
%
%  Written by Donggyu Lee
%  Date : 28th Dec 2016
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

var = std^2; % variance

k   = (1-mean)/mean; % auxililary variable

alpha = (k-(1+k)^2*var)/((1+k)^2*var*(1+k)); % first shape parameter
beta  = k*alpha;                             % second shape parameter

out = betapdf(x, alpha, beta); % calculate density of beta distribution by using matlab built-in function




