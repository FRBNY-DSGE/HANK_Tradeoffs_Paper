function out = gamma_pdf(x, mean, std)


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

theta = var/mean;   % scale parameter (shape?)
k     = mean/theta; % shape parameter (location?)



out = gampdf(x, k, theta); % calculate density of gamma distribution by using matlab built-in function




