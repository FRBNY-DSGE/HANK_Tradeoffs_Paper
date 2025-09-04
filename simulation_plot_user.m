%%%%%%%%%% Main script to make all plots:

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% simul_plot.m
%
% Description: Given simulation results, file recreates
% (1) aggregate and distributional variable volatility vs phi_pi values
% values ('standard deviation plots' in paper|'phi-pi' plots in our code)
% (2) disaggregation of aggregate and distributional volatility vs phi_pi
% parameter by supply and demand shocks
% (3) aggregate and distributional variable volatility vs inflation
% volatility ('frontier plots' in the paper)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all;
close all;
clc;

addpath(genpath('auxfiles'))
addpath(genpath('helper_functions'))
addpath(genpath('plotting_functions'))
addpath(genpath('paper_figures_replication'))

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% (1) General Settings:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%(1.1) Folders and save files
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Sub-root folder for plots
plt_folder = 'paper_figures_replication/Simulations/';

% Root folder for simulation results
results_folder = 'paper_simulations';

% Supply shock series string (default is 'shock_series_combined_hk_5')
simul_str = 'shock_series_combined_hk_5';


%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%(1.2) Plots
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pol1 = 'Taylor';



% Plot toggles
plot_main_figures = false;
plot_appendix_frontier = false;
plot_appendix_volatility = true;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%(1.3) Kappa Changeables
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

plt_labels = ["Baseline", "High Kappa"];

% Iterables for simulation results under both policy rules (for kappa
% baseline and high kappa results)
baseline_kappa_val = "0_053";
high_kappa_val = "0_1";

% Kappa parameter value presets for each policy rule
kappa_val1 = baseline_kappa_val;
kappa_val2 = high_kappa_val;


load('simul_Taylor_environment.mat');

% Iterables for simulation results under Taylor policy rule (by phi_pi
% values)
phipi_Taylor = ["1_5", "1_8533", "2", "2_5", "3", "3_5", "4"];
phipi_Taylor_vals = [1.5, 1.8533, 2.0, 2.5, 3.0, 3.5, 4.0];
baseline_phipi_T = 1.8533;
%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%(1.1) Initialize Settings
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

phipi1 = phipi_Taylor;
phipi1_vals = phipi_Taylor_vals;
baseline_phi_ind(1) = find(phipi_Taylor_vals == baseline_phipi_T);

baseline_phi_ind(2) = find(phipi_Taylor_vals == baseline_phipi_T);

%%%%%%%%%%%%%%%%%%%%%%%%%
% (2) Call Plotting Functions
%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%
% Plot Figures 9-12 
%%%%%%%%%%
if plot_main_figures
    simul_ELB = 'true';
    plot_supply_demand = true;
    simul_plot;
end


if plot_appendix_frontier
    %%%%%%%%%%
    % Plot Figure A4:
    %%%%%%%%%%
    simul_ELB = "true";
    plot_supply_demand = false;
    simul_plot_phiu;
    
    
    %%%%%%%%%
    % Plot Figure A5 and A7:
    %%%%%%%%%
    plot_supply_demand = false;
    simul_ELB = 'false';
    simul_plot;
    
    
    %%%%%%%%%
    % Plot Figure A6
    %%%%%%%%%
    simul_str = 'shock_series_rand_new';
    simul_ELB = 'true';
    simul_plot;
end
  
%%%%%%%%%
% Plot Figures A8/9
%%%%%%%%%
if plot_appendix_volatility
    simul_ELB = "true";
    simul_str = 'shock_series_combined_hk_5';
    simul_str_2 = "shock_series_combined_hk_5_TFP_only";
    simul_str_3 = "shock_series_combined_hk_5_MK_only";
    simul_plot_tfp_mk;
end

