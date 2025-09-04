clear all;
close all;


% Add tools directory:
addpath(genpath('auxfiles'))
addpath(genpath('helper_functions'))
addpath(genpath('plotting_functions'))

%%%%%%%%%%%%%%%%%%%%%%% USER %%%%%%%%%%%%%%%%%%%%%%%
%Plotting options:

spec_settings.str_addl = ''; %additional string to put at the end of the figure name

% Calculate comparison as high phi_pi?
spec_settings.calculate_high_phi_pi = true;


%No need to change:
spec_settings.high_phi_pi = 3.0;
spec_settings.calc_consumption_equivalence = true; %Calculate CE values
spec_settings.fix_hh_bins = false;
spec_settings.calculate_high_kappa = false;
spec_settings.high_kappa = 0.1;

%Plot saving:
spec_settings.pltfolder = "paper_figures_replication/"; %figures replication directory
spec_settings.CEpltfolder = "Distributional/"; %folder to put CE figs in

if spec_settings.calculate_high_phi_pi
    spec_settings.IRFpltfolder = "IRFs/Taylor_vs_High_Phi/"; %Ending of folder to put shock IRFs in
else
    spec_settings.IRFpltfolder = "IRFs/";
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



spec_settings.fig_names = {'FC','RP','Z','LT','LP','MP','IT','PM','WM'};
spec_settings.shock_names = {'Fixed Cost', 'Risk Premia', 'Productivity', 'Lump Sum Transfer','Discount Factor', 'Monetary Policy', 'Investment Technology', 'Price Markup','Wage Markup'};
spec_settings.fig_label = {'FixedCost', 'RiskPremia', 'Productivity', 'LumpSumTransfer','DiscountFactor', 'MonetaryPolicy', 'InvestmentTechnology', 'PriceMarkup','WageMarkup'};

spec_settings.shock_inds = [2 3 5 6 8]; %This subsets to only Risk Premia, Productivity, Discount Factor, Monetary Policy, and Price Markup shocks that are used in the paper. The user can change this to observe results from other shocks.


%Makes Taylor vs high phi IRFs and Consumption Equivalence Plots:
for lllllll = spec_settings.shock_inds %Change to 1:9 to generate results for all shocks
    dyn_analysis_Taylor;
    plot_agg_IRFs;

    if spec_settings.calc_consumption_equivalence && lllllll~=2
        plot_consumption_equiv_bars;
    end
    clearvars -except spec_settings
end