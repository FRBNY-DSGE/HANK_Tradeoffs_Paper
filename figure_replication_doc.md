
## Replicating Main Figures

### Section 3:
The main file is <code>main_shock_analysis.m</code>, which calls <code>dyn_analysis_Taylor</code> and <code>dyn_analysis_AIT</code> to conduct the primary analysis. All IRF plotting is done in <code>plot_agg_IRFs</code> and <code>plot_C_dist_IRFs</code>, and the consumption equivalence plots are made in <code>plot_consumption_equiv_bars</code>.


The following settings in <code>main_shock_analysis.m</code> determine the plots that will be made and how they will be saved.


<ul>
<li> <code>spec_settings.calculate_high_phi_pi</code>: Determines if the baseline IRF will be compared to a version with a higher monetary policy response to inflation. Defaults to true.
<li> <code>spec_settings.high_phi_pi</code>: Determines the value of $\phi_{\pi}$ in the "high $\phi_{\pi}$" scenario. Defaults to 3.0.
<li> <code>spec_settings.calc_consumption_equivalence</code>: Determines if consumption equivalences are calculated and plotted. Defaults to true.
<li> <code>spec_settings.plot_dist_C_IRFs</code>: Determines if plots with C-10/C-Q3/C-90 IRFs are generated for each shock. Defaults to true.
</ul>

The three main shocks considered in the paper are Discount Factor (demand), Monetary Policy, and Price Markup (supply) shocks. In order to make IRFs for these shocks comparing our baseline, estimated parameters against a stronger monetary policy response to inflation (high $\phi_{\pi}$), set <code>spec_settings.calculate_high_phi_pi = true</code>.

IRF plots are saved as "ShockName.png". Consumption Equivalence plots will save as "ShockName_CE.png", and Consumption Equivalence plots showing the difference between high and low $\phi_{\pi}$ are saved as "ShockName_CE_diff.png". After running the spec file as organized, the following figures are created:

<ol>
<li> <b> Monetary policy shocks: </b> </li>
<ul>
<li> <em> Figure 2 </em>: <em> Responses to monetary policy shocks: aggregate variables </em> : <code>paper_figures_replication/IRFs/Taylor_High_Phi/MonetaryPolicy.png</code>
<li> <em> Figure 3 </em>: <em> The distributional consequences of monetary policy shocks </em> : <code>paper_figures_replication/Distributional/MonetaryPolicy_CE.png</code>
</ul>
<li> <b> Supply shocks: </b> </li>
<ul>
<li> <em> Figure 4 </em>: <em> Responses to markup shocks: aggregate variables </em> : <code>paper_figures_replication/IRFs/Taylor_High_Phi/PriceMarkup.png</code>
<li> <em> Figure 5 </em>: <em> The distributional consequences of markup shocks </em> : <code>paper_figures_replication/Distributional/PriceMarkup_CE.png</code>
<li> <em> Figure 6 </em>: <em> The distributional consequences of a stronger policy response to markup shocks </em> : <code>paper_figures_replication/Distributional/PriceMarkup_CE_diff.png</code>
</ul>
<li> <b> Demand shocks: </b> </li>
<ul>
<li> <em> Figure 7 </em>: <em> Responses to demand (discount rate) shocks: aggregate variables </em> : <code>paper_figures_replication/IRFs/Taylor_High_Phi/DiscountFactor.png</code>
<li> <em> Figure 8 </em>: <em> The distributional consequences of demand (discount rate) shocks </em> : <code>paper_figures_replication/Distributional/DiscountFactor_CE.png</code>
</ul>
</ol>

<br />

### Section 4:

#### Simulating Figures 9-12

To generate the results for these figures, we run <code> simulate_spec.m </code>
twice. For the first-run, apply the following settings at the top of the file:

<ul>
  <li> <code>spec_settings.policy_rule = 'Taylor';</code> </li>
  <li> <code>spec_settings.simul_ELB = true;</code> </li>
  <li> <code>spec_settings.simul_disagg = true;</code> </li>
  <li> <code>spec_settings.use_covidGR = true;</code> </li>

  <li> <code>spec_settings.high_kappa = false;</code> </li>
</ul>

Then, assign the output folder to contain simulation results on the line <code>spec_settings.sim_save_folder</code>. The simulation code should produce a set of simulation files differentiated by choices of $\phi_{\pi}$ as specified in the vector <code> phi_pi_vals </code>. Moreover, simulation files will be differentiated by shock types {all, supply, demand}.

Then, for the second-run, apply the following settings to obtain the high-kappa
(high PC slope) results:

<ul>
  <li> <code>spec_settings.policy_rule = 'Taylor';</code> </li>
  <li> <code>spec_settings.simul_ELB = true;</code> </li>
  <li> <code>spec_settings.simul_disagg = false;</code> </li>
  <li> <code>spec_settings.use_covidGR = true;</code> </li>

  <li> <code>spec_settings.high_kappa = true;</code> </li>
</ul>

</br>

#### Plotting Figures 9-12

In <code> simulate_plot_user </code>, first specify the input directory <code>results_folder</code> for the simulation results and output directory <code>plt_folder</code> for the plots of interest.

Next, to create plots from the main section of the paper, set <code> plot_appendix = false </code>.

In the "Kappa Changeables" section, the file knows to read in high-kappa results to create
Figure 12.

The file then calls <code> simul_plot.m </code> to create the necessary figures:
<ul>
<li>
<em> Figure 9: Frontier curves- consumption vs. inflation volatility by wealth: </em>
<code>paper_figures_replication/Simulations/ELB=true/TaylorFrontierPlots/</code> and
<code>paper_figures_replication/Simulations/ELB=true/FrontierPlots/SupplyDemand/</code>
</li>
<li>
<em> Figure 10: Standard deviation of consumption by wealth: </em>
<code> paper_figures_replication/Simulations/ELB=true/Supply_Demand_Disagg/Taylor/</code>
</li>
<li>
<em> Figure 11: Standard deviation of aggregate variables: </em>
<code> paper_figures_replication/Simulations/ELB=true/Supply_Demand_Disagg/Taylor/</code>
</li>
<li>
<em> Figure 12: Frontier curves- Baseline vs. high PC slope: </em>
<code> paper_figures_replication/Simulations/ELB=true/FrontierPlots/HighKappa/</code>
</li>
</ul>


## Replicating Appendix Figures

### Additional impulse responses:

#### <b> Figure A1 </b>: <em> Responses to productivity shocks: aggregate variables</em>
To generate this figure, in <code> main_shock_analysis.m </code>, follow the aforementioned settings outlined in Section 3. It is saved in <code> paper_figures_replication/IRFs/Taylor_vs_High_Phi/Productivity.png</code>

#### <b> Figure A2 </b>: <em> Responses to risk premium shocks: aggregate variables </em>
To generate this figure, in <code> main_shock_analysis.m </code>, follow the aforementioned settings outlined in Section 3. It is saved in <code> paper_figures_replication/IRFs/Taylor_vs_High_Phi/RiskPremia.png</code>

### Additional distributional responses:

#### <b> Figure A3 </b>: <em> The distributional consequences of productivity shocks </em>
This figure is generated as a byproduct of generating figure A1. It is saved in <code> paper_figures_replication/Distributional/Productivity_CE.png</code>

### Additional frontier curves:

#### Simulating Figures A4-7

##### <b> Figure A4 </b>: <em> Frontier curves for different values of the response to unemployment $\phi_u$ </em>
To generate this result, run the file <code>simulate_spec_phiu.m</code>. This is identical to main spec file <code>simulate_spec.m</code>, but it loops over relevant $\phi_u$ values. Settings will not need to be changed in this file other than changing the desired $\phi_u$ values on line 76. The simulation code should produce simulation files differentiated by $\phi_{\pi}$ and $\phi_{u}$.

##### <b> Figure A5 </b>: <em> Frontier curves: consumption vs inflation volatility by wealth—simulations without imposing the ZLB </em>
To generate this result in <code> simulate_spec.m </code>, apply the following settings at the top of the file:

<ul>
  <li> <code>spec_settings.policy_rule = 'Taylor';</code> </li>
  <li> <code>spec_settings.simul_ELB = false;</code> </li>
  <li> <code>spec_settings.simul_disagg = false;</code> </li>
  <li> <code>spec_settings.use_covidGR = true;</code> </li>

  <li> <code>spec_settings.high_kappa = false;</code> </li>
</ul>

The simulation code should produce simulation files with <code> ELB=false </code> in the save-string, differentiated by $\phi_{\pi}$.

##### <b> Figure A6 </b>: <em> Frontier curves: consumption vs inflation volatility by wealth—simulations without Covid or GR episodes  </em>
To generate this result in <code> simulate_spec.m </code>, apply the following settings at the top of the file:

<ul>
  <li> <code>spec_settings.policy_rule = 'Taylor';</code> </li>
  <li> <code>spec_settings.simul_ELB = true;</code> </li>
  <li> <code>spec_settings.simul_disagg = false;</code> </li>
  <li> <code>spec_settings.use_covidGR = false;</code> </li>

  <li> <code>spec_settings.high_kappa = false;</code> </li>
</ul>

The simulation code should produce simulation files with <code>shock_series_nocovidGR</code> in the save-string, differentiated by $\phi_{\pi}$.

##### <b> Figure A7 </b>: <em> Frontier curves: baseline vs high PC slope-simulations without imposing the ZLB</em>
To generate this result in <code> simulate_spec.m </code>, apply the following settings at the top of the file:

<ul>
  <li> <code>spec_settings.policy_rule = 'Taylor';</code> </li>
  <li> <code>spec_settings.simul_ELB = false;</code> </li>
  <li> <code>spec_settings.simul_disagg = false;</code> </li>
  <li> <code>spec_settings.use_covidGR = true;</code> </li>

  <li> <code>spec_settings.high_kappa = true;</code> </li>
</ul>

The simulation code should produce simulation files with <code> kappa=0_1 </code> and <code> ELB=false </code> in the save-string, differentiated by $\phi_{\pi}$.

Note that to plot either A5 or A7, both simulations must be run.

#### Plotting Figures A4-7
Set the toggles <code>plot_main_figures = true</code> and <code>plot_appendix_frontier = true</code> in <code> simulation_plot_user </code>. Accordingly, the following figures are created:

<ul>
<li>
<em> Figure A4: Frontier curves for different values of the response to unemployment</em>:
<code>paper_figures_replication/Simulations/ELB=true/Phi_u/FrontierPlots/</code>
</li>
<li>
<em> Figure A5: Frontier curves: consumption vs inflation volatility by wealth—simulations without imposing the ZLB </em>:
<code> paper_figures_replication/Simulations/ELB=false/TaylorFrontierPlots/</code>
</li>
<li>
<em> Figure A6: Frontier curves: consumption vs inflation volatility by wealth—simulations without Covid or GR episodes</em>:
<code> paper_figures_replication/Simulations/ELB=true/nocovidGR/TaylorFrontierPlots/</code>
</li>
<li>
<em> Figure A7: Frontier curves: baseline vs high PC slope-simulations without imposing the ZLB</em>
<code> paper_figures_replication/Simulations/ELB=false/FrontierPlots/HighKappa/</code>
</li>
</ul>

### Additional volatility plots:

#### Simulating Figures A8 and A9
To generate these results, run <code> simulate_spec_TFP_MK.m </code>. To create the simulation results with productivity (TFP) shocks only, set the <code> spec_settings.use_TFP_or_MK </code> toggle to "TFP." Likewise, the same is done to generate simulation results with price markup ("MK") shocks. All other toggles remain unchanged.
The simulation code should produce simulation files differentiated by the save-strings
<code> shock_series_combined_hk5_{TFP/MK}_only</code>.

#### Plotting Figures A8 and A9
To plot these figures, in <code> simulation_plot_user </code>, set the toggle <code> plot_appendix_volatility = true</code>. Accordingly, the following figures are created:

<ul>
<li>
<em> Figure A8: Standard deviation of consumption by wealth: Markup vs productivity shocks</em>:
<code>paper_figures_replication/Simulations/ELB=true/TFP_MK_Disagg/</code>
</li>
<li>
<em> Figure A9: Standard deviation of aggregate variables: Markup vs productivity shocks</em>:
<code> paper_figures_replication/Simulations/ELB=true/TFP_MK_Disagg/</code>
</li>
</ul>



