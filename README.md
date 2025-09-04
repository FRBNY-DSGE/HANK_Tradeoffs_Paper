# HANK_Tradeoffs_Paper
This repository contains code used for production of results in the 2025 paper "Tradeoffs for the Poor, Divine Coincidence for the Rich" by Marco Del Negro, Ibrahima Diagne, Keshav Dogra, Pranay Gundam, Donggyu Lee, and Brian Pacula. If you have any questions or inquiries about the codebase, please reach out to Ibrahima.Diagne@ny.frb.org.

## Software Requirement
- MATLAB (Parallel Computing Toolbox is recommended)

## Package Structure
- <b> auxfiles/: </b> This folder contains MATLAB data files required to calculate IRFs and run simulations, including
  - data_est_covid: Contains data on aggregate variables from 1992Q1-2023Q4. In units of log deviations. Used for model estimation and shock filtering.
  - Jacob_base_{Taylor/AIT}: Contains model Jacobian matrices under AIT/Taylor regimes and during ELB periods
  - MCMC_posterior_mode/MCMC_posterior_precovid_ugap_v4: Contains environmental variables and posterior parameter estimation results (using RWMH)
  - shock_series_combined_hk_5: Contains 60,000 period sequence of supply and demand shocks, including filtered shocks from 2021Q2-2022Q3 (COVID shocks) and 2008Q1-2009Q2 (Great Recession shocks)
  - shock_series_combined_noCovidGR: Contains 60,000 period sequence of supply and demand shocks, not including COVID nor Great Recession shocks
  - shock_series_combined_hk_5_MK_Only: Contains 60,000 period sequence of price markup shocks, including filtered shocks from 2021Q2-2022Q3 (COVID shocks) and 2008Q1-2009Q2 (Great Recession shocks)
  - shock_series_combined_hk_5_TFP_Only: Contains 60,000 period sequence of productivity (TFP) shocks, including filtered shocks from 2021Q2-2022Q3 (COVID shocks) and 2008Q1-2009Q2 (Great Recession shocks)
  - ss_results: Contains steady state values for parameters, state variables, and control variables
- <b> paper_figures/: </b> This folder stores aggregate IRFs, consumption equivalence (CE), and simulation figures from the paper
- <b> paper_simulations/: </b> This folder stores .mat files containing simulation output files from <code> simulate_spec.m </code>
- <b> helper_functions/: </b> This folder contains secondary functions used to run simulations, calculate IRFs, and plot results
- <b> plotting_functions/: </b> This folder contains files used to create IRFs/CE and simulation figures
- <b> shockfiles/: </b> This folder stores helper functions and files used to generate and trim replicated shock series
- <b> shocks_replication/: </b> This folder should store the output folders containing the custom/replicated shock series
- <b> paper_simulations_replication/: </b> This folder should store replicated simulation results from the paper
- <b> paper_figures_replication/: </b> This folder should store replicated figures from the paper

## Additional Documentation
- The file <code> figure_replication_doc.md </code> contains instructions for replicating both main figures and additional appendix figures from the paper.
- The file <code> shock_generation_doc.md </code> contains instructions for replicating the shock series used to generate simulations in the paper.

# Disclaimer
Copyright Federal Reserve Bank of New York. You may reproduce, use, modify, make derivative works of, and distribute and this code in whole or in part so long as you keep this notice in the documentation associated with any distributed works. Neither the name of the Federal Reserve Bank of New York (FRBNY) nor the names of any of the authors may be used to endorse or promote works derived from this code without prior written permission. Portions of the code attributed to third parties are subject to applicable third party licenses and rights. By your use of this code you accept this license and any applicable third party license.

THIS CODE IS PROVIDED ON AN "AS IS" BASIS, WITHOUT ANY WARRANTIES OR CONDITIONS OF ANY KIND, EITHER EXPRESS OR IMPLIED, INCLUDING WITHOUT LIMITATION ANY WARRANTIES OR CONDITIONS OF TITLE, NON-INFRINGEMENT, MERCHANTABILITY OR FITNESS FOR A PARTICULAR PURPOSE, EXCEPT TO THE EXTENT THAT THESE DISCLAIMERS ARE HELD TO BE LEGALLY INVALID. FRBNY IS NOT, UNDER ANY CIRCUMSTANCES, LIABLE TO YOU FOR DAMAGES OF ANY KIND ARISING OUT OF OR IN CONNECTION WITH USE OF OR INABILITY TO USE THE CODE, INCLUDING, BUT NOT LIMITED TO DIRECT, INDIRECT, INCIDENTAL, CONSEQUENTIAL, PUNITIVE, SPECIAL OR EXEMPLARY DAMAGES, WHETHER BASED ON BREACH OF CONTRACT, BREACH OF WARRANTY, TORT OR OTHER LEGAL OR EQUITABLE THEORY, EVEN IF FRBNY HAS BEEN ADVISED OF THE POSSIBILITY OF SUCH DAMAGES OR LOSS AND REGARDLESS OF WHETHER SUCH DAMAGES OR LOSS IS FORESEEABLE.
