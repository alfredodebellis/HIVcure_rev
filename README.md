# HIVcure

## Overview 
GitHub repository for the paper "Model-based evaluation of the impact of a potential HIV cure on HIV transmission dynamics". All details are described in the manuscript:
> De Bellis, A., Willemsen, M., Guzzetta, G., van Sighem, A., Romijnders, K.A.G.J., Reiss, P., Schim van der Loeff, M.S., van de Wijgert, J.H.H.M., Nijhuis, M., Kretzschmar, M., & Rozhnova, G. (2024). Model-based evaluation of the impact of a potential HIV cure on HIV transmission dynamics in a low-incidence, high ART setting. 

## Survey
The folder [SURVEY]("SURVEY") contains the original questions asked within survey on sexual behavior data.

## Data
The folder [DATA/EPIDEMIOLOGICAL]("DATA/EPIDEMIOLOGICAL") contains epidemiological data and  PrEP data.
The folder [DATA/SURVEY_BEHAVIOR]("DATA/SURVEY_BEHAVIOR") contains aggregated survey data on sexual contact rates. There are 3 files with 529 rows containing the individual values of the variable $c$ (contact rate adjusted by condom use and steady partners, as described in Section 2.1 of the Supplementary Material).
- `total_curr.xlsx` refers to the scenario with no cure
- `total_sc1.xlsx` refers to the scenario with remission cure
- `total_sc2.xlsx` refers to the scenario with eradication cure

## Running the code
### Fitting 
To fit the model, run the following scripts in the suggested order:

1)	fit_contactrates.nb to fit Weibulls to pre-processed aggregated survey data for all scenarios
2)	diagnostic_delays.R to calibrate diagnostic delays as described in Section 2.2 of the Supplementary Material
3)	ABC.R to fit the baseline model before cure with Approximate Bayesian Computation. ABC_noimp.R and ABC_fixed_imp.R to fit the models for the sensitivity analyses without immigration and fixed proportion of undiagnosed to diagnosed imported individuals, respectively (see the paper for a detailed explanation of various assumptions on importation)
Initial conditions for the model variables in the fitting process are defined in lines 423 - 448 of ABC.R.

The results from steps 1) and 2) are already included in the DATA/SURVEY_BEHAVIOR and DATA/EPIDEMIOLOGICAL subfolders, respectively.

To avoid high computational costs of ABC procedures, we have provided pre-saved workspaces, that are "fit_threshold15.Rdata" (baseline), "fit_threshold15_noimp.Rdata" and "fit_threshold15_fixedimp.Rdata" (sensitivity analyses) that can be dowloaded from Zenodo, DOI: 10.5281/zenodo.14851476 (link: https://doi.org/10.5281/zenodo.14851476). After downloading, fit_threshold15.Rdata (or similarly the other two sensitivity analyses workspaces), place it in the RESULTS folder. This file contains the output from the Approximate Bayesian Computation (ABC) fit.
You can skip step 3) and proceed directly to simulating scenarios using the provided workspace (see instructions below).

### Simulating
To generate scenarios output files, run the following scripts:

-	scenario1.R to simulate remission
-	scenario1_noimp.R to simulate remission with no importation
-	scenario1_fixedimp.R to simulate remission with fixed proportion of undiagnosed to diagnosed imported individuals
-	scenario1_behav.R to simulate remission with behavioral changes
-	scenario1_erlang.R to simulate remission with Erlang distributed rebound times
-	scenario_mix.R to simulate remission with possibility of re-infection
-	scenario2.R to simulate eradication
-	scenario2_noimp.R to simulate eradication with no importation
-	scenario2_fixedimp.R to simulate eradication with fixed proportion of undiagnosed to diagnosed imported individuals
-	scenario2_behav.R to simulate eradication with behavioral changes
-	incr_lambda_1.R to simulate remission with HIV incidence increasing after 2022
-	incr_lambda_2.R to simulate eradication with HIV incidence increasing after 2022
-	check_proportions_groups.R to simulate the model without cure and check if prevalence by risk group remains constant

Computation for the ABC requires about 3-4 hours on a MacBook Pro (Apple M1 Pro).
Computation for each cure scenarios requires 1-2 hours.

### Plotting 
To produce final figures, run the following scripts:

-	plot_CDFcontactrates.R to plot CDF
-	plot_ABC.R to plot the fit results
-	plot_ABC_noimp.R to plot the fit results relative to ABC_noimp.R
-	plot_ABC_fixed.R to plot the fit results relative to ABC_fixedimp.R
-	plot_proportions_groups.R to plot the check results relative to check_proportions_groups.R
-	plot_scenario1.R to plot remission dynamics
-	plot_scenario1_noimp.R to plot remission dynamics simulated with scenario1_noimp.R
-	plot_scenario1_fixedimp.R.R to plot remission dynamics simulated with scenario1_fixedimp.R
-	plot_scenario1_behav.R to plot remission with behavioral changes dynamics
-	plot_scenarioerlang.R.R to plot remission dynamics with Erlang distributed rebound times
-	plot_scenariomix.R.R to plot remission dynamics simulated with scenario_mix.R
-	plot_scenario2.R to plot eradication dynamics
-	plot_scenario1_noimp.R to plot eradication dynamics simulated with scenario2_noimp.R
-	plot_scenario1_fixedimp.R.R to plot eradication dynamics simulated with scenario2_fixedimp.R
-	plot_scenario2_behav.R to plot eradication with behavioral changes dynamics
-	plot_SA_scenario1.R to plot sensitivity analyses for remission
-	plot_SA_scenario2.R to plot sensitivity analyses for eradication (after "plot_SA_scenario1.R")
-	plot_incr_1.R to plot sensitivity analyses of increasing incidence for remission
-	plot_incr_2.R to plot sensitivity analyses of increasing incidence for eradication

## Dependencies and hardware requirements 
This code is was developed on a MacBook Pro (14-inch, 2021), macOS: Sonoma 14.3.1, chip: Apple M1 Pro, memory: 32GB. 
Our study requires only a standard computer with enough RAM to support the in-memory operations.

- [Mathematica 13.3.1.0](https://www.wolfram.com/mathematica/)
- [R version 4.1.2](https://www.r-project.org/) with the following packages:
  - [deSolve](https://cran.r-project.org/web/packages/deSolve/index.html) vesion 1.30
  - [lhs](https://cran.r-project.org/web/packages/lhs/index.html) version 1.1.6
  - Additional libraries used for data processing and plotting: RColorBrewer, viridis, ggplot2, patchwork, ppcor, patchwork, readxl

No non-standard hardware is required.
Other than installation of these required software and packages, no installation is needed.

## Funding  
The authors gratefully acknowledge funding by the Aidsfonds Netherlands, grant number P53902, and funding by Aidsfonds & NWO, the SPIRAL project KICH2.V4P.AF23.001.
