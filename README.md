# Msc-Data-Science-Bayesian-Modelling-of-Heterogeneity-in-Mosquito-Control-Tool-Evaluation-Code

## Overview

This repository contains the R and Stan code used to fit and evaluate the Bayesian models developed for this project. The analysis investigates heterogeneity in experimental hut trial data for Anopheles arabiensis and Culex mosquitoes.

Four models of increasing complexity were considered:

Model 1: Baseline model with no random effects. 

Model 2: Includes day-level random effects on biting and mortality.

Model 3: Extends Model 2 by allowing the intervention parameters $\beta$ and $\kappa$ to vary between experimental days.

Model 4: Extends Model 2 by including volunteer-level variation in biting.

The repository contains code for fitting the models to observed data, conducting the simulation study, and reproducing the figures and numerical results reported in the project.

## Repository structure

The repository is organised as follows:
<img width="794" height="693" alt="image" src="https://github.com/user-attachments/assets/337c94be-bf64-47ae-8a09-6e66dbdacc21" />


ArabiensisNewResults and CulexNewResults contain the outputs from fitting Models 1–4 to the observed data. SimulationResults contains the saved results from the parameter-recovery simulations. FiguresAndTables contains the scripts used to generate the reported figures and tables. 

## Software requirements
The analysis was conducted in R using Stan through the rstan package.

The R packages required by the analysis include:
```
library(rstan)
library(loo)
library(ggplot2)
library(dplyr)
library(patchwork)
```
The analysis was conducted using:

- R 4.5.2
- RStan 2.32.7
- StanHeaders 2.32.10

The main R packages used to reproduce the reported analysis are:

- loo 2.9.0
- ggplot2 4.0.1
- dplyr 1.1.4
- tidyr 1.3.1
- patchwork 1.3.2

## Setting the working directory
Before running this analysis, download this repository and set the working directory in R to the main Models folder: setwd("path/to/Models"). You can ensure the working directory is correct using 
```
getwd()
list.files()
```
If it has been set correctly, you should see the working directory end in /Models. The files should include the folders Data, ArabiensisNewResults, CulexNewResults, SimulationResults, and FiguresAndTables, together with the model and analysis scripts e.g (RunModel1.R).

All scripts use paths relative to the Models directory. The repository structure should therefore be retained.

## Data
The analyis uses experimental hut trial data for _Anopholes arabiensis_ and _Culex_ mosquitoes. Separate datasets are used for the control and intervention arms.

The scripts expect the required datasets to be located in the Data folder with filenames:
```
arabiensis_cleanEH_BIT046_data0.csv
arabiensis_cleanEH_BIT046_data1.csv
culex_cleanEH_BIT046_data0.csv
culex_cleanEH_BIT046_data1.csv
```
The data0 csvs correspond to the control datasets and the data1 correspond to the interventions.
These datasets originate from [Fairbanks et al. (2026)](https://doi.org/10.64898/2026.08.24.746696).

The fitted Stan objects have been saved as .rds files for both the real-data and simulation analyses. Thus, the reported figures and table outputs remain reproducible using my code prior to the public release of the datasets.

## Reproduction workflow
<img width="800" height="960" alt="Blank diagram (1)" src="https://github.com/user-attachments/assets/42250cd4-ec08-4995-a9f9-81451b7d8820" />


The diagram shows the complete workflow to reproduce the results from my ERP.

The black arrow corresponds to steps taken for both species. The diagram highlights that the simulation study is only done for the _An. arabiensis_ species and must take place after running its model fit. The real-data analysis takes place after running the model fits for both the _An. arabiensis_ species and _Culex_ species.
## Reproducing the observed-data analysis
Models 1–4 are implemented in:
```
Model1.stan
Model2.stan
Model3.stan
Model4.stan
```
These can be run in R using rstan. The corresponding R scripts (RunModelX.R) prepare the data, compile the relevant Stan model and perform posterior sampling. The saved fit is saved as an .rds file and the posterior csv files are saved in the correct folder, e.g CulexNewResults/Model2 contains the Stan fit and posterior csvs for the Model 2 fit on the _Culex_ dataset.

The code supplied can be ran directly for the _An. arabiensis_ datasets. They must be changed manually for Culex. The R scripts are commented so it is clear which lines need to be changed. These lines change the datasets to the Culex datasets and change the species to Culex (to ensure outputs are saved to the correct folder.)

**Note:** For the Culex Model 4 fit max_treedepth must be changed to 15 and adapt_delta must be changed to 0.99. This is commented on the RunModel4.R code file.
The fitted Stan objects are saved within the corresponding model folders. For example:
```
ArabiensisNewResults/Model1/Model1_fit.rds
CulexNewResults/Model1/Model1_fit.rds
```
These fitted objects have also been supplied in the repository, allowing the reported results to be reproduced without refitting the models.

## Reproducing the simulation study
Parameter recovery was investigated using 100 simulated datasets for each model.

For each simulation, a posterior draw from the corresponding observed-data fit was selected and used as the true parameter values. The paramter values are taken from the saved csv files from the real-data fit to the _An .arabiensis_ data. Therefore, if you are running everything from scratch, ensure that you have run all 4 model fits for the _An. arabiensis_ datasets first. I have uploaded all posterior csvs to the correct folders already, so it is possible to reproduce the simulation study without rerunning the real-data model fits.

A new dataset with the same experimental structure was then generated and the corresponding model was refitted.

To ensure that the same posterior draws are used a seed was set to set.seed(123).

The resulting simultion outputs are stored in SimulationResults/ as:
```
Model1ParameterRecoveryResults.rds
Model2ParameterRecoveryResults.rds
Model3ParameterRecoveryResults.rds
Model4ParameterRecoveryResults.rds
```
These files can be used to reproduce the reported bias, RMSE, correlation and 95% credible interval coverage results without rerunning the full simulation study.

**Warning:** The simulations were computationally expensive, especially for Models 3 and 4. I used a computational shared facility to run the code, as it took more than 12 hours. Take this into account before running the code. The saved .rds files can be used to reproduce all reported figures, and be used any new analysis.

## Reproducing figures and tables
Scripts used to generate the reported analytical outputs are contained in FiguresAndTables/
They are separated by 2 folders: Real-data analysis code/ and Simulation analysis code/
The Real-data analysis code/ folder contains:
```
ModelParameterSummaries.R
    Calculates residual overdispersion and random-effect SD summaries

real_beta_kappa_posterior_means_and_credible_intervals_figures.R
    Produces posterior estimates and 95% credible intervals for beta and kappa

PPC_Counts_code.R
    Produces posterior predictive checks for total and fed mosquito counts

elpd_loo_values_code.R
    Calculates PSIS-LOO model comparison results
```
**Note:** The code analyses both the _Culex_ and _An. arabiensis_ files at the same time, so both model fits need to have occurred.

The Simulation analysis code/ folder contains:
```

TableA.1_ParameterRecovery.R
    Calculates parameter-recovery metrics reported in Appendix Table A.1

true_vs_simulated_beta_kappa_correlation_figures.R
    Produces true-versus-estimated parameter recovery plots for beta and kappa

Coverage_plot_code.R
    Produces the credible interval coverage figure
```
**Note:** Coverage_plot_code.R uses the table output from TableA.1_ParamterRecovery.R. Therefore, it must be run after this. I have uploaded the table, so to reproduce the plot, you can run them in either order. But, if you wish to run all the code from the beginning, without using any saved files, then ensure the order is followed correctly.

The scripts in FiguresAndTables/ use the supplied fitted-model and simulation .rds files as inputs.
## Computational Requirements
The Stan models were generally fitted using four chains with 6,000 iterations per chain, of which 3,000 were warm-up iterations. This produced 12,000 post-warm-up posterior draws.

The final Culex Model 4 fit used 8,000 iterations per chain with 4,000 warm-up iterations. For this fit, adapt_delta was set to 0.99 and max_treedepth to 15.

The simulation study involved refitting each model to 100 simulated datasets and is therefore considerably more computationally intensive than regenerating the reported figures and tables from the supplied results.
