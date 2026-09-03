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
The main folders are organised as follows:
## Repository structure

The repository is organised as follows:
<img width="823" height="708" alt="image" src="https://github.com/user-attachments/assets/6e457ad8-da64-4842-96a7-d418a79e3a9c" />
