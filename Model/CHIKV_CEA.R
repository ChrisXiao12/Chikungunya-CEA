rm(list = ls())
##*****************************************************************************************************
## This script attempts to model Chikungunya using a SVEIRD model
## Everyone in V stays in V until their immunity wanes at rate theta
## IXCHIQ = live-attenuated, Vimkunya = recombinant
## Theta is derived from extrapolating a linear decline of 1 month sero response to 0
## Individuals go to VE from V representing waning or from S representing those in which the 
## vaccine did not work
## These people progress at rate (1-phi)*psi*S
##*****************************************************************************************************

##****************************************************************************************
## Packages to load
##****************************************************************************************
packages_to_load <- c(
  "tidyverse", "mosaic", "ggplot2", "mc2d", "scales", "parallel", "dplyr", "tidyr", 
  "tibble", "readxl", "MASS", "minpack.lm", "boot", "flextable", "officer", "ggsci", 
  "future.apply", "akima", "patchwork", "DiagrammeR", "DiagrammeRsvg", "rsvg", "gt",
  "metR", "ggtext", "binom"
)

##****************************************************************************************
## Directories
##****************************************************************************************
home_dir <- "/Users/emmanuel/Dropbox/vaccines/BMJOpen/RR1"
model_dir <- file.path(home_dir, "Model")
data_dir <- file.path(home_dir, "Data")
figures_dir <- file.path(home_dir, "Figures")
tables_dir <- file.path(home_dir, "Tables")

##****************************************************************************************
## Simulation Settings
##****************************************************************************************
## cohort size
n_cohort <- 1000

## cycle length (in weeks)
cycle_length <- 1

## time horizon (years)
time_horizon_set <- unique(c(10, seq(from = 1, to = 20, by = 1)))

## willingness-to-pay thresholds
gdp_multipliers <- seq(from = 0.5, to = 4, by = 0.5)

## Vary xi, psi, phi, and pi
xi_set <- c(0.95, 1, 1.05)
zeta_set <- c(0.7, 0.8, 0.9, 1.0)
phi_set <- seq(from = 0.40, to = 0.75, by = 0.05)
psi_stated_set <- c(0.51, 0.55, 0.70)
pi_set <- seq(from = 0.05, to = 1.2, by = 0.05)
omicron_set <- c(0.50, 1.00, 1.50)
chi_set <- unique(c(1.25, seq(from = 0.75, to = 1.5, by = 0.25))) 
eta_set <- unique(c(0.5, seq(from = 0, to = 0.75, by = 0.05)))

## analytic perspectives (0 = healthcare sector perspective; 1 = societal perspective)
perspectives <- c(0, 1)

## scenarios to simulate
scenarios <- c("SOC", "IXCHIQ", "VIMKUNYA")

## scenario labels
scenarios_labs <- c(
  SOC = "No vaccination", 
  IXCHIQ = "Live-attenuated vaccine", 
  VIMKUNYA = "Recombinant vaccine"
)

## health states
state_names <- c("S", "V", "SV", "E", "I", "R", "C", "D")

## transitions for incident cases
trans_names <- c(
  "StoV", "StoSV", "StoE", "StoD",
  "VtoSV", "VtoD", 
  "SVtoE", "SVtoD",
  "EtoI", "EtoD", 
  "ItoR", "ItoC", "ItoD", 
  "RtoD", 
  "CtoR", "CtoD"
)

## number of PSA draws
n_psa_draws <- 1000

## number of bootstrap samples
nbootsamps <- 15

## seting seed
sim_seed <- 123

## number of workers for parallel computing
n_workers <- 4

## flow diagram
create_flow_diagram <- FALSE

## vaccine seroresponse
create_coeff_table <- TRUE
print_coef_plot <- TRUE
save_coeffs_plot <- TRUE

## use saved data
use_saved_cpi_adj <- FALSE
use_saved_seroresponse <- FALSE

## base case
use_saved_base_cea_data <- FALSE
create_base_cea_table_2 <- TRUE
create_base_cea_table_S1 <- TRUE

## owsa
use_saved_owsa_data <- FALSE
create_tornado_plots <- TRUE

## psa
use_saved_psa_params <- FALSE
use_saved_psa_data <- FALSE

create_psa_plots <- FALSE
use_saved_psa_boot_data <- FALSE

## sensitivity analyses
use_saved_sensi_cea_data <- FALSE
use_saved_phi_pi_sensi_cea_data <- FALSE
use_saved_zeta_sensi_cea_data <- FALSE
use_saved_sensi_time_horizon_cea_data <- FALSE
use_saved_sensi_time_eta_zeta_cea_data <- FALSE

##****************************************************************************************
## Functions
##****************************************************************************************
## load functions
source(file = file.path(model_dir, "Functions.R"))

## load packages
load_libraries(packages = packages_to_load)

##****************************************************************************************
## Model Flow Diagram
##****************************************************************************************
source(file = file.path(model_dir, "Flow_diagram.R"))

##****************************************************************************************
## Input Data
##****************************************************************************************
source(file = file.path(model_dir, "Data.R"))

##****************************************************************************************
## Analyses
##****************************************************************************************
#source(file = file.path(model_dir, "Analyses.R"))

##*****************************************************************************************************
## END OF SCRIPT
##***************************************************************************************************** 
  