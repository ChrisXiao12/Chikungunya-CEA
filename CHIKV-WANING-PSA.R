#----
#this script attempts to model the impact of vaccine waning
#The assumptions I make are as follows
#V no longer can go to E, everyone in V stays in V until their immunity wanes at rate theta
#Theta is derived from extrapolating a linear decline of 1 month sero response to 0
#S is vaccinated at rate psi mutliplied by 1 mo VE to simulate the fact that not all people respond
#----
#packages
library(tidyverse)
library(ggplot2)
library(mc2d)
library(scales)
#----
R0_draws <- rgamma(1000,shape = 54.8258, scale = 0.0620146)
Lambda_draws <- rgamma(1000,shape = 0.7314914, scale = 1.913898)
Gamma_draws <- rgamma(1000,shape = 115.7945, scale = 0.007111992)
Delta_draws <- rbeta(1000,61.40313,61341.73)
Phi_draws <- rbeta(1000,260.2326,8.602732)
Phi_VIM_draws <- rbeta(1000,1935.908,328.3118)
Psi_draws <- rgamma(1000,shape = 61.4656, scale = 0.0001708273)
Mu_draws <- rgamma(1000,shape = 61.4656, scale = 0.0000001130876)
Kappa_draws <- rbeta(1000,2.392294,122.3488)
Omega_draws <- rbeta(1000,95.95506,79.02587)
Beta_draws <- R0_draws * Gamma_draws
#transform to probabilities
pDelta_draws <- Delta_draws
pPhi_draws <- Phi_draws
pPhi_Vim_draws <- Phi_VIM_draws
pOmega_draws <- Omega_draws
pLambda_draws <- 1 - exp(-Lambda_draws)
pGamma_draws <- 1 - exp(-Gamma_draws)
pPsi_draws <- 1 - exp(-Psi_draws)
pMu_draws <- 1 - exp(-Mu_draws)
pKappa_draws <- 1 - exp(-Kappa_draws)
pBeta_draws <- 1 - exp(-Beta_draws)
discount_draws <- rep(0.03, 1000)
pop_draws <- rep(1000,1000)
fracinfected_draws <- rbeta(1000, 61.03465,30532.77)
fracsusceptible_draws <- 1 - fracinfected_draws
cyclelength_draws <- rep(1/52,1000)
#----
R_adverse_IXCHIQ <- rbeta(1000,60.27875,3112.287)
R_adverse_Vim <- rbeta(1000,10.81465, 405.1335)
U_S_draws <- rbeta(1000,11.75982,2.511806)
U_R_draws <- U_S_draws
Theta_I_draws <- rbeta(1000,61.44024,151269.2)
Theta_V_draws <- rbeta(1000,61.10256,11144.21)
#U_V_draws <- rbeta(1000,60.77709,4774.314)
#U_VIM_draws <- rbeta(1000,53.87676,3115.344)
U_E_draws_raw <- rbeta(1000,17.63476,7.557754)
U_I_draws_raw <- rbeta(1000,20.11337, 10.26937)
U_E_draws <- pmin(U_E_draws_raw, U_S_draws)
U_I_draws <- pmin(U_I_draws_raw, U_S_draws)
U_V_draws <- R_adverse_IXCHIQ * U_S_draws *cyclelength_draws - R_adverse_IXCHIQ * U_I_draws * cyclelength_draws
U_VIM_draws <- R_adverse_Vim * U_S_draws * cyclelength_draws - R_adverse_Vim * U_I_draws * cyclelength_draws
U_C_draws <- rbeta(1000,5.088638,5.212248)
U_D_draws <- rep(0,1000)
C_S_draws <- rep(0,1000)
Vax_Dose_draws_I <- rgamma(1000,shape = 104.7666, scale = 0.2243081)
Vax_Dose_draws_Vim <- rgamma(1000,shape = 104.7666, scale = 0.2243081)
Admin_draws <- rgamma(1000,shape = 61.4656, scale = 0.04522855)
Waste_draws <- rpert(1000, min = 0.05, mode = 0.1, max = 0.15)
#----
C_V_draws <- Vax_Dose_draws_I / (1 - Waste_draws) + Admin_draws
C_Vim_draws <- Vax_Dose_draws_Vim / (1 - Waste_draws) + Admin_draws
#----
D_hosp_draws <- rgamma(1000, shape = 61.57906, scale = 0.1761963)
Cost_hosp_draws <- rgamma(1000, shape = 61.4656, scale = 14.09113)
P_hosp_draws <- rbeta(1000, 60.43122,3291.949)
P_nhosp_draws <- 1 - P_hosp_draws
Outpatient_visit <- rgamma(1000, shape = 1.736605, scale = 8.176875)
Infectious_direct <- rgamma(1000, shape = 61.45182, scale = 0.5804548)
Chronic_direct <- rgamma(1000, shape = 61.492, scale = 1.874252)
Infectious_absenteeism_days <- rgamma(1000,shape = 1.356759, scale = 4.599195)
Infectious_absenteeism_freq <- rbeta(1000, 26.27682, 3.247697)
Infectious_caregiving_days <- rgamma(1000, shape = 11.8879, scale = 0.4054543)
Infectious_caregiving_freq <- rbeta(1000, 708.5952, 3228.045)
Chronic_absenteeism_days <- rgamma(1000, shape = 96.1739, scale = 0.1493129)
Chronic_absenteeism_freq <- rbeta(1000, 56.03872, 8.373601)
Chronic_caregiving_days <- rgamma(1000, shape = 95.20305, scale = 0.1197441)
Chronic_caregiving_freq <- rbeta(1000, 197.3991, 383.1864)
Absenteeism_cost <- rgamma(1000, shape = 61.4656, scale = 0.4822209)
#----
TC_indirect_infectious <- Absenteeism_cost * ((Infectious_absenteeism_days * Infectious_absenteeism_freq) + (Infectious_caregiving_days * Infectious_caregiving_freq))
TC_indirect_chronic <- Absenteeism_cost * ((Chronic_absenteeism_days * Chronic_absenteeism_freq) + (Chronic_caregiving_days * Chronic_caregiving_freq))
C_E_draws <- rep(0,1000)
#C_I_draws <- rgamma(1000,shape = 1.855343, scale = 126.8983)
C_I_draws <- (Cost_hosp_draws / D_hosp_draws * 7 * P_hosp_draws) + (Infectious_direct * P_nhosp_draws) + TC_indirect_infectious
C_R_draws <- rep(0,1000)
#C_C_draws <- rgamma(1000, shape = 61.4656, scale = 11.15486)
C_C_draws <- Chronic_direct + TC_indirect_chronic + 6 * Outpatient_visit
C_D_draws <- rep(0,1000)
