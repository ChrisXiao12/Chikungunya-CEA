#----
#this script attempts to model Chikungunya using a SVEIRD model
#Everyone in V stays in V until their immunity wanes at rate theta
#IXCHIQ = live-attenuated, Vimkunya = recombinant
#Theta is derived from extrapolating a linear decline of 1 month sero response to 0
#Individuals go to VE from V representing waning or from S representing those in which the vaccine did not work
#These people progress at rate (1-phi)*psi*S
#----
#packages
library(tidyverse)
library(ggplot2)
library(mosaic)
library(mc2d)
library(scales)
library(parallel)
library(dplyr)
library(tibble)
#----
#base_case parameter dataframes
liveattenuated_base_case <- list(
  pbeta = 0.939189937,
  plambda = 0.753403036,
  pgamma = 0.56112007,
  pdelta = 0.001,
  pphi = 0.991,
  ppsi = 0.010445067,
  pmu = 6.95098e-06,
  pkappa = 0.018995353,
  pomega = 0.521,
  discount = 0.03,
  ptheta = 0.000406315,
  popsize = 1000,
  fracsusceptible = 0.9980005,
  fracinfected = 0.0019995,
  cycle_length = 0.019230769,
  U_S = 0.824,
  U_E = 0.700,
  U_V = 0.00006,
  U_I = 0.662,
  U_R = 0.824,
  U_C = 0.494,
  U_D = 0.000,
  C_S = 0.00,
  C_V = 28.89,
  C_E = 0.00,
  C_I = 235.44,
  C_R = 0.00,
  C_C = 685.64,
  C_D = 0.00
)

recombinant_base_case <- list(
  pbeta = 0.939189937,
  plambda = 0.753403036,
  pgamma = 0.56112007,
  pdelta = 0.001,
  pphi = 0.978,
  ppsi = 0.010445067,
  pmu = 0.00000695098,
  pkappa = 0.018995353,
  pomega = 0.521,
  discount = 0.03,
  ptheta = 0.005453202,
  popsize = 1000,
  fracsusceptible = 0.9980005,
  fracinfected = 0.0019995,
  cycle_length = 0.019230769,
  U_S = 0.824,
  U_E = 0.700,
  U_V = 0.000081,
  U_I = 0.662,
  U_R = 0.824,
  U_C = 0.494,
  U_D = 0.000,
  C_S = 0.00,
  C_V = 28.89,
  C_E = 0.00,
  C_I = 235.44,
  C_R = 0.00,
  C_C = 685.64,
  C_D = 0.00
)
#----
#create variable draws for PSA
R0_draws <- rgamma(1000,shape = 54.8258, scale = 0.0620146)
Lambda_draws <- rgamma(1000,shape = 0.7314914, scale = 1.913898)
Gamma_draws <- rgamma(1000,shape = 115.7945, scale = 0.007111992)
Delta_draws <- rbeta(1000,61.40313,61341.73)
Phi_draws <- rbeta(1000,260.2326,8.602732)
Phi_R_draws <- rbeta(1000,1935.908,328.3118)
Psi_draws <- rgamma(1000,shape = 61.4656, scale = 0.0001708273)
Mu_draws <- rgamma(1000,shape = 61.4656, scale = 0.0000001130876)
Kappa_draws <- rbeta(1000,2.392294,122.3488)
Omega_draws <- rbeta(1000,95.95506,79.02587)
Beta_draws <- R0_draws * Gamma_draws
#transform to probabilities
pDelta_draws <- Delta_draws
pPhi_draws <- Phi_draws
pPhi_R_draws <- Phi_R_draws
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
R_adverse_L <- rbeta(1000,60.27875,3112.287)
R_adverse_R <- rbeta(1000,10.81465, 405.1335)
U_S_draws <- rbeta(1000,28691.07,6128.19)
U_R_draws <- U_S_draws
Theta_L_draws <- rbeta(1000,61.44024,151269.2)
Theta_R_draws <- rbeta(1000,61.10256,11144.21)
U_E_draws_raw <- rbeta(1000,23.90391,10.24453)
U_I_draws_raw <- rbeta(1000,138.2644, 70.59423)
U_E_draws <- pmin(U_E_draws_raw, U_S_draws) #force U_E and U_I to be less than U_S logically
U_I_draws <- pmin(U_I_draws_raw, U_S_draws)
U_V_draws <- R_adverse_L * U_S_draws *cyclelength_draws - R_adverse_L * U_I_draws * cyclelength_draws #disutility from adverse effects
U_VR_draws <- R_adverse_R * U_S_draws * cyclelength_draws - R_adverse_R * U_I_draws * cyclelength_draws
U_C_draws <- rbeta(1000,5.088638,5.212248)
U_D_draws <- rep(0,1000)
C_S_draws <- rep(0,1000)
Vax_Dose_draws_L <- rgamma(1000,shape = 104.7666, scale = 0.2243081)
Vax_Dose_draws_R <- rgamma(1000,shape = 104.7666, scale = 0.2243081)
Admin_draws <- rgamma(1000,shape = 61.4656, scale = 0.04522855)
Waste_draws <- rpert(1000, min = 0.05, mode = 0.1, max = 0.15)
C_V_draws <- Vax_Dose_draws_L / (1 - Waste_draws) + Admin_draws #adds vaccine waste and admin cost
C_VR_draws <- Vax_Dose_draws_R / (1 - Waste_draws) + Admin_draws
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
TC_indirect_infectious <- Absenteeism_cost * ((Infectious_absenteeism_days * Infectious_absenteeism_freq) + (Infectious_caregiving_days * Infectious_caregiving_freq))
TC_indirect_chronic <- Absenteeism_cost * ((Chronic_absenteeism_days * Chronic_absenteeism_freq) + (Chronic_caregiving_days * Chronic_caregiving_freq))
C_E_draws <- rep(0,1000)
C_I_draws <- (Cost_hosp_draws / D_hosp_draws * 7 * P_hosp_draws) + (Infectious_direct * P_nhosp_draws) + TC_indirect_infectious
C_R_draws <- rep(0,1000)
C_C_draws <- Chronic_direct + TC_indirect_chronic + 6 * Outpatient_visit
C_D_draws <- rep(0,1000)
#----
#create a dataframe of draws for PSA analysis
#live-attenuated
parametersdf <- data.frame(
  pbeta = pBeta_draws,
  plambda = pLambda_draws,
  pgamma = pGamma_draws,
  pdelta = pDelta_draws,
  pphi = pPhi_draws,
  ppsi = pPsi_draws,
  pmu = pMu_draws,
  pkappa = pKappa_draws,
  pomega = pOmega_draws,
  discount = discount_draws,
  popsize = pop_draws,
  fracsusceptible = fracsusceptible_draws,
  fracinfected = fracinfected_draws,
  cycle_length = cyclelength_draws,
  ptheta = Theta_L_draws,
  U_S = U_S_draws,
  U_V = U_V_draws,
  U_E = U_E_draws,
  U_I = U_I_draws,
  U_R = U_R_draws,
  U_C = U_C_draws,
  U_D = U_D_draws,
  C_S = C_S_draws,
  C_V = C_V_draws,
  C_E = C_E_draws,
  C_I = C_I_draws,
  C_R = C_R_draws,
  C_C = C_C_draws,
  C_D = C_D_draws
)
#Recombinant
Recomboparametersdf <- parametersdf
Recomboparametersdf$pphi <- pPhi_R_draws
Recomboparametersdf$U_V <- U_VR_draws
Recomboparametersdf$C_V <- C_VR_draws
Recomboparametersdf$ptheta <- Theta_R_draws
#----
#create a function which runs the SVEIRD model
#520 cycles with 1 week cycle lengths
run_SVEIRD5 <- function(params, return_trace = FALSE) {
  #load in the parameters
  pbeta <- params[["pbeta"]]
  plambda <- params[["plambda"]]
  pgamma <- params[["pgamma"]]
  pdelta <- params[["pdelta"]]
  pphi <- params[["pphi"]]
  ppsi <- params[["ppsi"]]
  pmu <- params[["pmu"]]
  pkappa <- params[["pkappa"]]
  pomega <- params[["pomega"]]
  discount <- params[["discount"]]
  ptheta <- params[["ptheta"]]
  popsize <- params[["popsize"]]
  fracsusceptible <- params[["fracsusceptible"]]
  fracinfected <- params[["fracinfected"]]
  cycle_length <- params[["cycle_length"]]

  U_S <- params[["U_S"]] * cycle_length #adjusted for length of state
  U_E <- params[["U_E"]] * cycle_length
  U_V <- params[["U_V"]] #already adjusted
  U_I <- params[["U_I"]] * cycle_length
  U_R <- params[["U_R"]] * cycle_length
  U_C <- params[["U_C"]] * cycle_length
  U_D <- 0

  C_S <- params[["C_S"]]
  C_E <- params[["C_E"]]
  C_V <- params[["C_V"]]
  C_I <- params[["C_I"]] #already adjusted in terms of week lengths
  C_R <- params[["C_R"]]
  C_C <- params[["C_C"]] * cycle_length #since cost was estimated in annual costs
  C_D <- 0
  names <- c("S", "E", "V", "VE", "I", "R", "C", "D", "N", "Check",
             "SV", "SE", "SVE", "SD", "VVE", "VD", "VEE", "VED", "EI", "ED", "IR", "IC", "ID", "RD", "CR", "CD")
  trace_v <- data.frame(matrix(0, nrow = 520, ncol = length(names)))
  colnames(trace_v) <- names

  #Create the trace matrices
  trace_v <- data.frame(matrix(0, nrow = 520, ncol = length(names)))
  colnames(trace_v) <- names

  # Initialize with only S and I non-zero; all else 0
  trace_v[1, ] <- list(
    fracsusceptible * popsize,  # S
    0,                          # E
    0,                          # V
    0,                          # VE
    fracinfected * popsize,     # I
    0,                          # R
    0,                          # C
    0,                          # D
    popsize,                    # N
    1000,                       # Check
    0, 0, 0, 0,                 # SV, SE, SVE, SD
    0, 0, 0, 0,                 # VVE, VD, VEE, VED
    0, 0, 0, 0, 0, 0, 0         # EI, ED, IR, IC, ID, RD, CR, CD
  )
  trace_nv <- trace_v
  #To model differential betas
  pbetaS <- pbeta
  pbetaSV <- pbeta
  #change to have pbetaS or pbetaSV be different where pbetaS is beta for
  #susceptible vaccine naive and pbetaSV is beta for susceptible vaccine exposed
  #loop through the equations which define state transitions
  for(i in 2:length(trace_v$S)){
    trace_v$SV[i] <- trace_v$S[i-1] * ppsi * pphi
    trace_v$SVE[i] <- trace_v$S[i-1] * ppsi * (1-pphi)
    trace_v$SE[i] <- trace_v$S[i-1] * pbetaS * trace_v$I[i-1] / trace_v$N[i-1]
    trace_v$SD[i] <- trace_v$S[i-1] * pmu
    trace_v$VVE[i] <- trace_v$V[i-1] * ptheta
    trace_v$VED[i] <- trace_v$VE[i-1] * pmu
    trace_v$VEE[i] <- trace_v$VE[i-1] * pbetaSV * trace_v$I[i-1] / trace_v$N[i-1]
    trace_v$VD[i] <- trace_v$V[i-1] * pmu
    trace_v$EI[i] <- trace_v$E[i-1] * plambda
    trace_v$ED[i] <- trace_v$E[i-1] * pmu
    trace_v$IR[i] <- trace_v$I[i-1] * (1 - pomega) * pgamma
    trace_v$IC[i] <- trace_v$I[i-1] * pomega * pgamma
    trace_v$ID[i] <- trace_v$I[i-1] * (pmu + pdelta)
    trace_v$RD[i] <- trace_v$R[i-1] * pmu
    trace_v$CR[i] <- trace_v$C[i-1] * pkappa
    trace_v$CD[i] <- trace_v$C[i-1] * pmu
    trace_v$S[i] <- trace_v$S[i-1] - trace_v$SV[i] - trace_v$SE[i] - trace_v$SD[i] -trace_v$SVE[i]
    trace_v$E[i] <- trace_v$E[i-1] - trace_v$EI[i] - trace_v$ED[i] + trace_v$SE[i] + trace_v$VEE[i]
    trace_v$V[i] <- trace_v$V[i-1] - trace_v$VVE[i] - trace_v$VD[i] + trace_v$SV[i]
    trace_v$VE[i] <- trace_v$VE[i-1] + trace_v$VVE[i] + trace_v$SVE[i] - trace_v$VED[i] - trace_v$VEE[i]
    trace_v$I[i] <- trace_v$I[i-1] - trace_v$IR[i] - trace_v$IC[i] - trace_v$ID[i] + trace_v$EI[i]
    trace_v$R[i] <- trace_v$R[i-1] - trace_v$RD[i] + trace_v$IR[i] + trace_v$CR[i]
    trace_v$C[i] <- trace_v$C[i-1] - trace_v$CD[i] - trace_v$CR[i] + trace_v$IC[i]
    trace_v$D[i] <- trace_v$D[i-1] + trace_v$CD[i] + trace_v$SD[i] + trace_v$VD[i] + trace_v$RD[i] + trace_v$ID[i] + trace_v$ED[i] + trace_v$VED[i]
    trace_v$N[i] <- trace_v$S[i] + trace_v$E[i] + trace_v$V[i] + trace_v$I[i] + trace_v$R[i] + trace_v$C[i] + trace_v$VE[i]
    trace_v$Check[i] <- trace_v$S[i] + trace_v$V[i] + trace_v$I[i] + trace_v$R[i] + trace_v$C[i] + trace_v$E[i] + trace_v$D[i] + trace_v$VE[i]
  }
  #repeat for the nonvaccinated
  for(i in 2:length(trace_nv$S)){
    trace_nv$SV[i] <- 0
    trace_nv$SVE[i] <- 0
    trace_nv$SE[i] <- trace_nv$S[i-1] * pbetaS * trace_nv$I[i-1] / trace_nv$N[i-1]
    trace_nv$SD[i] <- trace_nv$S[i-1] * pmu
    trace_nv$VVE[i] <- trace_nv$V[i-1] * ptheta
    trace_nv$VED[i] <- trace_nv$VE[i-1] * pmu
    trace_nv$VEE[i] <- trace_nv$VE[i-1] * pbetaSV * trace_nv$I[i-1] / trace_nv$N[i-1]
    trace_nv$VD[i] <- trace_nv$V[i-1] * pmu
    trace_nv$EI[i] <- trace_nv$E[i-1] * plambda
    trace_nv$ED[i] <- trace_nv$E[i-1] * pmu
    trace_nv$IR[i] <- trace_nv$I[i-1] * (1 - pomega) * pgamma
    trace_nv$IC[i] <- trace_nv$I[i-1] * pomega * pgamma
    trace_nv$ID[i] <- trace_nv$I[i-1] * (pmu + pdelta)
    trace_nv$RD[i] <- trace_nv$R[i-1] * pmu
    trace_nv$CR[i] <- trace_nv$C[i-1] * pkappa
    trace_nv$CD[i] <- trace_nv$C[i-1] * pmu
    trace_nv$S[i] <- trace_nv$S[i-1] - trace_nv$SV[i] - trace_nv$SE[i] - trace_nv$SD[i] - trace_nv$SVE[i]
    trace_nv$E[i] <- trace_nv$E[i-1] - trace_nv$EI[i] - trace_nv$ED[i] + trace_nv$SE[i] + trace_nv$VEE[i]
    trace_nv$V[i] <- trace_nv$V[i-1] - trace_nv$VVE[i] - trace_nv$VD[i] + trace_nv$SV[i]
    trace_nv$VE[i] <- trace_nv$VE[i-1] + trace_nv$VVE[i] + trace_nv$SVE[i] - trace_nv$VED[i] - trace_nv$VEE[i]
    trace_nv$I[i] <- trace_nv$I[i-1] - trace_nv$IR[i] - trace_nv$IC[i] - trace_nv$ID[i] + trace_nv$EI[i]
    trace_nv$R[i] <- trace_nv$R[i-1] - trace_nv$RD[i] + trace_nv$IR[i] + trace_nv$CR[i]
    trace_nv$C[i] <- trace_nv$C[i-1] - trace_nv$CD[i] - trace_nv$CR[i] + trace_nv$IC[i]
    trace_nv$D[i] <- trace_nv$D[i-1] + trace_nv$CD[i] + trace_nv$SD[i] + trace_nv$VD[i] + trace_nv$RD[i] + trace_nv$ID[i] + trace_nv$ED[i] + trace_nv$VED[i]
    trace_nv$N[i] <- trace_nv$S[i] + trace_nv$E[i] + trace_nv$V[i] + trace_nv$I[i] + trace_nv$R[i] + trace_nv$C[i]
    trace_nv$Check[i] <- trace_nv$S[i] + trace_nv$V[i] + trace_nv$I[i] + trace_nv$R[i] + trace_nv$C[i] + trace_nv$E[i] + trace_nv$D[i] + trace_nv$VE[i]
  }
  #create a vector of utilities and costs
  utility_vector <- c(U_S, U_E, U_S, U_S, U_I, U_R, U_C, U_D)
  cost_vector <- c(C_S, C_E, C_V, C_V, C_I, C_R, C_C, C_D)
  #create trace matrix over which to apply these vectors
  #Trace matrices include only the states and not the transition columns
  utility_trace_v <- trace_v[,1:8]
  utility_trace_v <- apply(trace_v[,1:8], 1, function(row) row * utility_vector)
  #transpose the trace matrix after
  utility_trace_v <- t(utility_trace_v)
  #need to subtract the disutilities from adverse events related to vaccination
  for(i in 1:nrow(utility_trace_v)) {
    utility_trace_v[i,3] <- utility_trace_v[i,3] - trace_v$SV[i] * U_V
    utility_trace_v[i,4] <- utility_trace_v[i,4] - trace_v$SVE[i] * U_V
  }
  utility_trace_nv <- trace_nv[,1:8]
  utility_trace_nv <- apply(trace_nv[,1:8], 1, function(row) row * utility_vector)
  utility_trace_nv <- t(utility_trace_nv)
  cost_trace_v <- trace_v[,1:8]
  cost_trace_v <- apply(cost_trace_v, 1, function(row) row * cost_vector)
  cost_trace_v <- t(cost_trace_v)
  #add in onetime vaccination costs (costs only come form those that newly transitioned to V)
  for(i in 1:nrow(cost_trace_v)) {
    cost_trace_v[i,3] <- trace_v$SV[i] * C_V
    cost_trace_v[i,4] <- trace_v$SVE[i] * C_V
  }
  trace_nv$V <- trace_nv$SV
  cost_trace_nv <- trace_nv[,1:8]
  cost_trace_nv <- apply(trace_nv[,1:8], 1, function(row) row * cost_vector)
  cost_trace_nv <- t(cost_trace_nv)
  discount_factors <- 1 / (1 + discount) ^ ((0:(520 - 1)) / 52)
  #apply discounting
  v_eff_d <- sum(rowSums(utility_trace_v) * discount_factors)
  nv_eff_d <- sum(rowSums(utility_trace_nv) * discount_factors)
  v_cost_d <- sum(rowSums(cost_trace_v) * discount_factors)
  nv_cost_d <- sum(rowSums(cost_trace_nv) * discount_factors)
  if (return_trace) {
    return(list(
      v_trace = trace_v,
      nv_trace = trace_nv,
      v_cost_d = v_cost_d,
      nv_cost_d = nv_cost_d,
      v_eff_d = v_eff_d,
      nv_eff_d = nv_eff_d
    ))
  }
  #return the discounted costs and effects for vaccination vs no vaccination
  return(c(v_cost_d = v_cost_d, nv_cost_d = nv_cost_d, v_eff_d = v_eff_d, nv_eff_d = nv_eff_d))
}
#----
#get base-case results and calculate NMB / INMB with a WTP of 7000
base_case_results_live <- run_SVEIRD5(liveattenuated_base_case)
base_case_results_recombinant <- run_SVEIRD5(recombinant_base_case)
NMB_live <- base_case_results_live[3] * 7000 - base_case_results_live[1]
NMB_novax <- base_case_results_live[4] * 7000 - base_case_results_live[2]
NMB_recombo <- base_case_results_recombinant[3] * 7000 - base_case_results_recombinant[1]
INMB_livevsnovax <- NMB_live - NMB_novax
INMB_recombovsnovax <- NMB_recombo - NMB_novax
INMB_livevsrecombo <- NMB_live - NMB_recombo
#----
#multicore processing for PSA
num_cores <- detectCores() - 1
cl <- makeCluster(num_cores)  # Create parallel cluster

clusterExport(cl, varlist = c("parametersdf", "run_SVEIRD5"))

clusterEvalQ(cl, { library(dplyr); library(tibble) })

#Run parallel execution
#Runs the SVEIRD function over each row of the PSA draws
results_list <- parLapply(cl, 1:1000, function(i) {
  # Load required libraries inside worker
  library(dplyr)
  library(tibble)

  run_SVEIRD5(parametersdf[i, ])
})

# Stop cluster
stopCluster(cl)

# Convert list to dataframe
resultsdf <- do.call(rbind, results_list)
resultsdf <- as.data.frame(resultsdf)
#----
#repeat this for the recombinant vaccine
num_cores <- detectCores() - 1
cl <- makeCluster(num_cores)  # Create parallel cluster

clusterExport(cl, varlist = c("Recomboparametersdf", "run_SVEIRD5"))

clusterEvalQ(cl, { library(dplyr); library(tibble) })

# Run parallel execution
Recomboresults_list <- parLapply(cl, 1:1000, function(i) {
  # Load required libraries inside worker
  library(dplyr)
  library(tibble)

  run_SVEIRD5(Recomboparametersdf[i, ])
})

# Stop cluster
stopCluster(cl)

# Convert list to dataframe
Recomboresultsdf <- do.call(rbind, Recomboresults_list)
Recomboresultsdf <- as.data.frame(Recomboresultsdf)
#----
#NMB calculations for PSA
#create a sequence of WTP values
wtp_values <- seq(1000,150000, by = 1000)
ceac_df <- data.frame(WTP = wtp_values, NoVax = NA, Live = NA, Recombo = NA)

nmb_list <- vector("list", length(wtp_values))
names(nmb_list) <- wtp_values
#for each WTP value calculate the NMB for each draw
#Calculate which strategy has the highest NMB
#Calculate how many times the strategy wins
for (i in seq_along(wtp_values)) {
  wtp <- wtp_values[i]
  nmbNoVax <- resultsdf$nv_eff_d * wtp - resultsdf$nv_cost_d
  nmbLive <- resultsdf$v_eff_d * wtp - resultsdf$v_cost_d
  nmbRecombo <- Recomboresultsdf$v_eff_d * wtp - Recomboresultsdf$v_cost_d
  nmbMatrix <- cbind(NoVax = nmbNoVax, Live = nmbLive, Recombo = nmbRecombo)
  nmb_list[[i]] <- nmbMatrix
  best_strategy <- apply(nmbMatrix, 1, function(x) names(which.max(x)))
  ceac_df$NoVax[i] <- mean(best_strategy == "NoVax")
  ceac_df$Live[i] <- mean(best_strategy == "Live")
  ceac_df$Recombo[i] <- mean(best_strategy == "Recombo")
}
#----
#conduct one-way deterministic sensitivty analysis and generate tornado plots
#create dataframes for IXCHIQ (live-attenuated) and VIMKUNYA (recombinant) containing
#base-case values, LB, and UB
df_Live <- data.frame(
  variable = c(
    "pbeta", "plambda", "pgamma", "pdelta", "pphi", "ppsi", "pmu", "pkappa", "pomega",
    "discount", "popsize", "fracsusceptible", "fracinfected", "cycle_length",
    "U_S", "U_V", "U_E", "U_I", "U_R", "U_C", "U_D",
    "C_S", "C_V", "C_E", "C_I", "C_R", "C_C", "C_D", "ptheta"
  ),
  base_case = c(
    0.939189937, 0.753403036, 0.56112007, 0.001, 0.991, 0.010445067, 6.95098e-06,
    0.018995353, 0.521, 0.03, 1000, 0.9980005, 0.0019995, 0.019230769,
    0.824, 0.00024, 0.700, 0.662, 0.824, 0.494, 0,
    0, 28.89, 0, 235.44, 0, 685.64, 0, 0.000406
  ),
  LB = c(
    0.813626024, 0.441964854, 0.503414696, 0.00075, 0.975, 0.007844073, 5.21324e-06,
    0.009543213, 0.445, 0.03, 1000, 0.9985, 0.0015, 0,
    0.820, 0.00018, 0.525, 0.593, 0.820, 0.176, 0,
    0, 22.09, 0, 41.75, 0, 514.23, 0, 0.000305
  ),
  UB = c(
    0.985004423, 0.999088118, 0.632120559, 0.00125, 0.998, 0.013039243, 8.68871e-06,
    0.055910442, 0.597, 0.03, 1000, 0.9975, 0.0025, 0,
    0.828, 0.00030, 0.828, 0.721, 0.828, 0.759, 0,
    0, 36.42, 0, 888.37, 0, 857.05, 0, 0.000508
  )
)

df_Recombo <- data.frame(
  variable = c(
    "pbeta", "plambda", "pgamma", "pdelta", "pphi", "ppsi", "pmu", "pkappa", "pomega",
    "discount", "popsize", "fracsusceptible", "fracinfected", "cycle_length",
    "U_S", "U_V", "U_E", "U_I", "U_R", "U_C", "U_D",
    "C_S", "C_V", "C_E", "C_I", "C_R", "C_C", "C_D", "ptheta"
  ),
  base_case = c(
    0.939189937, 0.753403036, 0.56112007, 0.001, 0.978, 0.010445067, 6.95098e-06,
    0.018995353, 0.521, 0.03, 1000, 0.9980005, 0.0019995, 0.019230769,
    0.824, 0.000331, 0.700, 0.662, 0.824, 0.494, 0,
    0, 28.89, 0, 235.44, 0, 685.64, 0, 0.005453
  ),
  LB = c(
    0.813626024, 0.441964854, 0.503414696, 0.00075, 0.972, 0.007844073, 5.21324e-06,
    0.009543213, 0.445, 0.03, 1000, 0.9985, 0.0015, 0,
    0.820, 0.000248, 0.525, 0.593, 0.820, 0.176, 0,
    0, 22.09, 0, 41.75, 0, 514.23, 0, 0.00409
  ),
  UB = c(
    0.985004423, 0.999088118, 0.632120559, 0.00125, 0.983, 0.013039243, 8.68871e-06,
    0.055910442, 0.597, 0.03, 1000, 0.9975, 0.0025, 0,
    0.828, 0.000414, 0.828, 0.721, 0.828, 0.759, 0,
    0, 36.42, 0, 888.37, 0, 857.05, 0, 0.006817
  )
)
wtp <- 7000
#----
#run the function SVEIRD over LB and UB
param_dfL <- df_Live
param_dfR <- df_Recombo
base_paramsL <- setNames(param_dfL$base_case, param_dfL$variable)
base_paramsR <- setNames(param_dfR$base_case, param_dfR$variable)
#Function to calculate the INMB
get_inmb <- function(params) {
  res <- run_SVEIRD5(params)
  v_nmb <- res["v_eff_d"] * wtp - res["v_cost_d"]
  nv_nmb <- res["nv_eff_d"] * wtp - res["nv_cost_d"]
  return(v_nmb - nv_nmb)
}
#initialize an empty datafrane to store results
DSA_resultsL <- data.frame(Parameter = character(), Bound = character(), INMB = numeric())
#loop for each variable calling get_inmb and storing the result, and then adding it to the dataframe
for (i in seq_len(nrow(param_dfL))) {
  param_name <- param_dfL$variable[i]
  params_lb <- base_paramsL
  params_lb[param_name] <- param_dfL$LB[i]
  inmb_lb <- get_inmb(as.list(params_lb))
  params_ub <- base_paramsL
  params_ub[param_name] <- param_dfL$UB[i]
  inmb_ub <- get_inmb(as.list(params_ub))
  DSA_resultsL <- rbind(
    DSA_resultsL,
    data.frame(Parameter = param_name, Bound = "Lower Bound", INMB = inmb_lb),
    data.frame(Parameter = param_name, Bound = "Upper Bound", INMB = inmb_ub)
  )
}
DSA_resultsR <- data.frame(Parameter = character(), Bound = character(), INMB = numeric())
for (i in seq_len(nrow(param_dfR))) {
  param_name <- param_dfR$variable[i]
  params_lb <- base_paramsR
  params_lb[param_name] <- param_dfR$LB[i]
  inmb_lb <- get_inmb(as.list(params_lb))
  params_ub <- base_paramsR
  params_ub[param_name] <- param_dfR$UB[i]
  inmb_ub <- get_inmb(as.list(params_ub))
  DSA_resultsR <- rbind(
    DSA_resultsR,
    data.frame(Parameter = param_name, Bound = "Lower Bound", INMB = inmb_lb),
    data.frame(Parameter = param_name, Bound = "Upper Bound", INMB = inmb_ub)
  )
}
DSA_results_Live <- DSA_resultsL
DSA_results_Recombo <- DSA_resultsR
#change format to wide
DSA_Live_wide <- DSA_results_Live %>% pivot_wider(names_from = Bound, values_from = INMB)
DSA_Recombo_wide <- DSA_results_Recombo %>% pivot_wider(names_from = Bound, values_from = INMB)
#calculate the range of absolute value of (UB - LB)
DSA_Live_wide <- DSA_Live_wide %>% mutate(range = abs(`Upper Bound` - `Lower Bound`))
DSA_Recombo_wide <- DSA_Recombo_wide %>% mutate(range = abs(`Upper Bound` - `Lower Bound`))
DSA_Live_wide <- DSA_Live_wide %>% arrange(desc(range))
#arrange by descending order for range and take only the values that are > 2.5% total INMB
#This happens to be only the 1-13th parameters
DSA_Live_wide <- DSA_Live_wide[1:13,]
DSA_Recombo_wide <- DSA_Recombo_wide %>% arrange(desc(range))
DSA_Recombo_wide <- DSA_Recombo_wide[1:13,]
BaseL <- 759068 #the base INMB
sorted_L <- DSA_Live_wide %>%
  arrange(desc(range)) %>%
  pull(Parameter) #pull the parameter name
plot_dataL <- DSA_Live_wide %>% #turn this make into long format for ggplot
  pivot_longer(cols = c("Lower Bound", "Upper Bound"),
               names_to = "Bound",
               values_to = "Value") %>%
  mutate(
    xmin = pmin(Value, BaseL),
    xmax = pmax(Value, BaseL),
    Parameter = factor(Parameter, levels = rev(sorted_L))
  )
BaseR <- 694106
sorted_R <- DSA_Recombo_wide %>%
  arrange(desc(range)) %>%
  pull(Parameter)
plot_dataR <- DSA_Recombo_wide %>%
  pivot_longer(cols = c("Lower Bound", "Upper Bound"),
               names_to = "Bound",
               values_to = "Value") %>%
  mutate(
    xmin = pmin(Value, BaseR),
    xmax = pmax(Value, BaseR),
    Parameter = factor(Parameter, levels = rev(sorted_R))
  )
#change the labels for the selected parameters
param_labels <- c(
  pkappa = "Chronic Disease Recovery Probability",
  U_C = "Utility Chronic Disease",
  C_C = "Cost Chronic Disease",
  C_I = "Cost Infection",
  ppsi = "Vaccination Rate",
  pomega = "Probability of Chronic Disease",
  plambda = "Infectious Rate",
  pbeta = "Transmission Rate",
  pgamma = "Probability Recover Infectious",
  fracinfected = "Initial Infection Proportion",
  U_I = "Utility Infectious CHIKV",
  U_E = "Utility Exposed",
  U_R = "Utility Recovered",
  U_S = "Utility Susceptible"
)
#----
#tornado
#Live-attenuated
ggplot(plot_dataL) +
  geom_segment(aes(x = xmin, xend = xmax, y = Parameter, yend = Parameter, color = Bound),
               size = 6) +
  geom_vline(xintercept = BaseL, linetype = "dashed", color = "black") +
  scale_color_manual(values = c("Lower Bound" = "#3182bd", "Upper Bound" = "#6baed6")) +
  scale_x_continuous(
    limits = c(0,1.5e6),
    breaks = seq(0, 1.5e6, by = 1e5),
    labels = function(x) sprintf("%.1f", x / 1e5)
  ) +
  scale_y_discrete(labels = param_labels) +
  labs(
    x = "INMB (Hundred Thousands USD)",
    y = "Parameter"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line.x = element_line(color = "black", size = 1),
    axis.line.y = element_line(color = "black", size = 1),
    axis.ticks.x = element_line(color = "black"),
    plot.title = element_text(hjust = 0.5, face = "bold"),
    axis.text.x = element_text(angle = 0, hjust = 0.5, size = 8),
    axis.text.y = element_text(size = 8),
    legend.title = element_blank(),
    legend.position = "bottom"
  )
#recombinant
ggplot(plot_dataR) +
  geom_segment(aes(x = xmin, xend = xmax, y = Parameter, yend = Parameter, color = Bound),
               size = 6) +
  geom_vline(xintercept = BaseR, linetype = "dashed", color = "black") +
  scale_color_manual(values = c("Lower Bound" = "#3182bd", "Upper Bound" = "#6baed6")) +
  scale_x_continuous(
    limits = c(0,1.5e6),
    breaks = seq(0, 1.5e6, by = 1e5),
    labels = function(x) sprintf("%.1f", x / 1e5)
  ) +
  scale_y_discrete(labels = param_labels) +
  labs(
    x = "INMB (Hundred Thousands USD)",
    y = "Parameter"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line.x = element_line(color = "black", size = 1),
    axis.line.y = element_line(color = "black", size = 1),
    axis.ticks.x = element_line(color = "black"),
    plot.title = element_text(hjust = 0.5, face = "bold"),
    axis.text.x = element_text(angle = 0, hjust = 0.5, size = 8),
    axis.text.y = element_text(size = 8),
    legend.title = element_blank(),
    legend.position = "bottom"
  )

#----
#bootstrapping
#adapted from Data Science in R: A Gentle Introduction, James G. Scott, August 2021, https://bookdown.org/jgscott/DSGI/
#take the nmbs from the WTP = 7000 PSA results
nmb7000 <- nmb_list$'7000'
nmb7000 <- as.data.frame(nmb7000)
nmb7000 <- nmb7000 %>% mutate(INMB_LNV = Live - NoVax, INMB_RNV = Recombo - NoVax, INMBLR = Live - Recombo)
#create 10,000 bootstrap samples from the PSA results
bootLCosts <- do(10000)*mean(~v_cost_d, data = mosaic::resample(resultsdf))
bootLEff <- do(10000)*mean(~v_eff_d, data = mosaic::resample(resultsdf))
bootNVCosts <- do(10000)*mean(~nv_cost_d, data = mosaic::resample(resultsdf))
bootNVEff <- do(10000)*mean(~nv_eff_d, data = mosaic::resample(resultsdf))
bootRCosts <- do(10000) * mean(~v_cost_d, data = mosaic::resample(Recomboresultsdf))
bootREff <- do(10000) * mean(~v_eff_d, data = mosaic::resample(Recomboresultsdf))
bootNMBNV <- do(10000)*mean(~NoVax, data = mosaic::resample(nmb7000))
bootNMBR <- do(10000)*mean(~Recombo, data = mosaic::resample(nmb7000))
bootNMBL <- do(10000)*mean(~Live, data = mosaic::resample(nmb7000))
bootINMB_LNV <- do(10000)*mean(~INMB_LNV, data = mosaic::resample(nmb7000))
bootINMB_RNV <- do(10000)*mean(~INMB_RNV, data = mosaic::resample(nmb7000))
bootINMBLR <- do(10000) * mean(~INMBLR, data = mosaic::resample(nmb7000))
#compute confidence intervals of the 10,000 bootstrap samples
CI_NV_costs <- confint(bootNVCosts, level = 0.95)
CI_NV_effects <- confint(bootNVEff, level = 0.95)
CI_R_costs <- confint(bootRCosts, level = 0.95)
CI_R_effects <- confint(bootREff, level = 0.95)
CI_L_costs <- confint(bootLCosts, level = 0.95)
CI_L_effects <- confint(bootLEff, level = 0.95)
CI_NMB_NV <- confint(bootNMBNV, level = 0.95)
CI_NMB_R <- confint(bootNMBR, level = 0.95)
CI_NMB_L <- confint(bootNMBL, level = 0.95)
CI_INMB_LNV <- confint(bootINMB_LNV, level = 0.95)
CI_INMB_RNV <- confint(bootINMB_RNV, level = 0.95)
CI_INMBLR <- confint(bootINMBLR, level = 0.95)

