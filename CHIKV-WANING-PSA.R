#----
#this script attempts to model the impact of vaccine waning
#V no longer can go to E, everyone in V stays in V until their immunity wanes at rate theta
#Theta is derived from extrapolating a linear decline of 1 month sero response to 0
#Individuasl go to VE from V representing waning or from S representing those in which the vaccine did not work
#These people progress at rate (1-phi)*psi*S
#----
#packages
library(tidyverse)
library(ggplot2)
library(mc2d)
library(scales)
library(parallel)
library(dplyr)
library(tibble)
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
#----
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
  ptheta = Theta_I_draws,
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
Vimkunyaparametersdf <- parametersdf
Vimkunyaparametersdf$pphi <- pPhi_Vim_draws
Vimkunyaparametersdf$U_V <- U_VIM_draws
Vimkunyaparametersdf$C_V <- C_Vim_draws
Vimkunyaparametersdf$ptheta <- Theta_V_draws
#----
run_SVEIRD5 <- function(params, return_trace = FALSE) {
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

  U_S <- params[["U_S"]] * cycle_length
  U_E <- params[["U_E"]] * cycle_length
  U_V <- params[["U_V"]]
  U_I <- params[["U_I"]] * cycle_length
  U_R <- params[["U_R"]] * cycle_length
  U_C <- params[["U_C"]] * cycle_length
  U_D <- 0

  C_S <- params[["C_S"]]
  C_E <- params[["C_E"]]
  C_V <- params[["C_V"]]
  C_I <- params[["C_I"]]
  C_R <- params[["C_R"]]
  C_C <- params[["C_C"]]
  C_D <- 0
  names <- c("S", "E", "V", "VE", "I", "R", "C", "D", "N", "Check",
             "SV", "SE", "SVE", "SD", "VVE", "VD", "VEE", "VED", "EI", "ED", "IR", "IC", "ID", "RD", "CR", "CD")
  trace_v <- data.frame(matrix(0, nrow = 520, ncol = length(names)))
  colnames(trace_v) <- names

  #initialize trace
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


  for(i in 2:length(trace_v$S)){
    trace_v$SV[i] <- trace_v$S[i-1] * ppsi * pphi
    trace_v$SVE[i] <- trace_v$S[i-1] * ppsi * (1-pphi)
    trace_v$SE[i] <- trace_v$S[i-1] * pbeta * trace_v$I[i-1] / trace_v$N[i-1]
    trace_v$SD[i] <- trace_v$S[i-1] * pmu
    trace_v$VVE[i] <- trace_v$V[i-1] * ptheta
    trace_v$VED[i] <- trace_v$VE[i-1] * pmu
    trace_v$VEE[i] <- trace_v$VE[i-1] * pbeta * trace_v$I[i-1] / trace_v$N[i-1]
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

  for(i in 2:length(trace_nv$S)){
    trace_nv$SV[i] <- 0
    trace_nv$SVE[i] <- 0
    trace_nv$SE[i] <- trace_nv$S[i-1] * pbeta * trace_nv$I[i-1] / trace_nv$N[i-1]
    trace_nv$SD[i] <- trace_nv$S[i-1] * pmu
    trace_nv$VVE[i] <- trace_nv$V[i-1] * ptheta
    trace_nv$VED[i] <- trace_nv$VE[i-1] * pmu
    trace_nv$VEE[i] <- trace_nv$VE[i-1] * pbeta * trace_nv$I[i-1] / trace_nv$N[i-1]
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

  utility_vector <- c(U_S, U_E, U_S, U_S, U_I, U_R, U_C, U_D)
  cost_vector <- c(C_S, C_E, C_V, C_V, C_I, C_R, C_C, C_D)
  utility_trace_v <- trace_v[,1:8]
  utility_trace_v <- apply(trace_v[,1:8], 1, function(row) row * utility_vector)
  utility_trace_v <- t(utility_trace_v)
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
  for(i in 1:nrow(cost_trace_v)) {
    cost_trace_v[i,3] <- trace_v$SV[i] * C_V
    cost_trace_v[i,4] <- trace_v$SVE[i] * C_V
  }
  trace_nv$V <- trace_nv$SV
  cost_trace_nv <- trace_nv[,1:8]
  cost_trace_nv <- apply(trace_nv[,1:8], 1, function(row) row * cost_vector)
  cost_trace_nv <- t(cost_trace_nv)
  discount_factors <- 1 / (1 + discount) ^ ((0:(520 - 1)) / 52)
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
  return(c(v_cost_d = v_cost_d, nv_cost_d = nv_cost_d, v_eff_d = v_eff_d, nv_eff_d = nv_eff_d))
}
#----
#multicore processing
num_cores <- detectCores() - 1
cl <- makeCluster(num_cores)  # Create parallel cluster

clusterExport(cl, varlist = c("parametersdf", "run_SVEIRD5"))

clusterEvalQ(cl, { library(dplyr); library(tibble) })

# Run parallel execution
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
num_cores <- detectCores() - 1
cl <- makeCluster(num_cores)  # Create parallel cluster

clusterExport(cl, varlist = c("Vimkunyaparametersdf", "run_SVEIRD5"))

clusterEvalQ(cl, { library(dplyr); library(tibble) })

# Run parallel execution
Vimresults_list <- parLapply(cl, 1:1000, function(i) {
  # Load required libraries inside worker
  library(dplyr)
  library(tibble)

  run_SVEIRD5(Vimkunyaparametersdf[i, ])
})

# Stop cluster
stopCluster(cl)

# Convert list to dataframe
Vimresultsdf <- do.call(rbind, Vimresults_list)
Vimresultsdf <- as.data.frame(Vimresultsdf)
#----
#NMB calculations
wtp_values <- seq(1000,150000, by = 1000)
ceac_df <- data.frame(WTP = wtp_values, NoVax = NA, IXCHIQ = NA, Vimkunya = NA)

nmb_list <- vector("list", length(wtp_values))
names(nmb_list) <- wtp_values

for (i in seq_along(wtp_values)) {
  wtp <- wtp_values[i]
  nmbNoVax <- resultsdf$nv_eff_d * wtp - resultsdf$nv_cost_d
  nmbIXCHIQ <- resultsdf$v_eff_d * wtp - resultsdf$v_cost_d
  nmbVIM <- Vimresultsdf$v_eff_d * wtp - Vimresultsdf$v_cost_d
  nmbMatrix <- cbind(NoVax = nmbNoVax, IXCHIQ = nmbIXCHIQ, Vimkunya = nmbVIM)
  nmb_list[[i]] <- nmbMatrix
  best_strategy <- apply(nmbMatrix, 1, function(x) names(which.max(x)))
  ceac_df$NoVax[i] <- mean(best_strategy == "NoVax")
  ceac_df$IXCHIQ[i] <- mean(best_strategy == "IXCHIQ")
  ceac_df$Vimkunya[i] <- mean(best_strategy == "Vimkunya")
}
#----
#EVPI
EVPI <- numeric(length(wtp_values))
for (i in seq_along(wtp_values)) {
  nmbMatrix <- nmb_list[[i]]
  E_max_NMB <- mean(apply(nmbMatrix, 1, max))
  max_E_NMB <- max(colMeans(nmbMatrix))
  EVPI[i] <- E_max_NMB - max_E_NMB
}
EVPI_df <- data.frame(WTP = wtp_values, EV_PerfectInfo = EVPI)
#----
#plots
ceac_renamed <- ceac_df %>%
  rename(`No Vaccination` = NoVax) %>% rename(`Vaccination with live-attenuated CHIKV vaccine` = IXCHIQ) %>% rename(`Vaccination with recombinant CHIKV vaccine` = Vimkunya)
ceac_long <- pivot_longer(ceac_renamed, cols = -WTP, names_to = "Strategy", values_to = "Probability")
n_sim <- 1000
ceac_long <- ceac_long %>%
  mutate(
    SE = sqrt(Probability * (1-Probability) / n_sim),
    Lower = pmax(Probability - 1.96* SE,0),
    Upper = Probability + 1.96 * SE
  )
#----
#standalone PSA
library(scales)  # for comma()

ggplot(ceac_long, aes(x = WTP, y = Probability, color = Strategy, fill = Strategy)) +
  geom_smooth(se = FALSE, method = "loess", span = 0.85, size = 1.2) +
  geom_ribbon(aes(ymin = Lower, ymax = Upper), alpha = 0.15, color = NA) +
  scale_x_continuous(labels = comma) +  # thousands separator
  labs(
    x = "Willingness to Pay (USD)",
    y = "Probability Cost-Effective"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    legend.position = c(0.75, 0.75),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line.x = element_line(color = "black", size = 1),
    axis.line.y = element_line(color = "black", size = 1),
    axis.ticks.x = element_line(color = "black"),
    axis.ticks.y = element_line(color = "black")
  )
#----
#overlayed PSA and EVPI
ceac_df$EVPI_scaled <- EVPI_df$EV_PerfectInfo / 50000
ggplot(ceac_long, aes(x = WTP, y = Probability, color = Strategy)) +
  geom_line(size = 1.2) +
  geom_ribbon(aes(ymin = Lower, ymax = Upper, fill = Strategy),
              alpha = 0.3, color = NA) +
  geom_line(data = ceac_df,
            aes(x = WTP, y = EVPI_scaled, color = "EVPI"),
            size = 1.2) +
  geom_segment(aes(x = min(ceac_df$WTP), xend = max(ceac_df$WTP),
                   y = 0.5, yend = 0.5),
               linetype = "dashed", color = "grey", size = 1) +
  geom_segment(aes(x = 7000, xend = 7000,
                   y = 0, yend = 1),
               linetype = "dashed", color = "grey", size = 1) +
  scale_y_continuous(
    name = "Probability of being cost-effective",
    limits = c(0, 1),
    expand = c(0, 0),
    sec.axis = sec_axis(trans = function(x) { x * 50000 }, name = "EVPI (USD)", labels = scales::comma)
  ) +
  scale_x_continuous(
    limits = c(0, max(ceac_df$WTP)),
    expand = c(0, 0),
    labels = scales::dollar_format()
  ) +
  scale_color_manual(
    values = c(
      "EVPI" = "black",
      "No Vaccination" = "red",
      "Vaccination with live-attenuated CHIKV vaccine" = "blue",
      "Vaccination with recombinant CHIKV vaccine" = "green"
    ),
    breaks = c(
      "No Vaccination",
      "Vaccination with live-attenuated CHIKV vaccine",
      "Vaccination with recombinant CHIKV vaccine",
      "EVPI"
    )
  ) +
  scale_fill_manual(
    values = c(
      "No Vaccination" = "red",
      "Vaccination with live-attenuated CHIKV vaccine" = "blue",
      "Vaccination with recombinant CHIKV vaccine" = "green"
    ),
    guide = "none"
  ) +
  labs(
    x = "Willingness-to-Pay Threshold",
    color = "Legend",
    fill = "Strategy"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    legend.position = c(0.5, 0.75),
    axis.title.y.right = element_text(color = "black"),
    axis.text.y.right = element_text(color = "black"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.ticks.x = element_line(color = "grey"),
    axis.ticks.y = element_line(color = "grey"),
    axis.line.x = element_line(color = "grey", size = 1),
    axis.line.y = element_line(color = "grey", size = 1)
  )



#----
#scatter
resultsdf <- resultsdf %>% mutate(IC = v_cost_d - nv_cost_d)
resultsdf <- resultsdf %>% mutate(IE = v_eff_d - nv_eff_d)
ggplot(resultsdf, aes(x = IE, y = IC)) +
  geom_point(alpha = 0.5, color = "blue") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "red") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  labs(
    title = "Cost-Effectiveness Plane",
    x = "Incremental Effectiveness",
    y = "Incremental Cost"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 12),
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 10)
  ) +
  scale_y_continuous(labels = scales::comma)

