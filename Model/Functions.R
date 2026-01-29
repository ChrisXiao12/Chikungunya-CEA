##****************************************************************************************
## FUNCTIONS
##****************************************************************************************

##*********************************************************************
## Utilities
##*********************************************************************
source(file = file.path(model_dir, "Utilities_Functions.R"))

##*********************************************************************
## Function to update the parameters for FX and Inflation
##*********************************************************************
params_adj_function <- function(
  pars, 
  fx_adjust, 
  cpi_adjust
){
  val <- pars
  ## fx adjustment
  val[names(fx_adjust)] <- as.numeric(val[names(fx_adjust)]) * as.numeric(fx_adjust[names(fx_adjust)]) 
  ## inflation adjustment
  val[names(cpi_adjust)] <- as.numeric(val[names(cpi_adjust)]) * as.numeric(cpi_adjust[names(cpi_adjust)])
  ## output
  return(val)  
}

##*********************************************************************
## Function for seroconversion rate parameters
##*********************************************************************
vac_seroresponse_func <- function(
  data = vac_seroresp_data,
  fit_model = c("expo", "poly")[1],
  facet_labels = c(
    "IXCHIQ" = "(a) Live-attenuated vaccine",
    "VIMKUNYA" = "(b) Recombinant vaccine"
  ),
  create_coeff_table = FALSE,
  print_coef_plot = FALSE,
  save_coeffs_plot = FALSE,
  figures_dir = figures_dir,
  tables_dir = tables_dir
){

  ## model fit
  df_vaccines <- data
  fit_model <- fit_model

  results <- fit_seroprotection(
    data = df_vaccines,
    model = fit_model,
    degree = 3,
    B = 1000,
    xlab = "Time since vaccination (Days)",
    ylab = expression(
      "Seroprotection/seroresponse Rate Differences," *
      epsilon[list(plain(v),plain(t))] *
      " (%)"
    ),
    title = NULL,
    facet_labs = facet_labels,
    strip_txt_size = 12
  )

  ## print plot
  if(isTRUE(print_coef_plot)){
    print(results$plot)
  }

  ## save coefficients plot
  if(isTRUE(save_coeffs_plot)){
    ## save the plot as a PNG file
    ggsave(
      filename = "seroprotection.png",
      plot = results$plot,
      path = figures_dir,
      width = 9, height = 5,
      units = c("in", "cm", "mm", "px")[1], 
      dpi = 300, bg = "white"
    )
  }

  ## create table of coefficients
  if(isTRUE(create_coeff_table)){
    ## map coefficient names
    if(fit_model == "expo"){
      latex_map <- c(
        "a" = "\\alpha",
        "b" = "\\theta"
      )
    }else{
      latex_map <- c(
        "alpha" = "\\alpha",
        "beta1" = "\\beta_1",
        "beta2" = "\\beta_2",
        "beta3" = "\\beta_3"
      )
    }
    
    coeffs_df <- as.data.frame(results$coeffs_df) %>%
      mutate(
        Coefficient = latex_map[coefficient],
        Median_fmt = sprintf("%.6f", median),
        Lower95_fmt = sprintf("%.6f", lower),
        Upper95_fmt = sprintf("%.6f", upper)
      ) %>%
      dplyr::rename(Vaccine = vaccine) %>%
      dplyr::select(
        Vaccine,
        Coefficient,
        Median = Median_fmt,
        `Lower 95% CI` = Lower95_fmt,
        `Upper 95% CI` = Upper95_fmt
      )
    
    ## create flextable
    ft <- flextable(coeffs_df) %>%
      theme_vanilla() %>%
      fontsize(size = 11, part = "all") %>%
      font(fontname = "Times New Roman", part = "all") %>%
      autofit() %>%      
      ## remove all borders
      border_remove() %>%
      ## add top and bottom table borders
      hline_top(border = fp_border(color = "black", width = 1), part = "all") %>%
      hline_bottom(border = fp_border(color = "black", width = 1), part = "all") %>%      
      ## add header underline only
      hline_bottom(border = fp_border(color = "black", width = 1), part = "header") %>%      
      ## remove all horizontal lines in the body
      hline(border = fp_border(color = "transparent"), part = "body") %>%
      ## alignment
      align(j = 1:5, align = "left", part = "all") %>%
      bold(part = "header")
    
    ## save to docx
    doc <- read_docx()
    doc <- body_add_flextable(doc, ft)
    print(doc, target = file.path(tables_dir, "seroprotection_coefficients.docx"))

  }
  
  ## output
  return(results)
}

##############################################
## Plot observed + bootstrap predicted
##############################################
## exponential fit
fit_exp <- function(y = value, time = days) {
  df <- data.frame(y = y, time = time)
  a0 <- max(y)-min(y)
  b0 <- 0.0001
  nlsLM(
    y ~ a * exp(-b * time),
    data = df,
    start = list(a = a0, b = b0),
    control = nls.lm.control(maxiter = 500)
  )
}

## polynomial fit
fit_poly <- function(y = value, time = days, degree = 2){
  df <- data.frame(y = y, time = time)
  lm(y ~ poly(time, degree, raw = TRUE), data = df)
}

## waning rate function: derivative of the polynomial
waning_function <- function(time = 100, coeffs, degree = 3){
  sum(sapply(1:degree, function(k) {k * coefs[k+1] * time^(k-1)}))
} 

## utility functions
logit <- function(p) log(p / (1 - p))
invlogit <- function(x) exp(x) / (1 + exp(x))

estimate_rho <- function(x) {
  x_centered <- x - mean(x)
  num <- sum(x_centered[-length(x)] * x_centered[-1])
  den <- sum(x_centered^2)
  num / den
}

## estimate the coefficients of the fit to the seroprotection data
fit_seroprotection <- function(
  data, 
  model = c("expo", "poly"),
  degree = 3, 
  B = 1000, 
  seed = sim_seed,
  xlab = "Time since vaccination (Days)", 
  ylab = "Seroprotection rate (%)", 
  title = NULL,
  facet_labs = facet_labels,
  strip_txt_size = 10
) {

  set.seed(seed)

  coeffs_df <- NULL
  pred_df <- NULL
  waning_df <- NULL
  
  for(v in unique(data$vaccine)) {

    df_vac <- data %>% filter(vaccine == v)
    days <- df_vac$days
    value <- df_vac$value
    lower <- df_vac$lower
    upper <- df_vac$upper
    
    ## logit transform
    mu <- logit(value)
    sd <- (logit(upper) - logit(lower)) / (2 * 1.96)
    
    ## AR(1) correlation
    rho <- estimate_rho(mu)
    Sigma <- outer(
      1:length(days), 
      1:length(days), 
      function(i, j){
        sd[i] * sd[j] * rho^(abs(i - j))
      }
    )
    
    ## Fit base model
    if(model == "poly"){
      fit_mid <- lm(value ~ poly(days, degree, raw = TRUE))
      pars_names <- c("alpha", paste0("beta", seq(1, degree, 1)))
    } else {
      fit_mid <- fit_exp(y = value, time = days)
      pars_names <- names(coef(fit_mid))
    }

    ## Initialize matrices
    pred_mat_vac <- matrix(NA, nrow = B, ncol = length(days))
    coeffs_mat_vac <- matrix(NA, nrow = B, ncol = length(pars_names))
    colnames(coeffs_mat_vac) <- pars_names

    ## Bootstrapping
    for(i in 1:B) {

      ## simulate data
      y_sim_logit <- mvrnorm(1, mu = mu, Sigma = Sigma)
      y_sim <- invlogit(y_sim_logit)

      ## fit simulated model
      if(model == "poly"){
        fit_sim <- try(lm(y_sim ~ poly(days, degree, raw = TRUE)), silent = TRUE)
      } else {
        fit_sim <- try(fit_exp(y = y_sim, time = days), silent = TRUE)
      }

      if(!inherits(fit_sim, "try-error")) {
        pred_mat_vac[i, ] <- predict(fit_sim, newdata = data.frame(days = days))
        coeffs_mat_vac[i, ] <- coef(fit_sim)
      }
    }

    ## remove failed fits
    good_rows <- rowSums(is.na(coeffs_mat_vac)) == 0 & rowSums(is.na(pred_mat_vac)) == 0
    coeffs_mat_vac <- coeffs_mat_vac[good_rows, , drop = FALSE]
    pred_mat_vac   <- pred_mat_vac[good_rows, , drop = FALSE]

    ## rename polynomial coefficients → alpha, beta1, beta2, …
    if(model == "poly"){
      orig_names <- colnames(coeffs_mat_vac)
      new_names <- c("alpha", paste0("beta", 1:(length(orig_names)-1)))
      colnames(coeffs_mat_vac) <- new_names
    }

    ## summarize coefficients
    coeffs_df_vac <- data.frame(
      vaccine = v,
      coefficient = colnames(coeffs_mat_vac),
      median = apply(coeffs_mat_vac, 2, median),
      lower = apply(coeffs_mat_vac, 2, quantile, 0.025),
      upper = apply(coeffs_mat_vac, 2, quantile, 0.975)
    )
    coeffs_df <- rbind(coeffs_df, coeffs_df_vac)
    
    ## summarize predictions
    pred_df_vac <- data.frame(
      vaccine = v,
      days = days,
      median = apply(pred_mat_vac, 2, median),
      lower = apply(pred_mat_vac, 2, quantile, 0.025),
      upper = apply(pred_mat_vac, 2, quantile, 0.975)
    )
    pred_df <- rbind(pred_df, pred_df_vac)

    ## waning (derivative)
    waning_rate_mat_vac <- t(apply(coeffs_mat_vac, 1, function(coef_vec){
      if(model == "poly"){
        sapply(days, function(d){
          sum(sapply(1:degree, function(k) k * coef_vec[k+1] * d^(k-1)))
        })
      } else {
        rep(coef_vec["b"], length(days))
      }
    }))

    waning_df_vac <- data.frame(
      vaccine = v,
      days = days,
      median = apply(waning_rate_mat_vac, 2, median),
      lower = apply(waning_rate_mat_vac, 2, quantile, 0.025),
      upper = apply(waning_rate_mat_vac, 2, quantile, 0.975)
    )
    waning_df <- rbind(waning_df, waning_df_vac)
  }
  
  ## plot object
  seroprotection_p <- ggplot() +
    geom_ribbon(
      data = pred_df, 
      aes(x = days, ymin = lower, ymax = upper, fill = "95% Prediction Interval"), 
      alpha = 0.3
    ) +
    geom_line(
      data = pred_df, 
      aes(x = days, y = median, color = "Median fitted curve"), 
      linewidth = 1.2
    ) +
    geom_point(
      data = data, 
      aes(x = days, y = value, color = "Observed data", shape = "Observed data"), 
      size = 3
    ) +
    geom_errorbar(
      data = data, 
      aes(x = days, ymin = lower, ymax = upper, color = "Observed data"), 
      width = 20, 
      linewidth = 0.9
    ) +
    facet_wrap(. ~ vaccine, ncol = 2, scales = "free", labeller = labeller(vaccine = facet_labs)) +
    scale_color_manual(name = NULL, values = c("Median fitted curve" = "#0072B2", "Observed data" = "#000000")) +
    scale_fill_manual(values = c("95% Prediction Interval" = "#56B4E9")) +
    scale_shape_manual(values = c("Observed data" = 16)) +
    guides(fill = "none", shape = "none") +
    scale_y_continuous(limits = c(0, 1.1), breaks = seq(0, 1, 0.25), labels = seq(0, 1, 0.25) * 100) +
    labs(x = xlab, y = ylab, title = title) +
    theme_bw(base_size = 14) +
    theme(
      legend.position = "bottom",
      legend.background = element_blank(),
      plot.title = element_text(face = "bold", size = 16),
      axis.title.y = element_text(size = 12),
      axis.title.x = element_text(size = 12),
      strip.text = element_text(face = "bold", size = strip_txt_size),
      strip.background = element_rect(fill = NA, color = NA)
    )
  
  list(
    coeffs_df = coeffs_df,
    pred_df = pred_df,
    waning_df = waning_df,
    plot = seroprotection_p
  )
}

##*********************************************************************
## Function to generate the initial population distribution
##*********************************************************************
inipop_func <- function(
  states = state_names, 
  params = params_value
) {
  with(as.list(params), {
    incidence_prob <- 1 - exp(-( (total_inf_cases / base_population) / 365) * cycle_length_days)
    val <- matrix(
      data = 0, 
      nrow = 1, 
      ncol = length(states), 
      dimnames = list(NULL, states)
    )
    val[,"S"] <- (1 - incidence_prob) * n_cohort
    val[,"I"] <- incidence_prob * n_cohort
    return(val)
  })
}


##*********************************************************************
## Markov model
##*********************************************************************
markov_model_func <- function(
  scenario = scenarios[1],
  params = params_value,
  mu,
  xi_factor = NULL,
  zeta_factor = NULL, 
  phi_factor = NULL,
  psi_stated_factor = NULL,
  pi_factor = NULL,
  chi_factor = NULL,
  eta_factor = NULL,
  omicron_factor = NULL,
  init_states = ini_pop_vec,
  inc_names = trans_names,
  n_cycles = NULL
){
  

  ## cycle length
  if(is.null(n_cycles)){
    n_cycles <- as.numeric(params[["time_horizon_years"]]) * 52
  } else {
    n_cycles <- n_cycles
  }
  
  ## prepare numeric indices and constants
  states_names <- colnames(init_states)
  n_states <- length(states_names)
  n_inc <- length(inc_names)
  alive_idx <- which(states_names != "D")
  D_idx <- which(states_names == "D")
  S_idx <- which(states_names == "S")
  V_idx <- which(states_names == "V")
  SV_idx <- which(states_names == "SV")
  E_idx <- which(states_names == "E")
  I_idx <- which(states_names == "I")
  R_idx <- which(states_names == "R")
  C_idx <- which(states_names == "C")

  ## turn params to a simple named vector for fast access
  p <- as.list(params)
  
  ## Precompute cycle-based transition scalars (these don't change inside loop)
  cycle_length_days <- as.numeric(p$cycle_length_days)

  ## xi
  if(is.null(xi_factor)){
    xi <- as.numeric(p$xi)
  } else {
    xi <- xi_factor
  }
  if(scenario == "SOC") {xi <- 1} 

  ## betas, accounting for behavioural response
  beta_S <- 1 - exp(-((1 / as.numeric(p$gamma_days)) * as.numeric(p$R0)) * cycle_length_days)
  beta_SV <- beta_S * xi

  ## lambda
  lambda <- 1 - exp(-(1 / as.numeric(p$lambda_days)) * cycle_length_days)

  ## kappa
  kappa <- 1 - exp(-(1 / as.numeric(p$kappa_days)) * cycle_length_days)
  
  ## delta
  delta_hosp <- 1 - exp(-as.numeric(p$delta_hosp_10years) / (10 * 52))
  delta_pop <- 1 - exp(-(as.numeric(p$deaths_CHIKV_pop) / as.numeric(p$CHIKV_pop)) / 52)
  delta <- as.numeric(p$prop_hosp) * delta_hosp + (1 - as.numeric(p$prop_hosp)) * delta_pop

  ## pi
  if (scenario == "VIMKUNYA") {
    pi <- if (is.null(pi_factor)) p$pi_VIMKUNYA_vs_IXCHIQ else pi_factor
  } else if (scenario == "IXCHIQ") {
    pi <- 1
  } else {
    pi <- 0
  }

  ## phi: Realization factor for vaccination
  if(is.null(phi_factor)){
    phi <- as.numeric(p$phi)
  } else {
    phi <- phi_factor
  }

  ## omicron
  if(is.null(omicron_factor)){
    omicron <- as.numeric(p$omicron)
  } else {
    omicron <- omicron_factor
  }

  ## psi
  if(is.null(psi_stated_factor)){
    psi_stated <- as.numeric(p$psi)
  } else {
    psi_stated <- psi_stated_factor
  }
  psi_eff <- max(c(0, min(c(1, phi * max(c(0, min(c(1, omicron * psi_stated))))))))
  psi <- max(c(0, min(c(1, (1 - (1 - psi_eff)^(1 / 3)) * pi))))

  ## epsilon
  epsilon_IXCHIQ_180_days <- as.numeric(p$seroresponse_treatment_IXCHIQ_180_days) - as.numeric(p$seroresponse_placebo_IXCHIQ_180_days)
  if (scenario == "IXCHIQ") {
    epsilon <- 1 - exp(-(epsilon_IXCHIQ_180_days / 180) * cycle_length_days)
  } else if (scenario == "VIMKUNYA") {
    epsilon <- 1 - exp(-(as.numeric(p$seroprotection_VIMKUNYA_183_days) / 183) * cycle_length_days)
  } else {
    epsilon <- 0
  }
  
  if(is.null(zeta_factor)){
    epsilon <- epsilon
  } else {
    epsilon <- epsilon * zeta_factor
  }

  ## theta
  theta_IXCHIQ <- 1 - exp(-(as.numeric(p$theta_IXCHIQ_days)) * cycle_length_days)
  theta_VIMKUNYA <- 1 - exp(-(as.numeric(p$theta_VIMKUNYA_days)) * cycle_length_days)  
  theta <- if (scenario == "IXCHIQ") theta_IXCHIQ else if (scenario == "VIMKUNYA") theta_VIMKUNYA else 0

  ## calibration of CCA and recovery removal rates
  rho_year <- as.numeric(p$rho_year)

  gamma_d <- 1 / as.numeric(p$gamma_days)
  gamma_CCA_d <- (rho_year * gamma_d) / (1 - exp(-gamma_d * 365))
  gamma_rec_d <- gamma_d - gamma_CCA_d

  gamma_CCA  <- 1 - exp(-7 * gamma_CCA_d)
  gamma_rec  <- 1 - exp(-7 * gamma_rec_d)

  ## base transition matrix M with zeros then fill only non-zero off-diagonal probs
  M_base <- matrix(0, n_states, n_states)
  rownames(M_base) <- colnames(M_base) <- states_names

  ## fill fixed off-diagonal rates that do not depend on I or mu
  ## S -> V and SV and D (D component uses mu[t] so set partial here)
  M_base[S_idx, V_idx] <- epsilon * psi
  M_base[S_idx, SV_idx] <- (1 - epsilon) * psi

  ## V -> SV
  M_base[V_idx, SV_idx] <- theta

  ## E -> I
  M_base[E_idx, I_idx] <- lambda

  ## I -> R and C (I->D depends on mu + delta)
  M_base[I_idx, R_idx] <- gamma_rec
  M_base[I_idx, C_idx] <- gamma_CCA

  ## C -> R
  M_base[C_idx, R_idx] <- kappa

  ## D absorbing
  M_base[D_idx, D_idx] <- 1.0

  ## rows that have D but mu[t] controls magnitude: we'll add mu[t] dynamically
  ## precompute which rows need to have mu[t] added to column D
  rows_with_mu_to_D <- c(S_idx, V_idx, SV_idx, E_idx, R_idx, C_idx)

  ## prepare pop matrix (keep numeric)
  pop <- matrix(
    data = 0, 
    nrow = n_cycles + 1, 
    ncol = n_states + n_inc, 
    dimnames = list(NULL, c(states_names, inc_names))
  )
  pop[1, states_names] <- as.numeric(init_states[1, states_names]) 

  ## convenience: store columns indices for incidence columns
  inc_col_idx <- seq(
    from = (n_states + 1), 
    to = (n_states + n_inc), 
    by = 1
  )
  names(inc_col_idx) <- inc_names

  ## main loop: only update necessary rows and compute next population via matrix multiplication
  for (t in seq_len(n_cycles)) {
    ## current population
    cur_pop_states <- pop[t, 1:n_states]
    I <- cur_pop_states[I_idx]
    N <- sum(cur_pop_states[alive_idx])
    
    ## avoid division by zero
    if (!is.na(N) && N > 0) {
      force_infection_S <- beta_S * (I / N) 
      force_infection_SV <- beta_SV * (I / N) 
    } else {
      force_infection_S <- 0
      force_infection_SV <- 0
    }

    ## copy M_base into Mt but only rows that are changed will be edited
    Mt <- M_base 

    ## add mu[t] to D probabilities for specified rows
    mu_t <- mu[t]
    Mt[rows_with_mu_to_D, D_idx] <- mu_t

    ## I -> D is mu + delta
    Mt[I_idx, D_idx] <- mu_t + delta

    ## S->E and SV->E depend on I/N
    Mt[S_idx, E_idx] <- force_infection_S
    Mt[SV_idx, E_idx] <- force_infection_SV

    ## now, set diagonal elements for rows other than D as 1 - sum(other columns)
    ## for rows S, V, SV, E, I, R, C compute diag:
    rows_to_normalize <- c(S_idx, V_idx, SV_idx, E_idx, I_idx, R_idx, C_idx)
    
    ## compute row sums of off-diagonal (exclude diagonal)
    ## off_diag_sum = sum of row except diagonal
    off_sum <- rowSums(Mt[rows_to_normalize, , drop = FALSE]) - diag(Mt)[rows_to_normalize]
    diag_vals <- 1 - off_sum
    
    ## ensure non-negative
    diag_vals[diag_vals < 0] <- 0
    Mt[cbind(rows_to_normalize, rows_to_normalize)] <- diag_vals

    ## ensure D row stays absorbing
    Mt[D_idx, ] <- 0
    Mt[D_idx, D_idx] <- 1

    ## now multiply: next_state = cur_pop_states %*% Mt
    next_states <- as.vector(cur_pop_states %*% Mt)

    ## incident flows (use the Mt entries for the relevant transitions)
    ## StoV StoSV StoE StoD
    pop[t + 1, "StoV"]  <- cur_pop_states[S_idx] * Mt[S_idx, V_idx]
    pop[t + 1, "StoSV"] <- cur_pop_states[S_idx] * Mt[S_idx, SV_idx]
    pop[t + 1, "StoE"]  <- cur_pop_states[S_idx] * Mt[S_idx, E_idx]
    pop[t + 1, "StoD"]  <- cur_pop_states[S_idx] * Mt[S_idx, D_idx]

    pop[t + 1, "VtoSV"] <- cur_pop_states[V_idx] * Mt[V_idx, SV_idx]
    pop[t + 1, "VtoD"]  <- cur_pop_states[V_idx] * Mt[V_idx, D_idx]

    pop[t + 1, "SVtoE"] <- cur_pop_states[SV_idx] * Mt[SV_idx, E_idx]
    pop[t + 1, "SVtoD"] <- cur_pop_states[SV_idx] * Mt[SV_idx, D_idx]

    pop[t + 1, "EtoI"] <- cur_pop_states[E_idx] * Mt[E_idx, I_idx]
    pop[t + 1, "EtoD"] <- cur_pop_states[E_idx] * Mt[E_idx, D_idx]

    pop[t + 1, "ItoR"] <- cur_pop_states[I_idx] * Mt[I_idx, R_idx]
    pop[t + 1, "ItoC"] <- cur_pop_states[I_idx] * Mt[I_idx, C_idx]
    pop[t + 1, "ItoD"] <- cur_pop_states[I_idx] * Mt[I_idx, D_idx]

    pop[t + 1, "RtoD"] <- cur_pop_states[R_idx] * Mt[R_idx, D_idx]

    pop[t + 1, "CtoR"] <- cur_pop_states[C_idx] * Mt[C_idx, R_idx]
    pop[t + 1, "CtoD"] <- cur_pop_states[C_idx] * Mt[C_idx, D_idx]

    ## fill next state columns
    pop[t + 1, 1:n_states] <- next_states

    ## re-add deaths to S as in original (I think you intended roll-back of deceased into S)
    ## Note: original code added all death incident counts into S (maybe representing replacement births)
    death_inc_sum <- sum(pop[t + 1, c("StoD", "VtoD", "SVtoD", "EtoD", "ItoD", "RtoD", "CtoD")])
    pop[t + 1, S_idx] <- pop[t + 1, S_idx] + death_inc_sum
  }

  ## return as data.frame and add cycle and scenario columns
  out <- as.data.frame(pop)
  out$cycle <- seq(0, n_cycles)
  out$scenario <- scenario
  return(out)
}

##*********************************************************************
## Function to calculate productivity cost due to premature mortality
##*********************************************************************
pv_human_capital_weekly <- function(
  age_val,
  retirement_age = 65,
  monthly_wage,
  le_data = le_death_data,
  discount_rate = 0.03,  
  wage_growth = 0,
  employment_rate = 1,
  adjust_for_employment = TRUE,
  years_rounding = c("floor", "ceiling", "none")[3]
) {
  
  ## year
  years_rounding <- match.arg(years_rounding)

  ## convert inputs to weekly equivalents
  weeks_per_year <- 52
  weekly_wage <- (monthly_wage * 12) / weeks_per_year
  weekly_r <- (1 + discount_rate)^(1 / weeks_per_year) - 1
  weekly_g <- (1 + wage_growth)^(1 / weeks_per_year) - 1

  ## remaining years of potential earning
  le_at_mean_age <- le_data %>% 
    dplyr::filter(age == round(age_val, 0)) %>% 
    dplyr::select(life_expectancy) %>%
    unlist() %>% as.vector() %>% as.numeric()

  remaining_years_raw <- max(0, min(retirement_age - age_val, le_at_mean_age))

  if (remaining_years_raw <= 0){
    return(0)
  }

  ## convert remaining years to weeks
  total_weeks_raw <- remaining_years_raw * weeks_per_year

  if (years_rounding == "floor") {
    N <- floor(total_weeks_raw)
    frac <- 0
  } else if (years_rounding == "ceiling") {
    N <- ceiling(total_weeks_raw)
    frac <- 0
  } else { 
    N <- floor(total_weeks_raw)
    frac <- total_weeks_raw - N
  }

  ## PV of full weeks as a growing annuity
  if (abs(weekly_r - weekly_g) > 1e-12) {
    ratio <- (1 + weekly_g) / (1 + weekly_r)
    PV_full <- weekly_wage * (1/(1 + weekly_r)) * (1 - ratio^N) / (1 - ratio)
  } else {
    # equal growth & discount → level discounted payments
    PV_full <- N * weekly_wage / (1 + weekly_r)
  }

  ## present value of fractional week costs
  if (frac > 0) {
    payment_next <- weekly_wage * (1 + weekly_g)^N / (1 + weekly_r)^(N + 1)
    PV_frac <- frac * payment_next
  } else {
    PV_frac <- 0
  }

  PV <- PV_full + PV_frac

  ## employment adjustment
  if (adjust_for_employment) {
    PV <- PV * employment_rate
  }

  ## output
  return(as.vector(PV))
}


##******************************************************************************************
## cost-effectiveness analysis (CEA)
##******************************************************************************************
##*********************************************************************
## cea function
##*********************************************************************
cea_func <- function(
  wtp = wtp,
  scenarios = scenarios, 
  ref_scenario = scenarios[1],
  params = params_value, 
  xi_factor = NULL,
  zeta_factor = NULL,
  phi_factor = NULL,
  psi_stated_factor = NULL,
  pi_factor = NULL,
  chi_factor = NULL,
  eta_factor = NULL,
  omicron_factor = NULL,
  init_states = ini_pop_vec,
  inc_names = trans_names,
  time_horizon_years = time_horizon_set[1],
  perspective_val = perspectives[2],
  le_death_data = le_death_data, 
  fx_adjust = fx_adjust, 
  cpi_adjust = cpi_adjust,
  bound_val = c("base", "low", "high")[1],
  par_name = "base",
  base_par_val = NA,
  bound_par_val = NA
){

  ## time horizon
  if(is.null(time_horizon_years)){
    time_horizon_val <- as.numeric(params[["time_horizon_years"]])
  } else {
    time_horizon_val <- time_horizon_years
  }
  n_cycles <- time_horizon_val * 52

  ## fx and cpi adjustment
  sim_params <- params_adj_function(
    pars = params, 
    fx_adjust = fx_adjust, 
    cpi_adjust = cpi_adjust
  )
    
  ## cea analyses
  with(
    as.list(sim_params), {
      ## willingness-to-pay threshold
      Base_WTP <- GDP_Per_Capita
      if(is.null(wtp)){
        wtp <- Base_WTP
      } else {
        wtp <- wtp
      }
      
      ## update age
      age_vals <- Mean_age + (seq(0, n_cycles, 1) - 1) / 52

      ## background mortality
      mu <- unlist(
        lapply(
          X = age_vals,
          FUN = function(x){
            le_death_data %>% 
            dplyr::filter(age == round(x, 0)) %>% 
            mutate(mu = 1 - exp(-((-log(1 - as.numeric(death_prob) / 1000)) / 52))) %>%
            dplyr::select(mu) %>% unlist() %>% as.vector() %>% as.numeric()
          }
        )
      )

      ## analytic perspective
      societal_flag <- perspective_val

      ## population trace
      ## list 
      trace_list <- lapply(
        X = as.vector(scenarios),
        FUN = markov_model_func,
        params = sim_params, 
        mu = mu,
        xi_factor = xi_factor,
        zeta_factor = zeta_factor,
        phi_factor = phi_factor,
        psi_stated_factor = psi_stated_factor,     
        pi_factor = pi_factor,
        chi_factor = chi_factor,
        eta_factor = eta_factor,
        omicron_factor = omicron_factor,
        init_states = init_states,
        inc_names = inc_names,
        n_cycles = n_cycles 
      )
      names(trace_list) <- as.vector(scenarios)

      ## data frame
      trace_df <- dplyr::bind_rows(trace_list) %>%
        dplyr::mutate(
          bound = bound_val,
          parameter = par_name,
          base_param_val = base_par_val,
          bound_param_val = bound_par_val,
          perspective = perspective_val         
        )

      ## utility and cost parameters
      ## utility weights
      qol_E <- qol_S * qol_E_factor
      qol_I <- qol_S * qol_I_factor
      qol_C <- qol_S * qol_C_factor
      qol_R <- qol_S
      qol_brasil_usa_multiplier <- qol_I / qol_I_acute_usa
      qol_loss_severe_vse_IXCHIQ <- prop_severe_vse_IXCHIQ * qol_brasil_usa_multiplier * qol_loss_severe_vae_usa
      qol_loss_serious_vse_IXCHIQ <- prop_serious_vse_IXCHIQ * qol_brasil_usa_multiplier * qol_loss_serious_vae_usa
      qol_loss_vse_IXCHIQ <- qol_loss_severe_vse_IXCHIQ + qol_loss_serious_vse_IXCHIQ
      qol_loss_severe_vse_VIMKUNYA <- prop_severe_vse_VIMKUNYA * qol_brasil_usa_multiplier * qol_loss_severe_vae_usa
      qol_loss_serious_vse_VIMKUNYA <- prop_serious_vse_VIMKUNYA * qol_brasil_usa_multiplier * qol_loss_serious_vae_usa
      qol_loss_vse_VIMKUNYA <- qol_loss_severe_vse_VIMKUNYA + qol_loss_serious_vse_VIMKUNYA

      ## unit costs
      ## direct costs
      prop_acute <- n_reported_cases_acute / sum(c(n_reported_cases_acute, n_reported_cases_postacute, n_reported_cases_chronic))
      prop_postacute <- n_reported_cases_postacute / sum(c(n_reported_cases_acute, n_reported_cases_postacute, n_reported_cases_chronic))
      prop_chronic <- n_reported_cases_chronic / sum(c(n_reported_cases_acute, n_reported_cases_postacute, n_reported_cases_chronic))

      C_hosp_per_stay <- (C_inpatient_stay_private * (1 - pub_share_ip_admissions)) + (C_inpatient_stay_public * pub_share_ip_admissions)
      C_direct_acute <- n_outpatient_visit_acute * C_outpatient_visit + mean(c(C_dipyrone, C_acetaminophen)) * minor_moderate_share_acute + mean(c(C_tramadol, C_codeine, C_oxycodone)) * (1 - minor_moderate_share_acute)
      C_direct_post_acute <- n_outpatient_visit_postacute * C_outpatient_visit + arthritis_share_postacute * C_prednisone + (1 - arthritis_share_postacute) * mean(c(C_ibuprofen, C_amitriptyline, C_gabapentin))
      C_direct_chronic <- n_outpatient_visit_chronic * C_outpatient_visit + C_hydroxychloroquine * mild_illness_share_chronic + sum(c(C_methotrexate, C_folic_acid)) * (1 - mild_illness_share_chronic)  
      
      ## indirect costs (absenteeism + caregiving costs)
      ## net income
      C_INSS_month <- C_Gross_Inc_month * C_INSS_Rate
      C_SimpDed_year <- min(C_Gross_Inc_month * C_Simp_Ded_Rate * 12, C_Simp_Ded_Cap_year)
      C_TaxableInc_year <- C_Gross_Inc_month * 12 - C_INSS_month * 12 - C_SimpDed_year
      TaxRate <- ifelse(
        C_TaxableInc_year <= 26963.20, 
        0, 
        ifelse(
          C_TaxableInc_year <= 33919.80, 
          min(0.08 * C_TaxableInc_year, 2022.24), 
          ifelse(
            C_TaxableInc_year <= 45012.60, 
            min(0.15 * C_TaxableInc_year, 4566.23), 
            ifelse(
              C_TaxableInc_year <= 55976.16, 
              min(0.32 * C_TaxableInc_year, 7942.17), 
              min(0.28 * C_TaxableInc_year, 10740.98)
            )
          )
        )
      )
      C_Parcela_year <- C_TaxableInc_year * TaxRate
      C_NetInc_year <- C_Gross_Inc_month * 12 - C_INSS_month * 12 - C_Parcela_year
      
      ## absenteeism costs
      C_Abs_day <- C_NetInc_year / 365
      Abs_Days_chronic <- Work_Days_Lost - Abs_Days_Inf
      abs_costs_inf <- C_Abs_day * (Abs_Days_Inf * Abs_Freq_Inf)
      abs_costs_chronic <- C_Abs_day * (Abs_Days_chronic * Abs_Freq_Chronic)

      ## caregiving costs
      C_Inf_Caregiving_day <- C_Gross_Inc_month / 30
      caregiving_costs_inf <- C_Inf_Caregiving_day * (Inf_Caregiving_Days_Inf * Inf_Caregiving_Freq_Inf)
      caregiving_costs_chronic <- C_Inf_Caregiving_day * (Inf_Caregiving_Days_Chronic * Inf_Caregiving_Freq_Chronic)

      ## indirect costs (absenteeism + caregiving costs)
      C_indirect_inf <- abs_costs_inf + caregiving_costs_inf
      C_indirect_chronic <- abs_costs_chronic + caregiving_costs_chronic

      ## total costs
      C_S_total <- 0.0
      C_V_NoVac_total <- 0.0

      ## vaccination costs
      ## price discount factor
      if(is.null(eta_factor)){
        vaccine_price_discount <- eta 
      } else{
        vaccine_price_discount <- eta_factor
      }

      ## IXCHIQ
      C_IXCHIQ_Dose <- C_IXCHIQ_Dose_public * (1 - vaccine_price_discount) * pub_share_ip_admissions + 
      C_IXCHIQ_Dose_private * (1 - vaccine_price_discount) * (1 - pub_share_ip_admissions)
      C_V_IXCHIQ_total <- C_IXCHIQ_Dose * (1 + Wastage_Rate) + C_Vaccine_Admin

      ## VIMKUNYA
      if(is.null(chi_factor)) {
        VIMKUNYA_cost_multiplier <-  C_VIMKUNYA_Dose_usa / C_IXCHIQ_Dose_usa
      } else {
        VIMKUNYA_cost_multiplier <- chi_factor
      }

      C_VIMKUNYA_Dose_public <- C_IXCHIQ_Dose_public * (1 - vaccine_price_discount) * VIMKUNYA_cost_multiplier
      C_VIMKUNYA_Dose_private <- C_IXCHIQ_Dose_private * (1 - vaccine_price_discount) * VIMKUNYA_cost_multiplier
      C_VIMKUNYA_Dose <- C_VIMKUNYA_Dose_public * pub_share_ip_admissions + C_VIMKUNYA_Dose_private * (1 - pub_share_ip_admissions)
      C_V_VIMKUNYA_total <- C_VIMKUNYA_Dose * (1 + Wastage_Rate) + C_Vaccine_Admin

      ## total costs
      C_E_total <- 0.0
      C_I_total <- C_hosp_per_stay * prop_hosp + C_direct_acute * (1 - prop_hosp) + C_indirect_inf * societal_flag
      C_R_total <- 0.0
      C_C_total <- (prop_postacute / (prop_postacute + prop_chronic)) * C_direct_post_acute + (prop_chronic/(prop_postacute + prop_chronic)) * C_direct_chronic + C_indirect_chronic * societal_flag
      C_D_total <- C_D
        
      ## cost of premature death: C_D_hc and C_D_total at each cycle (for societal perspective)
      C_D_hc <- unlist(
        lapply(
          X = age_vals, 
          FUN = pv_human_capital_weekly,
          retirement_age = 65,
          monthly_wage = C_Gross_Inc_month,
          le_data = le_death_data,
          discount_rate = disc_rate,  
          wage_growth = 0,
          employment_rate = 1,
          adjust_for_employment = TRUE,
          years_rounding = c("floor", "ceiling", "none")[3]
        )
      )
        
      if(societal_flag == 0){
        C_D_indirect <- 0    
      } else {
        C_D_indirect <- C_D_hc
      }

      xi_val <- ifelse(is.null(xi_factor), xi, xi_factor)
      zeta_val <- ifelse(is.null(zeta_factor), zeta, zeta_factor)
      phi_val <- ifelse(is.null(phi_factor), phi, phi_factor)
      psi_stated_val <- ifelse(is.null(psi_stated_factor), psi_stated, psi_stated_factor)
      pi_val <- ifelse(is.null(pi_factor), pi, pi_factor)
      chi_val <- VIMKUNYA_cost_multiplier
      eta_val <- vaccine_price_discount
      omicron_val <- ifelse(is.null(omicron_factor), omicron, omicron_factor)
      beta_S <-  1 - exp(-((1 / gamma_days) * R0) * cycle_length_days)
      beta_SV <- beta_S * xi_val
      perspective_val <- perspective_val
      time_horizon_val <- n_cycles

      ## unit costs
      unit_costs <- numeric(length(colnames(init_states)))
      names(unit_costs) <- colnames(init_states)

      trans_update <- c("VtoD", "SVtoD", "EtoI", "EtoD", "ItoC", "ItoD", "RtoD", "CtoD")
      trans_update_vals <- c(C_D_total, C_D_total, C_I_total, C_D_total, C_C_total, C_D_total, C_D_total, C_D_total)

      one_time_costs_SOC <- numeric(length(inc_names))
      names(one_time_costs_SOC) <- inc_names  
      one_time_costs_SOC[trans_update] <- trans_update_vals  
      one_time_costs_SOC[c("StoV", "StoSV")] <- c(C_V_NoVac_total, C_V_NoVac_total)    
      
      one_time_costs_IXCHIQ <- numeric(length(inc_names))
      names(one_time_costs_IXCHIQ) <- inc_names  
      one_time_costs_IXCHIQ[trans_update] <- trans_update_vals  
      one_time_costs_IXCHIQ[c("StoV", "StoSV")] <- c(C_V_IXCHIQ_total, C_V_IXCHIQ_total)  
        
      one_time_costs_VIMKUNYA <- numeric(length(inc_names))
      names(one_time_costs_VIMKUNYA) <- inc_names  
      one_time_costs_VIMKUNYA[trans_update] <- trans_update_vals  
      one_time_costs_VIMKUNYA[c("StoV", "StoSV")] <- c(C_V_VIMKUNYA_total, C_V_VIMKUNYA_total)  
        
      costs_list <- list(
        SOC = c(unit_costs, one_time_costs_SOC),
        IXCHIQ = c(unit_costs, one_time_costs_IXCHIQ), 
        VIMKUNYA = c(unit_costs, one_time_costs_VIMKUNYA) 
      )
      
      ## total costs  
      ## list
      total_costs_list <- lapply(
        X = as.vector(scenarios),
        FUN = function(x){
          costs <- costs_list[[x]]
          states <- names(costs)

          state_costs <- matrix(
            data = costs, 
            nrow = length(states), 
            ncol = 1
          )
          rownames(state_costs) <- states
          
          pop <- trace_list[[x]][, states]
          state_pop <- as.matrix(sapply(pop, as.numeric))
          tval <- state_pop %*% state_costs
          
          death_inc_states <- c("StoD", "VtoD", "SVtoD", "EtoD", "ItoD", "RtoD", "CtoD")
          inc_state_pop <- state_pop[, death_inc_states]
          ttval <- matrix(rowSums(inc_state_pop * C_D_indirect), ncol = 1)
          val <- tval + ttval 
          return(val)
        }
      )  
      names(total_costs_list) <- as.vector(scenarios) 

      
      ## data frame
      total_costs_usd_df <- total_costs_list %>%
      dplyr::bind_cols() %>%
      as.data.frame() %>%
      `colnames<-`(as.vector(scenarios)) %>%
      dplyr::mutate(
        cycle = seq_len(nrow(.)) - 1,
        bound = bound_val,
        parameter = par_name,
        base_param_val = base_par_val,
        bound_param_val = bound_par_val,
        xi = xi_val,
        zeta = zeta_val,
        phi = phi_val,
        psi_stated = psi_stated_val,
        pi = pi_val,
        chi = chi_val,
        eta = eta_val,
        omicron = omicron_val,
        beta_S = beta_S,
        beta_SV = beta_SV,
        perspective = perspective_val,
        time_horizon_years = time_horizon_val       
      )
      
      total_costs_brl_df <- total_costs_usd_df %>%
        dplyr::mutate(
          across(
            all_of(scenarios),
            ~ .x / brl_usd_fx_2024
          )
        )      
      
      ## hrqols
      unit_hrqol <- c(qol_S, qol_S, qol_S, qol_E, qol_I, qol_R, qol_C, qol_D)
      names(unit_hrqol) <- colnames(init_states)
      
      one_time_hrqol_SOC <- numeric(length(inc_names))
      names(one_time_hrqol_SOC) <- inc_names  
        
      one_time_hrqol_IXCHIQ <- numeric(length(inc_names))
      names(one_time_hrqol_IXCHIQ) <- inc_names  
      one_time_hrqol_IXCHIQ[c("StoV", "StoSV")] <- -c(qol_loss_vse_IXCHIQ, qol_loss_vse_IXCHIQ)
          
      one_time_hrqol_VIMKUNYA <- numeric(length(inc_names))
      names(one_time_hrqol_VIMKUNYA) <- inc_names  
      one_time_hrqol_VIMKUNYA[c("StoV", "StoSV")] <- -c(qol_loss_vse_VIMKUNYA, qol_loss_vse_VIMKUNYA)
        
      hrqols_list <- list(
        SOC = c(unit_hrqol, one_time_hrqol_SOC),
        IXCHIQ = c(unit_hrqol, one_time_hrqol_IXCHIQ),
        VIMKUNYA = c(unit_hrqol, one_time_hrqol_VIMKUNYA)
      )
          
      ## total qalys
      ## list
      total_qalys_list <- lapply(
        X = as.vector(scenarios),
        FUN = function(x) {
          utils <- hrqols_list[[x]]        
          states <- names(utils)
          state_utils <- matrix(data = as.vector(as.numeric(utils)), nrow = length(states), ncol = 1)
          pop <- trace_list[[x]][, names(utils)]
          state_pop <- as.matrix(sapply(pop, as.numeric))
          val <- state_pop %*% state_utils    
          return(val)
        }
      )  
      names(total_qalys_list) <- as.vector(scenarios)  

      ## data frame
      total_qalys_df <- total_qalys_list %>%
      dplyr::bind_cols() %>%
      as.data.frame() %>%
      `colnames<-`(as.vector(scenarios)) %>%
      dplyr::mutate(
        cycle = seq_len(nrow(.)) - 1,
        bound = bound_val,
        parameter = par_name,
        base_param_val = base_par_val,
        bound_param_val = bound_par_val,
        xi = xi_val,
        zeta = zeta_val,
        phi = phi_val,
        psi_stated = psi_stated_val,
        pi = pi_val,
        chi = chi_val,
        eta = eta_val,
        omicron = omicron_val,
        beta_S = beta_S,
        beta_SV = beta_SV,
        perspective = perspective_val,
        time_horizon_years = time_horizon_val  
      )

      ## apply discounting  
      ## discount factors
      cycle_disc_rate <- (1 + disc_rate)^(1 / 52) - 1
      disc_factors <- 1 / (1 + cycle_disc_rate)^(seq(from = 1, to = (n_cycles + 1), by = 1) -1)

      ## costs
      ## list
      disc_total_costs_list <- lapply(
        X = as.vector(scenarios),
        FUN = function(x){
          val <- total_costs_list[[x]]
          states <- setdiff(colnames(val), c("scenario", "cycle"))
          val[, states] <- disc_factors * val[, states]
          return(val)
        }
      )
      names(disc_total_costs_list) <- as.vector(scenarios)  

      ## data frame
      disc_total_costs_usd_df <- disc_total_costs_list %>%
        dplyr::bind_cols() %>%
        as.data.frame() %>%
        `colnames<-`(as.vector(scenarios)) %>%
        dplyr::mutate(
          cycle = seq_len(nrow(.)) - 1,
          bound = bound_val,
          parameter = par_name,
          base_param_val = base_par_val,
          bound_param_val = bound_par_val,
          xi = xi_val,
          zeta = zeta_val,
          phi = phi_val,
          psi_stated = psi_stated_val,
          pi = pi_val,
          chi = chi_val,
          eta = eta_val,
          omicron = omicron_val,
          beta_S = beta_S,
          beta_SV = beta_SV,
          perspective = perspective_val,
          time_horizon_years = time_horizon_val         
        )

      disc_total_costs_brl_df <- disc_total_costs_usd_df %>%
        mutate(
          across(
            all_of(as.vector(scenarios)),
            ~ .x / brl_usd_fx_2024
          )
        )     
      
      ## qalys
      ## list
      disc_total_qalys_list <- lapply(
        X = as.vector(scenarios),
        FUN = function(x){
          val <- disc_factors * total_qalys_list[[x]]
        }
      )
      names(disc_total_qalys_list) <- as.vector(scenarios)  

      ## data frame
      disc_total_qalys_df <- disc_total_qalys_list %>%
      dplyr::bind_cols() %>%
      as.data.frame() %>%
      `colnames<-`(as.vector(scenarios)) %>%
      dplyr::mutate(
        cycle = seq_len(nrow(.)) - 1,
        bound = bound_val,
        parameter = par_name,
        base_param_val = base_par_val,
        bound_param_val = bound_par_val,
        xi = xi_val,
        zeta = zeta_val,
        phi = phi_val,
        psi_stated = psi_stated_val,
        pi = pi_val,
        chi = chi_val,
        eta = eta_val,
        omicron = omicron_val,
        beta_S = beta_S,
        beta_SV = beta_SV,
        perspective = perspective_val,
        time_horizon_years = time_horizon_val      
      )
    
      ## cea  
      ## discounted costs
      d_costs <- do.call(
        rbind,
        lapply(
          X = as.vector(scenarios),
          FUN = function(x){
            sum(disc_total_costs_list[[x]])
          }
        )
      )
      colnames(d_costs) <- "costs"

      ## discounted qalys
      d_qalys <- do.call(
        rbind,
        lapply(
          X = as.vector(scenarios),
          FUN = function(x){
            sum(disc_total_qalys_list[[x]])
          }
        )
      )
      colnames(d_qalys) <- "qalys" 
      
      ## cea output
      ## brazilian real
      cea_out_brl_df <- cbind(d_costs, d_qalys) %>%
        as.data.frame() %>%
        mutate(
          alternative = as.vector(scenarios),
          wtp = wtp,
          nmb = qalys * wtp - costs,
          nhb = qalys - costs / wtp,
          inc_costs = costs - costs[which(alternative == ref_scenario)],
          inc_qalys = qalys - qalys[which(alternative == ref_scenario)],
          icer = inc_costs / inc_qalys,        
          inmb = inc_qalys * wtp - inc_costs,
          inhb = inc_qalys - inc_costs / wtp,
          bound = bound_val,
          parameter = par_name,
          base_param_val = base_par_val,
          bound_param_val = bound_par_val,
          xi = xi_val,
          zeta = zeta_val,
          phi = phi_val,
          psi_stated = psi_stated_val,
          pi = pi_val,
          chi = chi_val,
          eta = eta_val,
          omicron = omicron_val,
          beta_S = beta_S,
          beta_SV = beta_SV,
          perspective = perspective_val,
          time_horizon_years = time_horizon_val    
        ) %>%
        dplyr::select(
          alternative, wtp, costs, qalys, nmb, nhb, 
          inc_costs, inc_qalys, icer, inmb, inhb,
          bound, parameter, base_param_val, bound_param_val,
          xi, zeta, pi, zeta, psi_stated, pi, chi, eta, omicron,
          beta_S, beta_SV,
          perspective, time_horizon_years
        )
      
      ## usd
      cea_out_usd_df <- cea_out_brl_df %>%
        mutate(
          across(
            c(wtp, costs, nmb, inc_costs, icer, inmb),
            ~ .x * (1 / brl_usd_fx_2024)
          )
        ) 
      
      ## output
      return(
        list(
          ## trace
          trace = trace_df,
          
          ## undiscounted costs an qalys
          costs_brl = total_costs_usd_df,
          costs_usd = total_costs_brl_df,
          qalys = total_qalys_df,

          ## discounted costs and qalys
          disc_costs_brl = disc_total_costs_usd_df,
          disc_costs_usd = disc_total_costs_brl_df,
          disc_qalys = disc_total_qalys_df,

          ## cea
          cea_out_brl = cea_out_brl_df,
          cea_out_usd = cea_out_usd_df 
        )
      )
    }
  )
}

##*********************************************************************
## Base case CEA function for all perspectives
##*********************************************************************
base_case_cea_func <- function(
  wtp = NULL,
  scenarios = scenarios, 
  ref_scenario = scenarios[1],
  params = params_value,
  xi_factor = NULL,
  zeta_factor = NULL,
  phi_factor = NULL,
  psi_stated_factor = NULL,
  pi_factor = NULL,
  chi_factor = NULL,
  eta_factor = NULL,
  omicron_factor = NULL,
  init_states = ini_pop_vec,
  inc_names = trans_names,
  time_horizon_years = time_horizon_set[1],
  perspectives = c(0, 1),
  le_death_data = le_death_data, 
  fx_adjust = fx_adjust, 
  cpi_adjust = cpi_adjust,
  bound_val = c("base", "low", "high")[1],
  par_name = "base",
  base_par_val = NA,
  bound_par_val = NA
){
  ## run by perspective
  outer_results <- lapply(
    X = seq_along(perspectives),
    FUN = function(i){
      cea_func(
        wtp = wtp,
        scenarios = scenarios, 
        ref_scenario = ref_scenario,
        params = params, 
        xi_factor = xi_factor,
        zeta_factor = zeta_factor,
        phi_factor = phi_factor,
        psi_stated_factor = psi_stated_factor,      
        pi_factor = pi_factor,
        chi_factor = chi_factor,
        eta_factor = eta_factor,
        omicron_factor = omicron_factor,
        init_states = init_states,
        inc_names = inc_names,
        time_horizon_years = time_horizon_years,
        perspective_val = perspectives[i],
        le_death_data = le_death_data, 
        fx_adjust = fx_adjust, 
        cpi_adjust = cpi_adjust,
        bound_val = bound_val,
        par_name = par_name,
        base_par_val = base_par_val,
        bound_par_val = bound_par_val      
      )
    }
  )

  ## combine outputs by element name (safe, no duplication)
  out <- purrr::map(
    names(outer_results[[1]]),
    function(nm) {
      dplyr::bind_rows(
        purrr::map(outer_results, ~ .x[[nm]])
      )
    }
  )
  names(out) <- names(outer_results[[1]])

  ## output
  return(out)
}

##******************************************************************************
## Function for summary results table of base case analyses
##******************************************************************************
cea_table_func <- function(
  df,
  perspective_labs = c(`0` = "Healthcare", `1` = "Societal"),
  journal = "bmj",
  digits = digits_vals,
  scale = scale_vals,
  linebreak = "\n",
  big.mark = big.mark_vals,
  decimal.mark = decimal.mark_vals,
  ci_sep = NULL,
  perspective_var_lab = "Perspective",
  reference_lab = "Comparator",
  alternative_lab = "Alternative",
  Costs_lab = "Costs (US$)",
  QALYs_lab = "QALYs",
  NMB_lab = "NMB (US$ million)",
  Inc_Costs_lab = "Costs (US$)",
  Inc_QALYs_lab = "QALYs",
  INMB_lab = "INMB (US$ million)",
  ICER_lab = "ICER (US$/QALY)",
  dominance_status_lab = "Dominance Status",
  tab_caption = "Table 2. Base Case Cost-effectiveness Analysis Results",
  tab_dir = tables_dir,
  fname = "basecase_cea_results_table.docx"
) {

  ## extended dominance
  label_dominance <- function(df) {
    df %>%
      arrange(costs, qalys) %>%
      ## weak (simple) dominance
      mutate(
        weakly_dominated =
          (costs >= lag(costs) & qalys <= lag(qalys)) &
          (costs > lag(costs) | qalys < lag(qalys))
      ) %>%      
      ## ICERs only make sense for non–weakly dominated strategies
      mutate(
        icer = if_else(
          is.na(weakly_dominated) | !weakly_dominated,
          (costs - lag(costs)) / (qalys - lag(qalys)),
          NA_real_
        )
      ) %>%     
      ## extended dominance (only among non–weakly dominated)
      mutate(
        extended_dominated =
          !weakly_dominated &
          !is.na(lag(icer)) &
          lag(icer) < icer
      ) %>%     
      ## final dominance label
      mutate(
        dominance_status = case_when(
          weakly_dominated ~ "weakly_dominated",
          extended_dominated ~ "extended_dominated",
          TRUE ~ "efficient"
        )
      )
  }
  
  ## apply dominance rules
  combo_df <- df %>%
    group_by(perspective) %>%
    label_dominance() %>%
    mutate(strategy_order = row_number()) %>%
    ungroup()

  ## pairwise scenario comparisons
  incremental <- tidyr::expand_grid(
    reference = unique(combo_df$strategy_order),
    target = unique(combo_df$strategy_order),
    perspective = unique(combo_df$perspective),
    alternative = unique(combo_df$alternative)
  ) %>%
    dplyr::left_join(
      combo_df %>% dplyr::rename_with(~ paste0(.x, "_t")),
      by = c("target" = "strategy_order_t", "perspective" = "perspective_t")
    ) %>%
    dplyr::left_join(
      combo_df %>% dplyr::rename_with(~ paste0(.x, "_r")),
      by = c("reference" = "strategy_order_r", "perspective" = "perspective_r")
    ) %>%
    dplyr::mutate(
      inc_costs = costs_t - costs_r,
      inc_qalys = qalys_t - qalys_r,
      inmb = nmb_t - nmb_r,
      inhb = nhb_t - nhb_r,
      icer = dplyr::if_else(
        inc_qalys == 0 & inc_costs == 0,
        NA_real_,
        inc_costs / inc_qalys
      )
    ) %>%
    dplyr::select(
      alternative_t, alternative_r,
      perspective,
      costs_t, qalys_t, nmb_t,
      inc_costs, inc_qalys, inmb, icer,
      dominance_status_t
    ) %>%
    dplyr::rename(
      alternative = alternative_t,
      reference = alternative_r,
      costs = costs_t,
      qalys = qalys_t,
      nmb = nmb_t,
      dominance_status = dominance_status_t
    ) %>%
    dplyr::group_by(perspective, alternative) %>%
    dplyr::filter(
      (alternative == reference & dplyr::row_number() == 1) |
      (alternative != reference & inmb >= 0)
    ) %>%
    unique()

  ## formatting
  df_fmt <- incremental %>%
    mutate(
      perspective = factor(
        perspective,
        levels = c(0, 1),
        labels = perspective_labs
      ),
      Costs = fmt_ci_journal(
        mean = costs, 
        lower = NULL, 
        upper = NULL,
        journal = journal,
        digits = digits$Costs,
        scale = scale$Costs,
        linebreak = linebreak,
        big.mark = big.mark$Costs,
        decimal.mark = decimal.mark$Costs,
        ci_sep = ci_sep
      ),
      QALYs = fmt_ci_journal(
        mean = qalys, 
        lower = NULL, 
        upper = NULL,
        journal = journal,
        digits = digits$QALYs,
        scale = scale$QALYs,
        linebreak = linebreak,
        big.mark = big.mark$QALYs,
        decimal.mark = decimal.mark$QALYs,
        ci_sep = ci_sep
      ),
      NMB = fmt_ci_journal(
        mean = nmb, 
        lower = NULL, 
        upper = NULL,
        journal = journal,
        digits = digits$NMB,
        scale = scale$NMB,
        linebreak = linebreak,
        big.mark = big.mark$NMB,
        decimal.mark = decimal.mark$NMB,
        ci_sep = ci_sep
      ),
      Inc_Costs = fmt_ci_journal(
        mean = inc_costs, 
        lower = NULL, 
        upper = NULL,
        journal = journal,
        digits = digits$Inc_Costs,
        scale = scale$Inc_Costs,
        linebreak = linebreak,
        big.mark = big.mark$Inc_Costs,
        decimal.mark = decimal.mark$Inc_Costs,
        ci_sep = ci_sep
      ),
      Inc_QALYs = fmt_ci_journal(
        mean = inc_qalys, 
        lower = NULL, 
        upper = NULL,
        journal = journal,
        digits = digits$Inc_QALYs,
        scale = scale$Inc_QALYs,
        linebreak = linebreak,
        big.mark = big.mark$Inc_QALYs,
        decimal.mark = decimal.mark$Inc_QALYs,
        ci_sep = ci_sep
      ),
      INMB = fmt_ci_journal(
        mean = inmb, 
        lower = NULL, 
        upper = NULL,
        journal = journal,
        digits = digits$INMB,
        scale = scale$INMB,
        linebreak = linebreak,
        big.mark = big.mark$INMB,
        decimal.mark = decimal.mark$INMB,
        ci_sep = ci_sep
      ),
      ICER = case_when(
        is.na(icer) ~ "··",
        icer < 0 ~ "Dominant",
        TRUE ~ fmt_ci_journal(
          mean = icer, 
          lower = NULL, 
          upper = NULL,
          journal = journal,
          digits = digits$ICER,
          scale = scale$ICER,
          linebreak = linebreak,
          big.mark = big.mark$ICER,
          decimal.mark = decimal.mark$ICER,
          ci_sep = ci_sep
        )
      )
    ) %>%
    dplyr::select(
      perspective,
      reference,
      alternative,
      Costs,
      QALYs,
      NMB,
      Inc_Costs,
      Inc_QALYs,
      INMB,
      ICER,
      dominance_status
    )
    
  ## header rows
  df_with_headers <- df_fmt %>%
    dplyr::group_by(perspective, reference) %>%
    dplyr::reframe(
      alternative = c(paste0("Relative to ", unique(reference)), alternative),
      dplyr::across(-alternative, ~ c("", .x))
    ) %>%
    dplyr::ungroup()

  ## flextable
  ft <- flextable::flextable(df_with_headers) %>%
    flextable::theme_vanilla() %>%
    flextable::merge_v(j = 1) %>%
    flextable::bold(i = ~ grepl("Relative to", alternative), j = 3) %>%
    flextable::italic(i = ~ grepl("Relative to", alternative), j = 3) %>%
    flextable::set_header_labels(
      perspective = perspective_var_lab,
      reference = reference_lab,
      alternative = alternative_lab,
      Costs = Costs_lab,
      QALYs = QALYs_lab,
      NMB = NMB_lab ,
      Inc_Costs = Inc_Costs_lab,
      Inc_QALYs = Inc_QALYs_lab,
      INMB = INMB_lab,
      ICER = ICER_lab,
      dominance_status = dominance_status_lab 
    ) %>%
    flextable::add_header_row(
      values = c(
        "", "", "", "Total Effects", "", "",
        "Incremental Effects", "", "", "", ""
      ),
      colwidths = rep(1, 11)
    ) %>%
    flextable::autofit()
  
  ## save to word
  doc <- read_docx() %>%
    body_add_par(
      x = .,
      value = tab_caption, 
      style = "heading 1", 
      pos = "after"
    ) %>%
    body_add_flextable(
      ft
    ) %>%
    body_add_par(
      x = .,
      value = "", 
      style = "Normal", 
      pos = "after"
    )
  
  ## save word file
  print(
    x = doc, 
    target = file.path(tab_dir, fname)
  )
  
  return(ft)
}


##********************************************************************************************
## Functions for one-way sensitivity analysis (OWSA)
##********************************************************************************************

##*********************************************************************
## extract base outcome (ICER or INMB or INHB)
##*********************************************************************
extract_outcome <- function(
  X = 1,
  currencies_outcomes,
  result 
) {
  currency <- currencies_outcomes$currency[X]
  outcome <- currencies_outcomes$outcome[X]
  curr_result <- result[[paste0("cea_out_", currency)]]
  if (outcome == "icer") {
    val <- tryCatch({
      curr_result %>% dplyr::select(alternative, icer)
    }, error = function(e) NA)
  } else {
    if (outcome == "inmb") {
      val <- tryCatch({
        curr_result %>% dplyr::select(alternative, inmb)
      }, error = function(e) NA)
    } else {
      val <- tryCatch({
        curr_result %>% dplyr::select(alternative, inhb)
      }, error = function(e) NA)
    }
  }
  return(val)
}

##*********************************************************************
## inner OWSA function: iterates overs the lower and upper bound values
## of each parameter
##*********************************************************************
inner_owsa_func <- function(
  param = names(params_value)[1],
  base_params = params_value, 
  low_params = params_lower, 
  up_params = params_upper,
  xi_factor = NULL,
  zeta_factor = NULL,
  phi_factor = NULL,
  psi_stated_factor = NULL,
  pi_factor = NULL,
  chi_factor = NULL,
  eta_factor = NULL,
  omicron_factor = NULL,
  base_data = base_result,
  base_val_outcome, 
  wtp = NULL,
  scenarios = scenarios, 
  ref_scenario = scenarios[1],
  init_states = ini_pop_vec,
  inc_names = trans_names,
  time_horizon_years = time_horizon_set[1],
  perspectives = c(0, 1),
  le_death_data = le_death_data, 
  fx_adjust = fx_adjust, 
  cpi_adjust = cpi_adjust
) {

  ## initialize parameters
  base_val_param <- base_params[[param]]
  low_val_param  <- low_params[[param]]
  high_val_param <- up_params[[param]]
  
  ## safety check
  if (!is.numeric(base_val_param) || is.na(low_val_param) || is.na(high_val_param)){
    return(NULL)
  }
  if (low_val_param == high_val_param){
    return(NULL)
  }

  ## base value
  results_base <- base_data 

  ## lower bound
  results_low <- base_case_cea_func(
    wtp = wtp,
    scenarios = scenarios, 
    ref_scenario = ref_scenario,
    params = replace(base_params, param, low_val_param),
    xi_factor = xi_factor,
    zeta_factor = zeta_factor,
    phi_factor = phi_factor,
    psi_stated_factor = psi_stated_factor,
    pi_factor = pi_factor,
    chi_factor = chi_factor,
    eta_factor = eta_factor,
    omicron_factor = omicron_factor,
    init_states = init_states,
    inc_names = inc_names,
    time_horizon_years = time_horizon_years,
    perspectives = perspectives,
    le_death_data = le_death_data, 
    fx_adjust = fx_adjust, 
    cpi_adjust = cpi_adjust,
    bound_val = "low",
    par_name = param,
    base_par_val = as.vector(base_val_param),
    bound_par_val = as.vector(low_val_param)
  )
  
  ## upper bound
  results_high <- base_case_cea_func(
    wtp = wtp,
    scenarios = scenarios, 
    ref_scenario = ref_scenario,
    params = replace(base_params, param, high_val_param),
    xi_factor = xi_factor,
    zeta_factor = zeta_factor,
    phi_factor = phi_factor,
    psi_stated_factor = psi_stated_factor,
    pi_factor = pi_factor,
    chi_factor = chi_factor,
    eta_factor = eta_factor,
    omicron_factor = omicron_factor,
    init_states = init_states,
    inc_names = inc_names,
    time_horizon_years = time_horizon_years,
    perspectives = perspectives,
    le_death_data = le_death_data, 
    fx_adjust = fx_adjust, 
    cpi_adjust = cpi_adjust,
    bound_val = "high",
    par_name = param,
    base_par_val = as.vector(base_val_param),
    bound_par_val = as.vector(high_val_param)
  )
  
  ## append dataframes from all lists
  val <- Map(
    function(b, l, h) {
      out <- rbind(b, l, h)
      rownames(out) <- NULL
      out
    },
    results_base,
    results_low,
    results_high
  )
  return(val)
}


##*********************************************************************
## function for one-way sensitivity analysis
##*********************************************************************
run_owsa <- function(
  base_params = params_value, 
  low_params = params_lower, 
  up_params = params_upper,
  xi_factor = NULL,
  zeta_factor = NULL,
  phi_factor = NULL,
  psi_stated_factor = NULL,
  pi_factor = NULL,
  chi_factor = NULL,
  eta_factor = NULL,
  omicron_factor = NULL,
  currencies = c("brl", "usd"),
  outcomes = c("icer", "inmb", "inhb"),
  wtp = NULL,
  scenarios = scenarios, 
  ref_scenario = scenarios[1],
  init_states = ini_pop_vec,
  inc_names = trans_names,
  time_horizon_years = time_horizon_set[1],
  perspectives = c(0, 1),
  le_death_data = le_death_data, 
  fx_adjust = fx_adjust, 
  cpi_adjust = cpi_adjust
) {

  ##--------------------------------------------------
  ## Parameters present in all bounds
  ##--------------------------------------------------
  parameter_names <- Reduce(
    union,
    list(names(base_params), names(low_params), names(up_params))
  )

  ##--------------------------------------------------
  ## Base case (run once)
  ##--------------------------------------------------
  base_result <- base_case_cea_func(
    wtp = wtp,
    scenarios = scenarios, 
    ref_scenario = ref_scenario,
    params = base_params, 
    xi_factor = xi_factor,
    zeta_factor = zeta_factor,
    phi_factor = phi_factor,
    psi_stated_factor = psi_stated_factor,  
    pi_factor = pi_factor,
    chi_factor = chi_factor,
    eta_factor = eta_factor,
    omicron_factor = omicron_factor,
    init_states = init_states,
    inc_names = inc_names,
    time_horizon_years = time_horizon_years,
    perspectives = perspectives,
    le_death_data = le_death_data, 
    fx_adjust = fx_adjust, 
    cpi_adjust = cpi_adjust,
    bound_val = "base",
    par_name = "base",
    base_par_val = NA,
    bound_par_val = NA
  )

  ##--------------------------------------------------
  ## Run OWSA per parameter
  ##--------------------------------------------------  
  val <- list(
    trace = NULL,
    costs_brl = NULL,
    costs_usd = NULL,
    qalys = NULL,
    disc_costs_brl = NULL,
    disc_costs_usd = NULL,
    disc_qalys = NULL,
    cea_out_brl = NULL,
    cea_out_usd = NULL
  )
  
  for(p in parameter_names){

    df <- inner_owsa_func(
      param = p,
      base_params = base_params, 
      low_params = low_params, 
      up_params = up_params,
      xi_factor = xi_factor,
      zeta_factor = zeta_factor,
      phi_factor = phi_factor,
      psi_stated_factor = psi_stated_factor,
      pi_factor = pi_factor,
      chi_factor = chi_factor,
      eta_factor = eta_factor,
      omicron_factor = omicron_factor,
      base_data = base_result,
      wtp = wtp,
      scenarios = scenarios, 
      ref_scenario = ref_scenario,
      init_states = init_states,
      inc_names = inc_names,
      time_horizon_years = time_horizon_years,
      perspectives = perspectives,
      le_death_data = le_death_data, 
      fx_adjust = fx_adjust, 
      cpi_adjust = cpi_adjust
    )

    #print(head(df$trace))
    val$trace <- rbind(val$trace, df$trace)
    val$costs_brl <- rbind(val$costs_brl, df$costs_brl)
    val$costs_usd <- rbind(val$costs_usd, df$costs_usd)
    val$qalys <- rbind(val$qalys, df$qalys)
    val$disc_costs_brl <- rbind(val$disc_costs_brl, df$disc_costs_brl)
    val$disc_costs_usd <- rbind(val$disc_costs_usd, df$disc_costs_usd)
    val$disc_qalys <- rbind(val$disc_qalys, df$disc_qalys)
    val$cea_out_brl <- rbind(val$cea_out_brl, df$cea_out_brl)
    val$cea_out_usd <- rbind(val$cea_out_usd, df$cea_out_usd)
  }

  ##--------------------------------------------------
  ## Output
  ##--------------------------------------------------
  return(val)
}

##******************************************
## Function to create data for the OWSA tornado diagram
##******************************************
owsa_tornado_plot_data_func <- function(
  alt_val = scenarios[2],
  perspective_val = c(0,1)[1],
  curr_val = c("brl", "usd")[2],
  out_val = c("icer", "inmb", "inhb")[1],
  owsa_output_data = owsa_output_data,
  params_include = owsa_params,
  params_label = owsa_params_labels,
  diff_tolerance = 0.001,
  bar_width = 0.75
){

  ## safety checks
  if (is.null(owsa_output_data)) {
    stop("owsa_output cannot be empty")
  }

  ##******************************************
  ## Data for the Tornado diagram
  ##******************************************
  ## base case data
  base_data <- owsa_output_data[[paste0("cea_out_", curr_val)]] %>%
    dplyr::filter(
      perspective == perspective_val,
      alternative == alt_val,
      bound == "base"
    ) %>%
    dplyr::select(
      all_of(
        c(
          "alternative",
          "perspective",
          out_val
        )
      )
    ) %>%
    dplyr::rename(
      base = all_of(out_val)
    ) %>%
    unique()  

  ## owsa data
  owsa_data <- owsa_output_data[[paste0("cea_out_", curr_val)]] %>%
    dplyr::filter(
      perspective == perspective_val,
      alternative == alt_val,
      bound %in% c("base", "low", "high"),
      parameter %in% params_include
    ) %>%
    dplyr::select(
      all_of(
        c(
          "alternative",
          "perspective",
          "wtp",
          out_val,
          "parameter",
          "base_param_val",
          "bound"
        )
      )
    )

  ## calculate the difference between upper and lower bounds
  owsa_df <- owsa_data %>%
    pivot_wider(
      names_from  = bound,
      values_from = all_of(out_val)
    ) %>%
    left_join(base_data, by = c("alternative", "perspective")) %>%
    mutate(
      diff_low = low  - base,
      diff_high = high - base,
      range = abs(high - low)
    ) %>%
    dplyr::filter(range > 0)
      
  ## base case value of the outcome
  base_value <- owsa_df %>%
    dplyr::select(base) %>%
    unlist() %>% as.vector() %>% unique()

  ## max difference
  max_diff <- max(owsa_df$range, na.rm = TRUE)

  ## check
  if (!is.finite(max_diff)) {
    warning(
      "No meaningful OWSA variation for this outcome; tornado plot omitted."
    )
    return(NULL)
  } else {
    owsa_df <- owsa_df %>% 
      dplyr::filter(range >= (0 * max_diff))

    if (nrow(owsa_df) == 0) {
      stop("No parameters pass the inclusion / threshold filter.")
    }  
  }

  ## get order of parameters according to size of intervals
  params_order_vals <- owsa_df %>% 
    dplyr::arrange(range) %>%
    mutate(parameter = factor(x = parameter, levels = parameter))
      
  ## extract ordered parameter levels
  params_order <- levels(params_order_vals$parameter)

  ## labels mapping (parameter → label)
  pars_label <- tibble::tibble(
    parameter = params_order,
    labels = params_label[params_order]
  )
 
  ## get data frame in shape for ggplot and geom_rect
  plot_data <- owsa_df %>% 
    pivot_longer(
      cols = c(low, high),
      names_to = "bound",
      values_to = "value"
    ) %>%
    mutate(
      parameter = factor(parameter, levels = params_order),
      ymin = pmin(value, base),
      ymax = pmax(value, base),
      xmin = as.numeric(parameter) - bar_width/2,
      xmax = as.numeric(parameter) + bar_width/2,
      bound = factor(
        bound, 
        levels = c("low", "high"), 
        labels = c("Lower bound", "Upper bound")
      ),
      annotate_label = dplyr::case_when(
        out_val == "icer" ~ paste0("ICER = US$", formatC(base_value, format = "f", big.mark = ",", digits = 0), " per DALY averted"),
        out_val == "inmb" ~ paste0("INMB = US$", formatC(base_value, format = "f", big.mark = ",", digits = 0)),
        out_val == "inhb" ~ paste0("INHB = ", formatC(base_value, format = "f", big.mark = ",", digits = 0), " QALYs gained"),
        TRUE ~ NA_character_
      )
    ) %>%
    left_join(pars_label, by = "parameter")
  
  ## output
  return(plot_data)
  
}

##******************************************
## Tornado plot function
##******************************************
owsa_tornado_plot_func <- function(
  alternatives = scenarios[2:3],
  perspective_val = c(0,1)[1],
  curr_val = c("brl", "usd")[2],
  out_val = c("icer", "inmb", "inhb")[1],
  ##
  owsa_output_data = owsa_output_data,
  params_include = owsa_params,
  params_label = owsa_params_labels,
  diff_tolerance = 0.001,
  bar_width = 0.75,
  ##
  y_axis_rescale_by = 1000,
  y_min = 10000,
  y_max = 50000,
  y_by = 10000,
  y_exp_min = 0.01,
  y_exp_max = 0.01,
  x_exp_min = 0.01,
  x_exp_max = 0.01,
  ##
  xlab = "ICER, 2024 US$ per DALY averted",
  ylab = "Model input parameters",
  title = NULL,
  strip_text = facet_strip_text, 
  strip_txt_size = 12,
  ##
  annotate_txt_size = 4,
  x_offset_val = 0.2,
  hjust_val_min = 1,
  hjust_val_max = 0,
  ##
  x_axis_txt_size = 12,
  x_axis_title_size = 12,
  y_axis_txt_size = 12,
  y_axis_title_size = 12,
  legend_txt_size = 12, 
  ##
  file_name = "tornado_owsa_IXCHIQ_societal.png",
  file_path = figures_dir,
  fig_scale = 1,
  fig_width = 8, 
  fig_height = 6,
  units = c("in", "cm", "mm", "px")[1],
  dpi = 300,
  limitsize = TRUE, 
  bg = "white",
  ##
  journal = c("lancet", "nejm", "jama", "bmj", "nature", "default")[4]
){
  
  ## Combined plots data
  ini_plot_data <- do.call(
    rbind,
    lapply(
      X = alternatives,
      FUN = owsa_tornado_plot_data_func,
      perspective_val = perspective_val,
      curr_val = curr_val,
      out_val = out_val,
      owsa_output_data = owsa_output_data,
      params_include = params_include,
      params_label = params_label,
      diff_tolerance = diff_tolerance,
      bar_width = bar_width
    )
  )

  ## max difference
  max_diff <- max(ini_plot_data$range, na.rm = TRUE)

  temp_plot_data <- ini_plot_data %>%
    dplyr::filter(range > (diff_tolerance * max_diff))
  
  parameter_names <- unique(temp_plot_data$parameter)
  parameter_labels <- unique(temp_plot_data$labels)

  t_plot_data <- ini_plot_data %>%
    dplyr::filter(parameter %in% parameter_names)

  ## Fix factor ordering across all alternatives
  ordered_params <- t_plot_data %>%
    dplyr::arrange(desc(-range), parameter, labels) %>%  
    dplyr::distinct(parameter, labels) %>%
    mutate(parameter_index = row_number())

  plot_data <- t_plot_data %>%
    left_join(ordered_params, by = c("parameter", "labels")) %>%
    mutate(
      xmin = parameter_index - bar_width / 2,
      xmax = parameter_index + bar_width / 2
    )
  
  ## annotations of the base-case vertical lines
  annot_df <- expand.grid(
    facet_alternative = unique(plot_data$alternative), 
    line_alternative = alternatives,
    stringsAsFactors = FALSE
  ) %>%
  ## join base outcomes and labels for each line alternative
  left_join(
    plot_data %>%
    dplyr::distinct(alternative, base, annotate_label),
    by = c("line_alternative" = "alternative")
  ) %>%
  dplyr::mutate(
    y_val = base / y_axis_rescale_by,
    line_color = ifelse(
      facet_alternative == line_alternative, 
      "black", 
      "black"
    ),
    ## label positions: offset left/right depending on main or secondary line
    x_offset = ifelse(
      y_val == min(y_val),
      min(plot_data$parameter_index) + x_offset_val,
      min(plot_data$parameter_index) + x_offset_val
    ),
    hjust_val = ifelse(
      y_val == min(y_val), 
      hjust_val_min, 
      hjust_val_max
    ),
    label_text = annotate_label
  )

  plot_data %>%
    dplyr::filter(parameter == "lambda_days") %>%
    head() %>% print()

  ## plot
  p_owsa <- ggplot(data = plot_data) + 
  facet_wrap(
    .~ alternative, 
    scales = "free_x",
    labeller = labeller(alternative = strip_text)
  ) +
  geom_rect(
    aes(
      ymax = ymax / y_axis_rescale_by, 
      ymin = ymin / y_axis_rescale_by, 
      xmax = xmax, 
      xmin = xmin, 
      fill = bound
    )
  ) +
  geom_hline(
    data = annot_df,
    aes(yintercept = y_val, color = line_color, group = facet_alternative),
    linetype = "dashed",
    show.legend = FALSE
  ) +
  geom_text(
    data = annot_df,
    aes(x = x_offset, y = y_val, label = label_text, hjust = hjust_val, color = line_color),
    inherit.aes = FALSE,
    size = annotate_txt_size
  ) +
  scale_color_identity() +   
  labs(
    x = xlab, 
    y = ylab
  ) + 
  theme_minimal() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line.x = element_line(color = "black", size = 0.5),
    axis.line.y = element_line(color = "black", size = 0.5),
    axis.ticks.y = element_line(color = "black"),
    axis.ticks.x = element_line(color = "black"),
    axis.text.x = element_text(size = x_axis_txt_size, hjust = 0.5),
    axis.text.y = element_text(size = y_axis_txt_size),
    axis.title.x = element_text(size = x_axis_title_size, face = "bold"),
    axis.title.y = element_text(size = y_axis_title_size, face = "bold"),
    legend.title = element_blank(),
    legend.text = element_text(size = legend_txt_size), 
    legend.position = "bottom", 
    legend.background = element_rect(fill = "white", color = NA), 
    legend.key.size = unit(0.5, "cm"),
    legend.key.height = unit(0.5, 'cm'),
    legend.key.width = unit(0.5, 'cm'), 
    plot.margin = unit(c(0.3, 0.3, 0.2, 0.2),"cm"),
    legend.margin = margin(t = 1, r = 1, b = 1, l = 1), 
    strip.text = element_text(size = strip_txt_size, face = "bold"),
    strip.background = element_blank()  
  ) + 
  scale_x_continuous(
    breaks = ordered_params$parameter_index,
    labels = ordered_params$labels,
    expand = c(x_exp_min, x_exp_max) 
  ) +
  scale_y_continuous(
    labels = scales::comma, 
    limits = c(y_min, y_max) / y_axis_rescale_by, 
    breaks = seq(y_min, y_max, y_by) / y_axis_rescale_by, 
    expand = c(y_exp_min, y_exp_max)  
  ) +
  get_journal_fill_scale(journal) +
  coord_flip()

  ggsave(
    filename = file_name, 
    plot = p_owsa, 
    path = file_path,
    scale = fig_scale, 
    width = fig_width, 
    height = fig_height, 
    units = units, 
    dpi = dpi, 
    limitsize = limitsize, 
    bg = bg
  )

  ## output
  return(p_owsa) 
}

##******************************************
## Function to draw from appropriate distribution
##******************************************
draw_params <- function(
  dist, 
  mean, 
  lower, 
  upper, 
  n_sim = 1000
){
 
  dist_clean <- trimws(tolower(gsub("-", "", dist)))
  dist_clean[dist_clean == "" | is.na(dist_clean)] <- "na"
  
  ## beta distribution
  if (dist_clean == "beta") {
    ## transform bounds to 95% CI quantiles of beta;
    ## use method of moments if approximate (bounded [0,1])
    if (mean == 0 || mean == 1) {
      return(rep(mean, n_sim))
    } else{
      epsilon <- 1e-6
      lower <- max(lower, epsilon)
      upper <- min(upper, 1 - epsilon)
      if (mean < 0 || mean > 1 || lower < 0 || upper > 1) {
        stop("Beta distribution requires mean, lower, and upper in (0,1)")
      } else {
        sd <- (upper - lower) / (2 * 1.96)
        alpha <- ((1 - mean) / sd^2 - 1 / mean) * mean^2
        beta <- alpha * (1 / mean - 1)
        return(rbeta(n_sim, shape1 = alpha, shape2 = beta))
      }  
    }
  
  ## Chi-square distribution
  } else if (dist_clean == "chisq" || dist_clean == "chisquare") {
    if (mean <= 0) stop("Chi-square distribution requires mean > 0 (mean = df)")
    df <- mean
    return(rchisq(n_sim, df = df))
    
  } else {
    ## lognormal distribution
    if (dist_clean == "lognormal") {
      if (mean <= 0 || lower <= 0 || upper <= 0) {
        return(rep(mean, n_sim))
      }
      sdlog <- (log(upper) - log(lower)) / (2 * 1.96)
      meanlog <- log(mean) - 0.5 * sdlog^2
      return(rlnorm(n_sim, meanlog = meanlog, sdlog = sdlog))
            
    } else {
      ## gamma distribution
      if (dist_clean == "gamma") {
        if (mean == 0 || upper == lower) {
          return(rep(mean, n_sim))
        }
        sd <- (upper - lower) / (2 * 1.96)
        shape <- (mean / sd)^2
        rate  <- shape / mean
        if (shape <= 0 || rate <= 0) {
          return(rep(mean, n_sim))
        } else {
          return(rgamma(n_sim, shape = shape, rate = rate))
        }
        
      } else {
        ## normal distribution
        if (dist_clean == "normal") {
          sd <- (upper - lower) / (2 * 1.96)
          return(rnorm(n_sim, mean = mean, sd = sd))
          
        } else {
          ## PERT distribution
          if (dist_clean == "pert") {
            if (lower > mean || mean > upper) {
              stop("Require: min ≤ mode ≤ max")
            } else {
              ## transform pert to beta
              alpha <- (mean - lower) * 4 / (upper - lower) + 1
              beta  <- (upper - mean) * 4 / (upper - lower) + 1
              return(rbeta(n_sim, alpha, beta) * (upper - lower) + lower)
            }
            
          } else {
            ## Triangular distribution
            if (dist_clean == "triangular") {
              if (lower > mean || mean > upper) {
                stop("Require: min ≤ mode ≤ max")
              } else {
                if (lower == upper) {
                  stop("Lower and upper bounds must not be equal")
                } else {
                  ## transform triangular to beta
                  alpha <- 2 * (mean - lower) / (upper - lower)
                  beta  <- 2 * (upper - mean) / (upper - lower)                  
                  if (alpha <= 0 || beta <= 0) {
                    stop("Invalid alpha or beta calculated")
                  }
                  return(rbeta(n_sim, alpha, beta) * (upper - lower) + lower)
                }
              }
              
            } else {
                if (dist_clean == "na") {
                return(rep(mean, n_sim))           
              } else {
                stop(paste("Unsupported distribution:", dist))
              }
            }
          }
        }
      }
    } 
  }
}

##******************************************
## Wrapper function for draws
##******************************************
draw_psa_params <- function(dist, mean, lower, upper, n_sim = 1000) {
  do.call(
    cbind,
    mapply(
      function(d, m, lo, up) {
        draw_params(dist = d, mean = m, lower = lo, upper = up, n_sim = n_sim)
      }, 
      dist, mean, lower, upper, SIMPLIFY = FALSE
    )
  )
}

##******************************************
## Function to run CEA PSA
##******************************************
## function to run model for a single iteration x wtp
run_single <- function(
  i, 
  wtp_val = NULL, 
  scenarios = scenarios,
  ref_scenario = scenarios[1],
  params = psa_params[1, ],
  xi_factor = NULL,
  zeta_factor = NULL,
  phi_factor = NULL,
  psi_stated_factor = NULL,
  pi_factor = NULL,
  chi_factor = NULL,
  eta_factor = NULL,
  omicron_factor = NULL,
  init_states = ini_pop_vec,
  inc_names = trans_names,
  time_horizon_years = time_horizon_set[1],
  perspectives = c(0, 1),
  le_death_data = le_death_data,
  fx_adjust = fx_adjust,
  cpi_adjust = cpi_adjust,
  bound_val = "base",
  par_name = "base",
  base_par_val = NA,
  bound_par_val = NA
) {

  ## run cea function
  res <- base_case_cea_func(
    wtp = wtp_val,
    scenarios = scenarios,
    ref_scenario = ref_scenario,
    params = params,
    xi_factor = xi_factor,
    zeta_factor = zeta_factor,
    phi_factor = phi_factor,
    psi_stated_factor = psi_stated_factor,
    pi_factor = pi_factor,
    chi_factor = chi_factor,
    eta_factor = eta_factor,
    omicron_factor = omicron_factor,
    init_states = init_states,
    inc_names = inc_names,
    time_horizon_years = time_horizon_years,
    perspectives = perspectives,
    le_death_data = le_death_data,
    fx_adjust = fx_adjust,
    cpi_adjust = cpi_adjust,
    bound_val = bound_val,
    par_name = par_name,
    base_par_val = base_par_val,
    bound_par_val = bound_par_val
  )

  ## fx factor
  fx_factor <- (1 / as.numeric(params[["brl_usd_fx_2024"]]))

  ## add identifiers
  add_id <- function(df, id, wtp, fx_val) {
    df$iteration <- id
    df$wtp <- wtp * fx_val
    return(df)
  }

  ## output
  val <- list(
    trace = add_id(df = res$trace, id = i, wtp = wtp_val, fx_val = 1),
    costs_brl = add_id(df = res$costs_brl, id = i, wtp = wtp_val, fx_val = 1),
    costs_usd = add_id(df = res$costs_usd, id = i, wtp = wtp_val, fx_val = fx_factor),
    qalys = add_id(df = res$qalys, id = i, wtp = wtp_val, fx_val = 1),
    disc_costs_brl = add_id(df = res$disc_costs_brl, id = i, wtp = wtp_val, fx_val = 1),
    disc_costs_usd = add_id(df = res$disc_costs_usd, id = i, wtp = wtp_val, fx_val = fx_factor),
    disc_qalys = add_id(df = res$disc_qalys, id = i, wtp = wtp_val, fx_val = 1),
    cea_out_brl = add_id(df = res$cea_out_brl, id = i, wtp = wtp_val, fx_val = 1),
    cea_out_usd = add_id(df = res$cea_out_usd, id = i, wtp = wtp_val, fx_val = fx_factor)
  )

  return(val)
}

##******************************
## psa cea function
##******************************
cea_psa_func <- function(
  gdp_factors = gdp_multipliers,
  scenarios = scenarios,
  ref_scenario = scenarios[1],
  psa_params = psa_params,
  xi_factor = NULL,
  zeta_factor = NULL,
  phi_factor = NULL,
  psi_stated_factor = NULL,
  pi_factor = NULL,
  chi_factor = NULL,
  eta_factor = NULL,
  omicron_factor = NULL,
  init_states = ini_pop_vec,
  inc_names = trans_names,
  time_horizon_years = time_horizon_set[1],
  perspectives = c(0, 1),
  le_death_data = le_death_data,
  fx_adjust = fx_adjust,
  cpi_adjust = cpi_adjust,
  bound_val = "base",
  par_name = "base",
  base_par_val = NA,
  bound_par_val = NA
){
  ## lapply over PSA draws
  outer_results <- future.apply::future_lapply(
    X = seq_len(nrow(psa_params)),
    FUN = function(i){
      pars <- as.numeric(psa_params[i, ])
      names(pars) <- colnames(psa_params)
      gdp <- as.numeric(pars[["GDP_Per_Capita"]])
      wtp_vec <- gdp * gdp_factors

      ## lapply over wtp values
      inner_results <- lapply(
        X = wtp_vec,
        FUN = function(x){
          run_single(
            i = i, 
            wtp_val = x, 
            scenarios = as.vector(scenarios),
            ref_scenario = ref_scenario,
            params = pars,
            xi_factor = xi_factor,
            zeta_factor = zeta_factor,
            phi_factor = phi_factor,
            psi_stated_factor = psi_stated_factor,          
            pi_factor = pi_factor,
            chi_factor = chi_factor,
            eta_factor = eta_factor,
            omicron_factor = omicron_factor,
            init_states = init_states,
            inc_names = inc_names,
            time_horizon_years = time_horizon_years,
            perspectives = perspectives,
            le_death_data = le_death_data,
            fx_adjust = fx_adjust,
            cpi_adjust = cpi_adjust,  
            bound_val = bound_val,
            par_name = par_name,
            base_par_val = base_par_val,
            bound_par_val = bound_par_val
          )
        }
      )
      return(inner_results)
    }
  )

  ## flatten nested list
  flat <- do.call(c, outer_results)  

  ## combine by element name
  out <- list(
    trace = dplyr::bind_rows(lapply(flat, `[[`, "trace")),
    costs_brl = dplyr::bind_rows(lapply(flat, `[[`, "costs_brl")),
    costs_usd = dplyr::bind_rows(lapply(flat, `[[`, "costs_usd")),
    qalys = dplyr::bind_rows(lapply(flat, `[[`, "qalys")),
    disc_costs_brl = dplyr::bind_rows(lapply(flat, `[[`, "disc_costs_brl")),
    disc_costs_usd = dplyr::bind_rows(lapply(flat, `[[`, "disc_costs_usd")),
    disc_qalys = dplyr::bind_rows(lapply(flat, `[[`, "disc_qalys")),
    cea_out_brl = dplyr::bind_rows(lapply(flat, `[[`, "cea_out_brl")),
    cea_out_usd = dplyr::bind_rows(lapply(flat, `[[`, "cea_out_usd"))
  )

  ## output
  return(out)
}

##******************************************
## function to prepare CEAC plot data
##******************************************
inner_ceac_data_func <- function(
  currency_val = c("brl", "usd")[1],
  data = psa_data
) {

  ## prepare the data
  val <- data[[paste0("cea_out_", currency_val)]] %>%
    dplyr::select(
      iteration, 
      wtp, 
      perspective, 
      alternative, 
      nmb
    ) %>%
    ## keep psa structure
    pivot_wider(
      names_from = alternative,
      values_from = nmb
    ) %>%
    ## best strategy per PSA draw
    mutate(
      best_strategy = pmap_chr(
        list(SOC, IXCHIQ, VIMKUNYA),
        function(a, b, c) {
          x <- c(SOC = a, IXCHIQ = b, VIMKUNYA = c)
          if (all(is.na(x))) {
            NA_character_
          } else {
            names(which.max(x))
          }
        }
      )
    ) %>%    
    ## ceac by wtp & perspective
    dplyr::group_by(wtp, perspective) %>%
    dplyr::summarise(
      n = n(),
      ##
      SOC_k = sum(best_strategy == "SOC", na.rm = TRUE),
      IXCHIQ_k = sum(best_strategy == "IXCHIQ", na.rm = TRUE),
      VIMKUNYA_k = sum(best_strategy == "VIMKUNYA", na.rm = TRUE),
      ##
      SOC_prob = SOC_k / n,
      IXCHIQ_prob = IXCHIQ_k / n,
      VIMKUNYA_prob = VIMKUNYA_k / n,
      SOC_ci = list(binom.confint(SOC_k, n, method = "wilson")),
      IXCHIQ_ci = list(binom.confint(IXCHIQ_k, n, method = "wilson")),
      VIMKUNYA_ci = list(binom.confint(VIMKUNYA_k, n, method = "wilson")),
      .groups = "drop"
    ) %>%   
    ## unnest CIs
    tidyr::unnest_wider(SOC_ci, names_sep = "_") %>%
    tidyr::unnest_wider(IXCHIQ_ci, names_sep = "_") %>%
    tidyr::unnest_wider(VIMKUNYA_ci, names_sep = "_") %>%
    dplyr::rename(
      SOC_lower = SOC_ci_lower,
      SOC_upper = SOC_ci_upper,
      IXCHIQ_lower = IXCHIQ_ci_lower,
      IXCHIQ_upper = IXCHIQ_ci_upper,
      VIMKUNYA_lower = VIMKUNYA_ci_lower,
      VIMKUNYA_upper = VIMKUNYA_ci_upper
    ) %>%
    mutate(currency = currency_val)

  ## long format for plotting
  val_long <- val %>%
    pivot_longer(
      cols = matches("SOC|IXCHIQ|VIMKUNYA"),
      names_to = c("alternative", ".value"),
      names_pattern = "(SOC|IXCHIQ|VIMKUNYA)_(prob|lower|upper)"
    ) %>%
    dplyr::filter(!is.na(alternative)) %>%
    dplyr::arrange(perspective, wtp, alternative) %>%
    mutate(currency = currency_val)

  return(list(wide = val, long = val_long))
}


## function for the ceac plot data; provides long and wide form data
ceac_data_func <- function(
  currencies = c("brl", "usd"),
  data = psa_data
) {
  lst <- lapply(
    X = currencies, 
    FUN = function(i) {
      inner_ceac_data_func(currency = i, data = data)
    }
  )
  final <- Map(
    function(...) do.call(rbind, list(...)), lst[[1]], lst[[2]]
  )
  return(final)
}

##******************************************
## function to plot psa results
##******************************************
psa_plot_func <- function(
  data = psa_plot_data$long,
  x = "wtp", y = "prob", color = "alternative", fill = "alternative",
  ymin = "lower", ymax = "upper",
  ##
  alt_labs = c(
    "IXCHIQ" = "Live-attenuated vaccine",
    "VIMKUNYA" = "Recombinant vaccine"
  ),
  perspective_labs = c(
    "0" = "(a) health Care Sector Perspective",
    "1" = "(b) Societal Perspective"
  ),
  x_lab = "Willingness-to-pay (WTP) threshold, 2024 US$'000",
  y_lab = "Probability of being cost-effective",
  ##
  color_lab = "Vaccination strategy",
  fill_lab = "Vaccination strategy",
  ##
  x_axis_rescale_by = 1000,
  x_min = 25,
  x_max = 225,
  x_by = 50,
  x_exp_min = 0.025,
  x_exp_max = 0.03,
  ##
  y_min = 0,
  y_max = 1,
  y_by = 0.2,
  y_exp_min = 0.01,
  y_exp_max = 0.01,
  ##
  strip_text = facet_strip_text,
  strip_txt_size = 12,
  ##
  annotate_txt_size = 4,
  x_offset_val = -1.0,
  hjust_val_min = 1.1,
  hjust_val_max = -0.1,
  ##
  x_axis_txt_size = 10,
  x_axis_title_size = 12,
  y_axis_txt_size = 10,
  y_axis_title_size = 12,
  legend_txt_size = 12, 
  journal = "bmj",
  ##
  save_plot = FALSE,
  file_name = "ceac.png", 
  file_path = figures_dir,
  fig_scale = 1, 
  fig_width = 7, 
  fig_height = 5, 
  units = c("in", "cm", "mm", "px")[1], 
  dpi = 300, 
  limitsize = TRUE,
  bg = "white"
){

  ## tidy-eval symbols
  x_sym <- sym(x)
  y_sym <- sym(y)
  ymin_sym <- sym(ymin)
  ymax_sym <- sym(ymax)
  col_sym <- sym(color)
  fill_sym <- sym(fill)

  ## plot data
  plot_data <- data %>%
    filter(.data[[color]] %in% names(alt_labs)) %>%
    mutate(
      alternative = factor(
        .data[[color]],
        levels = names(alt_labs),
        labels = alt_labs
      ),
      perspective = factor(
        perspective,
        levels = names(perspective_labs),
        labels = perspective_labs
      )
    )
  
    ## produce the plot
    psa_plot <- ggplot(
      plot_data,
      aes(
        x = !!x_sym / x_axis_rescale_by,
        y = !!y_sym,
        color = alternative,
        fill  = alternative
      )
    ) +
      geom_ribbon(
        aes(
          ymin = !!ymin_sym, 
          ymax = !!ymax_sym
        ),
        alpha = 0.15,
        color = NA
      ) +
      geom_line(linewidth = 1) +
      facet_grid(~ perspective) +
      labs(
        x = x_lab,
        y = y_lab,
        color = color_lab,
        fill = fill_lab
      ) +
      theme_minimal() +
      theme(
        panel.grid = element_blank(),
        axis.line = element_line(color = "black", linewidth = 0.5),
        axis.ticks = element_line(color = "black"),
        axis.text.x = element_text(size = x_axis_txt_size),
        axis.text.y = element_text(size = y_axis_txt_size),
        axis.title.x = element_text(size = x_axis_title_size, face = "bold"),
        axis.title.y = element_text(size = y_axis_title_size, face = "bold"),
        legend.text = element_text(size = legend_txt_size),
        legend.position = "bottom",
        legend.background = element_rect(fill = "white", color = NA),
        legend.key.size = unit(0.5, "cm"),
        plot.margin = unit(c(0.3, 0.3, 0.2, 0.2), "cm"),
        strip.text = element_text(size = strip_txt_size, face = "bold"),
        strip.background = element_blank()
      ) +
      scale_x_continuous(
        labels = scales::comma,
        limits = c(x_min, x_max) / x_axis_rescale_by,
        breaks = seq(x_min, x_max, x_by) / x_axis_rescale_by,
        expand = expansion(mult = c(x_exp_min, x_exp_max))
      ) +
      scale_y_continuous(
        limits = c(y_min, y_max),
        breaks = seq(y_min, y_max, y_by),
        expand = expansion(mult = c(y_exp_min, y_exp_max))
      ) +
      get_journal_color_scale(journal) +
      get_journal_fill_scale(journal)
 
  ## save
  if(isTRUE(save_plot)){
    ggsave(
      filename = file_name, 
      plot = psa_plot, 
      path = file_path,
      scale = fig_scale, 
      width = fig_width, 
      height = fig_height, 
      units = units, 
      dpi = dpi, 
      limitsize = limitsize, 
      bg = bg
    )
  }
  
  return(psa_plot)
}

# ##******************************************
# ## Function to calculate bootstrap means
# ##******************************************
# inner_boot_mean_ci_func <- function(
#   currency_val = c("brl", "usd")[2],
#   data = psa_data,
#   nbootsamps = 1000,
#   params = params_value,
#   scenarios = scenarios
# ) {
# 
#   ## wtp
#   gdp <- as.numeric(params[["GDP_Per_Capita"]])
#   if(currency_val == "usd"){
#     wtp_val <- gdp * (1 / as.numeric(params[["brl_usd_fx_2024"]]))
#   } else {
#     wtp_val <- gdp
#   }
#  
#   ## psa data for the selected currency and WTP
#   boot_data <- data[[paste0("cea_out_", currency_val)]] %>%
#     dplyr::filter(wtp == wtp_val)
# 
#   boot_results <- replicate(
#     nbootsamps, 
#     {
#       ## resample within each alternative
#       samp <- boot_data %>%
#         dplyr::select(perspective, alternative, costs, qalys) %>%
#         group_by(perspective, alternative) %>%
#         mosaic::resample() %>%
#         dplyr::summarise(
#           mean_costs = if (n() == 0) NA_real_ else mean(costs),
#           mean_qalys = if (n() == 0) NA_real_ else mean(qalys),
#           .groups = "drop"
#         )
# 
#       ## compute incremental outcomes sequentially
#       samp <- samp %>%
#         dplyr::filter(!is.na(mean_costs), !is.na(mean_qalys)) %>%
#         mutate(
#           nmb = mean_qalys * wtp_val - mean_costs,
#           nhb = mean_qalys - mean_costs / wtp_val
#         ) %>%
#         group_by(perspective) %>%
#         arrange(nmb) %>%
#         mutate(strategy_order = row_number())
# 
#       # Create all pairwise comparisons
#       incremental <- samp %>%
#       group_by(perspective) %>%
#       dplyr::filter(n() > 1) %>% 
#       expand_grid(
#         reference = strategy_order,
#         target = strategy_order
#       ) %>%
#       #filter(reference != target) %>%
#       # join target and reference rows
#       left_join(
#         samp %>% rename_with(~ paste0(.x, "_t"), -perspective), 
#         by = c("perspective", "target" = "strategy_order_t")
#       ) %>%
#       left_join(
#         samp %>% rename_with(~ paste0(.x, "_r"), -perspective),
#         by = c("perspective", "reference" = "strategy_order_r")
#       ) %>%
#       mutate(
#         inc_costs = mean_costs_t - mean_costs_r,
#         inc_qalys = mean_qalys_t - mean_qalys_r,
#         inmb = nmb_t - nmb_r,
#         inhb = nhb_t - nhb_r,
#         icer = ifelse(
#           inc_qalys <= 0,
#           NA_real_,
#           inc_costs / inc_qalys
#         )
#       ) %>%
#       dplyr::select(
#         perspective, alternative_t, alternative_r, 
#         mean_costs_t, mean_qalys_t,
#         nmb_t, nhb_t,
#         inc_costs, inc_qalys, 
#         inmb, inhb, icer
#       ) %>%
#       dplyr::rename(
#         alternative = alternative_t,
#         reference = alternative_r,
#         costs = mean_costs_t, 
#         qalys = mean_qalys_t,
#         nmb = nmb_t,
#         nhb = nhb_t
#       ) #%>%
#       #dplyr::filter(reference < alternative)
#       
#       incremental
#     }, 
#     simplify = FALSE
#   ) %>%
#   bind_rows(.id = "bootstrap")
# 
#   ## summarise bootstrap distribution: mean, SE, percentile CI
#   ## helper function to summarise mean, SE, and 95% CI
#   summarise_boot <- function(x) {
#     x <- na.omit(x)
#   
#     if (length(x) == 0) {
#       return(tibble::tibble(
#         statistic = c("mean", "se", "lower", "upper"),
#         value = NA_real_
#       ))
#     }
#   
#     tibble::tibble(
#       statistic = c("mean", "se", "lower", "upper"),
#       value = c(
#         mean(x),
#         sd(x),
#         quantile(x, 0.025, names = FALSE),
#         quantile(x, 0.975, names = FALSE)
#       )
#     )
#   }
#   
#   ## columns to summarise
#   cols_to_summarise <- c(
#     "costs", "qalys", 
#     "nmb", "nhb",
#     "inc_costs", "inc_qalys", 
#     "icer", "inmb", "inhb"
#   )
#   
#   ## summarize
#   summary_boot <- boot_results %>%
#     group_by(perspective, alternative, reference) %>%
#     reframe(
#       across(
#         all_of(cols_to_summarise),
#         ~ list(summarise_boot(.x))
#       )
#     ) %>%
#     ungroup()
#   
#   summary_boot_long <- summary_boot %>%
#     tidyr::pivot_longer(
#       cols = all_of(cols_to_summarise),
#       names_to = "measure",
#       values_to = "stats"
#     ) %>%
#     tidyr::unnest(stats)
# 
#     if(nrow(summary_boot_long) == 0){
#       out_boot <- NULL
#     } else {
#       summary_boot_wide <- summary_boot_long %>%
#         tidyr::pivot_wider(
#           names_from = statistic,
#           values_from = value
#         )
#       
#     
#       out_boot <- summary_boot_wide %>%
#         arrange(perspective, reference, alternative) %>%
#         mutate(currency = currency_val)
#     }
#   
#   return(out_boot)
# }
# 
# ##******************************************
# ## Main function for the CEA bootstrap CIs
# ##******************************************
# boot_mean_ci_func <- function(
#   currencies = c("brl", "usd"),
#   data = psa_data,
#   nbootsamps = 1000,
#   params = params_value,
#   scenarios = scenarios
# ) {
#   do.call(
#     rbind,
#     lapply(
#       X = currencies,
#       FUN = inner_boot_mean_ci_func,
#       data = data,
#       params = params,
#       scenarios = scenarios
#     ) 
#   )
# }




##******************************************
## Inner function: Bayesian CEA summaries
##******************************************
inner_bayes_ci_func <- function(
  currency_val = c("brl", "usd")[2],
  data = psa_data,
  params = params_value,
  scenarios = scenarios
) {

  ## willingness-to-pay
  gdp <- as.numeric(params[["GDP_Per_Capita"]])
  if (currency_val == "usd") {
    wtp_val <- gdp * (1 / as.numeric(params[["brl_usd_fx_2024"]]))
  } else {
    wtp_val <- gdp
  }

  ## psa data for currency + wtp
  psa_data_wtp <- data[[paste0("cea_out_", currency_val)]] %>%
    dplyr::filter(wtp == wtp_val)

  if (nrow(psa_data_wtp) == 0) {
    return(NULL)
  }

  ## compute NMB / NHB per draw
  psa_base <- psa_data_wtp %>%
    dplyr::mutate(
      nmb = qalys * wtp_val - costs,
      nhb = qalys - costs / wtp_val
    )

  ## sequential ordering per draw
  psa_ordered <- psa_base %>%
    dplyr::group_by(iteration, perspective) %>%
    dplyr::arrange(nmb, .by_group = TRUE) %>%
    dplyr::mutate(strategy_order = dplyr::row_number()) %>%
    dplyr::filter(dplyr::n() > 1) %>%
    dplyr::ungroup()

  ## pairwise incremental outcomes
  incremental_psa <- psa_ordered %>%
    dplyr::inner_join(
      psa_ordered,
      by = c("iteration", "perspective"),
      suffix = c("_t", "_r")
    ) %>%
    dplyr::mutate(
      inc_costs = costs_t - costs_r,
      inc_qalys = qalys_t - qalys_r,
      inmb      = nmb_t - nmb_r,
      inhb      = nhb_t - nhb_r,
      icer      = dplyr::if_else(
        inc_qalys <= 0,
        NA_real_,
        inc_costs / inc_qalys
      )
    ) %>%
    dplyr::select(
      iteration,
      perspective,
      alternative = alternative_t,
      reference   = alternative_r,
      costs       = costs_t,
      qalys       = qalys_t,
      nmb         = nmb_t,
      nhb         = nhb_t,
      inc_costs,
      inc_qalys,
      inmb,
      inhb,
      icer
    )
  
  ## bayesian posterior summariser
  summarise_bayes <- function(x) {
    x <- na.omit(x)
    if (length(x) == 0) {
      return(
        tibble::tibble(
          mean  = NA_real_,
          lower = NA_real_,
          upper = NA_real_
        )
      )
    }

    tibble::tibble(
      mean  = mean(x),
      lower = quantile(x, 0.025, names = FALSE),
      upper = quantile(x, 0.975, names = FALSE)
    )
  }

  ## columns to summarise
  cols_to_summarise <- c(
    "costs", "qalys", 
    "nmb", "nhb",
    "inc_costs", "inc_qalys", 
    "icer", "inmb", "inhb"
  )

  ## summarise across psa draws
  bayes_summary <- incremental_psa %>%
    group_by(perspective, alternative, reference) %>%
    reframe(
      across(
        all_of(cols_to_summarise),
        ~ list(summarise_bayes(.x))
      )
    ) %>%
    ungroup()
  
  bayes_summary_long <- bayes_summary %>%
    tidyr::pivot_longer(
      cols = all_of(cols_to_summarise),
      names_to = "measure",
      values_to = "stats"
    ) %>%
    tidyr::unnest(stats)
  
  ## final formatting
  out_bayes <- bayes_summary_long %>%
    dplyr::arrange(perspective, reference, alternative) %>%
    dplyr::mutate(
      currency = currency_val,
      interval_type = "bayesian_credible"
    )

  return(out_bayes)
}


##******************************************
## Main function: Bayesian CEA CIs
##******************************************
bayes_mean_ci_func <- function(
  currencies = c("brl", "usd"),
  data = psa_data,
  params = params_value,
  scenarios = scenarios
) {

  do.call(
    rbind,
    lapply(
      X = currencies,
      FUN = inner_bayes_ci_func,
      data = data,
      params = params,
      scenarios = scenarios
    )
  )
}


##******************************************
## Function to create a table of CEA results from the multivariate analysis
##******************************************
mult_cea_table_func <- function(
  psa_boot_data = cea_psa_boot_data,
  currency_val = c("brl", "usd")[1],
  scenarios = scenarios, 
  scenarios_labs = scenarios_labs,
  journal = "bmj",
  digits = digits_vals,
  scale = scale_vals,
  linebreak = "\n",
  big.mark = big.mark_vals,
  decimal.mark = decimal.mark_vals,
  ci_sep = NULL,
  table_caption = "Table 3. Cost-effectiveness Analysis Results from Multivariate Sensitivity Analyses",
  Costs_lab = "Costs (US$)\n Est. (95% CI)",
  QALYs_lab = "QALYs\n Est. (95% CI)",
  NMB_lab = "NMB (US$ million)\n Est. (95% CI)",
  INMB_lab = "INMB (US$ million)\n Est. (95% CI)",
  tables_dir = tables_dir, 
  table_fname = "multivariate_cea_results_table.docx"
){

  # ---- reshape long to wide ----
  table_data <- psa_boot_data %>%
    dplyr::filter(currency == currency_val) %>%
    dplyr::select(
      perspective, alternative, reference, measure, mean, lower, upper
    ) %>%
    tidyr::pivot_wider(
      names_from = measure,
      values_from = c(mean, lower, upper),
      names_sep = "_"
    ) %>%
    dplyr::rename(
      costs_mean = mean_costs,
      costs_lower = lower_costs,
      costs_upper = upper_costs,
      qalys_mean = mean_qalys,
      qalys_lower = lower_qalys,
      qalys_upper = upper_qalys,
      nmb_mean = mean_nmb,
      nmb_lower = lower_nmb,
      nmb_upper = upper_nmb,
      inmb_mean = mean_inmb,
      inmb_lower = lower_inmb,
      inmb_upper = upper_inmb
    ) %>%
    mutate(
      reference = factor(reference, levels = scenarios, labels = scenarios_labs),
      alternative = factor(alternative, levels = scenarios, labels = scenarios_labs)
    ) %>%
    dplyr::arrange(perspective, reference, alternative)

  # ---- format values ----
  table_data <- table_data %>%
    mutate(
      Costs = fmt_ci_journal(
        mean = costs_mean, 
        lower = costs_lower, 
        upper = costs_upper,
        journal = journal,
        digits = digits$Costs,
        scale = scale$Costs,
        linebreak = linebreak,
        big.mark = big.mark$Costs,
        decimal.mark = decimal.mark$Costs,
        ci_sep = ci_sep
      ),
      QALYs = fmt_ci_journal(
        mean = qalys_mean, 
        lower = qalys_lower, 
        upper = qalys_upper,
        journal = journal,
        digits = digits$QALYs,
        scale = scale$QALYs,
        linebreak = linebreak,
        big.mark = big.mark$QALYs,
        decimal.mark = decimal.mark$QALYs,
        ci_sep = ci_sep
      ),
      NMB = fmt_ci_journal(
        mean = nmb_mean, 
        lower = nmb_lower, 
        upper = nmb_upper,
        journal = journal,
        digits = digits$NMB,
        scale = scale$NMB,
        linebreak = linebreak,
        big.mark = big.mark$NMB,
        decimal.mark = decimal.mark$NMB,
        ci_sep = ci_sep
      ),
      INMB = fmt_ci_journal(
        mean = inmb_mean, 
        lower = inmb_lower, 
        upper = inmb_upper,
        journal = journal,
        digits = digits$INMB,
        scale = scale$INMB,
        linebreak = linebreak,
        big.mark = big.mark$INMB,
        decimal.mark = decimal.mark$INMB,
        ci_sep = ci_sep
      )
    ) %>%
    dplyr::select(perspective, alternative, reference, Costs, QALYs, NMB, INMB) %>%
    mutate(
      across(
        c(Costs, QALYs, NMB, INMB),
        ~ ifelse(. %in% c("0\n(0, 0)", "0\n(0,0)", "0.00\n(0.00, 0.00)"), "–", .)
      )
    )

  # ---- create flextable ----
  ft <- flextable(table_data) %>%
    set_header_labels(
      perspective = "Perspective",
      alternative = "Alternative",
      reference = "Reference",
      Costs = Costs_lab,
      QALYs = QALYs_lab,
      NMB = NMB_lab,
      INMB = INMB_lab
    ) %>%
    align(align = "center", part = "header") %>%
    theme_vanilla() %>%
    autofit() %>%
    bold(part = "header")

  # ---- save docx ----
  doc <- read_docx() %>%
    body_add_par(value = table_caption, style = "heading 1") %>%
    body_add_flextable(ft) %>%
    body_add_par(value = "", style = "Normal")

  print(doc, target = file.path(tables_dir, table_fname))
}


##************************************************************
## Sensitivity Analyses
##************************************************************
sensi_func <- function(
  xi_vals = NULL,
  zeta_vals = NULL,
  phi_vals = NULL,
  psi_stated_vals = NULL,
  pi_vals = NULL,
  chi_vals = NULL,
  eta_vals = NULL,
  omicron_vals = NULL,
  perspectives = perspectives[2],
  time_horizon_years = time_horizon_set[1],
  wtp = NULL,
  scenarios = scenarios,
  ref_scenario = scenarios[1],
  params, 
  init_states,
  inc_names,
  le_death_data,
  fx_adjust,
  cpi_adjust
) {

  ## create combinations
  iter_combo <- expand.grid(
    xi = if (is.null(xi_vals)) NA_real_ else xi_vals,
    zeta = if (is.null(zeta_vals)) NA_real_ else zeta_vals,
    phi = if (is.null(phi_vals)) NA_real_ else phi_vals,
    psi_stated = if (is.null(psi_stated_vals)) NA_real_ else psi_stated_vals,
    pi = if (is.null(pi_vals)) NA_real_ else pi_vals,
    chi = if (is.null(chi_vals)) NA_real_ else chi_vals,
    eta = if (is.null(eta_vals)) NA_real_ else eta_vals,
    omicron = if (is.null(omicron_vals)) NA_real_ else omicron_vals,
    perspective = perspectives,
    time_horizon_years = as.numeric(time_horizon_years),
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
  )
 
  beta_S_base <- 1 - exp(
    -((1 / as.numeric(params["gamma_days"])) *
        as.numeric(params["R0"])) *
      as.numeric(params["cycle_length_days"])
  )
  
  cea_vals_list <- lapply(
    X = seq_len(nrow(iter_combo)),
    FUN = function(i) {

      xi_val <- iter_combo$xi[i]
      zeta_val <- iter_combo$zeta[i]
      phi_val <- iter_combo$phi[i]
      psi_stated_val <- iter_combo$psi_stated[i]
      pi_val <- iter_combo$pi[i]
      chi_val <- iter_combo$chi[i]
      eta_val <- iter_combo$eta[i]      
      omicron_val <- iter_combo$omicron[i]   
      
      beta_S  <- beta_S_base
      beta_SV <- beta_S_base * xi_val

      perspective_val <- iter_combo$perspective[i]
      time_horizon_val <- iter_combo$time_horizon_years[i]

      ## cea analysis
      cea_vals <- base_case_cea_func(
        wtp = wtp,
        scenarios = scenarios,
        ref_scenario = ref_scenario,
        params = params,
        xi_factor = if (is.na(xi_val)) NULL else xi_val,
        zeta_factor = if (is.na(zeta_val)) NULL else zeta_val,
        phi_factor = if (is.na(phi_val)) NULL else phi_val,
        psi_stated_factor = if (is.na(psi_stated_val)) NULL else psi_stated_val,      
        pi_factor = if (is.na(pi_val)) NULL else pi_val,
        chi_factor = if (is.na(chi_val)) NULL else chi_val,
        eta_factor = if (is.na(eta_val)) NULL else eta_val,
        omicron_factor = if (is.na(omicron_val)) NULL else omicron_val,
        init_states = init_states,
        inc_names = inc_names,
        time_horizon_years = time_horizon_val,
        perspectives = perspective_val,
        le_death_data = le_death_data,
        fx_adjust = fx_adjust,
        cpi_adjust = cpi_adjust,
        bound_val = "base",
        par_name = "base",
        base_par_val = NA,
        bound_par_val = NA
      )

      ## attach metadata (base R faster than mutate)
      val <- lapply(
        X = cea_vals, 
        FUN = function(obj) {

          if (inherits(obj, "data.frame")) {

            obj$xi <- xi_val
            obj$zeta <- zeta_val
            obj$phi <- phi_val
            obj$psi_stated <- psi_stated_val
            obj$pi <- pi_val
            obj$chi <- chi_val
            obj$eta <- eta_val
            obj$omicron <- omicron_val

            obj$beta_S <- beta_S
            obj$beta_SV <- beta_SV

            obj$perspective <- perspective_val
            obj$time_horizon_years <- time_horizon_val
          }

          return(obj)

        }
      )

      return(val)
    }
  )
  
  ## combine results by name
  output_names <- names(cea_vals_list[[1]])
  combo_cea_vals <- lapply(
    X = output_names, 
    FUN = function(nm) {
      objs <- lapply(
        X = cea_vals_list, 
        FUN = `[[`, 
        nm
      )
      if (all(vapply(objs, inherits, logical(1), "data.frame"))) {
        do.call(rbind, objs)
      } else {
        objs
      }
    }
  )
  names(combo_cea_vals) <- output_names
  combo_cea_vals
}

##************************************************************
## behavioral feedback function
##************************************************************
behavior_table_func <- function(
  data = sensi_cea_data,
  currency = "usd",
  tables_dir = tables_dir,
  table_fname = "sensi_table.docx"  
) {

  ## incidence data
  sensi_inc <- data$trace %>%
    dplyr::select(cycle, perspective, scenario, beta_S, xi, beta_SV, EtoI) %>%
    group_by(perspective, scenario, xi) %>%
    dplyr::arrange(cycle, .by_group = TRUE) %>%
    dplyr::summarise(
      mean_beta_S  = mean(beta_S, na.rm = TRUE),
      mean_beta_SV = mean(beta_SV, na.rm = TRUE),
      total_EtoI = round(sum(EtoI, na.rm = TRUE), 6),
      .groups = "drop"
    ) %>%
    dplyr::rename(
      beta_S = mean_beta_S, 
      beta_SV = mean_beta_SV
    ) %>%
    dplyr::arrange(xi, perspective, scenario)  
  
  ## cea data
  sensi_cea <- data[[paste0("cea_out_", currency)]] %>%
    dplyr::rename(scenario = alternative) %>%
    dplyr::select(
      perspective, scenario, 
      beta_S, xi, beta_SV, costs, qalys, 
      nmb, inc_costs, inc_qalys, inmb, icer
    ) %>%
    mutate(
      nmb = nmb * 1e-6,
      inmb = inmb * 1e-6
    ) %>%
    group_by(perspective, scenario, xi) %>%
    dplyr::arrange(xi, perspective, scenario)  
  
  ## sensitivity table data
  sensi_table <- merge(
    x = sensi_inc, 
    y = sensi_cea, 
    by = c("perspective", "scenario", "beta_S", "xi", "beta_SV")
  )

  ## convert NMB to millions with formatting
  sensi_table_ft <- sensi_table %>%
    mutate(
      beta_S = sprintf("%.3f", beta_S),
      xi = sprintf("%.3f", xi),
      beta_SV = sprintf("%.3f", beta_SV),
      total_EtoI = sprintf("%s", format(round(total_EtoI, 0), big.mark = ",", nsmall = 0)),
      nmb = sprintf("%s", format(round(nmb, 0), big.mark = ",", nsmall = 0)),
      inmb = sprintf("%s", format(round(inmb, 2), big.mark = ",", nsmall = 2)),
      icer = sprintf("%s", format(round(icer, 0), big.mark = ",", nsmall = 0)),
    ) %>%
    dplyr::select(perspective, scenario, beta_S, xi, beta_SV, total_EtoI, nmb, inmb, icer)
  
  # Create flextable with plain text headers
  ft <- flextable(sensi_table_ft) %>%
    set_header_labels(
      pertspective = "",
      scenario = "",
      beta_S = "β_S",
      xi = "ξ",
      beta_SV = "β_SV",
      total_EtoI = "Cumulative Incidence",
      nmb = "NMB (US$ M)",
      inmb = "INMB (US$ M)",
      icer = "ICER (US$/QALY)"
    ) %>%
    bold(part = "header") %>%
    align(j = c("beta_S", "xi", "beta_SV", "total_EtoI", "nmb", "inmb", "icer"),
          align = "center", part = "all") %>%
    font(fontname = "Times New Roman", part = "all") %>%
    fontsize(size = 12, part = "all") %>%
    autofit()
  
  ## Add table title
  ft <- add_header_lines(ft, values = "Behavioural Assumption About Vaccination") %>%
    bold(i = 1, part = "header") %>%
    fontsize(i = 1, size = 12, part = "header")
  
  # Save to Word doc
  doc <- read_docx() %>% body_add_flextable(ft)
  print(
    doc, 
    target = file.path(tables_dir, table_fname)
  )
}

##************************************************************
## pi, phi and omicron sensitivity analysis function
##************************************************************
phi_pi_omicron_heatmap_func <- function(
  data = phi_pi_sensi_cea_data$cea_out_usd, 
  ref_strategy = "SOC",
  perspective_val = 1,
  outcome_var = "inmb",
  outcome_scaling_factor = 1e-6,
  perspective_var = "perspective",
  alt_var = "alternative",
  phi_var = "phi", 
  pi_var = "pi", 
  omicron_var = "omicron",
  x_facet_var = "alternative",
  y_facet_var = "omicron",
  x_facet_labs = c(
    "IXCHIQ" = "Live-attenuated vaccine", 
    "VIMKUNYA" = "Recombinant vaccine"
  ),
  y_facet_labs = c(
    "0.5" = "&omicron; = 0.50", 
    "1" = "&omicron; = 1.00",
    "1.5" = "&omicron; = 1.50"
  ),
  scale_lab = "INMB (2024 US$ million): ",
  midpoint_val = 2.5,
  title = NULL,
  x_lab = expression("Vaccination acceptance realization factor, " * phi),
  y_lab = expression("Relative vaccine acceptance parameter, " * pi[v]),
  legend_position = "bottom",
  file_name = "phi_pi_sensi_cea_plot.png", 
  file_path = figures_dir,
  fig_scale = 1, 
  fig_width = 7, 
  fig_height = 6, 
  units = c("in", "cm", "mm", "px")[1],
  dpi = 300,
  limitsize = TRUE, 
  bg = "white"
) {
  
  ## prepare plot data
  plot_data <- data %>%
    dplyr::select(
      !!sym(perspective_var),
      !!sym(alt_var),
      !!sym(phi_var),
      !!sym(pi_var),
      !!sym(omicron_var),
      !!sym(outcome_var)
    ) %>%
    dplyr::filter(
      !(.data[[alt_var]] == ref_strategy),
      !is.na(.data[[outcome_var]]),
      !is.nan(.data[[outcome_var]]) ,
      .data[[perspective_var]] == perspective_val
    ) %>%
    distinct() %>%
    mutate(
      !!sym(outcome_var) := outcome_scaling_factor * !!sym(outcome_var),
      !!sym(x_facet_var) := factor(!!sym(x_facet_var), levels = names(x_facet_labs)),
      !!sym(y_facet_var) := factor(!!sym(y_facet_var), levels = names(y_facet_labs))    
    )
  
  ## plot
  p <- ggplot(
    plot_data,
    aes(
      x = .data[[phi_var]],
      y = .data[[pi_var]],
      fill = .data[[outcome_var]]
    )
  ) +
  geom_tile(color = "white", linewidth = 0) + 
  facet_grid(
    rows = vars(!!sym(x_facet_var)),
    cols = vars(!!sym(y_facet_var)),
    labeller = labeller(
      !!x_facet_var := x_facet_labs,
      !!y_facet_var := y_facet_labs
    )
  ) + 
  scale_fill_gradient2(
    name = scale_lab,
    high = "#CD202C",
    mid = "#B4CFEE",
    low = "#2A6EBB",
    midpoint = midpoint_val,
    labels = scales::comma
  ) +
  labs(
    x = x_lab,
    y = y_lab,
    title = title
  ) +
  scale_x_continuous(
    limits = c(0.40, 0.75),
    breaks = seq(0.40, 0.80, 0.1),
    expand = c(0, 0)
  ) +
  scale_y_continuous(
    limits = c(0.0, 1.20),
    breaks = seq(0.0, 1.20, 0.20),
    expand = c(0, 0)
  ) +
  theme_minimal(base_size = 11) +
  theme(
    panel.grid = element_blank(),
    panel.spacing.y = unit(6, "mm"),
    axis.line = element_line(color = "black", linewidth = 0.5),
    axis.ticks = element_line(color = "black"),
    # strip.text = element_text(face = "bold", size = 11),
    strip.text = ggtext::element_markdown(face = "bold", size = 10),
    plot.title = element_text(
      face = "bold",
      size = 11,
      hjust = 0.5
    ),
    legend.position = legend_position,
    legend.title = element_text(size = 10), 
    legend.text = element_text(size = 8),
    legend.key.height = unit(2, "mm")
  )

  ggsave(
    filename = file_name, 
    plot = p, 
    path = file_path,
    scale = fig_scale, 
    width = fig_width, 
    height = fig_height, 
    units = units, 
    dpi = dpi, 
    limitsize = limitsize, 
    bg = bg
  )

  ## output
  return(p) 
  
}

## combine perspectives
combo_phi_pi_omicron_heatmap_func <- function(
  data = phi_pi_sensi_cea_data$cea_out_usd, 
  perspectives = perspectives,
  ref_strategy = "SOC",
  outcome_var = "inmb",
  outcome_scaling_factor = 1e-6,
  perspective_var = "perspective",
  alt_var = "alternative",
  phi_var = "phi", 
  pi_var = "pi", 
  omicron_var = "omicron",
  x_facet_var = "alternative",
  y_facet_var = "omicron",
  x_facet_labs = c(
    "IXCHIQ" = "Live-attenuated vaccine", 
    "VIMKUNYA" = "Recombinant vaccine"
  ),
  y_facet_labs = c(
    "0.5" = "&omicron; = 0.50", 
    "1" = "&omicron; = 1.00",
    "1.5" = "&omicron; = 1.50"
  ),
  perspective_labs = c(
    "0" = "(a) Health Care Sector Perspective", 
    "1" = "(b) Societal Perspective"
  ),
  scale_lab = "INMB (2024 US$ million): ",
  midpoint_val = 2.5,
  title = NULL,
  x_lab = expression("Vaccination acceptance realization factor, " *phi),
  y_lab = expression("Relative vaccine acceptance parameter, " * pi[v]),
  legend_position = "bottom",
  file_name = "phi_pi_sensi_heatmap.png", 
  file_path = figures_dir,
  fig_scale = 1, 
  fig_width = 7, 
  fig_height = 6,
  combo_fig_scale = 0.8,  
  combo_fig_width = 13, 
  combo_fig_height = 13, 
  units = c("in", "cm", "mm", "px")[1],
  dpi = 300,
  limitsize = TRUE, 
  bg = "white"
){

  p_list <- lapply(
    X = perspectives,
    FUN = function(i){
      ## perspective-specific parameters
      if(i == 0){
        legend_pos <- "none"
        perspective_tag <- "hc"
      } else {
        legend_pos <- legend_position
        perspective_tag <- "soc"
      }

      ## plot
      phi_pi_omicron_heatmap_func(
        data = data,
        perspective_val = i,
        ref_strategy = ref_strategy,
        outcome_var = outcome_var,
        outcome_scaling_factor = outcome_scaling_factor,
        perspective_var = perspective_var,
        alt_var = alt_var,
        phi_var = phi_var, 
        pi_var = pi_var, 
        omicron_var = omicron_var,
        x_facet_var = x_facet_var,
        y_facet_var = y_facet_var,
        x_facet_labs = x_facet_labs,
        y_facet_labs = y_facet_labs,
        scale_lab = scale_lab,
        midpoint_val = midpoint_val,
        title = perspective_labs[[as.character(i)]],
        x_lab = x_lab,
        y_lab = y_lab,
        legend_position = legend_pos,
        file_name = paste0("phi_pi_sensi_heatmap_", perspective_tag, ".png"), 
        file_path = file_path,
        fig_scale = fig_scale, 
        fig_width = fig_width , 
        fig_height = fig_height, 
        units = units,
        dpi = dpi,
        limitsize = limitsize, 
        bg = bg
      )
    }
  )
  names(p_list) <- c("hc", "soc")

  ## stacked vertically
  combined_plot <- p_list$hc / p_list$soc +
    plot_annotation(
      theme = theme(
        plot.title = element_text(size = 12, face = "bold")
      )
    )

  ## display
  print(combined_plot)

  ## save the plot
  ggsave(
    filename = file_name, 
    plot = combined_plot, 
    path = figures_dir,
    scale = combo_fig_scale, 
    width = combo_fig_width, 
    height = combo_fig_height, 
    units = units, 
    dpi = dpi, 
    limitsize = limitsize, 
    bg = bg
  )
    
  invisible(combined_plot)

}

##************************************************************
## zeta sensitivity analysis function
##************************************************************
plot_zeta_outcome_heatmap <- function(
  data = zeta_sensi_cea_data$cea_out_usd, 
  ref_strategy = "SOC",
  alternative_var = "alternative",
  outcome_var = "inmb",
  outcome_scaling_factor = 1e-6,
  zeta_var = "zeta", 
  x_facet_var = "alternative",
  y_facet_var = "perspective",
  x_facet_labs = c(
    "IXCHIQ" = "Live-attenuated vaccine", 
    "VIMKUNYA" = "Recombinant vaccine"
  ),
  y_facet_labs = c(
    "0" = "Healthcare sector perspective", 
    "1" = "Societal perspective"
  ),
  y_lab = "INMB, 2024 US$ million",
  x_lab = expression("Immunogenicity attenuation factor, " * zeta),
  title = NULL,
  legend_position = "bottom",
  file_name = "zeta_sensi_cea_plot.png", 
  file_path = figures_dir,
  fig_scale = 1, 
  fig_width = 7, 
  fig_height = 6, 
  units = c("in", "cm", "mm", "px")[1],
  dpi = 300,
  limitsize = TRUE, 
  bg = "white",
  journal = "bmj"
) {
  
  ## prepare plot data
  plot_data <- data %>%
    dplyr::select(
      all_of(unique(c(alternative_var, x_facet_var, y_facet_var, zeta_var, outcome_var)))
    ) %>%
    dplyr::filter(
      !(.data[[alternative_var]] == ref_strategy),
      !is.na(.data[[outcome_var]]),
      !is.nan(.data[[outcome_var]])
    ) %>%
    distinct() %>%
    mutate(
      !!sym(outcome_var) := outcome_scaling_factor * !!sym(outcome_var),
      !!sym(x_facet_var) := factor(!!sym(x_facet_var), levels = names(x_facet_labs), labels = x_facet_labs),
      !!sym(y_facet_var) := factor(!!sym(y_facet_var), levels = names(y_facet_labs), labels = y_facet_labs)    
    )
  
  ## plot
  p <- ggplot(
    plot_data,
    aes(
      x = .data[[zeta_var]],
      y = .data[[outcome_var]],
      group = .data[[x_facet_var]],
      color = .data[[x_facet_var]],
    )
  ) +
  geom_line() + 
  facet_grid(
    cols = vars(!!sym(y_facet_var)),
    labeller = labeller(
      !!y_facet_var := y_facet_labs
    )
  ) + 
  get_journal_color_scale(journal) + 
  labs(
    x = x_lab,
    y = y_lab,
    title = title
  ) +
  scale_x_continuous(
    limits = c(0.70, 1.00),
    breaks = seq(0.70, 1.00, 0.1),
    expand = c(0.01, 0.01)
  ) +
  scale_y_continuous(
    limits = c(0, 6),
    breaks = seq(0, 6, 1),
    labels = scales::label_number(
      big.mark = ",",
      decimal.mark = ".",
      accuracy = 0.01
    ),
    expand = c(0.01, 0.01)
  ) +
  theme_minimal(base_size = 11) +
  theme(
    panel.grid = element_blank(),
    panel.spacing.y = unit(6, "mm"),
    axis.line = element_line(color = "black", linewidth = 0.5),
    axis.ticks = element_line(color = "black"),
    # strip.text = element_text(face = "bold", size = 11),
    strip.text = ggtext::element_markdown(face = "bold", size = 11),
    plot.title = element_text(
      face = "bold",
      size = 11,
      hjust = 0.5
    ),
    legend.position = legend_position,
    legend.title = element_blank(), 
    legend.text = element_text(size = 10),
    legend.key.height = unit(2, "mm")
  )

  ggsave(
    filename = file_name, 
    plot = p, 
    path = file_path,
    scale = fig_scale, 
    width = fig_width, 
    height = fig_height, 
    units = units, 
    dpi = dpi, 
    limitsize = limitsize, 
    bg = bg
  )

  ## output
  return(p) 
  
}

##************************************************************
## time-horizon sensitivity plot function
##************************************************************
one_way_time_horizon_plot_func <- function(
  data = sensi_time_horizon_cea_data$cea_out_usd,
  outcome_col = "inmb",
  x_var = "time_horizon_years",
  ref_strategy = "SOC",
  facet_var = "perspective",
  group_var = "alternative",
  facet_labs = c(
    "0" = "(a) Health Care Sector Perspective", 
    "1" = "(b) Societal Perspective"
  ),
  group_labs = c(
    "IXCHIQ" = "Live-attenuated vaccine", 
    "VIMKUNYA" = "Recombinant vaccine"
  ),
  y_lab = "INMB, 2024 US$ million",
  x_lab = "Time horizon, years",
  title = NULL,
  legend_position = "bottom",
  file_name = "time_horizon_sensi_cea_plot.png", 
  file_path = figures_dir,
  fig_scale = 1, 
  fig_width = 7, 
  fig_height = 4, 
  units = c("in", "cm", "mm", "px")[1],
  dpi = 300,
  limitsize = TRUE, 
  bg = "white",
  journal = "bmj"
){

  ## prepare plot data
  plot_data <- data %>%
    dplyr::select(
      !!sym(facet_var),
      !!sym(x_var),
      !!sym(group_var),
      !!sym(outcome_col)
    ) %>%
    dplyr::filter(
      !(alternative == ref_strategy),
      !is.na(.data[[outcome_col]]),
      !is.nan(.data[[outcome_col]])
    ) %>%
    distinct() %>%
    mutate(
      !!sym(facet_var) := factor(!!sym(facet_var), levels = names(facet_labs)),
      !!sym(group_var) := factor(!!sym(group_var), levels = names(group_labs), labels = group_labs)    
    )
  
  ## plot
  p <- ggplot(
    plot_data,
    aes(
      x = .data[[x_var]],
      y = .data[[outcome_col]] * 1e-6,
      group = .data[[group_var]],
      color = .data[[group_var]],
    )
  ) +
  geom_line() + 
  facet_grid(
    cols = vars(!!sym(facet_var)),
    labeller = labeller(
      !!facet_var := facet_labs
    )
  ) + 
  get_journal_color_scale(journal) + 
  labs(
    x = x_lab,
    y = y_lab,
    title = title
  ) +
  scale_x_continuous(
    limits = c(1, 20),
    breaks = c(1, seq(0, 20, 5)),
    expand = c(0.01, 0.01)
  ) +
  scale_y_continuous(
    limits = c(0, 6),
    breaks = seq(0, 6, 1),
    labels = scales::number_format(accuracy = 0.01),
    expand = c(0.01, 0.01)
  ) +
  theme_minimal(base_size = 11) +
  theme(
    panel.grid = element_blank(),
    panel.spacing.y = unit(6, "mm"),
    axis.line = element_line(color = "black", linewidth = 0.5),
    axis.ticks = element_line(color = "black"),
    strip.text = ggtext::element_markdown(face = "bold", size = 11),
    plot.title = element_text(
      face = "bold",
      size = 11,
      hjust = 0.5
    ),
    legend.position = legend_position,
    legend.title = element_blank(), 
    legend.text = element_text(size = 10),
    legend.key.height = unit(2, "mm")
  )

  ggsave(
    filename = file_name, 
    plot = p, 
    path = file_path,
    scale = fig_scale, 
    width = fig_width, 
    height = fig_height, 
    units = units, 
    dpi = dpi, 
    limitsize = limitsize, 
    bg = bg
  )

  ## output
  return(p) 

}

##************************************************************
## heatmap plot function
##************************************************************
sensi_time_eta_zeta_heatmap_func <- function(
  data = sensi_time_eta_zeta_cea_data$cea_out_usd,
  ref_strategy = "SOC",
  perspective_val = 0,
  alt_var = "alternative",
  perspective_var = "perspective",
  x_var = "time_horizon_years",
  y_var = "eta",
  z_var = "inmb",
  z_scaling_factor = 1e-6,
  x_facet_var = "alternative",
  y_facet_var = "zeta",
  ##
  x_facet_labs = c(
    "IXCHIQ" = "Live-attenuated vaccine", 
    "VIMKUNYA" = "Recombinant vaccine"
  ),
  y_facet_labs = c(
    "0.7" = "&zeta; = 0.70",
    "0.8" = "&zeta; = 0.80", 
    "0.9" = "&zeta; = 0.90", 
    "1" = "&zeta; = 1.00"
  ),
  x_lab = "Time horizon, years",
  y_lab = expression("Vaccine price discount, " *eta),
  scale_lab = "INMB (2024 US$ million): ",
  title = "(a) Health Care Sector Perspective",
  midpoint_val = 2, 
  legend_position = "bottom",
  file_name = "time_horizon_sensi_heatmap_plot.png", 
  file_path = figures_dir,
  fig_scale = 1, 
  fig_width = 7, 
  fig_height = 4, 
  units = c("in", "cm", "mm", "px")[1],
  dpi = 300,
  limitsize = TRUE, 
  bg = "white",
  journal = "bmj"
){

  if(is.null(x_facet_var)){
    x_facet_var <- alt_var
  } else {
    x_facet_var <- x_facet_var
  }

  ## prepare plot data
  plot_data <- data %>%
    dplyr::select(
      perspective,
      !!sym(alt_var),
      !!sym(perspective_var),
      !!sym(y_facet_var),
      !!sym(x_var),
      !!sym(y_var),
      !!sym(z_var)
    ) %>%
    dplyr::filter(
      !(!!sym(alt_var) == ref_strategy),
      !is.na(.data[[z_var]]),
      !is.nan(.data[[z_var]]) ,
      !!sym(perspective_var) == perspective_val
    ) %>%
    distinct() %>%
    mutate(
      !!sym(z_var) := z_scaling_factor * !!sym(z_var),
      !!sym(x_facet_var) := factor(!!sym(x_facet_var), levels = names(x_facet_labs)),    
      !!sym(y_facet_var) := factor(!!sym(y_facet_var), levels = names(y_facet_labs))
    )
  
  ## plot
  p <- ggplot(
    plot_data,
    aes(
      x = .data[[x_var]],
      y = .data[[y_var]],
      fill = .data[[z_var]]
    )
  ) +
  geom_tile(color = "white", linewidth = 0) + 
  facet_grid(
    rows = vars(!!sym(x_facet_var)),
    cols = vars(!!sym(y_facet_var)),
    labeller = labeller(
      !!x_facet_var := x_facet_labs,
      !!y_facet_var := y_facet_labs
    )
  ) + 
  scale_fill_gradient2(
    name = scale_lab,
    high = "#CD202C",
    mid = "#B4CFEE",
    low = "#2A6EBB",
    midpoint = midpoint_val,
    labels = scales::comma
  ) +
  labs(
    x = x_lab,
    y = y_lab,
    title = title
  ) +
  # scale_x_continuous(
  #   limits = c(0.40, 0.75),
  #   breaks = seq(0.40, 0.80, 0.1),
  #   expand = c(0, 0)
  # ) +
  # scale_y_continuous(
  #   limits = c(0.0, 1.20),
  #   breaks = seq(0.0, 1.20, 0.20),
  #   expand = c(0, 0)
  # ) +
  theme_minimal(base_size = 11) +
  theme(
    panel.grid = element_blank(),
    panel.spacing.y = unit(6, "mm"),
    axis.line = element_line(color = "black", linewidth = 0.5),
    axis.ticks = element_line(color = "black"),
    # strip.text = element_text(face = "bold", size = 11),
    strip.text = ggtext::element_markdown(face = "bold", size = 10),
    plot.title = element_text(
      face = "bold",
      size = 11,
      hjust = 0.5
    ),
    legend.position = legend_position,
    legend.title = element_text(size = 10), 
    legend.text = element_text(size = 8),
    legend.key.height = unit(2, "mm")
  )

  ggsave(
    filename = file_name, 
    plot = p, 
    path = file_path,
    scale = fig_scale, 
    width = fig_width, 
    height = fig_height, 
    units = units, 
    dpi = dpi, 
    limitsize = limitsize, 
    bg = bg
  )

  ## output
  return(p) 
}


## Combined heatmap plot
combo_heatmap_func <- function(
  data = sensi_time_eta_zeta_cea_data$cea_out_usd,
  ref_strategy = "SOC",
  perspectives = perspectives,
  alt_var = "alternative",
  perspective_var = "perspective",
  x_var = "time_horizon_years",
  y_var = "eta",
  z_var = "inmb",
  z_scaling_factor = 1e-6,
  x_facet_var = "alternative",
  y_facet_var = "zeta",
  perspective_labs = c(
    "0" = "(a) Health Care Sector Perspective", 
    "1" = "(b) Societal Perspective"
  ),
  x_facet_labs = c(
    "IXCHIQ" = "Live-attenuated vaccine", 
    "VIMKUNYA" = "Recombinant vaccine"
  ),
  y_facet_labs = c(
    "0.7" = "&zeta; = 0.70",
    "0.8" = "&zeta; = 0.80", 
    "0.9" = "&zeta; = 0.90", 
    "1" = "&zeta; = 1.00"
  ),
  x_lab = "Time horizon, years",
  y_lab = expression("Vaccine price discount, " *eta),
  scale_lab = "INMB (2024 US$ million): ",
  midpoint_val = 2.5, 
  legend_position = "bottom",
  file_name = "time_horizon_sensi_heatmap_plot.png", 
  file_path = figures_dir,
  fig_scale = 1, 
  fig_width = 7, 
  fig_height = 4, 
  combo_fig_width = 11, 
  combo_fig_height = 9, 
  units = c("in", "cm", "mm", "px")[1],
  dpi = 300,
  limitsize = TRUE, 
  bg = "white",
  journal = "bmj"  
){

  p_list <- lapply(
    X = perspectives,
    FUN = function(i){
      ## perspective-specific variables
      if(i == 0) {
        iter_x_lab <- NULL
        legend_pos <- "none"
        perspective_tag <- "hc"
      } else {
        legend_pos <- "bottom"
        perspective_tag <- "soc"
      }
      
      ## plot objects for each perspective
      sensi_time_eta_zeta_heatmap_func(
        data = data,
        ref_strategy = ref_strategy,
        perspective_val = i,
        alt_var = alt_var,
        perspective_var = perspective_var,
        x_var = x_var,
        y_var = y_var,
        z_var = z_var,
        z_scaling_factor = z_scaling_factor,
        x_facet_var = x_facet_var,
        y_facet_var = y_facet_var,
        ##
        x_facet_labs = x_facet_labs,
        y_facet_labs = y_facet_labs,
        x_lab = x_lab,
        y_lab = y_lab,
        scale_lab = scale_lab,
        title = perspective_labs[as.character(i)],
        midpoint_val = midpoint_val, 
        legend_position = legend_pos,
        file_name = paste0("time_horizon_sensi_heatmap_plot_", perspective_tag, ".png"), 
        file_path = file_path,
        fig_scale = fig_scale, 
        fig_width = fig_width, 
        fig_height = fig_height, 
        units = units,
        dpi = dpi,
        limitsize = limitsize, 
        bg = bg,
        journal = journal
      )
    }
  )

  names(p_list) <- c("hc", "soc")  

  ## stacked vertically
  combined_plot <- p_list$hc / p_list$soc +
    plot_annotation(
      theme = theme(
        plot.title = element_text(size = 12, face = "bold")
      )
    )

  ## display
  print(combined_plot)

  ## save the plot
  ggsave(
    filename = file_name, 
    plot = combined_plot, 
    path = figures_dir,
    scale = fig_scale, 
    width = combo_fig_width, 
    height = combo_fig_height, 
    units = units, 
    dpi = dpi, 
    limitsize = limitsize, 
    bg = bg
  )
    
  invisible(combined_plot)

}



##*****************************************************************************************************
## END OF SCRIPT
##***************************************************************************************************** 
