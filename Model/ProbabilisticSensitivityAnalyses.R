##*****************************************************************************************************
## PROBABILISTIC SENSITIVITY ANALYSIS (PSA)
##***************************************************************************************************** 

##************************************************************
## sampled parameters
##************************************************************
psa_params_data <- file.path(data_dir, "psa_params.rds")
if(isTRUE(use_saved_psa_params) & file.exists(psa_params_data)){
  psa_params <- readRDS(file = psa_params_data)
}else{
  psa_params <- draw_psa_params(
    dist = params_dist,
    mean = params_value,
    lower = params_lower,
    upper = params_upper,
    n_sim = n_psa_draws 
  )  
  saveRDS(object = psa_params, file = psa_params_data) 
}

##************************************************************
## run psa cea
##************************************************************
psa_data_file <- file.path(data_dir, "psa_data.rds")
if(isTRUE(use_saved_psa_data) & file.exists(psa_data_file)){
  psa_data  <- readRDS(file = psa_data_file)
} else {
  # plan(multisession, workers = parallel::detectCores() - 1)
  old_plan <- plan()
  on.exit(plan(old_plan), add = TRUE)
  plan(multisession, workers = n_workers)

  psa_data <- cea_psa_func(
    gdp_factors = gdp_multipliers,
    scenarios = scenarios,
    ref_scenario = scenarios[1],
    psa_params = psa_params,
    xi_factor = NULL,
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
  )
  saveRDS(object = psa_data, file = psa_data_file) 
}

##************************************************************
## plot PSA results
##************************************************************
if(isTRUE(create_psa_plots)){
  psa_data_file <- file.path(data_dir, "psa_data.rds")
  ceac_plot_file <- file.path(figures_dir, "ceac_plot.png")
  if(!exists("psa_data")){
    if(file.exists(psa_data_file)){
      psa_data <- readRDS(file = psa_data_file)
    } else { 
      old_plan <- plan()
      on.exit(plan(old_plan), add = TRUE)
      plan(multisession, workers = n_workers)

      psa_data <- cea_psa_func(
        gdp_factors = gdp_multipliers,
        scenarios = scenarios,
        ref_scenario = scenarios[1],
        psa_params = psa_params,
        xi_factor = NULL,
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
      )
      saveRDS(object = psa_data, file = psa_data_file)     
    }
  }

  ## plot data
  psa_plot_data <- ceac_data_func(
    currencies = c("brl", "usd"),
    data = psa_data 
  )

  ## plot
  psa_plot_func(
    data = psa_plot_data$long %>% dplyr::filter(currency == "usd"),
    x = "wtp", y = "prob", 
    color = "alternative", 
    fill = "alternative",
    ymin = "lower", 
    ymax = "upper",
    ##
    alt_labs = c(
      "IXCHIQ" = "Live-attenuated vaccine",
      "VIMKUNYA" = "Recombinant vaccine"
    ),  
    perspective_labs = c(
      "0" = "(a) Health Care Sector Perspective",
      "1" = "(b) Societal Perspective"
    ),
    x_lab = "WTP threshold, 2024 US$'000",
    y_lab = "Probability of being cost-effective",
    ##
    color_lab = NULL,
    fill_lab = NULL,
    ##
    x_axis_rescale_by = 1000,
    x_min = 5000,
    x_max = 42000,
    x_by = 5000,
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
    save_plot = TRUE,
    file_name = "ceac.png", 
    file_path = figures_dir,
    fig_scale = 1, 
    fig_width = 9, 
    fig_height = 5, 
    units = c("in", "cm", "mm", "px")[1], 
    dpi = 300, 
    limitsize = TRUE,
    bg = "white"
  )  

}

##************************************************************
## bootstrap confidence intervals
##************************************************************
cea_psa_boot_file <- file.path(data_dir, "cea_psa_boot_data.rds")
if(isTRUE(use_saved_psa_boot_data) & file.exists(cea_psa_boot_file)){
  cea_psa_boot_data  <- readRDS(file = cea_psa_boot_file)
} else {
  cea_psa_boot_data <- bayes_mean_ci_func(
    currencies = c("brl", "usd"),
    data = psa_data,
    params = params_value,
    scenarios = scenarios
  )   
  saveRDS(object = cea_psa_boot_data, file = cea_psa_boot_file)   
}

colnames(cea_psa_boot_data)
head(cea_psa_boot_data)

tail(cea_psa_boot_data)

##******************************************
## create a table of CEA results from the multivariate analysis
##******************************************
scale_vals <- list(Costs = 1, QALYs = 1, NMB = 1e6, INMB = 1)
digits_vals <- list(Costs = 0, QALYs = 0, NMB = 0, INMB = 0)
big.mark_vals <- list(Costs = NULL, QALYs = NULL, NMB = NULL, INMB = NULL)
decimal.mark_vals <- list(Costs = NULL, QALYs = NULL, NMB = NULL, INMB = NULL)

mult_cea_table_func(
  psa_boot_data = cea_psa_boot_data,
  currency_val = c("brl", "usd")[2],
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
  NMB_lab = "NMB (US$ M)\n Est. (95% CI)",
  INMB_lab = "INMB (US$)\n Est. (95% CI)",
  tables_dir = tables_dir, 
  table_fname = "multivariate_cea_results_table.docx"
)


## appendix table
mult_cea_table_func(
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
  table_caption = "Table S9. Cost-effectiveness Analysis Results from Multivariate Sensitivity Analyses",
  Costs_lab = "Costs (R$)\n Est. (95% CI)",
  QALYs_lab = "QALYs\n Est. (95% CI)",
  NMB_lab = "NMB (R$ M)\n Est. (95% CI)",
  INMB_lab = "INMB (R$)\n Est. (95% CI)",
  tables_dir = tables_dir, 
  table_fname = "multivariate_cea_results_table_appendix.docx"
)

##*****************************************************************************************************
## END OF SCRIPT
##***************************************************************************************************** 