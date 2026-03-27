##****************************************************************************************
## ANALYSES
##****************************************************************************************

##************************************************************
## Base case CEA
##************************************************************
base_cea_data_file <- file.path(data_dir, "base_cea_data.rds")
if(isTRUE(use_saved_base_cea_data) & file.exists(base_cea_data_file)){
  base_case_cea <- readRDS(base_cea_data_file)
} else {
  base_case_cea_data <- base_case_cea_func(
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
    n_cycles = time_horizon_set[1] * 52,
    perspectives = perspectives,
    le_death_data = le_death_data, 
    fx_adjust = fx_adjust, 
    cpi_adjust = cpi_adjust,
    bound_val = "base",
    par_name = "base",
    base_par_val = NA,
    bound_par_val = NA
  )
  saveRDS(object = base_case_cea_data, file = base_cea_data_file)
}

## formatting
scale_vals <- list(Costs = 1, QALYs = 1, NMB = 1e6, Inc_Costs = 1, Inc_QALYs = 1, INMB = 1e6, ICER = 1)
digits_vals <- list(Costs = 0, QALYs = 0, NMB = 0, Inc_Costs = 0, Inc_QALYs = 0, INMB = 0, ICER = 0)
big.mark_vals <- list(Costs = NULL, QALYs = NULL, NMB = NULL, INMB = NULL, ICER = NULL)
decimal.mark_vals <- list(Costs = NULL, QALYs = NULL, NMB = NULL, Inc_Costs = NULL, Inc_QALYs = NULL, INMB = NULL, ICER = NULL)

## Main text, Table 2
if(isTRUE(create_base_cea_table_2)){
  cea_table_func(
    df = base_case_cea_data$cea_out_usd,
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
  ) 
}

## Supplemental Appendix, Table S1
if(isTRUE(create_base_cea_table_S1)){
  cea_table_func(
    df = base_case_cea_data$cea_out_brl,
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
    tab_caption = "Table S1. Base Case Cost-effectiveness Analysis Results (Brazilian Real)",
    tab_dir = tables_dir,
    fname = "basecase_cea_results_table_brazilian.docx"
  )
}

##************************************************************
## OWSA
##************************************************************
##----------------------------------------
## data
##----------------------------------------
owsa_data_file <- file.path(data_dir, "owsa_data.rds")
if(isTRUE(use_saved_owsa_data) & file.exists(owsa_data_file)){
  owsa_data <- readRDS(owsa_data_file)
} else {
  owsa_data <- run_owsa(
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
    n_cycles = time_horizon * 52,
    perspectives = c(0, 1),
    le_death_data = le_death_data, 
    fx_adjust = fx_adjust, 
    cpi_adjust = cpi_adjust
  )     
  saveRDS(object = owsa_data, file = owsa_data_file)
}

##----------------------------------------
## produce tornado plots, by perspective
##----------------------------------------
if(isTRUE(create_tornado_plots)){
  
  owsa_params <- names(params_label)
  owsa_params_labels <- params_label
  
  facet_strip_text <- c(
    "IXCHIQ" = "(a) Live-attenuated vaccine vs No vaccination",
    "VIMKUNYA" = "(b) Recombinant vaccine vs No vaccination"
  )
  
  ## Main text, Figure 2
  ## Produce tornado plots -- healthcare sector
  owsa_tornado_plot_func(
    alternatives = scenarios[2:3],
    perspective_val = c(0,1)[1],
    curr_val = c("brl", "usd")[2],
    out_val = c("icer", "inmb", "inhb")[2],
    ##
    owsa_output_data = owsa_data,
    params_include = owsa_params,
    params_label = owsa_params_labels,
    diff_tolerance = 0.01,
    bar_width = 0.75,
    ##
    y_axis_rescale_by = 1000000,
    y_min = -2000000,
    y_max = 12000000,
    y_by = 1000000,
    y_exp_min = 0.01,
    y_exp_max = 0.01,
    x_exp_min = 0.025,
    x_exp_max = 0.03,
    ##
    ylab = "INMB, 2024 million US$ per 1,000 population",
    xlab = "Model input parameters",
    title = NULL,
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
    ##
    file_name = "tornado_owsa_healthcare.png",
    file_path = figures_dir,
    fig_scale = 1,
    fig_width = 17,
    fig_height = 7,
    units = c("in", "cm", "mm", "px")[1],
    dpi = 300,
    limitsize = TRUE,
    bg = "white",
    ##
    journal = "bmj"
  )
  
  ## Main text, Figure 2
  ## Produce tornado plots -- societal perspective
  owsa_tornado_plot_func(
    alternatives = scenarios[2:3],
    perspective_val = c(0,1)[2],
    curr_val = c("brl", "usd")[2],
    out_val = c("icer", "inmb", "inhb")[2],
    ##
    owsa_output_data = owsa_data,
    params_include = owsa_params,
    params_label = owsa_params_labels,
    diff_tolerance = 0.01,
    bar_width = 0.75,
    ##
    y_axis_rescale_by = 1000000,
    y_min = -2000000,
    y_max = 12000000,
    y_by = 1000000,
    y_exp_min = 0.01,
    y_exp_max = 0.01,
    x_exp_min = 0.025,
    x_exp_max = 0.03,
    ##
    ylab = "INMB, 2024 million US$ per 1,000 population",
    xlab = "Model input parameters",
    title = NULL,
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
    ##
    file_name = "tornado_owsa_societal.png",
    file_path = figures_dir,
    fig_scale = 1,
    fig_width = 17,
    fig_height = 7,
    units = c("in", "cm", "mm", "px")[1],
    dpi = 300,
    limitsize = TRUE,
    bg = "white",
    ##
    journal = "bmj"
  )
}

##************************************************************
## Sensitivity Analyses
##************************************************************



##************************************************************
## Probabilistic Sensitivity Analyses (PSA)
##************************************************************




##*****************************************************************************************************
## END OF SCRIPT
##***************************************************************************************************** 
  