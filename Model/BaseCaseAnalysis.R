##*****************************************************************************************************
## BASE CASE ANALYSIS
##*****************************************************************************************************
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
    time_horizon_years = time_horizon_set[1],
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

##*****************************************************************************************************
## END OF SCRIPT
##***************************************************************************************************** 

