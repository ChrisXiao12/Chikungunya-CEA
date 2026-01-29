##*****************************************************************************************************
## SENSITIVITY ANALYSES
##*****************************************************************************************************   
sensi_cea_data_file <- file.path(data_dir, "sensi_cea_data.rds")
if(isTRUE(use_saved_sensi_cea_data) & file.exists(sensi_cea_data_file)){
  sensi_cea_data <- readRDS(file = sensi_cea_data_file)
} else {
  sensi_cea_data <- sensi_func(
    xi_vals = xi_set,
    zeta_vals = NULL,
    phi_vals = NULL,
    psi_stated_vals = NULL,
    pi_vals = NULL,
    chi_vals = NULL,
    eta_vals = NULL,
    omicron_vals = NULL,
    time_horizons = time_horizon_set[1],
    wtp = NULL,
    scenarios = scenarios, 
    ref_scenario = scenarios[1],
    params = params_value, 
    init_states = ini_pop_vec,
    inc_names = trans_names,
    perspectives = perspectives,
    le_death_data = le_death_data, 
    fx_adjust = fx_adjust, 
    cpi_adjust = cpi_adjust
  )
  saveRDS(object = sensi_cea_data, file = sensi_cea_data_file)   
}


## create table
behavior_table_func(
  data = sensi_cea_data,
  currency = "usd",
  tables_dir = tables_dir,
  table_fname = "sensi_table.docx"  
)


##************************************************************
## Sensitivity Analyses -- Differential Vaccine Uptake and realization factor for vaccine uptake
##************************************************************
phi_pi_sensi_cea_data_file <- file.path(data_dir, "xi_pi_sensi_cea_data.rds")
if(isTRUE(use_saved_phi_pi_sensi_cea_data) & file.exists(phi_pi_sensi_cea_data_file)){
  phi_pi_sensi_cea_data <- readRDS(file = phi_pi_sensi_cea_data_file)
} else {
  phi_pi_sensi_cea_data <- sensi_func(
    xi_vals = NULL,
    zeta_vals = NULL,
    phi_vals = phi_set,
    psi_stated_vals = NULL,
    pi_vals = pi_set,
    chi_vals = NULL,
    eta_vals = NULL,
    omicron_vals = omicron_set,
    time_horizons = time_horizon_set[1],
    wtp = NULL,
    scenarios = scenarios, 
    ref_scenario = scenarios[1],
    params = params_value, 
    init_states = ini_pop_vec,
    inc_names = trans_names,
    perspectives = perspectives,
    le_death_data = le_death_data, 
    fx_adjust = fx_adjust, 
    cpi_adjust = cpi_adjust
  )
  saveRDS(object = phi_pi_sensi_cea_data, file = phi_pi_sensi_cea_data_file)   
}

## appendix figure
combo_phi_pi_omicron_heatmap_func(
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
  y_lab = expression("Relative vaccine uptake parameter, " * pi[v]),
  legend_position = "bottom",
  file_name = "phi_pi_sensi_heatmap.png", 
  file_path = figures_dir,
  fig_scale = 1, 
  fig_width = 11, 
  fig_height = 6,
  combo_fig_scale = 0.8, 
  combo_fig_width = 11, 
  combo_fig_height = 12, 
  units = c("in", "cm", "mm", "px")[1],
  dpi = 300,
  limitsize = TRUE, 
  bg = "white"
)

##************************************************************
## Sensitivity Analyses on zeta (Immunogenicity atenuation factor)
##************************************************************
zeta_sensi_cea_data_file <- file.path(data_dir, "zeta_sensi_cea_data.rds")
if(isTRUE(use_saved_zeta_sensi_cea_data) & file.exists(zeta_sensi_cea_data_file)){
  zeta_sensi_cea_data <- readRDS(file = zeta_sensi_cea_data_file)
} else {
  zeta_sensi_cea_data <- sensi_func(
    xi_vals = NULL,
    zeta_vals = zeta_set,
    phi_vals = NULL,
    psi_stated_vals = NULL,
    pi_vals = NULL,
    chi_vals = NULL,
    eta_vals = NULL,
    omicron_vals = NULL,
    time_horizons = time_horizon_set[1],
    wtp = NULL,
    scenarios = scenarios, 
    ref_scenario = scenarios[1],
    params = params_value, 
    init_states = ini_pop_vec,
    inc_names = trans_names,
    perspectives = perspectives,
    le_death_data = le_death_data, 
    fx_adjust = fx_adjust, 
    cpi_adjust = cpi_adjust
  )
  saveRDS(object = zeta_sensi_cea_data, file = zeta_sensi_cea_data_file)   
}

## plot
plot_zeta_outcome_heatmap(
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
    "0" = "(a) Health Care Sector Perspective", 
    "1" = "(b) Societal Perspective"
  ),
  y_lab = "INMB, 2024 US$ million",
  x_lab = expression("Immunogenicity attenuation factor, " * zeta),
  title = NULL,
  legend_position = "bottom",
  file_name = "zeta_sensi_cea_plot.png", 
  file_path = figures_dir,
  fig_scale = 1, 
  fig_width = 7, 
  fig_height = 3, 
  units = c("in", "cm", "mm", "px")[1],
  dpi = 300,
  limitsize = TRUE, 
  bg = "white",
  journal = "bmj"
) 


##************************************************************
## Sensitivity Analyses -- Time horizon only
##************************************************************
sensi_time_horizon_cea_data_file <- file.path(data_dir, "sensi_time_horizon_cea_data.rds")
if(isTRUE(use_saved_sensi_time_horizon_cea_data) & file.exists(sensi_time_horizon_cea_data_file)){
  sensi_time_horizon_cea_data <- readRDS(file = sensi_time_horizon_cea_data_file)
} else {
  sensi_time_horizon_cea_data <- sensi_func(
    xi_vals = NULL,
    zeta_vals = NULL,
    phi_vals = NULL,
    psi_stated_vals = NULL,
    pi_vals = NULL,
    chi_vals = NULL,
    eta_vals = NULL,
    omicron_vals = NULL,
    time_horizons = time_horizon_set,
    wtp = NULL,
    scenarios = scenarios, 
    ref_scenario = scenarios[1],
    params = params_value, 
    init_states = ini_pop_vec,
    inc_names = trans_names,
    perspectives = perspectives,
    le_death_data = le_death_data, 
    fx_adjust = fx_adjust, 
    cpi_adjust = cpi_adjust
  )
  saveRDS(object = sensi_time_horizon_cea_data, file = sensi_time_horizon_cea_data_file)   
}


## plot
one_way_time_horizon_plot_func(
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
)

##************************************************************
## Sensitivity Analyses -- Time horizon, vaccine price 
## discount (eta), and zeta
##************************************************************

##**************************************
## Data
##**************************************
sensi_time_eta_zeta_cea_data_file <- file.path(data_dir, "sensi_time_eta_zeta_cea_data.rds")
if(isTRUE(use_saved_sensi_time_eta_zeta_cea_data) & file.exists(sensi_time_eta_zeta_cea_data_file)){
  sensi_time_eta_zeta_cea_data <- readRDS(file = sensi_time_eta_zeta_cea_data_file)
} else {
  sensi_time_eta_zeta_cea_data <- sensi_func(
    xi_vals = NULL,
    zeta_vals = zeta_set,
    phi_vals = NULL,
    psi_stated_vals = NULL,
    pi_vals = NULL,
    chi_vals = NULL,
    eta_vals = eta_set,
    omicron_vals = NULL,
    time_horizons = time_horizon_set,
    wtp = NULL,
    scenarios = scenarios, 
    ref_scenario = scenarios[1],
    params = params_value, 
    init_states = ini_pop_vec,
    inc_names = trans_names,
    perspectives = perspectives,
    le_death_data = le_death_data, 
    fx_adjust = fx_adjust, 
    cpi_adjust = cpi_adjust
  )
  saveRDS(object = sensi_time_eta_zeta_cea_data, file = sensi_time_eta_zeta_cea_data_file)   
}

##**************************************
## plot
##**************************************
combo_heatmap_func(
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
  combo_fig_height = 10, 
  units = c("in", "cm", "mm", "px")[1],
  dpi = 300,
  limitsize = TRUE, 
  bg = "white",
  journal = "bmj"  
)

##*****************************************************************************************************
## END OF SCRIPT
##***************************************************************************************************** 
  

