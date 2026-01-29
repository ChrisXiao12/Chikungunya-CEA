##*****************************************************************************************************
## INPUT PARAMETERS
##*****************************************************************************************************

## parameters file
params_fname <- "parameters.xlsx"

## read parameters file
params_data <- read_excel(
  path = file.path(data_dir, params_fname), 
  sheet = "parameters"
)

## categories
params_categories <- as.vector(unique(na.omit(params_data$Category)))

## remove calculations for now
base_params <- params_data %>% 
  filter(!(Category %in% c(NA, "Calculations")) & !(Parameter %in% c(NA)))

## assumed distribution
params_dist <- base_params$Distribution
names(params_dist) <- base_params$Parameter

## currency
params_currency <- base_params$Currency
names(params_currency) <- base_params$Parameter

## year
params_year <- base_params$Year
names(params_year) <- base_params$Parameter

## label
params_label <- base_params$Definition
names(params_label) <- base_params$Parameter

## Category
params_category <- base_params$Category
names(params_category) <- base_params$Parameter

## values
params_value <- as.numeric(base_params$Value)
names(params_value) <- base_params$Parameter

## lower bound
params_lower <- as.numeric(base_params$Lower)
names(params_lower) <- base_params$Parameter

## upper bound
params_upper <- as.numeric(base_params$Upper)
names(params_upper) <- base_params$Parameter

##*************************************************************
## USD to BRL conversion
##*************************************************************
## fx coversion factors
usd_cost_param <- base_params %>% 
  dplyr::filter((Currency == "USD") & (Category == "Costs"))%>%
  dplyr::select(Parameter, Year)

brl_to_usd_adj_factor <- base_params %>% 
  dplyr::filter((Category == "CurrencyAdjustments")) %>%
  mutate(brl_to_usd = as.numeric(Value)) %>%
  dplyr::select(Year, brl_to_usd) 

brl_cost_param <- merge(
  usd_cost_param,
  brl_to_usd_adj_factor,
  by = "Year",
  all.x = TRUE
) %>%
mutate(Value = coalesce(brl_to_usd, 1))

fx_adjust <- brl_cost_param$Value
names(fx_adjust) <- brl_cost_param$Parameter

rm(brl_to_usd_adj_factor)
rm(usd_cost_param)

##*************************************************************
## Inflation adjustment
##*************************************************************
## cpi
cpi_adj_data <- file.path(data_dir, "cpi_adjust.rds")
if(isTRUE(use_saved_cpi_adj) & file.exists(cpi_adj_data) ){
  cpi_adjust <- readRDS(file = cpi_adj_data)
} else {
  cost_params <- base_params %>% 
    dplyr::filter((Category == "Costs") & !(Year == "NA"))%>%
    dplyr::select(Parameter, Year)
  
  cpi_params <- base_params %>% 
    dplyr::filter(Category == "CPIAdjustments") %>%
    dplyr::select(Year, Value) %>%
    mutate(
      Value = as.numeric(Value),
      RefValue = Value[Year == 2024],
      inflation = RefValue / Value
    ) %>%
    dplyr::select(Year, inflation)
  
  cost_params <- merge(cost_params, cpi_params, "Year")  
  cpi_adjust <- cost_params$inflation
  names(cpi_adjust) <- cost_params$Parameter
  
  saveRDS(object = cpi_adjust, file = cpi_adj_data)
  rm(cost_params)
}

##*************************************************************
## Life expectancy and death probability
##*************************************************************
le_death_data <- read_excel(
  path = file.path(data_dir, params_fname), 
  sheet = "life_tables",
  range = "A7:G97"
) 

##*************************************************************
## Vaccine Seroprotection Data
## Sources: 
## IXCHIQ: 18-64 years; Figure 2 in McMahon, R., Toepfer, S., Sattler, N., et al. Antibody persistence and safety of a live-attenuated chikungunya virus vaccine up to 2 years after single-dose administration in adults in the USA: a single-arm, multicentre, phase 3b study. Lancet Infectious Diseases. 2024;24(12):1383–1392.	
## VIMKUNYA: Table 2 in Richardson JS, Anderson DM, Mendy J, Tindale LC, Muhammad S, Loreth T, Tredo SR, Warfield KL, Ramanathan R, Caso JT, Jenkins VA, Ajiboye P, Bedell L; EBSI-CV-317-004 Study Group. Chikungunya virus virus-like particle vaccine safety and immunogenicity in adolescents and adults in the USA: a phase 3, randomised, double-blind, placebo-controlled trial. Lancet. 2025 Apr 19;405(10487):1343-1352. doi: 10.1016/S0140-6736(25)00345-9. Epub 2025 Mar 27. PMID: 40158526.
##*************************************************************
seroresponse_params_data <- file.path(data_dir, "seroresponse_params.rds")
if(isTRUE(use_saved_seroresponse) & file.exists(seroresponse_params_data) ){
  seroresponse_params <- readRDS(file = seroresponse_params_data)
} else {
  vac_seroresp_data <- read_excel(
    path = file.path(data_dir, params_fname), 
    sheet = "seroresponse",
    range = "A1:H10"
  ) 
  vac_seroresp_data <- vac_seroresp_data %>% filter(days > 8)

  ## vaccine seroprotection parameters
  seroresponse_params <- vac_seroresponse_func(
    data = vac_seroresp_data,
    fit_model = c("expo", "poly")[1],
    facet_labels = c(
      "IXCHIQ" = "(a) Live-attenuated vaccine",
      "VIMKUNYA" = "(b) Recombinant vaccine"
    ),
    create_coeff_table = create_coeff_table,
    print_coef_plot = print_coef_plot,
    save_coeffs_plot = save_coeffs_plot,
    figures_dir = figures_dir,
    tables_dir = tables_dir
  )
  
  saveRDS(object = seroresponse_params, file = seroresponse_params_data)
  rm(vac_seroresp_data)
}

##************************************************************
## Update the parameters with the estimated coefficients
##************************************************************
## values
omega_value <- seroresponse_params$coeffs_df %>% 
  dplyr::filter(coefficient == "b") %>% 
  dplyr::select(median) %>% 
  unlist() %>% 
  as.vector()
params_value[c("theta_IXCHIQ_days", "theta_VIMKUNYA_days")] <- as.numeric(omega_value)

## lower bound
omega_lower <- seroresponse_params$coeffs_df %>% 
  dplyr::filter(coefficient == "b") %>% 
  dplyr::select(lower) %>% 
  unlist() %>% 
  as.vector()
params_lower[c("theta_IXCHIQ_days", "theta_VIMKUNYA_days")] <- as.numeric(omega_lower)

## upper bound
omega_upper <- seroresponse_params$coeffs_df %>% 
  dplyr::filter(coefficient == "b") %>% 
  dplyr::select(upper) %>% 
  unlist() %>% 
  as.vector()
params_upper[c("theta_IXCHIQ_days", "theta_VIMKUNYA_days")] <- as.numeric(omega_upper)

## free memory
rm(omega_value, omega_lower, omega_upper)

##*********************************************************************
## Initial population
##*********************************************************************
ini_pop_vec <- inipop_func(
  states = state_names, 
  params = params_value
)

##*****************************************************************************************************
## END OF SCRIPT
##***************************************************************************************************** 