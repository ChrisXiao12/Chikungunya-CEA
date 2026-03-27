# 
# 
# 
# #   flat_results <- do.call(c, results)
# #   
# #   owsa_df <- do.call(
# #     rbind, 
# #     Map(
# #       function(df, nm) {
# #         cbind(
# #           currency = nm, 
# #           df
# #         )
# #       }, 
# #       flat_results, 
# #       names(flat_results)
# #     )
# #   )
# # 
# #   rownames(owsa_df) <- NULL
# 
#   ## output value
#   # val <- list(
#   #   owsa_table = owsa_df,
#   #   base_result = base_result
#   # )
# 
# # ## grid of currencies and outcomes
# # currencies_outcomes <- expand.grid(
# #   currency = currencies,
# #   outcome  = outcomes
# # )
# # 
# # ## pull the outcome's value
# # base_val_outcome <- lapply(
# #   X = 1:nrow(currency_outcome),
# #   FUN = extract_outcome,
# #   currencies_outcomes = currencies_outcomes, 
# #   result = base_result 
# # )
# # names(base_val_outcome) <- currencies
# 
# 
# 
# 
#   
#   low_val_outcome <- lapply(
#     X = 1:nrow(currencies_outcomes),
#     FUN = extract_outcome,
#     currencies_outcomes = currencies_outcomes,
#     result = low_result
#   )
#   names(low_val_outcome) <- currencies_outc
# 
# 
#   high_val_outcome <- lapply(
#     X = currencies,
#     FUN = extract_outcome,
#     result = high_result, 
#     outcome = outcome
#   )
#   names(high_val_outcome) <- currencies
# 
# 
#   ## output value
#   val <- lapply(
#     X = currencies,
#     FUN = function(x){
#       val <- data.frame(
#         outcome = outcome,
#         parameter = param,
#         base_par_val = base_val_param,
#         low_par_val = low_val_param,
#         high_par_val = high_val_param,
#         alternative = base_val_outcome[[x]] %>% dplyr::select(alternative),
#         base_outcome_val = base_val_outcome[[x]] %>% dplyr::select(all_of(outcome)) %>% unlist() %>% as.vector(),
#         low_outcome_val = low_val_outcome[[x]] %>% dplyr::select(all_of(outcome)) %>% unlist() %>% as.vector(),
#         high_outcome_val = high_val_outcome[[x]] %>% dplyr::select(all_of(outcome)) %>% unlist() %>% as.vector()
#       )
#       return(val)
#     } 
#   )
#   names(val) <- currencies


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
#     wtp_val <- gdp * (1 / as.numeric(params[["brl_usd_fx_2025"]]))
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

