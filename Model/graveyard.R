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