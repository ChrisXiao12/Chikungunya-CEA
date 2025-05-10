#RDSA - redoing my DSA in R vs excel
df_IXCHIQ <- data.frame(
  variable = c(
    "pbeta", "plambda", "pgamma", "pdelta", "pphi", "ppsi", "pmu", "pkappa", "pomega",
    "discount", "popsize", "fracsusceptible", "fracinfected", "cycle_length",
    "U_S", "U_V", "U_E", "U_I", "U_R", "U_C", "U_D",
    "C_S", "C_V", "C_E", "C_I", "C_R", "C_C", "C_D"
  ),
  base_case = c(
    0.939189937, 0.753403036, 0.56112007, 0.001, 0.989, 0.010445067, 6.95098e-06,
    0.018995353, 0.521, 0.03, 1000, 0.9980005, 0.0019995, 0.019230769,
    0.824, 0.00024, 0.700, 0.662, 0.816, 0.494, 0,
    0, 313.89, 0, 235.44, 0, 685.64, 0
  ),
  LB = c(
    0.813626024, 0.441964854, 0.503414696, 0.00075, 0.967, 0.007844073, 5.21324e-06,
    0.009543213, 0.445, 0.03, 1000, 0.9985, 0.0015, 0,
    0.618, 0.00018, 0.525, 0.496, 0.612, 0.176, 0,
    0, 2.09, 0, 41.75, 0, 514.23, 0
  ),
  UB = c(
    0.985004423, 0.999088118, 0.632120559, 0.00125, 0.998, 0.013039243, 8.68871e-06,
    0.055910442, 0.597, 0.03, 1000, 0.9975, 0.0025, 0,
    1.000, 0.00030, 0.876, 0.827, 1.000, 0.759, 0,
    0, 532.89, 0, 888.37, 0, 857.05, 0
  )
)

df_VIMKUNYA <- data.frame(
  variable = c(
    "pbeta", "plambda", "pgamma", "pdelta", "pphi", "ppsi", "pmu", "pkappa", "pomega",
    "discount", "popsize", "fracsusceptible", "fracinfected", "cycle_length",
    "U_S", "U_V", "U_E", "U_I", "U_R", "U_C", "U_D",
    "C_S", "C_V", "C_E", "C_I", "C_R", "C_C", "C_D"
  ),
  base_case = c(
    0.939189937, 0.753403036, 0.56112007, 0.001, 0.978, 0.010445067, 6.95098e-06,
    0.018995353, 0.521, 0.03, 1000, 0.9980005, 0.0019995, 0.019230769,
    0.824, 0.000331, 0.700, 0.662, 0.816, 0.494, 0,
    0, 313.89, 0, 235.44, 0, 685.64, 0
  ),
  LB = c(
    0.813626024, 0.441964854, 0.503414696, 0.00075, 0.962, 0.007844073, 5.21324e-06,
    0.009543213, 0.445, 0.03, 1000, 0.9985, 0.0015, 0,
    0.618, 0.000248, 0.525, 0.496, 0.612, 0.176, 0,
    0, 2.09, 0, 41.75, 0, 514.23, 0
  ),
  UB = c(
    0.985004423, 0.999088118, 0.632120559, 0.00125, 0.987, 0.013039243, 8.68871e-06,
    0.055910442, 0.597, 0.03, 1000, 0.9975, 0.0025, 0,
    1.000, 0.000414, 0.876, 0.827, 1.000, 0.759, 0,
    0, 532.89, 0, 888.37, 0, 857.05, 0
  )
)

#----

wtp <- 6808
param_df <- df_IXCHIQ
base_params <- setNames(param_df$base_case, param_df$variable)
get_inmb <- function(params) {
  res <- run_SVEIRD4(params)
  v_nmb <- res["v_eff_d"] * wtp - res["v_cost_d"]
  nv_nmb <- res["nv_eff_d"] * wtp - res["nv_cost_d"]
  return(v_nmb - nv_nmb)
}
DSA_results <- data.frame(Parameter = character(), Bound = character(), INMB = numeric())
for (i in seq_len(nrow(param_df))) {
  param_name <- param_df$variable[i]
  params_lb <- base_params
  params_lb[param_name] <- param_df$LB[i]
  inmb_lb <- get_inmb(as.list(params_lb))
  params_ub <- base_params
  params_ub[param_name] <- param_df$UB[i]
  inmb_ub <- get_inmb(as.list(params_ub))
  DSA_results <- rbind(
    DSA_results,
    data.frame(Parameter = param_name, Bound = "Lower Bound", INMB = inmb_lb),
    data.frame(Parameter = param_name, Bound = "Upper Bound", INMB = inmb_ub)
  )
}
DSA_results_IXCHIQ <- DSA_results
DSA_results_VIMKUNYA <- DSA_results
#----
DSA_IXCHIQ_wide <- DSA_results_IXCHIQ %>% pivot_wider(names_from = Bound, values_from = INMB)
DSA_VIMKUNYA_wide <- DSA_results_VIMKUNYA %>% pivot_wider(names_from = Bound, values_from = INMB)
DSA_IXCHIQ_wide <- DSA_IXCHIQ_wide %>% mutate(range = abs(`Upper Bound` - `Lower Bound`))
DSA_VIMKUNYA_wide <- DSA_VIMKUNYA_wide %>% mutate(range = abs(`Upper Bound` - `Lower Bound`))
DSA_IXCHIQ_wide <- DSA_IXCHIQ_wide %>% arrange(desc(range))
DSA_IXCHIQ_wide <- DSA_IXCHIQ_wide[1:14,]
DSA_VIMKUNYA_wide <- DSA_VIMKUNYA_wide %>% arrange(desc(range))
DSA_VIMKUNYA_wide <- DSA_VIMKUNYA_wide[1:14,]
#----
BaseI <- 7372824
sorted_I <- DSA_IXCHIQ_wide %>%
  arrange(desc(range)) %>%
  pull(Parameter)
plot_data <- DSA_IXCHIQ_wide %>%
  pivot_longer(cols = c("Lower Bound", "Upper Bound"),
               names_to = "Bound",
               values_to = "Value") %>%
  mutate(
    xmin = pmin(Value, BaseI),
    xmax = pmax(Value, BaseI),
    Parameter = factor(Parameter, levels = rev(sorted_I))
  )
BaseV <- 7282347
sorted_V <- DSA_VIMKUNYA_wide %>%
  arrange(desc(range)) %>%
  pull(Parameter)
plot_data <- DSA_VIMKUNYA_wide %>%
  pivot_longer(cols = c("Lower Bound", "Upper Bound"),
               names_to = "Bound",
               values_to = "Value") %>%
  mutate(
    xmin = pmin(Value, BaseV),
    xmax = pmax(Value, BaseV),
    Parameter = factor(Parameter, levels = rev(sorted_V))
  )
#----
param_labels <- c(
  pkappa = "Chronic Disease Recovery Probability",
  U_S = "Utility Susceptible",
  U_R = "Utility Recovered",
  U_C = "Utility Chronic Disease",
  C_C = "Cost Chronic Disease",
  C_I = "Cost Infection",
  C_V = "Cost Vaccination",
  ppsi = "Vaccination Rate",
  pomega = "Probability of Chronic Disease",
  pphi = "Vaccine Efficacy",
  plambda = "Infectious Rate",
  pbeta = "Transmission Rate",
  pgamma = "Probability Recover Infectious",
  fracinfected = "Initial Infection Proportion"
)

ggplot(plot_data) +
  geom_segment(aes(x = xmin, xend = xmax, y = Parameter, yend = Parameter, color = Bound),
               size = 6) +
  geom_vline(xintercept = BaseV, linetype = "dashed", color = "black") +
  scale_color_manual(values = c("Lower Bound" = "#3182bd", "Upper Bound" = "#6baed6")) +
  scale_x_continuous(
    limits = c(0,15e6),
    breaks = seq(0, 15e6, by = 1e6),
    labels = function(x) sprintf("%.1f", x / 1e6)
  ) +
  scale_y_discrete(labels = param_labels) +
  labs(
    x = "INMB (Millions USD)",
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


