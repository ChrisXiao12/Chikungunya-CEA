library(ggplot2)
library(dplyr)
library(tidyr)

INMB_data <- data.frame(
  Parameter = c(
    "Vaccine Cost", "Utility Chronic disease", "Chronic disease probability",
    "Chronic disease cost", "Vaccine efficacy", "Disease transmission rate",
    "Vaccination rate", "Infectious rate", "Discount rate", "Initial fraction infected", "Cost of infection"
  ),
  LB = c(
    7576495.45, 7785054.91, 6322707.79, 5677170.28, 7190980.18,
    8002100.57, 5766324.27, 8700089.99, 7616724.41, 7871432.84,
    7246273.73
  ),
  UB = c(
    7229776.79, 7030501.06, 8422940.84, 9068478.35, 7446153.54,
    6917088.94, 8698970.22, 6799094.94, 7144998.84, 6985259.79,
    7695026.02
  ),
  Base = rep(7372824.32, 11)
)

INMB_data <- INMB_data %>%
  mutate(Range = abs(UB - LB))

sorted_order <- INMB_data %>%
  arrange(desc(Range)) %>%
  pull(Parameter)

plot_data <- INMB_data %>%
  pivot_longer(cols = c("LB", "UB"),
               names_to = "Bound",
               values_to = "Value") %>%
  mutate(
    xmin = pmin(Value, Base),
    xmax = pmax(Value, Base),
    Bound = ifelse(Bound == "LB", "Lower Bound", "Upper Bound"),
    Parameter = factor(Parameter, levels = rev(sorted_order))
  )


#----
#fixed
ggplot(plot_data) +
  geom_segment(aes(x = xmin, xend = xmax, y = Parameter, yend = Parameter, color = Bound),
               size = 6) +
  geom_vline(xintercept = INMB_data$Base[1], linetype = "dashed", color = "black") +
  scale_color_manual(values = c("Lower Bound" = "#3182bd", "Upper Bound" = "#6baed6")) +
  scale_x_continuous(
    breaks = seq(6e6, 10e6, by = 1e6),
    labels = function(x) sprintf("%.1f", x / 1e6)
  ) +
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
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.title = element_blank(),
    legend.position = "bottom"
  )

