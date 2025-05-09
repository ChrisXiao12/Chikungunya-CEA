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
    7485090.54, 7689604.79, 6244899.91, 5607150.37, 7149146.36,
    7949792.12, 5685101.40, 8612160.22, 7523126.44, 7777544.95,
    7157323.66
  ),
  UB = c(
    7139951.79, 6944154.13, 8319795.07, 8957544.60, 7356436.26,
    6821642.22, 8609363.08, 6710511.94, 7057434.87, 6897870.98,
    7600662.02
  ),
  Base = rep(7282347.49, 11)
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
    Parameter = factor(Parameter, levels = rev(sorted_order))  # this is key!
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
