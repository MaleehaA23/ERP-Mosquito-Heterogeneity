library(ggplot2)

# Load parameter recovery results

recovery <- read.csv(
  "FiguresAndTables/TableA.1_ParameterRecovery.csv"
)


# Create parameter type and net number

recovery$Type <- ifelse(
  grepl("beta", recovery$parameter),
  "Beta",
  "Kappa"
)

recovery$Net <- rep(1:4, 8)




recovery$model <- factor(
  recovery$model,
  levels = c("Model 1", "Model 2", "Model 3", "Model 4")
)

recovery$Type <- factor(
  recovery$Type,
  levels = c("Beta", "Kappa")
)


# Plot coverage

coverage_plot <- ggplot(
  recovery,
  aes(
    x = factor(Net),
    y = coverage * 100
  )
) +
  geom_col(
    fill = "grey40"
  ) +
  geom_hline(
    yintercept = 95,
    linetype = "dashed"
  ) +
  facet_grid(
    Type ~ model
  ) +
  scale_y_continuous(
    limits = c(0, 105),
    breaks = c(0, 25, 50, 75, 95)
  ) +
  labs(
    x = "Net",
    y = "95% credible interval coverage (%)"
  ) +
  theme_bw()

coverage_plot