library(ggplot2)
library(dplyr)
library(tidyr)

# load simulation results
model1 <- readRDS("SimulationResults/Model1ParameterRecoveryResults.rds")
model2 <- readRDS("SimulationResults/Model2ParameterRecoveryResults.rds")
model3 <- readRDS("SimulationResults/Model3ParameterRecoveryResults.rds")
model4 <- readRDS("SimulationResults/Model4ParameterRecoveryResults.rds")


# create data for model 1
beta_model1 <- data.frame(
  True_beta1 = model1$true_beta[, 1],
  Estimated_beta1 = model1$est_beta[, 1],
  True_beta2 = model1$true_beta[, 2],
  Estimated_beta2 = model1$est_beta[, 2],
  True_beta3 = model1$true_beta[, 3],
  Estimated_beta3 = model1$est_beta[, 3],
  True_beta4 = model1$true_beta[, 4],
  Estimated_beta4 = model1$est_beta[, 4],
  Model = "Model 1"
)


# create data for model 2
beta_model2 <- data.frame(
  True_beta1 = model2$true_beta[, 1],
  Estimated_beta1 = model2$est_beta[, 1],
  True_beta2 = model2$true_beta[, 2],
  Estimated_beta2 = model2$est_beta[, 2],
  True_beta3 = model2$true_beta[, 3],
  Estimated_beta3 = model2$est_beta[, 3],
  True_beta4 = model2$true_beta[, 4],
  Estimated_beta4 = model2$est_beta[, 4],
  Model = "Model 2"
)


# create data for model 3
beta_model3 <- data.frame(
  True_beta1 = model3$true_beta[, 1],
  Estimated_beta1 = model3$est_beta[, 1],
  True_beta2 = model3$true_beta[, 2],
  Estimated_beta2 = model3$est_beta[, 2],
  True_beta3 = model3$true_beta[, 3],
  Estimated_beta3 = model3$est_beta[, 3],
  True_beta4 = model3$true_beta[, 4],
  Estimated_beta4 = model3$est_beta[, 4],
  Model = "Model 3"
)


# create data for model 4
beta_model4 <- data.frame(
  True_beta1 = model4$true_beta[, 1],
  Estimated_beta1 = model4$est_beta[, 1],
  True_beta2 = model4$true_beta[, 2],
  Estimated_beta2 = model4$est_beta[, 2],
  True_beta3 = model4$true_beta[, 3],
  Estimated_beta3 = model4$est_beta[, 3],
  True_beta4 = model4$true_beta[, 4],
  Estimated_beta4 = model4$est_beta[, 4],
  Model = "Model 4"
)


# combine all models
beta_data <- rbind(
  beta_model1,
  beta_model2,
  beta_model3,
  beta_model4
)


# create data for beta 1
beta1 <- data.frame(
  Model = beta_data$Model,
  True = beta_data$True_beta1,
  Estimated = beta_data$Estimated_beta1,
  Parameter = "beta[1]"
)


# create data for beta 2
beta2 <- data.frame(
  Model = beta_data$Model,
  True = beta_data$True_beta2,
  Estimated = beta_data$Estimated_beta2,
  Parameter = "beta[2]"
)


# create data for beta 3
beta3 <- data.frame(
  Model = beta_data$Model,
  True = beta_data$True_beta3,
  Estimated = beta_data$Estimated_beta3,
  Parameter = "beta[3]"
)


# create data for beta 4
beta4 <- data.frame(
  Model = beta_data$Model,
  True = beta_data$True_beta4,
  Estimated = beta_data$Estimated_beta4,
  Parameter = "beta[4]"
)


# combine beta 1 to beta 4
beta_long <- rbind(
  beta1,
  beta2,
  beta3,
  beta4
)


# set the order of the models
beta_long$Model <- factor(
  beta_long$Model,
  levels = c(
    "Model 1",
    "Model 2",
    "Model 3",
    "Model 4"
  )
)


# set the order of the parameters
beta_long$Parameter <- factor(
  beta_long$Parameter,
  levels = c(
    "beta[1]",
    "beta[2]",
    "beta[3]",
    "beta[4]"
  )
)


# create the plot
ggplot(
  beta_long,
  aes(
    x = True,
    y = Estimated
  )
) +
  geom_point(
    alpha = 0.5,
    colour = "black",
    size = 1.5
  ) +
  geom_abline(
    intercept = 0,
    slope = 1,
    linetype = "dashed"
  ) +
  facet_grid(
    Model ~ Parameter,
    scales = "free"
  ) +
  labs(
    x = "True value",
    y = "Posterior mean"
  ) +
  theme_bw()


# Figure 4.6: Kappa parameter recovery

kappa_model1 <- data.frame(
  true_kappa1 = model1$true_kappa[,1],
  true_kappa2 = model1$true_kappa[,2],
  true_kappa3 = model1$true_kappa[,3],
  true_kappa4 = model1$true_kappa[,4],
  est_kappa1 = model1$est_kappa[,1],
  est_kappa2 = model1$est_kappa[,2],
  est_kappa3 = model1$est_kappa[,3],
  est_kappa4 = model1$est_kappa[,4],
  Model = "Model 1"
)

kappa_model2 <- data.frame(
  true_kappa1 = model2$true_kappa[,1],
  true_kappa2 = model2$true_kappa[,2],
  true_kappa3 = model2$true_kappa[,3],
  true_kappa4 = model2$true_kappa[,4],
  est_kappa1 = model2$est_kappa[,1],
  est_kappa2 = model2$est_kappa[,2],
  est_kappa3 = model2$est_kappa[,3],
  est_kappa4 = model2$est_kappa[,4],
  Model = "Model 2"
)

kappa_model3 <- data.frame(
  true_kappa1 = model3$true_kappa[,1],
  true_kappa2 = model3$true_kappa[,2],
  true_kappa3 = model3$true_kappa[,3],
  true_kappa4 = model3$true_kappa[,4],
  est_kappa1 = model3$est_kappa[,1],
  est_kappa2 = model3$est_kappa[,2],
  est_kappa3 = model3$est_kappa[,3],
  est_kappa4 = model3$est_kappa[,4],
  Model = "Model 3"
)

kappa_model4 <- data.frame(
  true_kappa1 = model4$true_kappa[,1],
  true_kappa2 = model4$true_kappa[,2],
  true_kappa3 = model4$true_kappa[,3],
  true_kappa4 = model4$true_kappa[,4],
  est_kappa1 = model4$est_kappa[,1],
  est_kappa2 = model4$est_kappa[,2],
  est_kappa3 = model4$est_kappa[,3],
  est_kappa4 = model4$est_kappa[,4],
  Model = "Model 4"
)

kappa_data <- rbind(
  kappa_model1,
  kappa_model2,
  kappa_model3,
  kappa_model4
)

kappa1 <- kappa_data %>%
  select(Model, true_kappa1, est_kappa1) %>%
  rename(True = true_kappa1, Estimated = est_kappa1) %>%
  mutate(Parameter = "kappa[1]")

kappa2 <- kappa_data %>%
  select(Model, true_kappa2, est_kappa2) %>%
  rename(True = true_kappa2, Estimated = est_kappa2) %>%
  mutate(Parameter = "kappa[2]")

kappa3 <- kappa_data %>%
  select(Model, true_kappa3, est_kappa3) %>%
  rename(True = true_kappa3, Estimated = est_kappa3) %>%
  mutate(Parameter = "kappa[3]")

kappa4 <- kappa_data %>%
  select(Model, true_kappa4, est_kappa4) %>%
  rename(True = true_kappa4, Estimated = est_kappa4) %>%
  mutate(Parameter = "kappa[4]")

kappa_long <- rbind(kappa1, kappa2, kappa3, kappa4)

kappa_long$Model <- factor(
  kappa_long$Model,
  levels = c("Model 1", "Model 2", "Model 3", "Model 4")
)

kappa_long$Parameter <- factor(
  kappa_long$Parameter,
  levels = c("kappa[1]", "kappa[2]", "kappa[3]", "kappa[4]")
)

ggplot(kappa_long, aes(x = True, y = Estimated)) +
  geom_point(alpha = 0.5, colour = "black", size = 1.5) +
  geom_abline(
    intercept = 0,
    slope = 1,
    linetype = "dashed"
  ) +
  facet_grid(
    Model ~ Parameter,
    scales = "free"
  ) +
  labs(
    x = "True value",
    y = "Posterior mean"
  ) +
  theme_bw()