library(dplyr)

model1 <- readRDS("SimulationResults/Model1ParameterRecoveryResults.rds")
model2 <- readRDS("SimulationResults/Model2ParameterRecoveryResults.rds")
model3 <- readRDS("SimulationResults/Model3ParameterRecoveryResults.rds")
model4 <- readRDS("SimulationResults/Model4ParameterRecoveryResults.rds")


# Model 1

model1_recovery <- data.frame(
  model = "Model 1",
  parameter = c(
    "beta[1]", "beta[2]", "beta[3]", "beta[4]",
    "kappa[1]", "kappa[2]", "kappa[3]", "kappa[4]"
  ),
  n = 100
)

model1_recovery$bias <- c(
  mean(model1$est_beta[,1] - model1$true_beta[,1]),
  mean(model1$est_beta[,2] - model1$true_beta[,2]),
  mean(model1$est_beta[,3] - model1$true_beta[,3]),
  mean(model1$est_beta[,4] - model1$true_beta[,4]),
  mean(model1$est_kappa[,1] - model1$true_kappa[,1]),
  mean(model1$est_kappa[,2] - model1$true_kappa[,2]),
  mean(model1$est_kappa[,3] - model1$true_kappa[,3]),
  mean(model1$est_kappa[,4] - model1$true_kappa[,4])
)

model1_recovery$rmse <- c(
  sqrt(mean((model1$est_beta[,1] - model1$true_beta[,1])^2)),
  sqrt(mean((model1$est_beta[,2] - model1$true_beta[,2])^2)),
  sqrt(mean((model1$est_beta[,3] - model1$true_beta[,3])^2)),
  sqrt(mean((model1$est_beta[,4] - model1$true_beta[,4])^2)),
  sqrt(mean((model1$est_kappa[,1] - model1$true_kappa[,1])^2)),
  sqrt(mean((model1$est_kappa[,2] - model1$true_kappa[,2])^2)),
  sqrt(mean((model1$est_kappa[,3] - model1$true_kappa[,3])^2)),
  sqrt(mean((model1$est_kappa[,4] - model1$true_kappa[,4])^2))
)

model1_recovery$correlation <- c(
  cor(model1$true_beta[,1], model1$est_beta[,1]),
  cor(model1$true_beta[,2], model1$est_beta[,2]),
  cor(model1$true_beta[,3], model1$est_beta[,3]),
  cor(model1$true_beta[,4], model1$est_beta[,4]),
  cor(model1$true_kappa[,1], model1$est_kappa[,1]),
  cor(model1$true_kappa[,2], model1$est_kappa[,2]),
  cor(model1$true_kappa[,3], model1$est_kappa[,3]),
  cor(model1$true_kappa[,4], model1$est_kappa[,4])
)

model1_recovery$coverage <- c(
  mean(model1$true_beta[,1] >= model1$lower_beta[,1] &
         model1$true_beta[,1] <= model1$upper_beta[,1]),
  mean(model1$true_beta[,2] >= model1$lower_beta[,2] &
         model1$true_beta[,2] <= model1$upper_beta[,2]),
  mean(model1$true_beta[,3] >= model1$lower_beta[,3] &
         model1$true_beta[,3] <= model1$upper_beta[,3]),
  mean(model1$true_beta[,4] >= model1$lower_beta[,4] &
         model1$true_beta[,4] <= model1$upper_beta[,4]),
  mean(model1$true_kappa[,1] >= model1$lower_kappa[,1] &
         model1$true_kappa[,1] <= model1$upper_kappa[,1]),
  mean(model1$true_kappa[,2] >= model1$lower_kappa[,2] &
         model1$true_kappa[,2] <= model1$upper_kappa[,2]),
  mean(model1$true_kappa[,3] >= model1$lower_kappa[,3] &
         model1$true_kappa[,3] <= model1$upper_kappa[,3]),
  mean(model1$true_kappa[,4] >= model1$lower_kappa[,4] &
         model1$true_kappa[,4] <= model1$upper_kappa[,4])
)

model1_recovery$mean_CI_width <- c(
  mean(model1$upper_beta[,1] - model1$lower_beta[,1]),
  mean(model1$upper_beta[,2] - model1$lower_beta[,2]),
  mean(model1$upper_beta[,3] - model1$lower_beta[,3]),
  mean(model1$upper_beta[,4] - model1$lower_beta[,4]),
  mean(model1$upper_kappa[,1] - model1$lower_kappa[,1]),
  mean(model1$upper_kappa[,2] - model1$lower_kappa[,2]),
  mean(model1$upper_kappa[,3] - model1$lower_kappa[,3]),
  mean(model1$upper_kappa[,4] - model1$lower_kappa[,4])
)


# Model 2

model2_recovery <- data.frame(
  model = "Model 2",
  parameter = c(
    "beta[1]", "beta[2]", "beta[3]", "beta[4]",
    "kappa[1]", "kappa[2]", "kappa[3]", "kappa[4]"
  ),
  n = 100
)

model2_recovery$bias <- c(
  mean(model2$est_beta[,1] - model2$true_beta[,1]),
  mean(model2$est_beta[,2] - model2$true_beta[,2]),
  mean(model2$est_beta[,3] - model2$true_beta[,3]),
  mean(model2$est_beta[,4] - model2$true_beta[,4]),
  mean(model2$est_kappa[,1] - model2$true_kappa[,1]),
  mean(model2$est_kappa[,2] - model2$true_kappa[,2]),
  mean(model2$est_kappa[,3] - model2$true_kappa[,3]),
  mean(model2$est_kappa[,4] - model2$true_kappa[,4])
)

model2_recovery$rmse <- c(
  sqrt(mean((model2$est_beta[,1] - model2$true_beta[,1])^2)),
  sqrt(mean((model2$est_beta[,2] - model2$true_beta[,2])^2)),
  sqrt(mean((model2$est_beta[,3] - model2$true_beta[,3])^2)),
  sqrt(mean((model2$est_beta[,4] - model2$true_beta[,4])^2)),
  sqrt(mean((model2$est_kappa[,1] - model2$true_kappa[,1])^2)),
  sqrt(mean((model2$est_kappa[,2] - model2$true_kappa[,2])^2)),
  sqrt(mean((model2$est_kappa[,3] - model2$true_kappa[,3])^2)),
  sqrt(mean((model2$est_kappa[,4] - model2$true_kappa[,4])^2))
)

model2_recovery$correlation <- c(
  cor(model2$true_beta[,1], model2$est_beta[,1]),
  cor(model2$true_beta[,2], model2$est_beta[,2]),
  cor(model2$true_beta[,3], model2$est_beta[,3]),
  cor(model2$true_beta[,4], model2$est_beta[,4]),
  cor(model2$true_kappa[,1], model2$est_kappa[,1]),
  cor(model2$true_kappa[,2], model2$est_kappa[,2]),
  cor(model2$true_kappa[,3], model2$est_kappa[,3]),
  cor(model2$true_kappa[,4], model2$est_kappa[,4])
)

model2_recovery$coverage <- c(
  mean(model2$true_beta[,1] >= model2$lower_beta[,1] &
         model2$true_beta[,1] <= model2$upper_beta[,1]),
  mean(model2$true_beta[,2] >= model2$lower_beta[,2] &
         model2$true_beta[,2] <= model2$upper_beta[,2]),
  mean(model2$true_beta[,3] >= model2$lower_beta[,3] &
         model2$true_beta[,3] <= model2$upper_beta[,3]),
  mean(model2$true_beta[,4] >= model2$lower_beta[,4] &
         model2$true_beta[,4] <= model2$upper_beta[,4]),
  mean(model2$true_kappa[,1] >= model2$lower_kappa[,1] &
         model2$true_kappa[,1] <= model2$upper_kappa[,1]),
  mean(model2$true_kappa[,2] >= model2$lower_kappa[,2] &
         model2$true_kappa[,2] <= model2$upper_kappa[,2]),
  mean(model2$true_kappa[,3] >= model2$lower_kappa[,3] &
         model2$true_kappa[,3] <= model2$upper_kappa[,3]),
  mean(model2$true_kappa[,4] >= model2$lower_kappa[,4] &
         model2$true_kappa[,4] <= model2$upper_kappa[,4])
)

model2_recovery$mean_CI_width <- c(
  mean(model2$upper_beta[,1] - model2$lower_beta[,1]),
  mean(model2$upper_beta[,2] - model2$lower_beta[,2]),
  mean(model2$upper_beta[,3] - model2$lower_beta[,3]),
  mean(model2$upper_beta[,4] - model2$lower_beta[,4]),
  mean(model2$upper_kappa[,1] - model2$lower_kappa[,1]),
  mean(model2$upper_kappa[,2] - model2$lower_kappa[,2]),
  mean(model2$upper_kappa[,3] - model2$lower_kappa[,3]),
  mean(model2$upper_kappa[,4] - model2$lower_kappa[,4])
)


# Model 3

model3_recovery <- data.frame(
  model = "Model 3",
  parameter = c(
    "beta[1]", "beta[2]", "beta[3]", "beta[4]",
    "kappa[1]", "kappa[2]", "kappa[3]", "kappa[4]"
  ),
  n = 100
)

model3_recovery$bias <- c(
  mean(model3$est_beta[,1] - model3$true_beta[,1]),
  mean(model3$est_beta[,2] - model3$true_beta[,2]),
  mean(model3$est_beta[,3] - model3$true_beta[,3]),
  mean(model3$est_beta[,4] - model3$true_beta[,4]),
  mean(model3$est_kappa[,1] - model3$true_kappa[,1]),
  mean(model3$est_kappa[,2] - model3$true_kappa[,2]),
  mean(model3$est_kappa[,3] - model3$true_kappa[,3]),
  mean(model3$est_kappa[,4] - model3$true_kappa[,4])
)

model3_recovery$rmse <- c(
  sqrt(mean((model3$est_beta[,1] - model3$true_beta[,1])^2)),
  sqrt(mean((model3$est_beta[,2] - model3$true_beta[,2])^2)),
  sqrt(mean((model3$est_beta[,3] - model3$true_beta[,3])^2)),
  sqrt(mean((model3$est_beta[,4] - model3$true_beta[,4])^2)),
  sqrt(mean((model3$est_kappa[,1] - model3$true_kappa[,1])^2)),
  sqrt(mean((model3$est_kappa[,2] - model3$true_kappa[,2])^2)),
  sqrt(mean((model3$est_kappa[,3] - model3$true_kappa[,3])^2)),
  sqrt(mean((model3$est_kappa[,4] - model3$true_kappa[,4])^2))
)

model3_recovery$correlation <- c(
  cor(model3$true_beta[,1], model3$est_beta[,1]),
  cor(model3$true_beta[,2], model3$est_beta[,2]),
  cor(model3$true_beta[,3], model3$est_beta[,3]),
  cor(model3$true_beta[,4], model3$est_beta[,4]),
  cor(model3$true_kappa[,1], model3$est_kappa[,1]),
  cor(model3$true_kappa[,2], model3$est_kappa[,2]),
  cor(model3$true_kappa[,3], model3$est_kappa[,3]),
  cor(model3$true_kappa[,4], model3$est_kappa[,4])
)

model3_recovery$coverage <- c(
  mean(model3$true_beta[,1] >= model3$lower_beta[,1] &
         model3$true_beta[,1] <= model3$upper_beta[,1]),
  mean(model3$true_beta[,2] >= model3$lower_beta[,2] &
         model3$true_beta[,2] <= model3$upper_beta[,2]),
  mean(model3$true_beta[,3] >= model3$lower_beta[,3] &
         model3$true_beta[,3] <= model3$upper_beta[,3]),
  mean(model3$true_beta[,4] >= model3$lower_beta[,4] &
         model3$true_beta[,4] <= model3$upper_beta[,4]),
  mean(model3$true_kappa[,1] >= model3$lower_kappa[,1] &
         model3$true_kappa[,1] <= model3$upper_kappa[,1]),
  mean(model3$true_kappa[,2] >= model3$lower_kappa[,2] &
         model3$true_kappa[,2] <= model3$upper_kappa[,2]),
  mean(model3$true_kappa[,3] >= model3$lower_kappa[,3] &
         model3$true_kappa[,3] <= model3$upper_kappa[,3]),
  mean(model3$true_kappa[,4] >= model3$lower_kappa[,4] &
         model3$true_kappa[,4] <= model3$upper_kappa[,4])
)

model3_recovery$mean_CI_width <- c(
  mean(model3$upper_beta[,1] - model3$lower_beta[,1]),
  mean(model3$upper_beta[,2] - model3$lower_beta[,2]),
  mean(model3$upper_beta[,3] - model3$lower_beta[,3]),
  mean(model3$upper_beta[,4] - model3$lower_beta[,4]),
  mean(model3$upper_kappa[,1] - model3$lower_kappa[,1]),
  mean(model3$upper_kappa[,2] - model3$lower_kappa[,2]),
  mean(model3$upper_kappa[,3] - model3$lower_kappa[,3]),
  mean(model3$upper_kappa[,4] - model3$lower_kappa[,4])
)


# Model 4

model4_recovery <- data.frame(
  model = "Model 4",
  parameter = c(
    "beta[1]", "beta[2]", "beta[3]", "beta[4]",
    "kappa[1]", "kappa[2]", "kappa[3]", "kappa[4]"
  ),
  n = 100
)

model4_recovery$bias <- c(
  mean(model4$est_beta[,1] - model4$true_beta[,1]),
  mean(model4$est_beta[,2] - model4$true_beta[,2]),
  mean(model4$est_beta[,3] - model4$true_beta[,3]),
  mean(model4$est_beta[,4] - model4$true_beta[,4]),
  mean(model4$est_kappa[,1] - model4$true_kappa[,1]),
  mean(model4$est_kappa[,2] - model4$true_kappa[,2]),
  mean(model4$est_kappa[,3] - model4$true_kappa[,3]),
  mean(model4$est_kappa[,4] - model4$true_kappa[,4])
)

model4_recovery$rmse <- c(
  sqrt(mean((model4$est_beta[,1] - model4$true_beta[,1])^2)),
  sqrt(mean((model4$est_beta[,2] - model4$true_beta[,2])^2)),
  sqrt(mean((model4$est_beta[,3] - model4$true_beta[,3])^2)),
  sqrt(mean((model4$est_beta[,4] - model4$true_beta[,4])^2)),
  sqrt(mean((model4$est_kappa[,1] - model4$true_kappa[,1])^2)),
  sqrt(mean((model4$est_kappa[,2] - model4$true_kappa[,2])^2)),
  sqrt(mean((model4$est_kappa[,3] - model4$true_kappa[,3])^2)),
  sqrt(mean((model4$est_kappa[,4] - model4$true_kappa[,4])^2))
)

model4_recovery$correlation <- c(
  cor(model4$true_beta[,1], model4$est_beta[,1]),
  cor(model4$true_beta[,2], model4$est_beta[,2]),
  cor(model4$true_beta[,3], model4$est_beta[,3]),
  cor(model4$true_beta[,4], model4$est_beta[,4]),
  cor(model4$true_kappa[,1], model4$est_kappa[,1]),
  cor(model4$true_kappa[,2], model4$est_kappa[,2]),
  cor(model4$true_kappa[,3], model4$est_kappa[,3]),
  cor(model4$true_kappa[,4], model4$est_kappa[,4])
)

model4_recovery$coverage <- c(
  mean(model4$true_beta[,1] >= model4$lower_beta[,1] &
         model4$true_beta[,1] <= model4$upper_beta[,1]),
  mean(model4$true_beta[,2] >= model4$lower_beta[,2] &
         model4$true_beta[,2] <= model4$upper_beta[,2]),
  mean(model4$true_beta[,3] >= model4$lower_beta[,3] &
         model4$true_beta[,3] <= model4$upper_beta[,3]),
  mean(model4$true_beta[,4] >= model4$lower_beta[,4] &
         model4$true_beta[,4] <= model4$upper_beta[,4]),
  mean(model4$true_kappa[,1] >= model4$lower_kappa[,1] &
         model4$true_kappa[,1] <= model4$upper_kappa[,1]),
  mean(model4$true_kappa[,2] >= model4$lower_kappa[,2] &
         model4$true_kappa[,2] <= model4$upper_kappa[,2]),
  mean(model4$true_kappa[,3] >= model4$lower_kappa[,3] &
         model4$true_kappa[,3] <= model4$upper_kappa[,3]),
  mean(model4$true_kappa[,4] >= model4$lower_kappa[,4] &
         model4$true_kappa[,4] <= model4$upper_kappa[,4])
)

model4_recovery$mean_CI_width <- c(
  mean(model4$upper_beta[,1] - model4$lower_beta[,1]),
  mean(model4$upper_beta[,2] - model4$lower_beta[,2]),
  mean(model4$upper_beta[,3] - model4$lower_beta[,3]),
  mean(model4$upper_beta[,4] - model4$lower_beta[,4]),
  mean(model4$upper_kappa[,1] - model4$lower_kappa[,1]),
  mean(model4$upper_kappa[,2] - model4$lower_kappa[,2]),
  mean(model4$upper_kappa[,3] - model4$lower_kappa[,3]),
  mean(model4$upper_kappa[,4] - model4$lower_kappa[,4])
)


# Combine models

beta_kappa_recovery <- rbind(
  model1_recovery,
  model2_recovery,
  model3_recovery,
  model4_recovery
)

write.csv(
  beta_kappa_recovery,
  "FiguresAndTables/TableA.1_ParameterRecovery.csv",
  row.names = FALSE
)
