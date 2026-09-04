library(rstan)
library(dplyr)
library(ggplot2)


# FIT FILES


fit_files <- list(
  Arabiensis = list(
    `Model 1` = "ArabiensisNewResults/Model1/Model1_fit.rds",
    `Model 2` = "ArabiensisNewResults/Model2/Model2_fit.rds",
    `Model 3` = "ArabiensisNewResults/Model3/Model3_fit.rds",
    `Model 4` = "ArabiensisNewResults/Model4/Model4_fit.rds"
  ),
  
  Culex = list(
    `Model 1` = "CulexNewResults/Model1/Model1_fit.rds",
    `Model 2` = "CulexNewResults/Model2/Model2_fit.rds",
    `Model 3` = "CulexNewResults/Model3/Model3_fit.rds",
    `Model 4` = "CulexNewResults/Model4/Model4_fit.rds"
  )
)


# EXTRACT BETA AND KAPPA SUMMARIES


plot_data <- data.frame()

for (species in names(fit_files)) {
  
  for (model in names(fit_files[[species]])) {
    
    fit <- readRDS(
      fit_files[[species]][[model]]
    )
    
    post <- rstan::extract(
      fit,
      pars = c("beta_out", "kappa_out")
    )
    
    
    # Beta
    
    for (n in 1:4) {
      
      beta_draws <- post$beta_out[, n]
      
      plot_data <- rbind(
        plot_data,
        data.frame(
          species = species,
          model_type = model,
          parameter = paste0("beta[", n, "]"),
          mean = mean(beta_draws),
          lower = quantile(beta_draws, 0.025),
          upper = quantile(beta_draws, 0.975)
        )
      )
    }
    
    
    # Kappa
    
    for (n in 1:4) {
      
      kappa_draws <- post$kappa_out[, n]
      
      plot_data <- rbind(
        plot_data,
        data.frame(
          species = species,
          model_type = model,
          parameter = paste0("kappa[", n, "]"),
          mean = mean(kappa_draws),
          lower = quantile(kappa_draws, 0.025),
          upper = quantile(kappa_draws, 0.975)
        )
      )
    }
  }
}


# Model order

plot_data$model_type <- factor(
  plot_data$model_type,
  levels = c(
    "Model 1",
    "Model 2",
    "Model 3",
    "Model 4"
  )
)



# BETA PLOT


beta_plot <- plot_data |>
  filter(grepl("^beta", parameter)) |>
  ggplot(
    aes(
      x = mean,
      y = model_type
    )
  ) +
  geom_errorbarh(
    aes(
      xmin = lower,
      xmax = upper
    ),
    height = 0.15
  ) +
  geom_point(size = 2) +
  facet_grid(
    parameter ~ species,
    scales = "free_x",
    labeller = labeller(
      parameter = label_parsed
    )
  ) +
  labs(
    x = expression(
      paste(
        "Posterior mean and 95% credible interval for ",
        beta
      )
    ),
    y = NULL
  ) +
  theme_bw()

beta_plot



# KAPPA PLOT


kappa_plot <- plot_data |>
  filter(grepl("^kappa", parameter)) |>
  ggplot(
    aes(
      x = mean,
      y = model_type
    )
  ) +
  geom_errorbarh(
    aes(
      xmin = lower,
      xmax = upper
    ),
    height = 0.15
  ) +
  geom_point(size = 2) +
  facet_grid(
    parameter ~ species,
    scales = "free_x",
    labeller = labeller(
      parameter = label_parsed
    )
  ) +
  labs(
    x = expression(
      paste(
        "Posterior mean and 95% credible interval for ",
        kappa
      )
    ),
    y = NULL
  ) +
  theme_bw()

kappa_plot
