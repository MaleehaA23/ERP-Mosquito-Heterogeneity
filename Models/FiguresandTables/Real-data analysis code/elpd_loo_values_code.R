library(rstan)
library(loo)

fit_files <- list(
  Culex = list(
    `Model 1` = "CulexNewResults/Model1/Model1_fit.rds",
    `Model 2` = "CulexNewResults/Model2/Model2_fit.rds",
    `Model 3` = "CulexNewResults/Model3/Model3_fit.rds",
    `Model 4` = "CulexNewResults/Model4/Model4_fit.rds"
  ),
  
  Arabiensis = list(
    `Model 1` = "ArabiensisNewResults/Model1/Model1_fit.rds",
    `Model 2` = "ArabiensisNewResults/Model2/Model2_fit.rds",
    `Model 3` = "ArabiensisNewResults/Model3/Model3_fit.rds",
    `Model 4` = "ArabiensisNewResults/Model4/Model4_fit.rds"
  )
)

for (species in names(fit_files)) {
  
  for (model in names(fit_files[[species]])) {
    
    fit <- readRDS(
      fit_files[[species]][[model]]
    )
    
    log_lik <- rstan::extract(
      fit,
      pars = "log_lik"
    )$log_lik
    
    loo_result <- loo::loo(log_lik)
    
    cat(
      "\n",
      species,
      "-",
      model,
      "\n"
    )
    
    print(loo_result)
  }
}