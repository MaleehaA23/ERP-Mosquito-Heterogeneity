library(rstan)


Data0 <- read.csv("Data/arabiensis_cleanEH_BIT046_data0.csv")
Data1 <- read.csv("Data/arabiensis_cleanEH_BIT046_data1.csv")

mod1 <- stan_model("Model1.stan")

## Posterior samples
alphaB_post <- read.csv("ArabiensisNewResults/Model1/Post_alphaB.csv")[,2]
alphaM_post <- read.csv("ArabiensisNewResults/Model1/Post_alphaM.csv")[,2]

beta1 <- read.csv("ArabiensisNewResults/Model1/Post_Beta_Net_1.csv")[,1]
beta2 <- read.csv("ArabiensisNewResults/Model1/Post_Beta_Net_2.csv")[,1]
beta3 <- read.csv("ArabiensisNewResults/Model1/Post_Beta_Net_3.csv")[,1]
beta4 <- read.csv("ArabiensisNewResults/Model1/Post_Beta_Net_4.csv")[,1]

kappa1 <- read.csv("ArabiensisNewResults/Model1/Post_Kappa_Net_1.csv")[,1]
kappa2 <- read.csv("ArabiensisNewResults/Model1/Post_Kappa_Net_2.csv")[,1]
kappa3 <- read.csv("ArabiensisNewResults/Model1/Post_Kappa_Net_3.csv")[,1]
kappa4 <- read.csv("ArabiensisNewResults/Model1/Post_Kappa_Net_4.csv")[,1]

phiC <- read.csv("ArabiensisNewResults/Model1/Post_phiC.csv")[,1]

## Read ALL intervention dispersions
phiI <- read.csv("ArabiensisNewResults/Model1/Post_phiI.csv")

nsim <- 100
lower_alphaB <- upper_alphaB <- numeric(nsim)
lower_alphaM <- upper_alphaM <- numeric(nsim)

lower_beta  <- upper_beta  <- matrix(NA, nsim, 4)
lower_kappa <- upper_kappa <- matrix(NA, nsim, 4)

lower_phiC <- upper_phiC <- numeric(nsim)
lower_phiI <- upper_phiI <- matrix(NA, nsim, 4)

# ----------------------------
# Parameter recovery
# ----------------------------



true_alphaB <- numeric(nsim)
est_alphaB  <- numeric(nsim)

true_alphaM <- numeric(nsim)
est_alphaM  <- numeric(nsim)

true_beta  <- matrix(NA, nsim, 4)
est_beta   <- matrix(NA, nsim, 4)

true_kappa <- matrix(NA, nsim, 4)
est_kappa  <- matrix(NA, nsim, 4)

true_phiC <- numeric(nsim)
est_phiC  <- numeric(nsim)

true_phiI <- matrix(NA, nsim, 4)
est_phiI  <- matrix(NA, nsim, 4)
set.seed(123)

for(sim in 1:nsim){
  
  cat("Simulation", sim, "of", nsim, "\n")
  draw <- sample(seq_along(alphaB_post), 1)
  
  ## True parameter values
  alphaB <- alphaB_post[draw]
  alphaM <- alphaM_post[draw]
  
  beta <- c(beta1[draw],
            beta2[draw],
            beta3[draw],
            beta4[draw])
  
  kappa <- c(kappa1[draw],
             kappa2[draw],
             kappa3[draw],
             kappa4[draw])
  
  phiC_true <- phiC[draw]
  phiI_true <- as.numeric(phiI[draw, ])
  
  ## Copy original data
  Data0_sim <- Data0
  Data1_sim <- Data1
  
  ## ------------------------
  ## CONTROL
  ## ------------------------
  
  Data0_sim$FA <- 0
  Data0_sim$FD <- 0
  Data0_sim$UA <- 0
  Data0_sim$UD <- 0
  
  for(i in 1:nrow(Data0_sim)){
    
    mu <- alphaB + alphaM
    q  <- alphaB/(alphaB + alphaM)
    
    total <- rnbinom(1,
                     mu = mu,
                     size = phiC_true)
    
    fed <- rbinom(1,
                  size = total,
                  prob = q)
    
    Data0_sim$FA[i] <- fed
    Data0_sim$FD[i] <- 0
    Data0_sim$UD[i] <- total - fed
    Data0_sim$UA[i] <- 0
  }
  

  ## INTERVENTION

  
  Data1_sim$FA <- 0
  Data1_sim$FD <- 0
  Data1_sim$UA <- 0
  Data1_sim$UD <- 0
  
  for(i in 1:nrow(Data1_sim)){
    
    net <- Data1_sim$Net[i]
    
    feed <- alphaB * beta[net]
    kill <- alphaB * kappa[net]
    natural <- alphaM
    
    mu <- feed + kill + natural
    p <- feed / mu
    
    total <- rnbinom(1,
                     mu = mu,
                     size = phiI_true[net])
    
    fed <- rbinom(1,
                  size = total,
                  prob = p)
    
    Data1_sim$FA[i] <- fed
    Data1_sim$FD[i] <- 0
    Data1_sim$UD[i] <- total - fed
    Data1_sim$UA[i] <- 0
  }
  
  ## Store true values
  true_alphaB[sim] <- alphaB
  true_alphaM[sim] <- alphaM
  true_beta[sim, ] <- beta
  true_kappa[sim, ] <- kappa
  true_phiC[sim] <- phiC_true
  true_phiI[sim, ] <- phiI_true
  
  ## ------------------------
  ## Stan data
  ## ------------------------
  
  weight_count <- c()
  
  for(t in unique(Data0_sim$Day)){
    weight_count <- c(
      weight_count,
      sum(Data0_sim$FA[Data0_sim$Day == t] +
            Data0_sim$FD[Data0_sim$Day == t] +
            Data0_sim$UD[Data0_sim$Day == t])
    )
  }
  
  mod1_data <- list(
    n0 = nrow(Data0_sim),
    n1 = nrow(Data1_sim),
    Nnets = length(unique(Data1_sim$Net)),
    Ndays = length(unique(Data1_sim$Day)),
    Net = Data1_sim$Net,
    Day0 = Data0_sim$Day,
    Day1 = Data1_sim$Day,
    fed0 = Data0_sim$FA + Data0_sim$FD,
    total0 = Data0_sim$FA + Data0_sim$FD + Data0_sim$UD,
    fed1 = Data1_sim$FA + Data1_sim$FD,
    total1 = Data1_sim$FA + Data1_sim$FD + Data1_sim$UD,
    priorscale = 1,
    hierarchy = 1,
    samples0 = rep(1, nrow(Data0_sim)),
    samples1 = rep(1, nrow(Data1_sim)),
    weight_count = weight_count
  )
  
  ## ------------------------
  ## Fit model
  ## ------------------------
  cat("  Fitting Stan model...\n")
  flush.console()
  
  mod1_fit <- sampling(
    mod1,
    data = mod1_data,
    iter = 6000,      
    chains = 4,       
    cores = 4
  )
  
  posterior <- as.matrix(mod1_fit)
  
  ## ------------------------
  ## Posterior means
  ## ------------------------
  
  est_alphaB[sim] <- mean(posterior[, "alpha_out"])
  est_alphaM[sim] <- mean(posterior[, "alpha2_out"])
  
  beta_out <- posterior[, grep("^beta_out\\[", colnames(posterior))]
  kappa_out <- posterior[, grep("^kappa_out\\[", colnames(posterior))]
  phiI_out <- posterior[, grep("^phi_valuesI\\[", colnames(posterior))]
  
  est_beta[sim, ] <- apply(beta_out, 2, mean)
  est_kappa[sim, ] <- apply(kappa_out, 2, mean)
  
  est_phiC[sim] <- mean(posterior[, "phi_valuesC"])
  est_phiI[sim, ] <- apply(phiI_out, 2, mean)
  

  lower_alphaB[sim] <- quantile(posterior[, "alpha_out"], 0.025)
  upper_alphaB[sim] <- quantile(posterior[, "alpha_out"], 0.975)
  
  lower_alphaM[sim] <- quantile(posterior[, "alpha2_out"], 0.025)
  upper_alphaM[sim] <- quantile(posterior[, "alpha2_out"], 0.975)
  
  lower_beta[sim, ]  <- apply(beta_out, 2, quantile, probs=0.025)
  upper_beta[sim, ]  <- apply(beta_out, 2, quantile, probs=0.975)
  
  lower_kappa[sim, ] <- apply(kappa_out, 2, quantile, probs=0.025)
  upper_kappa[sim, ] <- apply(kappa_out, 2, quantile, probs=0.975)
  
  lower_phiC[sim] <- quantile(posterior[, "phi_valuesC"], 0.025)
  upper_phiC[sim] <- quantile(posterior[, "phi_valuesC"], 0.975)
  
  lower_phiI[sim, ] <- apply(phiI_out, 2, quantile, probs=0.025)
  upper_phiI[sim, ] <- apply(phiI_out, 2, quantile, probs=0.975)
}

# Recovery correlations


cat("PARAMETER RECOVERY")

cat("alphaB:", cor(true_alphaB, est_alphaB, use = "complete.obs"), "\n")
cat("alphaM:", cor(true_alphaM, est_alphaM, use = "complete.obs"), "\n\n")

for(i in 1:4){
  cat("beta", i, ":",
      cor(true_beta[, i], est_beta[, i], use = "complete.obs"),
      "\n")
}

cat("\n")

for(i in 1:4){
  cat("kappa", i, ":",
      cor(true_kappa[, i], est_kappa[, i], use = "complete.obs"),
      "\n")
}

cat("\n")

cat("phiC:",
    cor(true_phiC, est_phiC, use = "complete.obs"),
    "\n\n")

for(i in 1:4){
  cat("phiI", i, ":",
      cor(true_phiI[, i], est_phiI[, i], use = "complete.obs"),
      "\n")
}
results <- list(
  true_alphaB = true_alphaB,
  est_alphaB = est_alphaB,
  true_alphaM = true_alphaM,
  est_alphaM = est_alphaM,
  true_beta = true_beta,
  est_beta = est_beta,
  true_kappa = true_kappa,
  est_kappa = est_kappa,
  true_phiC = true_phiC,
  est_phiC = est_phiC,
  true_phiI = true_phiI,
  est_phiI = est_phiI,
  lower_alphaB = lower_alphaB,
  upper_alphaB = upper_alphaB,
  lower_alphaM = lower_alphaM,
  upper_alphaM = upper_alphaM,
  
  lower_beta = lower_beta,
  upper_beta = upper_beta,
  
  lower_kappa = lower_kappa,
  upper_kappa = upper_kappa,
  
  lower_phiC = lower_phiC,
  upper_phiC = upper_phiC,
  
  lower_phiI = lower_phiI,
  upper_phiI = upper_phiI
)

dir.create("SimulationResults", showWarnings = FALSE)

saveRDS(
  results,
  "SimulationResults/Model1ParameterRecoveryResults.rds"
)