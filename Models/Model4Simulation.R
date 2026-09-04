library(rstan)

Data0 <- read.csv("Data/arabiensis_cleanEH_BIT046_data0.csv")
Data1 <- read.csv("Data/arabiensis_cleanEH_BIT046_data1.csv")

mod1 <- stan_model("Model4.stan")

# Volunteer indexing

all_volunteers <- sort(
  unique(c(Data0$Volunteer, Data1$Volunteer))
)

volunteer_lookup <- setNames(
  seq_along(all_volunteers),
  all_volunteers
)

Data0$VolunteerID <- unname(
  volunteer_lookup[as.character(Data0$Volunteer)]
)

Data1$VolunteerID <- unname(
  volunteer_lookup[as.character(Data1$Volunteer)]
)

Nvolunteers <- length(all_volunteers)

# Posterior samples

alphaB_post <- read.csv(
  "ArabiensisNewResults/Model4/Post_alphaB.csv"
)[,1]

alphaM_post <- read.csv(
  "ArabiensisNewResults/Model4/Post_alphaM.csv"
)[,1]

sigma_alpha_post <- read.csv(
  "ArabiensisNewResults/Model4/Post_sigma_alpha.csv"
)[,1]

sigma_alpha2_post <- read.csv(
  "ArabiensisNewResults/Model4/Post_sigma_alpha2.csv"
)[,1]

sigma_vol_post <- read.csv(
  "ArabiensisNewResults/Model4/Post_sigma_vol.csv"
)[,1]

alpha_day_post <- as.matrix(
  read.csv(
    "ArabiensisNewResults/Model4/Post_alpha_day.csv",
    check.names = FALSE
  )
)

alpha2_day_post <- as.matrix(
  read.csv(
    "ArabiensisNewResults/Model4/Post_alpha2_day.csv",
    check.names = FALSE
  )
)

volunteer_effect_post <- as.matrix(
  read.csv(
    "ArabiensisNewResults/Model4/Post_volunteer_effect.csv",
    check.names = FALSE
  )
)

beta_post <- as.matrix(
  read.csv(
    "ArabiensisNewResults/Model4/Post_beta.csv",
    check.names = FALSE
  )
)

kappa_post <- as.matrix(
  read.csv(
    "ArabiensisNewResults/Model4/Post_kappa.csv",
    check.names = FALSE
  )
)

phiC_post <- read.csv(
  "ArabiensisNewResults/Model4/Post_phiC.csv"
)[,1]

phiI_post <- as.matrix(
  read.csv(
    "ArabiensisNewResults/Model4/Post_phiI.csv",
    check.names = FALSE
  )
)

Nnets <- ncol(beta_post)
Ndays <- ncol(alpha_day_post)

nsim <- 100

# Storage

lower_alphaB <- upper_alphaB <- numeric(nsim)
lower_alphaM <- upper_alphaM <- numeric(nsim)

lower_beta <- upper_beta <- matrix(NA, nsim, Nnets)
lower_kappa <- upper_kappa <- matrix(NA, nsim, Nnets)

lower_phiC <- upper_phiC <- numeric(nsim)
lower_phiI <- upper_phiI <- matrix(NA, nsim, Nnets)

lower_sigma_alpha <- upper_sigma_alpha <- numeric(nsim)
lower_sigma_alpha2 <- upper_sigma_alpha2 <- numeric(nsim)
lower_sigma_vol <- upper_sigma_vol <- numeric(nsim)

true_alphaB <- numeric(nsim)
est_alphaB <- numeric(nsim)

true_alphaM <- numeric(nsim)
est_alphaM <- numeric(nsim)

true_sigma_alpha <- numeric(nsim)
est_sigma_alpha <- numeric(nsim)

true_sigma_alpha2 <- numeric(nsim)
est_sigma_alpha2 <- numeric(nsim)

true_sigma_vol <- numeric(nsim)
est_sigma_vol <- numeric(nsim)

true_beta <- matrix(NA, nsim, Nnets)
est_beta <- matrix(NA, nsim, Nnets)

true_kappa <- matrix(NA, nsim, Nnets)
est_kappa <- matrix(NA, nsim, Nnets)

true_phiC <- numeric(nsim)
est_phiC <- numeric(nsim)

true_phiI <- matrix(NA, nsim, Nnets)
est_phiI <- matrix(NA, nsim, Nnets)

set.seed(123)

for(sim in 1:nsim){
  
  draw <- sample(seq_along(alphaB_post), 1)
  
  # True parameter values
  
  alphaB <- alphaB_post[draw]
  alphaM <- alphaM_post[draw]
  
  sigma_alpha <- sigma_alpha_post[draw]
  sigma_alpha2 <- sigma_alpha2_post[draw]
  sigma_vol <- sigma_vol_post[draw]
  
  alpha_day <- as.numeric(
    alpha_day_post[draw, ]
  )
  
  alpha2_day <- as.numeric(
    alpha2_day_post[draw, ]
  )
  
  volunteer_effect <- as.numeric(
    volunteer_effect_post[draw, ]
  )
  
  beta <- as.numeric(
    beta_post[draw, ]
  )
  
  kappa <- as.numeric(
    kappa_post[draw, ]
  )
  
  phiC_true <- phiC_post[draw]
  
  phiI_true <- as.numeric(
    phiI_post[draw, ]
  )
  
  # Store true values
  
  true_alphaB[sim] <- alphaB
  true_alphaM[sim] <- alphaM
  
  true_sigma_alpha[sim] <- sigma_alpha
  true_sigma_alpha2[sim] <- sigma_alpha2
  true_sigma_vol[sim] <- sigma_vol
  
  true_beta[sim, ] <- beta
  true_kappa[sim, ] <- kappa
  
  true_phiC[sim] <- phiC_true
  true_phiI[sim, ] <- phiI_true
  
  # Copy original data
  
  Data0_sim <- Data0
  Data1_sim <- Data1
  
  # Control
  
  Data0_sim$FA <- 0
  Data0_sim$FD <- 0
  Data0_sim$UA <- 0
  Data0_sim$UD <- 0
  
  for(i in 1:nrow(Data0_sim)){
    
    day <- Data0_sim$Day[i]
    vol <- Data0_sim$VolunteerID[i]
    
    alphaB_day_vol <-
      alpha_day[day] *
      volunteer_effect[vol]
    
    alphaM_day <- alpha2_day[day]
    
    mu <- alphaB_day_vol + alphaM_day
    
    q <- alphaB_day_vol /
      (alphaB_day_vol + alphaM_day)
    
    total <- rnbinom(
      1,
      mu = mu,
      size = phiC_true
    )
    
    fed <- rbinom(
      1,
      size = total,
      prob = q
    )
    
    Data0_sim$FA[i] <- fed
    Data0_sim$FD[i] <- 0
    Data0_sim$UD[i] <- total - fed
    Data0_sim$UA[i] <- 0
  }
  
  # Intervention
  
  Data1_sim$FA <- 0
  Data1_sim$FD <- 0
  Data1_sim$UA <- 0
  Data1_sim$UD <- 0
  
  for(i in 1:nrow(Data1_sim)){
    
    net <- Data1_sim$Net[i]
    day <- Data1_sim$Day[i]
    vol <- Data1_sim$VolunteerID[i]
    
    alphaB_day_vol <-
      alpha_day[day] *
      volunteer_effect[vol]
    
    alphaM_day <- alpha2_day[day]
    
    feed <- alphaB_day_vol * beta[net]
    kill <- alphaB_day_vol * kappa[net]
    natural <- alphaM_day
    
    mu <- feed + kill + natural
    p <- feed / mu
    
    total <- rnbinom(
      1,
      mu = mu,
      size = phiI_true[net]
    )
    
    fed <- rbinom(
      1,
      size = total,
      prob = p
    )
    
    Data1_sim$FA[i] <- fed
    Data1_sim$FD[i] <- 0
    Data1_sim$UD[i] <- total - fed
    Data1_sim$UA[i] <- 0
  }
  
  # Weights
  
  weight_count <- sapply(
    seq_len(Ndays),
    function(t){
      sum(
        Data0_sim$FA[Data0_sim$Day == t] +
          Data0_sim$FD[Data0_sim$Day == t] +
          Data0_sim$UD[Data0_sim$Day == t]
      )
    }
  )
  
  # Stan data
  
  mod1_data <- list(
    n0 = nrow(Data0_sim),
    n1 = nrow(Data1_sim),
    Nnets = Nnets,
    Ndays = Ndays,
    Nvolunteers = Nvolunteers,
    Net = Data1_sim$Net,
    Day0 = Data0_sim$Day,
    Day1 = Data1_sim$Day,
    Volunteer0 = Data0_sim$VolunteerID,
    Volunteer1 = Data1_sim$VolunteerID,
    fed0 = Data0_sim$FA + Data0_sim$FD,
    total0 = Data0_sim$FA + Data0_sim$FD + Data0_sim$UD,
    fed1 = Data1_sim$FA + Data1_sim$FD,
    total1 = Data1_sim$FA + Data1_sim$FD + Data1_sim$UD,
    priorscale = 1,
    hierarchy = 1,
    samples0 = rep(1, nrow(Data0_sim)),
    samples1 = rep(1, nrow(Data1_sim)),
    weight_count = as.vector(weight_count)
  )
  
  # Fit model
  
  mod1_fit <- sampling(
    mod1,
    data = mod1_data,
    iter = 6000,
    chains = 4,
    cores = 4,
    control = list(
      max_treedepth = 12,
      adapt_delta = 0.95
    )
  )
  
  posterior <- as.matrix(mod1_fit)
  
  # Posterior means
  
  est_alphaB[sim] <- mean(
    posterior[, "alpha_out"]
  )
  
  est_alphaM[sim] <- mean(
    posterior[, "alpha2_out"]
  )
  
  est_sigma_alpha[sim] <- mean(
    posterior[, "sigma_alpha"]
  )
  
  est_sigma_alpha2[sim] <- mean(
    posterior[, "sigma_alpha2"]
  )
  
  est_sigma_vol[sim] <- mean(
    posterior[, "sigma_vol"]
  )
  
  beta_out <- posterior[
    ,
    grep("^beta_out\\[", colnames(posterior)),
    drop = FALSE
  ]
  
  kappa_out <- posterior[
    ,
    grep("^kappa_out\\[", colnames(posterior)),
    drop = FALSE
  ]
  
  phiI_out <- posterior[
    ,
    grep("^phi_valuesI\\[", colnames(posterior)),
    drop = FALSE
  ]
  
  est_beta[sim, ] <- apply(
    beta_out,
    2,
    mean
  )
  
  est_kappa[sim, ] <- apply(
    kappa_out,
    2,
    mean
  )
  
  est_phiC[sim] <- mean(
    posterior[, "phi_valuesC"]
  )
  
  est_phiI[sim, ] <- apply(
    phiI_out,
    2,
    mean
  )
  
  # 95% credible intervals
  
  lower_alphaB[sim] <- quantile(
    posterior[, "alpha_out"],
    0.025
  )
  
  upper_alphaB[sim] <- quantile(
    posterior[, "alpha_out"],
    0.975
  )
  
  lower_alphaM[sim] <- quantile(
    posterior[, "alpha2_out"],
    0.025
  )
  
  upper_alphaM[sim] <- quantile(
    posterior[, "alpha2_out"],
    0.975
  )
  
  lower_sigma_alpha[sim] <- quantile(
    posterior[, "sigma_alpha"],
    0.025
  )
  
  upper_sigma_alpha[sim] <- quantile(
    posterior[, "sigma_alpha"],
    0.975
  )
  
  lower_sigma_alpha2[sim] <- quantile(
    posterior[, "sigma_alpha2"],
    0.025
  )
  
  upper_sigma_alpha2[sim] <- quantile(
    posterior[, "sigma_alpha2"],
    0.975
  )
  
  lower_sigma_vol[sim] <- quantile(
    posterior[, "sigma_vol"],
    0.025
  )
  
  upper_sigma_vol[sim] <- quantile(
    posterior[, "sigma_vol"],
    0.975
  )
  
  lower_beta[sim, ] <- apply(
    beta_out,
    2,
    quantile,
    probs = 0.025
  )
  
  upper_beta[sim, ] <- apply(
    beta_out,
    2,
    quantile,
    probs = 0.975
  )
  
  lower_kappa[sim, ] <- apply(
    kappa_out,
    2,
    quantile,
    probs = 0.025
  )
  
  upper_kappa[sim, ] <- apply(
    kappa_out,
    2,
    quantile,
    probs = 0.975
  )
  
  lower_phiC[sim] <- quantile(
    posterior[, "phi_valuesC"],
    0.025
  )
  
  upper_phiC[sim] <- quantile(
    posterior[, "phi_valuesC"],
    0.975
  )
  
  lower_phiI[sim, ] <- apply(
    phiI_out,
    2,
    quantile,
    probs = 0.025
  )
  
  upper_phiI[sim, ] <- apply(
    phiI_out,
    2,
    quantile,
    probs = 0.975
  )
}

# Recovery correlations

cor(true_alphaB, est_alphaB, use = "complete.obs")
cor(true_alphaM, est_alphaM, use = "complete.obs")

cor(
  true_sigma_alpha,
  est_sigma_alpha,
  use = "complete.obs"
)

cor(
  true_sigma_alpha2,
  est_sigma_alpha2,
  use = "complete.obs"
)

cor(
  true_sigma_vol,
  est_sigma_vol,
  use = "complete.obs"
)

for(n in 1:Nnets){
  cor(
    true_beta[,n],
    est_beta[,n],
    use = "complete.obs"
  )
}

for(n in 1:Nnets){
  cor(
    true_kappa[,n],
    est_kappa[,n],
    use = "complete.obs"
  )
}

cor(
  true_phiC,
  est_phiC,
  use = "complete.obs"
)

for(n in 1:Nnets){
  cor(
    true_phiI[,n],
    est_phiI[,n],
    use = "complete.obs"
  )
}

# Save results

results <- list(
  true_alphaB = true_alphaB,
  est_alphaB = est_alphaB,
  lower_alphaB = lower_alphaB,
  upper_alphaB = upper_alphaB,
  
  true_alphaM = true_alphaM,
  est_alphaM = est_alphaM,
  lower_alphaM = lower_alphaM,
  upper_alphaM = upper_alphaM,
  
  true_sigma_alpha = true_sigma_alpha,
  est_sigma_alpha = est_sigma_alpha,
  lower_sigma_alpha = lower_sigma_alpha,
  upper_sigma_alpha = upper_sigma_alpha,
  
  true_sigma_alpha2 = true_sigma_alpha2,
  est_sigma_alpha2 = est_sigma_alpha2,
  lower_sigma_alpha2 = lower_sigma_alpha2,
  upper_sigma_alpha2 = upper_sigma_alpha2,
  
  true_sigma_vol = true_sigma_vol,
  est_sigma_vol = est_sigma_vol,
  lower_sigma_vol = lower_sigma_vol,
  upper_sigma_vol = upper_sigma_vol,
  
  true_beta = true_beta,
  est_beta = est_beta,
  lower_beta = lower_beta,
  upper_beta = upper_beta,
  
  true_kappa = true_kappa,
  est_kappa = est_kappa,
  lower_kappa = lower_kappa,
  upper_kappa = upper_kappa,
  
  true_phiC = true_phiC,
  est_phiC = est_phiC,
  lower_phiC = lower_phiC,
  upper_phiC = upper_phiC,
  
  true_phiI = true_phiI,
  est_phiI = est_phiI,
  lower_phiI = lower_phiI,
  upper_phiI = upper_phiI
)

dir.create(
  "SimulationResults",
  showWarnings = FALSE
)

saveRDS(
  results,
  "SimulationResults/Model4ParameterRecoveryResults.rds"
)