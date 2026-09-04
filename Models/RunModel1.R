#Adaptation of FitEH_minimal_large9.stan taken from Fairbanks et al. experimental hut model

library(rstan)

options(mc.cores = parallel::detectCores())

# Compile Stan model
mod1 <- stan_model("Model1.stan")

species <- "Arabiensis"  #Change to Culex for Culex datasets
model_name <- "Model1"

# Load control and intervention datasets. Either Arabiensis or Culex.
Data0 <- read.csv("Data/arabiensis_cleanEH_BIT046_data0.csv")
Data1 <- read.csv("Data/arabiensis_cleanEH_BIT046_data1.csv")


# Calculate total mosquito count for each day in the control dataset
weight_count <- c()

for(t in 1:length(unique(Data0$Day))){
  weight_count <- cbind(
    weight_count,
    sum(Data0$FA[which(Data0$Day == t)] +
          Data0$FD[which(Data0$Day == t)] +
          Data0$UD[which(Data0$Day == t)])
  )
}


# Prepare data for the Stan model
mod1_data <- list(
  n0     = dim(Data0)[1],
  n1     = dim(Data1)[1],
  Nnets  = length(unique(Data1$Net)),
  Ndays  = length(unique(Data1$Day)),
  Net    = Data1$Net,
  Day0   = Data0$Day,
  Day1   = Data1$Day,
  fed0   = Data0$FA + Data0$FD,
  total0 = Data0$FA + Data0$FD + Data0$UD,
  fed1   = Data1$FA + Data1$FD,
  total1 = Data1$FA + Data1$FD + Data1$UD,
  priorscale   = 1,
  hierarchy    = 1,
  samples0     = rep(1, dim(Data0)[1]),
  samples1     = rep(1, dim(Data1)[1]),
  weight_count = as.vector(weight_count)
)


# Fit the Bayesian model
mod1_fit <- sampling(
  mod1,
  data = mod1_data,
  iter = 6000,
  chains = 4,
  cores = 4,
  control = list(max_treedepth = 12)
)


# Convert posterior samples to a matrix
posterior1 <- as.matrix(mod1_fit)


# Extract by parameter name 
extract_cols <- function(post, pattern) {
  post[, grep(pattern, colnames(post)), drop = FALSE]
}


alphaB <- posterior1[, "alpha_out"]
alphaM <- posterior1[, "alpha2_out"]

# beta and kappa are vectors over nets only
beta_out  <- extract_cols(posterior1, "^beta_out\\[([0-9]+)\\]$")
kappa_out <- extract_cols(posterior1, "^kappa_out\\[([0-9]+)\\]$")

#create the output directory: for this it will be ArabiensisNewResults/Model1
output_dir <- file.path(
  paste0(species, "NewResults"),
  model_name
)

dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

saveRDS(
  mod1_fit,
  file.path(output_dir, "Model1_fit.rds")
)
# Save posterior draws
write.csv(alphaB,
          file.path(output_dir, "Post_alphaB.csv"))

write.csv(alphaM,
          file.path(output_dir, "Post_alphaM.csv"))

write.csv(beta_out,
          file.path(output_dir, "Post_beta.csv"),
          row.names = FALSE)

write.csv(kappa_out,
          file.path(output_dir, "Post_kappa.csv"),
          row.names = FALSE)

# Create labels for the nets in the dataset
Nnets <- ncol(beta_out)
Net_names <- paste("Net", seq_len(Nnets))

# Save beta and kappa separately for each net
for (i in seq_len(Nnets)) {
  
  write.csv(
    beta_out[, i],
    file.path(output_dir, paste0("Post_Beta_Net_", i, ".csv")),
    row.names = FALSE
  )
  
  write.csv(
    kappa_out[, i],
    file.path(output_dir, paste0("Post_Kappa_Net_", i, ".csv")),
    row.names = FALSE
  )
}

# Save dispersion parameter posteriors
phiC <- posterior1[, "phi_valuesC"]
phiI <- extract_cols(posterior1, "^phi_valuesI\\[")

write.csv(phiC,
          file.path(output_dir, "Post_phiC.csv"),
          row.names = FALSE)

write.csv(phiI,
          file.path(output_dir, "Post_phiI.csv"),
          row.names = FALSE)