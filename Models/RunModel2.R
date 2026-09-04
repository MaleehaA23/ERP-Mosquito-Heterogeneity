#Adaptation of FitEH_minimal_large9.stan taken from Fairbanks et al. experimental hut model
library(rstan)

options(mc.cores = parallel::detectCores())

mod1 <- stan_model("Model2.stan")

species <- "Arabiensis"  #Change to Culex for Culex datasets
model_name <- "Model2"

# Load control and intervention datasets. Either Arabiensis or Culex.
Data0 <- read.csv("Data/arabiensis_cleanEH_BIT046_data0.csv")
Data1 <- read.csv("Data/arabiensis_cleanEH_BIT046_data1.csv")

weight_count = c()
for(t in 1:length(unique(Data0$Day))){
  weight_count = cbind(weight_count,
                       sum(Data0$FA[which(Data0$Day == t)] +
                             Data0$FD[which(Data0$Day == t)] +
                             Data0$UD[which(Data0$Day == t)]))
}

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

mod1_fit <- sampling(mod1, data = mod1_data,
                     iter = 6000, chains = 4, cores = 4,
                     control = list(max_treedepth = 12))

posterior1 <- as.matrix(mod1_fit)


extract_cols <- function(post, pattern) {
  post[, grep(pattern, colnames(post)), drop = FALSE]
}


alphaB <- posterior1[, "alpha_out"]
alphaM <- posterior1[, "alpha2_out"]

sigma_alpha  <- posterior1[, "sigma_alpha"]  # Save the day effects parameters
sigma_alpha2 <- posterior1[, "sigma_alpha2"] # Save the day effects parameters

alpha_day <- extract_cols(posterior1, "^alpha_day_out\\[")  #Save the day effects parameters
alpha2_day <- extract_cols(posterior1, "^alpha2_day_out\\[") #Save the day effects parameters


beta_out  <- extract_cols(posterior1, "^beta_out\\[([0-9]+)\\]$")
kappa_out <- extract_cols(posterior1, "^kappa_out\\[([0-9]+)\\]$")

phiC <- posterior1[, "phi_valuesC"]
phiI <- extract_cols(posterior1, "^phi_valuesI\\[")

# Create labels for the nets in the dataset
Nnets <- ncol(beta_out)
Net_names <- paste("Net", seq_len(Nnets))

# Create output directory
# For this run: ArabiensisNewResults/Model2
output_dir <- file.path(
  paste0(species, "NewResults"),
  model_name
)

dir.create(
  output_dir,
  showWarnings = FALSE,
  recursive = TRUE
)
saveRDS(
  mod1_fit,
  file.path(output_dir, "Model2_fit.rds")
)

# Save posterior draws
write.csv(
  data.frame(alphaB = alphaB),
  file.path(output_dir, "Post_alphaB.csv"),
  row.names = FALSE
)

write.csv(
  data.frame(alphaM = alphaM),
  file.path(output_dir, "Post_alphaM.csv"),
  row.names = FALSE
)

write.csv(
  data.frame(sigma_alpha = sigma_alpha),
  file.path(output_dir, "Post_sigma_alpha.csv"),
  row.names = FALSE
)

write.csv(
  data.frame(sigma_alpha2 = sigma_alpha2),
  file.path(output_dir, "Post_sigma_alpha2.csv"),
  row.names = FALSE
)

write.csv(
  as.data.frame(alpha_day),
  file.path(output_dir, "Post_alpha_day.csv"),
  row.names = FALSE
)

write.csv(
  as.data.frame(alpha2_day),
  file.path(output_dir, "Post_alpha2_day.csv"),
  row.names = FALSE
)

write.csv(
  as.data.frame(beta_out),
  file.path(output_dir, "Post_beta.csv"),
  row.names = FALSE
)

write.csv(
  as.data.frame(kappa_out),
  file.path(output_dir, "Post_kappa.csv"),
  row.names = FALSE
)


# Save dispersion parameters
write.csv(
  data.frame(phiC = phiC),
  file.path(output_dir, "Post_phiC.csv"),
  row.names = FALSE
)

write.csv(
  as.data.frame(phiI),
  file.path(output_dir, "Post_phiI.csv"),
  row.names = FALSE
)

cat("\nPosterior samples saved to",
    paste0(output_dir, "/"), "\n")