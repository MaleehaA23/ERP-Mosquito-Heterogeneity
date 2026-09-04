# Adaptation of FitEH_minimal_large9.stan taken from Fairbanks et al. experimental hut model
#NOTE: For Culex, change lines 82 and 83 to 15 and 0.99 respectively.
library(rstan)

mod1 <- stan_model("Model4.stan")

species <- "Arabiensis" # Change to Culex for Culex datasets
model_name <- "Model4"

Data0 <- read.csv("Data/arabiensis_cleanEH_BIT046_data0.csv")
Data1 <- read.csv("Data/arabiensis_cleanEH_BIT046_data1.csv")

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

# Control weights

weight_count <- c()

for(t in 1:length(unique(Data0$Day))) {
  weight_count <- cbind(
    weight_count,
    sum(
      Data0$FA[Data0$Day == t] +
        Data0$FD[Data0$Day == t] +
        Data0$UD[Data0$Day == t]
    )
  )
}

# Stan data

mod1_data <- list(
  n0 = nrow(Data0),
  n1 = nrow(Data1),
  Nnets = length(unique(Data1$Net)),
  Ndays = length(unique(Data1$Day)),
  Nvolunteers = Nvolunteers,
  Net = Data1$Net,
  Day0 = Data0$Day,
  Day1 = Data1$Day,
  Volunteer0 = Data0$VolunteerID,
  Volunteer1 = Data1$VolunteerID,
  fed0 = Data0$FA + Data0$FD,
  total0 = Data0$FA + Data0$FD + Data0$UD,
  fed1 = Data1$FA + Data1$FD,
  total1 = Data1$FA + Data1$FD + Data1$UD,
  priorscale = 1,
  hierarchy = 1,
  samples0 = rep(1, nrow(Data0)),
  samples1 = rep(1, nrow(Data1)),
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
    max_treedepth = 12,   #for Culex change to 15
    adapt_delta = 0.95    #for Culex change to 0.99
  )
)

posterior1 <- as.matrix(mod1_fit)

# Extract parameters

extract_cols <- function(post, pattern) {
  post[
    ,
    grep(pattern, colnames(post)),
    drop = FALSE
  ]
}

alphaB <- posterior1[, "alpha_out"]
alphaM <- posterior1[, "alpha2_out"]

sigma_alpha <- posterior1[, "sigma_alpha"]
sigma_alpha2 <- posterior1[, "sigma_alpha2"]

alpha_day <- extract_cols(
  posterior1,
  "^alpha_day_out\\["
)

alpha2_day <- extract_cols(
  posterior1,
  "^alpha2_day_out\\["
)

sigma_vol <- posterior1[, "sigma_vol"]

volunteer_effect <- extract_cols(
  posterior1,
  "^volunteer_effect_out\\["
)

beta_out <- extract_cols(
  posterior1,
  "^beta_out\\["
)

kappa_out <- extract_cols(
  posterior1,
  "^kappa_out\\["
)

phiC <- posterior1[, "phi_valuesC"]

phiI <- extract_cols(
  posterior1,
  "^phi_valuesI\\["
)

# Save posterior samples

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
  file.path(output_dir, "Model4_fit.rds")
)
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
  data.frame(sigma_vol = sigma_vol),
  file.path(output_dir, "Post_sigma_vol.csv"),
  row.names = FALSE
)

write.csv(
  as.data.frame(volunteer_effect),
  file.path(output_dir, "Post_volunteer_effect.csv"),
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