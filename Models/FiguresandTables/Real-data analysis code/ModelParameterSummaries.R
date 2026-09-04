library(rstan)


# LOAD FITS

A1 <- readRDS("ArabiensisNewResults/Model1/Model1_fit.rds")
A2 <- readRDS("ArabiensisNewResults/Model2/Model2_fit.rds")
A3 <- readRDS("ArabiensisNewResults/Model3/Model3_fit.rds")
A4 <- readRDS("ArabiensisNewResults/Model4/Model4_fit.rds")

C1 <- readRDS("CulexNewResults/Model1/Model1_fit.rds")
C2 <- readRDS("CulexNewResults/Model2/Model2_fit.rds")
C3 <- readRDS("CulexNewResults/Model3/Model3_fit.rds")
C4 <- readRDS("CulexNewResults/Model4/Model4_fit.rds")


# Extract posterior draws

A1_post <- rstan::extract(A1)
A2_post <- rstan::extract(A2)
A3_post <- rstan::extract(A3)
A4_post <- rstan::extract(A4)

C1_post <- rstan::extract(C1)
C2_post <- rstan::extract(C2)
C3_post <- rstan::extract(C3)
C4_post <- rstan::extract(C4)


# PHI VALUES

# Control phi posterior summaries

mean(A1_post$phi_valuesC)
quantile(A1_post$phi_valuesC, c(0.025, 0.975))

mean(A2_post$phi_valuesC)
quantile(A2_post$phi_valuesC, c(0.025, 0.975))

mean(A3_post$phi_valuesC)
quantile(A3_post$phi_valuesC, c(0.025, 0.975))

mean(A4_post$phi_valuesC)
quantile(A4_post$phi_valuesC, c(0.025, 0.975))

mean(C1_post$phi_valuesC)
quantile(C1_post$phi_valuesC, c(0.025, 0.975))

mean(C2_post$phi_valuesC)
quantile(C2_post$phi_valuesC, c(0.025, 0.975))

mean(C3_post$phi_valuesC)
quantile(C3_post$phi_valuesC, c(0.025, 0.975))

mean(C4_post$phi_valuesC)
quantile(C4_post$phi_valuesC, c(0.025, 0.975))


# Intervention phi posterior means
# One value for each of the four nets

apply(A1_post$phi_valuesI, 2, mean)
apply(A2_post$phi_valuesI, 2, mean)
apply(A3_post$phi_valuesI, 2, mean)
apply(A4_post$phi_valuesI, 2, mean)

apply(C1_post$phi_valuesI, 2, mean)
apply(C2_post$phi_valuesI, 2, mean)
apply(C3_post$phi_valuesI, 2, mean)
apply(C4_post$phi_valuesI, 2, mean)


# 95% credible intervals for intervention phi

apply(
  A1_post$phi_valuesI,
  2,
  quantile,
  probs = c(0.025, 0.975)
)

apply(
  A2_post$phi_valuesI,
  2,
  quantile,
  probs = c(0.025, 0.975)
)

apply(
  A3_post$phi_valuesI,
  2,
  quantile,
  probs = c(0.025, 0.975)
)

apply(
  A4_post$phi_valuesI,
  2,
  quantile,
  probs = c(0.025, 0.975)
)

apply(
  C1_post$phi_valuesI,
  2,
  quantile,
  probs = c(0.025, 0.975)
)

apply(
  C2_post$phi_valuesI,
  2,
  quantile,
  probs = c(0.025, 0.975)
)

apply(
  C3_post$phi_valuesI,
  2,
  quantile,
  probs = c(0.025, 0.975)
)

apply(
  C4_post$phi_valuesI,
  2,
  quantile,
  probs = c(0.025, 0.975)
)


# RANDOM-EFFECT STANDARD DEVIATIONS


# Arabiensis Model 2

mean(A2_post$sigma_alpha)
quantile(A2_post$sigma_alpha, c(0.025, 0.975))

mean(A2_post$sigma_alpha2)
quantile(A2_post$sigma_alpha2, c(0.025, 0.975))

# Arabiensis Model 3

mean(A3_post$sigma_alpha)
quantile(A3_post$sigma_alpha, c(0.025, 0.975))

mean(A3_post$sigma_alpha2)
quantile(A3_post$sigma_alpha2, c(0.025, 0.975))

apply(A3_post$sigma_beta, 2, mean)

apply(
  A3_post$sigma_beta,
  2,
  quantile,
  probs = c(0.025, 0.975)
)

apply(A3_post$sigma_kappa, 2, mean)

apply(
  A3_post$sigma_kappa,
  2,
  quantile,
  probs = c(0.025, 0.975)
)


# Arabiensis Model 4

mean(A4_post$sigma_alpha)
quantile(A4_post$sigma_alpha, c(0.025, 0.975))

mean(A4_post$sigma_alpha2)
quantile(A4_post$sigma_alpha2, c(0.025, 0.975))

mean(A4_post$sigma_vol)
quantile(A4_post$sigma_vol, c(0.025, 0.975))


# Culex Model 2

mean(C2_post$sigma_alpha)
quantile(C2_post$sigma_alpha, c(0.025, 0.975))

mean(C2_post$sigma_alpha2)
quantile(C2_post$sigma_alpha2, c(0.025, 0.975))


# Culex Model 3

mean(C3_post$sigma_alpha)
quantile(C3_post$sigma_alpha, c(0.025, 0.975))

mean(C3_post$sigma_alpha2)
quantile(C3_post$sigma_alpha2, c(0.025, 0.975))

apply(C3_post$sigma_beta, 2, mean)

apply(
  C3_post$sigma_beta,
  2,
  quantile,
  probs = c(0.025, 0.975)
)

apply(C3_post$sigma_kappa, 2, mean)

apply(
  C3_post$sigma_kappa,
  2,
  quantile,
  probs = c(0.025, 0.975)
)

# Culex Model 4

mean(C4_post$sigma_alpha)
quantile(C4_post$sigma_alpha, c(0.025, 0.975))

mean(C4_post$sigma_alpha2)
quantile(C4_post$sigma_alpha2, c(0.025, 0.975))

mean(C4_post$sigma_vol)
quantile(C4_post$sigma_vol, c(0.025, 0.975))