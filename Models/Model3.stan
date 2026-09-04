// Single hut-type version of FitEH_minimal_large9.stan
// Model 3:hierarchical Bayesian model which expands Model 2 by varying β and κ by day.
data {
  int<lower=0> n0; // number of control observations
  int<lower=0> n1; // number of intervention observations

  int<lower=1> Nnets; // number of net types
  int<lower=1> Ndays; // number of experiment days

  int<lower=1, upper=Nnets> Net[n1];

  int<lower=1, upper=Ndays> Day0[n0];
  int<lower=1, upper=Ndays> Day1[n1];

  int<lower=0> fed0[n0];   // FA + FD control
  int<lower=0> total0[n0]; // FA + FD + UD control

  int<lower=0> fed1[n1];   // FA + FD intervention
  int<lower=0> total1[n1]; // FA + FD + UD intervention

  real<lower=0> priorscale;
  real<lower=0> hierarchy;

  real<lower=0> samples0[n0];
  real<lower=0> samples1[n1];

  real<lower=0> weight_count[Ndays];
}


parameters {
   // Baseline and day-effect parameters as in Model 2
  real alpha_value;
  real alpha2_value;

  real<lower=0> sigma_alpha;
  real<lower=0> sigma_alpha2;

  vector[Ndays] z_alpha;
  vector[Ndays] z_alpha2;

  vector[Nnets] beta_mean;              // Mean beta for each net
  vector<lower=0>[Nnets] sigma_beta; // Amount beta varies between days, separately for each net
  matrix[Nnets, Ndays] z_beta;    // Standardised day effects for beta

  vector[Nnets] kappa_mean;           // Mean kappa for each net
  vector<lower=0>[Nnets] sigma_kappa; // Amount kappa varies between days, separately for each net
  matrix[Nnets, Ndays] z_kappa;   // Standardised day effects for kappa

  real<lower=0> phi_valuesC;
  real<lower=0> phi_valuesI[Nnets];
}

transformed parameters {


  matrix[Nnets, Ndays] beta;  //b_n,d
  matrix[Nnets, Ndays] kappa; //k_n,d


  // Control
  vector[n0] log_totalC;
  vector[n0] log_qC;

  // Intervention
  vector[n1] log_totalI;
  vector[n1] log_pI;

  vector<lower=0>[n1] phiI;


  // construct hierarchical day-specific beta and kappa

  for (n in 1:Nnets) {

    for (d in 1:Ndays) {

      beta[n, d] =
        beta_mean[n]
        + sigma_beta[n] * z_beta[n, d];

      kappa[n, d] =
        kappa_mean[n]
        + sigma_kappa[n] * z_kappa[n, d];
    }
  }


  // Control

  for (i in 1:n0) {

    real log_alphaB_i =
      alpha_value
      + sigma_alpha * z_alpha[Day0[i]]
      + log(samples0[i]);

    real log_alphaM_i =
      alpha2_value
      + sigma_alpha2 * z_alpha2[Day0[i]]
      + log(samples0[i]);


    log_totalC[i] =
      log_sum_exp(
        log_alphaB_i,
        log_alphaM_i
      );


    log_qC[i] =
      log_alphaB_i
      - log_sum_exp(
          log_alphaB_i,
          log_alphaM_i
        );
  }


  // Intervention

  for (j in 1:n1) {

    real log_alphaB_j =
      alpha_value
      + sigma_alpha * z_alpha[Day1[j]]
      + log(samples1[j]);

    real log_alphaM_j =
      alpha2_value
      + sigma_alpha2 * z_alpha2[Day1[j]]
      + log(samples1[j]);


    
    
    real log_beta_j =
      beta[Net[j], Day1[j]];           // beta and kappa now depend on both net and day

    real log_kappa_j =
      kappa[Net[j], Day1[j]];


    // Feeding component:
    // alpha_B * beta
    real log_feed_j =
      log_alphaB_j
      + log_beta_j;


    // Killing component:
    // alpha_B * kappa
    real log_kill_j =
      log_alphaB_j
      + log_kappa_j;


    // Total intervention rate:
    //
    // alpha_B * beta
    // + alpha_B * kappa
    // + alpha_M
    log_totalI[j] =
      log_sum_exp(
        log_sum_exp(
          log_feed_j,
          log_kill_j
        ),
        log_alphaM_j
      );


    // P(fed | total)
    log_pI[j] =
      log_feed_j
      - log_totalI[j]
      + log(samples1[j]);


    phiI[j] =
      phi_valuesI[Net[j]];
  }
}

model {

  

  alpha_value  ~ normal(0, 1);
  alpha2_value ~ normal(0, 1);

  sigma_alpha  ~ cauchy(0, hierarchy);
  sigma_alpha2 ~ cauchy(0, hierarchy);

  z_alpha  ~ normal(0, 1);
  z_alpha2 ~ normal(0, 1);


 
  
  beta_mean ~ normal(-1, 1);

  sigma_beta ~ cauchy(0, hierarchy);   // priors for between-day variation in beta
  to_vector(z_beta) ~ normal(0, 1);    // Non-centred standardised effects



  kappa_mean ~ normal(-1, 1);

  sigma_kappa ~ cauchy(0, hierarchy);  // priors for between-day variation in kappa
  to_vector(z_kappa) ~ normal(0, 1);



  target +=
    2 * cauchy_lpdf(
      phi_valuesC | 0, priorscale
    );

  target +=
    2 * cauchy_lpdf(
      phi_valuesI | 0, priorscale
    );


  // control observations

  for (k0 in 1:n0) {

    real wC =
      1.0 * n1 / n0
      + log1p(
          weight_count[Day0[k0]]
        );


    target +=
      wC *
      neg_binomial_2_log_lpmf(
        total0[k0]
        | log_totalC[k0],
        phi_valuesC
      );


    if (total0[k0] > 0) {

      target +=
        wC *
        binomial_lpmf(
          fed0[k0]
          | total0[k0],
          exp(log_qC[k0])
        );
    }
  }


  // intervention observations

  for (k1 in 1:n1) {

    real wI =
      1
      + log1p(
          weight_count[Day1[k1]]
        );


    target +=
      wI *
      neg_binomial_2_log_lpmf(
        total1[k1]
        | log_totalI[k1],
        phiI[k1]
      );


    if (total1[k1] > 0) {

      target +=
        wI *
        binomial_lpmf(
          fed1[k1]
          | total1[k1],
          exp(log_pI[k1])
        );
    }
  }
}

generated quantities {

  vector[n0 + n1] log_lik;

  real alpha_out;
  real alpha2_out;



  vector[Nnets] beta_out;
  vector[Nnets] kappa_out;


  vector[Nnets] p_feed_out;

  vector[Ndays] alpha_day_out;
  vector[Ndays] alpha2_day_out;

  // day-specific beta and kappa
  matrix[Nnets, Ndays] beta_day_out;
  matrix[Nnets, Ndays] kappa_day_out;


  for (d in 1:Ndays) {

    alpha_day_out[d] =
      exp(
        alpha_value
        + sigma_alpha * z_alpha[d]
      );

    alpha2_day_out[d] =
      exp(
        alpha2_value
        + sigma_alpha2 * z_alpha2[d]
      );
  }


  alpha_out =
    exp(alpha_value);

  alpha2_out =
    exp(alpha2_value);



  for (n in 1:Nnets) {

    // Typical effect for each net
    beta_out[n] =
      exp(beta_mean[n]);

    kappa_out[n] =
      exp(kappa_mean[n]);


    p_feed_out[n] =
      exp(
        beta_mean[n]
        - log_sum_exp(
            beta_mean[n],
            kappa_mean[n]
          )
      );


    // Day-specific effects
    for (d in 1:Ndays) {

      beta_day_out[n, d] =
        exp(beta[n, d]);

      kappa_day_out[n, d] =
        exp(kappa[n, d]);
    }
  }



  for (i in 1:n0) {

    log_lik[i] =
      neg_binomial_2_log_lpmf(
        total0[i]
        | log_totalC[i],
        phi_valuesC
      );

    if (total0[i] > 0) {

      log_lik[i] +=
        binomial_lpmf(
          fed0[i]
          | total0[i],
          exp(log_qC[i])
        );
    }
  }

  for (i in 1:n1) {

    log_lik[n0 + i] =
      neg_binomial_2_log_lpmf(
        total1[i]
        | log_totalI[i],
        phiI[i]
      );

    if (total1[i] > 0) {

      log_lik[n0 + i] +=
        binomial_lpmf(
          fed1[i]
          | total1[i],
          exp(log_pI[i])
        );
    }
  }
}
