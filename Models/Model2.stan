// Single hut-type version of FitEH_minimal_large9.stan
// Model 2: Day effects model: adds day-level random effects to the biting and mortality rates
data {
  int<lower=0> n0; // number of control data points
  int<lower=0> n1; // number of net data points
  int<lower=0> Nnets; // number of net types
  int<lower=0> Ndays; // number of experiment days
  int<lower=1,upper=Nnets> Net[n1]; // net used for data point
  int<lower=0> Day0[n0]; // day of data point for control
  int<lower=0> Day1[n1]; // day of data point for nets
  int<lower=0> fed0[n0];    // FA + FD control
  int<lower=0> total0[n0];  // FA + FD + UD control
  int<lower=0> fed1[n1];    // FA + FD intervention
  int<lower=0> total1[n1];  // FA + FD + UD intervention
  real<lower=0> priorscale;
  real<lower=0> hierarchy;
  real<lower=0> samples0[n0];
  real<lower=0> samples1[n1];
  real<lower=0> weight_count[Ndays];
}


parameters {
  real alpha_value;
  real alpha2_value;

  real<lower=0> sigma_alpha;  #day-level for biting
  real<lower=0> sigma_alpha2; #day-level for mortality

  vector[Ndays] z_alpha;   #day-level for biting
  vector[Ndays] z_alpha2; #day-level for mortality

  vector[Nnets] beta;
  vector[Nnets] kappa;
 
  real<lower=0> phi_valuesC;
  real<lower=0> phi_valuesI[Nnets];
}

transformed parameters {
  // control
  vector[n0] log_totalC; 
  vector[n0] log_qC;     

  // intervention
  vector[n1] log_totalI; 
  vector[n1] log_pI;     
  vector<lower=0>[n1] phiI;



  for(i in 1:n0) {
    real log_alphaB_i =
    alpha_value +
    sigma_alpha * z_alpha[Day0[i]] +
    log(samples0[i]);

    real log_alphaM_i =
    alpha2_value +
    sigma_alpha2 * z_alpha2[Day0[i]] +
    log(samples0[i]);


    log_totalC[i] = log_sum_exp(log_alphaB_i, log_alphaM_i);


    log_qC[i] = log_alphaB_i - log_sum_exp(log_alphaB_i, log_alphaM_i);
  }

  for(j in 1:n1) {
    real log_alphaB_j =
    alpha_value +
    sigma_alpha * z_alpha[Day1[j]] +
    log(samples1[j]);

    real log_alphaM_j =
        alpha2_value +
        sigma_alpha2 * z_alpha2[Day1[j]] +
        log(samples1[j]);
    real log_beta_j   = beta[Net[j]];
    real log_kappa_j  = kappa[Net[j]];


    real log_feed_j  = log_alphaB_j + log_beta_j;


    real log_kill_j  = log_alphaB_j + log_kappa_j;

    
    log_totalI[j] = log_sum_exp(log_sum_exp(log_feed_j, log_kill_j), log_alphaM_j);

    
    log_pI[j] = log_feed_j - log_totalI[j] + log(samples1[j]);

    phiI[j] = phi_valuesI[Net[j]];
  }
}


model {
  
  alpha_value ~ normal(0,1);
  alpha2_value ~ normal(0,1);
  sigma_alpha  ~ cauchy(0, hierarchy); // priors - day-varying rates
  sigma_alpha2 ~ cauchy(0, hierarchy); // priors - day-varying rates

  z_alpha  ~ normal(0,1); // priors - day-varying rates
  z_alpha2 ~ normal(0,1); // priors - day-varying rates


  beta ~ normal(-1, 1);
  kappa ~ normal(-1, 1);

  target += 2*cauchy_lpdf(phi_valuesC | 0, priorscale);
  target += 2*cauchy_lpdf(phi_valuesI | 0, priorscale);


  for(k0 in 1:n0) {
    real wC = 1.0 * n1/n0 + log1p(weight_count[Day0[k0]]);


    target += wC * neg_binomial_2_log_lpmf(total0[k0] | log_totalC[k0], phi_valuesC);

    
    if (total0[k0] > 0)
      target += wC * binomial_lpmf(fed0[k0] | total0[k0], exp(log_qC[k0]));
  }


  for(k1 in 1:n1) {
    real wI = 1 + log1p(weight_count[Day1[k1]]);

    target += wI * neg_binomial_2_log_lpmf(total1[k1] | log_totalI[k1], phiI[k1]);

    if (total1[k1] > 0)
      target += wI * binomial_lpmf(fed1[k1] | total1[k1], exp(log_pI[k1]));
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
  for(d in 1:Ndays){
  alpha_day_out[d] =
    exp(alpha_value + sigma_alpha * z_alpha[d]);

  alpha2_day_out[d] =
    exp(alpha2_value + sigma_alpha2 * z_alpha2[d]);
}
  alpha_out  = exp(alpha_value);
  alpha2_out = exp(alpha2_value);

  for(n in 1:Nnets) {
    real log_b = beta[n];
    real log_k = kappa[n];
    beta_out[n]   = exp(log_b);
    kappa_out[n]  = exp(log_k);
    p_feed_out[n] = exp(log_b - log_sum_exp(log_b, log_k));
  }

  // control observations
  for(i in 1:n0) {

    log_lik[i] =
      neg_binomial_2_log_lpmf(total0[i] | log_totalC[i], phi_valuesC);

    if(total0[i] > 0)
      log_lik[i] +=
        binomial_lpmf(fed0[i] | total0[i], exp(log_qC[i]));
  }

  // intervention observations
  for(i in 1:n1) {

    log_lik[n0+i] =
      neg_binomial_2_log_lpmf(total1[i] | log_totalI[i], phiI[i]);

    if(total1[i] > 0)
      log_lik[n0+i] +=
        binomial_lpmf(fed1[i] | total1[i], exp(log_pI[i]));
  }
}

