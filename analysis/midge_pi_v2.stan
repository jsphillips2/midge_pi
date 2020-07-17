//=======================================================================================

// Effect of midge / sediment treatment on PI curves

//=======================================================================================




//=======================================================================================

data {
  // Declare variables
  int n_obs; // number of observations
  int xs[3]; // matrix of predictors
  int g_b[n_obs]; // model matrix for beta
  int g_a[n_obs]; // model matrix for alpha
  int g_r[n_obs]; // model matrix for rho
  real par[n_obs]; // PAR observations
  real rateo2[n_obs]; // DO flux observations
  int n_cores; // number of cores (experimental units)
  int core[n_obs]; // core (experimental unit)
  int n_sum; // number of levels over which to generate fitted curves
  int g_sum_b[n_sum]; // model matrix for beta (fitted curves)
  int g_sum_a[n_sum]; // model matrix for alpha (fitted curves)
  int g_sum_r[n_sum]; // model matrix for rho (fitted curves)
  real par_sum[n_sum]; // PAR for fitted curves
  
}


//=======================================================================================


transformed data {
  // Declare variables
  real mu; // mean observed DO flux
  real<lower=0> tau; // sd observed DO flux
  real<lower=0> eta; // maximum observed PAR
  real y[n_obs]; // z-scored DO flux as response variable
  real<lower=0> l[n_obs]; // scale observed PAR
  real<lower=0> l_sum[n_sum]; // scaled PAR for fitted curves
  
  // Calculate values for scaling
  mu = mean(rateo2);
  tau = sd(rateo2);
  eta = mean(par);
  
  // Scale variables
  for (n in 1:n_obs){
    y[n] = (rateo2[n] - mu)/tau;
    l[n] = par[n]/eta;
  }
  
  // Scale variables for summary
  for (n in 1:n_sum) {
   l_sum[n] = par_sum[n]/eta; 
  }
  
  
}


//=======================================================================================


parameters {
  
  // Declare variables
  real<lower=0> b[xs[1]]; // coefficients for beta
  real<lower=0> a[xs[2]]; // coefficients for alpha
  real<lower=0> r[xs[3]]; // coefficients for rho
  real ran[n_cores]; // core-level variation in intercept
  real<lower=0> sig_ran; // core-level standard deviation 
  real<lower=0> sig_res; // residual standard deviation
  
}


//=======================================================================================


transformed parameters {
  
  // Declare variables 
  real gpp_z[n_obs]; // GPP (on z-scored response scale)
  real nep_z[n_obs]; // NEP (on z-scored response scale)
  real y_pred[n_obs]; // predicted values (on z-scored response scale)
  
  // Calculate metablism
  for(n in 1:n_obs){
    gpp_z[n] = b[g_b[n]] * tanh((a[g_a[n]] / b[g_b[n]]) * l[n]);
    nep_z[n] = gpp_z[n] - r[g_r[n]];
    y_pred[n] = nep_z[n] - mu / tau + sig_ran * ran[core[n]];
  }
  
}


//=======================================================================================


model {
  
  // Priors
  b ~ exponential(0.5); 
  a ~ exponential(0.5);
  r ~ exponential(0.5);
  ran ~ normal(0, 1);
  sig_ran ~ gamma(1.5, 1.5 / 0.5);
  sig_res ~ gamma(1.5, 1.5 / 0.5);
  
  // Likelihood
  y ~ normal(y_pred, sig_res);
  
}


//=======================================================================================


generated quantities {
  
  // Declare variables 
  real gpp[n_obs]; // GPP
  real er[n_obs]; // ER
  real nep[n_obs]; // NEP
  real gpp_sum[n_sum]; // GPP (fitted curves)
  real er_sum[n_sum]; // ER (fitted curves)
  real nep_sum[n_sum]; // NEP (fitted curves)
  real sat[n_sum]; // half-saturation constant
  real log_lik [n_obs]; // pointwise log-likelihood
  real log_lik_sum; // total log-likelihood
  
  // Backscale estimates
  for(n in 1:n_obs){
    gpp[n] = tau * gpp_z[n];
    nep[n] = tau * nep_z[n];
    er[n] = gpp[n] - nep[n];
  }
  
  // Calculate summarized
   for (n in 1:n_sum) {
    gpp_sum[n] = tau * b[g_sum_b[n]] * tanh((a[g_sum_a[n]] / b[g_sum_b[n]]) * l_sum[n]);
    er_sum[n] = tau * r[g_sum_r[n]];
    nep_sum[n] = gpp_sum[n] - er_sum[n];
    sat[n] = eta * (b[g_sum_b[n]] / a[g_sum_a[n]]) * atanh(0.5); 
   }
  
  // Pointwise log-likelihood
    for (n in 1 : n_obs) {
      log_lik[n] = normal_lpdf(y[n] | y_pred[n], sig_res);
    }
    
  // Total log-likelihood
  log_lik_sum = sum(log_lik);

}


//=======================================================================================
