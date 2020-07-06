//=======================================================================================

// Effect of midges on PI curve

//=======================================================================================




//=======================================================================================

data {
  // Declare variables
  int n_obs;
  int xs[3];
  matrix[n_obs, xs[1]] x_b;
  matrix[n_obs, xs[2]] x_a;
  matrix[n_obs, xs[3]] x_r;
  real par[n_obs];
  real rateo2[n_obs];
  int n_sum;
  matrix[n_sum, xs[1]] x_sum_b;
  matrix[n_sum, xs[2]] x_sum_a;
  matrix[n_sum, xs[3]] x_sum_r;
  real par_sum[n_sum];
  
}


//=======================================================================================


transformed data {
  // Declare variables
  real mu;
  real<lower=0> tau;
  real<lower=0> eta;
  real y[n_obs];
  real<lower=0> l[n_obs];
  real<lower=0> l_sum[n_sum];
  
  // SCalculate variable scales
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
  vector[xs[1]] coef_b;
  vector[xs[2]] coef_a;
  vector[xs[3]] coef_r;
  real<lower=0> sig_res;
  
}


//=======================================================================================


transformed parameters {
  
  // Declare variables 
  real b[n_obs];
  real a[n_obs];
  real r[n_obs];
  real gpp_z[n_obs];
  real nep_z[n_obs];
  real y_pred[n_obs];
  
  // Calculate metablism
  for(n in 1:n_obs){
    b[n] = exp(dot_product(coef_b, x_b[n,]));
    a[n] = exp(dot_product(coef_a, x_a[n,]));
    r[n] = exp(dot_product(coef_r, x_r[n,]));
    gpp_z[n] = b[n] * tanh((a[n] / b[n]) * l[n]);
    nep_z[n] = b[n] * tanh((a[n] / b[n]) * l[n]) - r[n];
    y_pred[n] = nep_z[n] - mu/tau;
  }
  
}


//=======================================================================================


model {
  
  // Priors
  coef_b ~ normal(0, 1); 
  coef_a ~ normal(0, 1); 
  coef_r ~ normal(0, 1); 
  sig_res ~ exponential(2);
  
  // Likelihood
  y ~ normal(y_pred, sig_res);
}


//=======================================================================================


generated quantities {
  
  // Declare variables 
  real gpp[n_obs];
  real er[n_obs];
  real nep[n_obs];
  real gpp_sum[n_sum];
  real er_sum[n_sum];
  real nep_sum[n_sum];
  real sat[n_sum];
  real log_lik [n_obs];
  real log_lik_sum;
  
  // Backscale estimates
  for(n in 1:n_obs){
    gpp[n] = tau * gpp_z[n];
    er[n] = tau * r[n];
    nep[n] = tau * nep_z[n];
  }
  
  // Calculate summarized
  {
    real b_t;
    real a_t;
    real r_t;
    for (n in 1:n_sum) {
     b_t = exp(dot_product(coef_b, x_sum_b[n,]));
     a_t = exp(dot_product(coef_a, x_sum_a[n,]));
     r_t = exp(dot_product(coef_r, x_sum_r[n,]));
     gpp_sum[n] = tau * b_t * tanh((a_t / b_t) * l_sum[n]);
     er_sum[n] = tau * r_t;
     nep_sum[n] = gpp_sum[n] - er_sum[n];
     sat[n] = eta * (b_t / a_t) * atanh(0.5); 
    }
  }
  
  // Pointwise log-likelihood
    for (n in 1 : n_obs) {
      log_lik[n] = normal_lpdf(y[n] | y_pred[n], sig_res);
    }
    
  // Total log-likelihood
  log_lik_sum = sum(log_lik);

}


//=======================================================================================
