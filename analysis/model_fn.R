model_fn <- function(b_, r_, a_, dd_, dd_sum_) {
  
  # define model matrices
  x_b_ <- model.matrix(b_, data = dd_)
  x_a_ <- model.matrix(a_, data = dd_)
  x_r_ <- model.matrix(r_, data = dd_)
  
  # define column indeces
  xs_ <- c(ncol(x_b_), 
           ncol(x_a_),
           ncol(x_r_))
  
  # define number of observations
  n_obs_ <- nrow(dd_) 
  
  
  # define summarized model matrixs
  x_sum_b_ <- model.matrix(b_, data = dd_sum_)
  x_sum_a_ <- model.matrix(a_, data = dd_sum_)
  x_sum_r_ <- model.matrix(r_, data = dd_sum_)
  
  # define number of summarized values
  n_sum_ <- nrow(dd_sum_)
  
  # package dataa
  data_list_ <- list(n_obs = n_obs_,
                     xs = xs_,
                     x_b = x_b_,
                     x_a = x_a_,
                     x_r = x_r_,
                     par = dd$par,
                     rateo2 = dd$rateo2,
                     n_cores = max(dd$core),
                     core = dd$core,
                     n_sum = n_sum_,
                     x_sum_b = x_sum_b_,
                     x_sum_a = x_sum_a_,
                     x_sum_r = x_sum_r_,
                     par_sum = dd_sum_$par)

  # return
  return(data_list_)
  
}

