model_fn_v2 <- function(b_, r_, a_, dd_, dd_sum_) {
  
  # define model matrices
  x_b_ <- model.matrix(b_, data = dd_)
  x_a_ <- model.matrix(a_, data = dd_)
  x_r_ <- model.matrix(r_, data = dd_)
  
  # define column indeces
  xs_ <- c(ncol(x_b_), 
           ncol(x_a_),
           ncol(x_r_))
  
  # define parameter groupings
  g_b_ <- as.numeric(factor(apply(x_b_, 1, paste, collapse = "")))
  g_a_ <- as.numeric(factor(apply(x_a_, 1, paste, collapse = "")))
  g_r_ <- as.numeric(factor(apply(x_r_, 1, paste, collapse = "")))
  
  # define number of observations
  n_obs_ <- nrow(dd_) 
  
  
  # define summarized model matrixs
  x_sum_b_ <- model.matrix(b_, data = dd_sum_)
  x_sum_a_ <- model.matrix(a_, data = dd_sum_)
  x_sum_r_ <- model.matrix(r_, data = dd_sum_)
  
  # define parameter groupings
  g_sum_b_ <- as.numeric(factor(apply(x_sum_b_, 1, paste, collapse = ""), 
                                labels = c(1:xs_[1])))
  g_sum_a_ <- as.numeric(factor(apply(x_sum_a_, 1, paste, collapse = ""), 
                                labels = c(1:xs_[2])))
  g_sum_r_ <- as.numeric(factor(apply(x_sum_r_, 1, paste, collapse = ""), 
                                labels = c(1:xs_[3])))
  
  # define number of summarized values
  n_sum_ <- nrow(dd_sum_)
  
  # package dataa
  data_list_ <- list(n_obs = n_obs_,
                     xs = xs_,
                     g_b = g_b_,
                     g_a = g_a_,
                     g_r = g_r_,
                     par = dd$par,
                     rateo2 = dd$rateo2,
                     n_cores = max(dd$core),
                     core = dd$core,
                     n_sum = n_sum_,
                     g_sum_b = g_sum_b_,
                     g_sum_a = g_sum_a_,
                     g_sum_r = g_sum_r_,
                     par_sum = dd_sum_$par)

  # return
  return(data_list_)
  
}

