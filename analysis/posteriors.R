options(mc.cores = parallel::detectCores()-2)

coef_b <- lapply(model_names, function(x) {
  extract(models[[x]][["fit"]], pars = "coef_b") %>%
    lapply(as_tibble) %>%
    bind_rows() %>%
    gather() %>%
    mutate(model = x)
}
) %>%
  bind_rows()

coef_b %>%
  ggplot(aes(value))+
  facet_grid(model~key)+
  geom_vline(xintercept = 0,
             linetype = 2)+
  geom_histogram(data = coef_b %>% 
                   expand(model, key, value = rnorm(12000, 0, 1)),
                 bins = 1e3,
                 position = "identity"
                 )+
  geom_histogram(position = "identity",
                 bins = 1e3,
                 alpha = 0.5,
                 color= "dodgerblue")


coef_a <- lapply(model_names, function(x) {
  extract(models[[x]][["fit"]], pars = "coef_a") %>%
    lapply(as_tibble) %>%
    bind_rows() %>%
    gather() %>%
    mutate(model = x)
}
) %>%
  bind_rows()

coef_a %>%
  ggplot(aes(value))+
  facet_grid(model~key)+
  geom_vline(xintercept = 0,
             linetype = 2)+
  geom_histogram(data = coef_a %>% 
                   expand(model, key, value = rnorm(12000, 0, 1)),
                 bins = 1e3,
                 position = "identity"
  )+
  geom_histogram(position = "identity",
                 bins = 1e3,
                 alpha = 0.5,
                 color= "dodgerblue")

coef_r <- lapply(model_names, function(x) {
  extract(models[[x]][["fit"]], pars = "coef_r") %>%
    lapply(as_tibble) %>%
    bind_rows() %>%
    gather() %>%
    mutate(model = x)
}
) %>%
  bind_rows()

coef_r %>%
  ggplot(aes(value))+
  facet_grid(model~key)+
  geom_vline(xintercept = 0,
             linetype = 2)+
  geom_histogram(data = coef_r %>% 
                   expand(model, key, value = rnorm(12000, 0, 1)),
                 bins = 1e3,
                 position = "identity"
  )+
  geom_histogram(position = "identity",
                 bins = 1e3,
                 alpha = 0.5,
                 color= "dodgerblue")



sig_ran <- lapply(model_names, function(x) {
  extract(models[[x]][["fit"]], pars = "sig_ran") %>%
    lapply(as_tibble) %>%
    bind_rows() %>%
    gather() %>%
    mutate(model = x)
}
) %>%
  bind_rows()

sig_ran %>%
  ggplot(aes(value))+
  facet_wrap(~model)+
  geom_vline(xintercept = 0.5,
             linetype = 2)+
  geom_histogram(data = sig_ran %>% 
                   expand(model,value = rgamma(12000, 1.5, 3)),
                 bins = 1e3,
                 position = "identity"
  )+
  geom_histogram(position = "identity",
                 bins = 1e3,
                 alpha = 0.5,
                 color= "dodgerblue")


sig_res <- lapply(model_names, function(x) {
  extract(models[[x]][["fit"]], pars = "sig_res") %>%
    lapply(as_tibble) %>%
    bind_rows() %>%
    gather() %>%
    mutate(model = x)
}
) %>%
  bind_rows()

sig_res %>%
  ggplot(aes(value))+
  facet_wrap(~model)+
  geom_vline(xintercept = 0.5,
             linetype = 2)+
  geom_histogram(data = sig_res %>% 
                   expand(model,value = rgamma(12000, 1.5, 3)),
                 bins = 1e3,
                 position = "identity"
  )+
  geom_histogram(position = "identity",
                 bins = 1e3,
                 alpha = 0.5,
                 color= "dodgerblue")




b <- parallel::mclapply(model_names, function(x) {
  extract(models[[x]][["fit"]], pars = "b") %>%
    lapply(as_tibble) %>%
    bind_rows() %>%
    gather() %>%
    mutate(id = str_split(key, "V") %>% 
             map_int(~as.integer(.x[2])),
           model = x) %>%
    full_join(dd %>%
                mutate(id = row_number()) %>%
                select(id, time, midge)) 
 }
)

b %>%
  bind_rows() %>%
  ggplot(aes(value, fill = factor(midge)))+
  facet_grid(model~time)+
  geom_histogram(position = "identity",
                 bins = 1e3,
                 alpha = 0.5)+
  scale_x_continuous(limits = c(0, 10))




a <- parallel::mclapply(model_names, function(x) {
  extract(models[[x]][["fit"]], pars = "a") %>%
    lapply(as_tibble) %>%
    bind_rows() %>%
    gather() %>%
    mutate(id = str_split(key, "V") %>% 
             map_int(~as.integer(.x[2])),
           model = x) %>%
    full_join(dd %>%
                mutate(id = row_number()) %>%
                select(id, time, midge)) 
}
)

a %>%
  bind_rows() %>%
  ggplot(aes(value, fill = factor(midge)))+
  facet_grid(model~time)+
  geom_histogram(position = "identity",
                 bins = 1e3,
                 alpha = 0.5)+
  scale_x_continuous(limits = c(0, 10))


r <- parallel::mclapply(model_names, function(x) {
  extract(models[[x]][["fit"]], pars = "r") %>%
    lapply(as_tibble) %>%
    bind_rows() %>%
    gather() %>%
    mutate(id = str_split(key, "V") %>% 
             map_int(~as.integer(.x[2])),
           model = x) %>%
    full_join(dd %>%
                mutate(id = row_number()) %>%
                select(id, time, midge)) 
}
)

r %>%
  bind_rows() %>%
  ggplot(aes(value, fill = factor(midge)))+
  facet_grid(model~time)+
  geom_histogram(position = "identity",
                 bins = 1e3,
                 alpha = 0.5)



sat_post <- parallel::mclapply(model_names, function(x) {
  extract(models[[x]][["fit"]], pars = "sat") %>%
    lapply(as_tibble) %>%
    bind_rows() %>%
    gather() %>%
    mutate(id = str_split(key, "V") %>% 
             map_int(~as.integer(.x[2])),
           model = x) %>%
    full_join(dd_sum %>%
                mutate(id = row_number()) %>%
                select(id, time, midge)) %>%
    select(-id, - key) %>%
    unique()
}
)

sat_post %>%
  bind_rows() %>%
  ggplot(aes(value, fill = factor(midge)))+
  facet_grid(model~time)+
  geom_histogram(position = "identity",
                 bins = 1e3,
                 alpha = 0.5)+
  scale_x_continuous(limits = c(0,300))
