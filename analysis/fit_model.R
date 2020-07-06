#==========
#========== Preliminaries
#==========

# load packages
library(tidyverse)
library(rstan)
source("analysis/model_fn.R")

# stan settings
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores()-2)

# load data
data = read_csv("data/metabolism_2015.csv")

# process data
dd = data %>%
  filter(trial=="a") %>%
  mutate(
    rateo2 = 15*(do_f - do_i)/(100*duration), 
    day = as.numeric(date) - min(as.numeric(date)) + 3,
    time = as.numeric(as.factor(day))-1,
    temp = (temp_i + temp_f)/2
  ) %>%
  rename(par = light) %>%
  select(core, rack, date, day, time, par, midge, rateo2, temp) %>%
  na.omit() %>%
  filter(rateo2 > -3)

# data frame for summary
dd_sum <- dd %>%
  mutate(id = interaction(time, midge)) %>%
  split(.$id) %>%
  lapply(function(x) {
    x %>% 
      tidyr::expand(par = seq(min(par), max(par), length.out = 100)) %>%
      mutate(time = unique(x$time),
             midge = unique(x$midge))
  }) %>%
  bind_rows()

# set theme
theme_set(theme_bw() %+replace%
            theme(panel.grid = element_blank(),
                  strip.background = element_blank(),
                  legend.margin = margin(0,0,0,0),
                  strip.text = element_text(size=12),
                  legend.text = element_text(size=12),
                  legend.title = element_text(size=14),
                  axis.text = element_text(size=12, color="black"),
                  axis.title.y = element_text(size=14,angle = 90,
                                              margin=margin(0,15,0,0)),
                  axis.title.x = element_text(size=14,margin=margin(15,0,0,0)),
                  strip.text.x = element_text(margin=margin(0,0,10,0)),
                  strip.text.y = element_text(margin=margin(0,0,0,10), angle=270),
                  axis.title = element_text(size=14)))



#==========
#========== Fit model
#==========

# model form
forms <- c("both","midge_b","midge_a","none")
form <- forms[1]
if(form == "both") {b = formula(~ 1 + time * midge)
                    a = formula(~ 1 + time * midge)
                    r = formula(~ 1 + time * midge)}
if(form == "midge_b") {b = formula(~ 1 + time * midge)
                       a = formula(~ 1 + time)
                       r = formula(~ 1 + time * midge)}
if(form == "midge_a") {b = formula(~ 1 + time)
                       a = formula(~ 1 + time * midge)
                       r = formula(~ 1 + time * midge)}
if(form == "none") {b = formula(~ 1 + time)
                       a = formula(~ 1 + time)
                       r = formula(~ 1 + time * midge)}


# package data
data_list <- model_fn(b_ = b,
                      a_ = a,
                      r_ = r,
                      dd_ = dd,
                      dd_sum_ = dd_sum)

# MCMC specifications
chains <- 4
iter <- 2000
adapt_delta <- 0.9
max_treedepth <- 10

# fit model
fit <- stan(file = "analysis/midge_pi.stan", data = data_list, seed=2e3,
            chains = chains, iter = iter,
            control = list(adapt_delta = adapt_delta, max_treedepth = max_treedepth))


# summary of fit
fit_summary <- summary(fit, probs=c(0.16, 0.5, 0.84))$summary %>% 
  {as_tibble(.) %>%
      mutate(var = rownames(summary(fit)$summary))}

# export
write_rds(list(b = b,
			   a = a,
			   r = r, 
			   data_list = data_list,
               fit = fit, 
               fit_summary = fit_summary),
           paste0("analysis/model_fit/",form,".rds"))

