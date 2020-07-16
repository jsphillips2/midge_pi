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
sonde = read_csv("data/sonde_final.csv")

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



# import model fits
models <- lapply(list.files("analysis/model_fit",
                            full.names = T),
                 read_rds)
model_names <- list.files("analysis/model_fit")
names(models) <- model_names

# labels
midge_lab <- tibble(midge = c(0, 1),
                    Sediment = c("Sieved", "Intact") %>% 
                      factor(levels = c("Sieved", "Intact")))
day_lab <- tibble(time = c(0, 1),
                  Day= c("Day~3", "Day~15") %>% 
                    factor(levels = c("Day~3", "Day~15")))
model_lab <- tibble(model = c("both.rds",
                              "midge_b.rds",
                              "midge_a.rds",
                              "none.rds"),
                    Model = c("beta~and~alpha",
                              "beta",
                              "alpha",
                              "neither")%>% 
                      factor(levels = c("beta",
                                        "alpha",
                                        "beta~and~alpha",
                                        "neither")))

# set theme
theme_set(theme_bw() %+replace%
            theme(panel.grid = element_blank(),
                  strip.background = element_blank(),
                  plot.margin = margin(1,1,1,1),
                  legend.margin = margin(0,0,0,-4),
                  strip.text = element_text(size = 10),
                  strip.text.y = element_text(angle = -90, margin=margin(0,0,0,2)),
                  strip.text.x = element_text(margin=margin(0,0,2,0)),
                  legend.text = element_text(size = 8),
                  axis.text = element_text(size = 10, color="black"),
                  axis.title = element_text(size =10),
                  axis.title.y = element_text(angle = 90, margin=margin(0,10,0,0)),
                  axis.title.x = element_text(margin = margin(10,0,0,0)),
                  panel.spacing = unit(0.1, "lines"),
                  axis.ticks = element_line(size = 0.25)))





#==========
#========== Model fit
#==========

# extract model fit
nep_sum <- lapply(model_names, function(x) {
  {models[[x]][["fit_summary"]] %>% 
      filter(str_detect(.$var, "nep_sum")) %>%
      mutate(name = strsplit(var, "\\[|\\]|,") %>% map_chr(~.x[1]),
             id = strsplit(var, "\\[|\\]|,") %>% map_int(~as.integer(.x[2]))) %>%
      filter(name %in% c("nep_sum")) %>%
      select(id, mean, sd) %>%
      unique() %>%
      left_join(dd_sum %>%
                  mutate(id = row_number())) %>%
      select(-id) %>%
      unique() %>%
      mutate(model = x)}
}
) %>%
  bind_rows() %>%
  full_join(midge_lab) %>%
  full_join(day_lab) %>% 
  full_join(model_lab) 

# plot
p1 <- dd %>%
  full_join(midge_lab) %>%
  full_join(day_lab) %>%
  ggplot()+
  facet_grid(Model~Day, labeller=label_parsed)+
  geom_hline(yintercept = 0, size = 0.1)+
  geom_ribbon(data = nep_sum %>% filter(model == "both.rds"), 
              aes(par,
                  ymin = mean - sd,
                  ymax = mean + sd,
                  group = interaction(Model, Sediment)),
              fill = "goldenrod2",
              alpha = 0.2)+
  geom_ribbon(data = nep_sum %>% filter(model == "midge_a.rds"), 
              aes(par,
                  ymin = mean - sd,
                  ymax = mean + sd,
                  group = interaction(Model, Sediment)),
              fill = "dodgerblue",
              alpha = 0.3)+
  geom_point(aes(par, rateo2, fill = Sediment), 
             shape = 21,
             size = 1.3)+
  geom_ribbon(data = nep_sum %>% filter(model == "midge_b.rds"), 
              aes(par,
                  ymin = mean - sd,
                  ymax = mean + sd,
                  group = interaction(Model, Sediment)),
              fill = "firebrick",
              alpha = 0.3)+
  geom_line(data = nep_sum %>% filter(model != "none.rds"), 
            aes(par,
                mean, 
                linetype = Sediment,
                color = model),
            size = 0.4)+
  scale_fill_manual("",values = c("gray90","black"))+
  scale_linetype_manual("",values = c(2, 1), guide = F)+
  scale_color_manual(values = c("goldenrod2","dodgerblue","firebrick"),
                     guide = F)+
  scale_y_continuous(limits=c(-0.135,0.135),
                     breaks=c(-0.1,0,0.1),
    expression("NEP (g "*O[2]~m^{-2}~h^{-1}*")"))+
  scale_x_continuous(expression("PAR ("*mu*mol~photons~m^{-2}~s^{-1}*")"), 
    breaks=c(40,120, 200), limits=c(0, 240))+
  theme(legend.position = c(0.875, 0.75),
        legend.text = element_text(margin = margin(l = -6)),
        legend.key.size = unit(0.7, "lines"),
        legend.spacing.y = unit(0, "lines"),
        panel.border = element_rect(size = 0.25))+
  guides(fill = guide_legend(override.aes = list(size = 2)))

p1
# ggsave(file = "analysis/figures/fig_1.pdf",
#           width = 3.5, height = 4)




#==========
#========== Half saturation
#==========

# extract
sats <- lapply(model_names, function(x) {
  {models[[x]][["fit_summary"]] %>% 
      filter(str_detect(.$var, "sat")) %>%
      mutate(id = strsplit(var, "\\[|\\]|,") %>% map_int(~as.integer(.x[2]))) %>%
      select(id, mean, sd) %>%
      left_join(dd_sum %>%
                  select(time, midge) %>%
                  mutate(id = row_number())) %>%
      select(-id) %>%
      unique() %>%
      mutate(model = x)}
}
) %>%
  bind_rows()%>%
  full_join(midge_lab) %>%
  full_join(day_lab) %>% 
  full_join(model_lab) 

# plot
p2 <- sats %>%
  filter(model != "none.rds") %>%
  mutate(midge = midge + 
           ifelse(model == "both.rds",  - 0.15, 0) + 
           ifelse(model == "midge_b.rds",  + 0.15, 0)) %>%
  ggplot(aes(midge, mean, color = Model))+
  facet_wrap(~Day, nrow = 2, labeller=label_parsed)+
  geom_errorbar(aes(ymin = mean - sd, ymax = mean + sd), width = 0, size = 0.5)+
  geom_line(size = 0.5)+
  geom_point(aes(fill = Model), 
             size = 1.75, 
             shape = 21, 
             color = "black",
             stroke = 0.3)+
  scale_x_continuous("Sediment Treatment",
                     limits = c(-0.5, 1.5),
                     breaks = c(0, 1),
                     labels = c("Sieved","Intact"))+
  scale_y_continuous(expression("Half saturation ("*mu*mol~photons~m^{-2}~s^{-1}*")"), 
                     breaks=c(40,100, 160), limits=c(0,200))+
  scale_color_manual("",
                     values = c("firebrick","dodgerblue","goldenrod2"),
                     labels = c(bquote(beta),
                                bquote(alpha),
                                bquote(beta~and~alpha)))+
  scale_fill_manual("",
                     values = c("firebrick","dodgerblue","goldenrod2"),
                     labels = c(bquote(beta),
                                bquote(alpha),
                                bquote(beta~and~alpha)))+
  theme(legend.position = c(0.80, 0.89),
        legend.text = element_text(margin = margin(l = -6)),
        legend.key.size = unit(0.9, "lines"),
        legend.spacing.y = unit(0, "lines"),
        panel.border = element_rect(size = 0.25))+
  guides(color = guide_legend(override.aes = list(size = 2, linetype = 0)))

p2
# ggsave(file = "analysis/figures/fig_2.pdf",
#           width = 2.5, height = 4)




#==========
#========== Benthic light
#==========

s = 30

ben_par <- sonde %>%
  select(year, yday, hour, par0, ext, pcyv) %>%
  filter(year %in% c(2018)) %>%
  mutate(par = par0 * exp(-ext * 3.3)) %>%
  group_by(year, yday) %>%
  summarize(par = mean(par, na.rm = T),
            pcyv = mean(pcyv, na.rm = T)) %>%
  na.omit() %>%
  ungroup() %>%
  mutate(pcyv_z = (pcyv - mean(pcyv))/sd(pcyv),
         pcyv_z = pcyv_z + abs(min(pcyv_z)),
         pcyv_z = s * pcyv_z) %>%
  gather(var, val, pcyv_z, par)

mean_sat <- test <- extract(models[["midge_a.rds"]]$fit, pars = "sat")$sat %>%
  as_tibble() %>%
  mutate(step = row_number()) %>%
  gather(id, val, -step) %>%
  ungroup() %>%
  left_join(dd_sum %>%
              mutate(id = paste0("V",row_number())) %>%
              select(id, time, midge)) %>%
  group_by(time, midge, step) %>%
  summarize(val = unique(val)) %>%
  group_by(step, midge) %>%
  summarize(val = mean(val)) %>%
  group_by(midge) %>%
  summarize(sd = sd(val),
            sat = mean(val)) %>%
  expand(nesting(sat, sd, midge), yday = c(150,230)) %>%
  mutate(label = factor(midge, levels = c(0,1), 
                        labels = c("Sieved","Intact")))

p3 <- ben_par %>%
  ggplot(aes(x = yday))+
  geom_ribbon(data = mean_sat,
              aes(ymin = sat - sd,
                  ymax = sat + sd,
                  group = midge),
              alpha = 0.3, 
              fill = "dodgerblue")+
  geom_line(aes(y = val, color = var), size = 0.4)+
  scale_y_continuous(expression("PAR ("*mu*mol~photons~m^{-2}~s^{-1}*")"), 
                     breaks=c(25,75,125), limits=c(0,150),
                     sec.axis = sec_axis(name = "Phycocyanin index (SDs)",
                                         ~ . / s))+
  scale_x_continuous("Day of year",
                     limits = c(147, 230),
                     breaks = c(160, 190, 220))+
  scale_color_manual("",
                     values = c("black","cyan4"),
                     labels = c("Benthic PAR",
                                "Phycocyanin"))+
  geom_text(data = mean_sat,
            aes(x = 147, 
                y = sat, 
                label = label),
            angle = 90,
            size = 2.8)+
  theme(legend.position = c(0.5, 0.942),
        legend.direction = "horizontal",
        legend.text = element_text(margin = margin(l = -10)),
        legend.spacing.x = unit(1, "lines"),
        legend.spacing.y = unit(0, "lines"),
        axis.title.y.right = element_text(angle = -90, margin=margin(0,0,0,15)),
        panel.border = element_rect(size = 0.25))
p3

# ggsave(file = "analysis/figures/fig_3.pdf",
#           width = 3.5, height = 2.75)
  










test %>% filter(step == 1)
test %>%
  group_by(id) %>%
  summarize(l = length(val))

  select(-id) %>%
  unique() %>%
  group_by(midge) %>%
  summarize(sat = mean(val),
            sd = sd(val))



sat_est <- lapply(model_names, function(x) {
  t(as.matrix(extract(models[[x]][["fit"]], pars = "sat")[[1]])) %>%
    as_tibble() %>%
    mutate(id = row_number()) %>%
    left_join(dd_sum %>%
                mutate(id = row_number()) %>%
                select(id, time, midge)) %>%
    gather(var, val, -id, -time, -midge) %>%
    select(-id) %>%
    unique() %>%
    mutate(model = x)
  # {models[[x]][["fit"]] %>% 
  #     filter(str_detect(.$var, "sat")) %>%
  #     mutate(id = strsplit(var, "\\[|\\]|,") %>% map_int(~as.integer(.x[2]))) %>%
  #     select(id, mean, sd) %>%
  #     left_join(dd_sum %>%
  #                 select(time, midge) %>%
  #                 mutate(id = row_number())) %>%
  #     select(-id) %>%
  #     unique() %>%
  #     mutate(model = x)}
}
) %>%
  bind_rows()
sat_est %>%
  filter(model != "none.rds") %>%
  spread(midge, val) %>%
  mutate(effect = `1` - `0`) %>%
  group_by(model, time) %>%
  summarize(se = sd(effect),
            effect = mean(effect)) %>%
  ungroup() %>%
  ggplot(aes(model, effect, color = model))+
  facet_wrap(~time)+
  geom_point(size = 3)+
  geom_errorbar(aes(ymin = effect - se, 
                    ymax = effect + se), width = 0)+
  geom_hline(yintercept = 0)+
  scale_color_manual(values = c("goldenrod","dodgerblue","firebrick"),
                     guide = F)
  


library(loo)

# function for LOOIC
loo_fn <- function(x){
  
  log_lik <- extract_log_lik(x, merge_chains = FALSE)
  r_eff <- relative_eff(exp(log_lik))
  
  loo <- loo(log_lik, r_eff = r_eff, cores = 10)
  
  return(loo)
}



loo_list <- lapply(model_names, function(x) {
  models[[x]][["fit"]] %>% loo_fn
}
)

comp = {loo_compare(loo_list[[1]],
                    loo_list[[2]],
                    loo_list[[3]],
                    loo_list[[4]]) %>%
    as_tibble() %>%
    as.matrix() %>%
    as.tibble() %>%
    mutate(model = model_names[as.numeric(rownames(.))])} %>%
  select(model, looic) %>%
  mutate(loo_dev = round(looic - looic[1], 1))
  


