models <- lapply(list.files("analysis/model_fit",
                            full.names = T),
                 read_rds)
model_names <- list.files("analysis/model_fit")
names(models) <- model_names

sats <- lapply(model_names, function(x) {
  {models[[x]][["fit_summarry"]] %>% 
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
  bind_rows()


sats %>%
  filter(model != "none.rds") %>%
  mutate(midge = midge + 
                    ifelse(model == "both.rds",  - 0.15, 0) + 
                    ifelse(model == "midge_b.rds",  + 0.15, 0)) %>%
  ggplot(aes(midge, mean, color = factor(model)))+
  facet_wrap(~time, nrow = 2)+
  geom_point(size = 4)+
  geom_line(size = 0.5)+
  geom_errorbar(aes(ymin = mean - sd, ymax = mean + sd), width = 0, size = 0.3)+
  scale_x_continuous("Sediment Treatment",
                     limits = c(-0.5, 1.5),
                     breaks = c(0, 1),
                     labels = c("Sieved","Intact"))+
  scale_y_continuous("Half saturation")+
  scale_color_manual("",
                     values = c("goldenrod","dodgerblue","firebrick"))



nep_sum <- lapply(model_names, function(x) {
  {models[[x]][["fit_summarry"]] %>% 
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
  bind_rows()



dd %>%
  ggplot()+
  facet_grid(model~time)+
  geom_ribbon(data = nep_sum %>% filter(model == "both.rds"), 
              aes(par,
                  ymin = mean - sd,
                  ymax = mean + sd,
                  group = interaction(model, midge)),
              fill = "goldenrod",
              alpha = 0.2)+
  geom_ribbon(data = nep_sum %>% filter(model == "midge_a.rds"), 
              aes(par,
                  ymin = mean - sd,
                  ymax = mean + sd,
                  group = interaction(model, midge)),
              fill = "dodgerblue",
              alpha = 0.3)+
  geom_ribbon(data = nep_sum %>% filter(model == "midge_b.rds"), 
              aes(par,
                  ymin = mean - sd,
                  ymax = mean + sd,
                  group = interaction(model, midge)),
              fill = "firebrick",
              alpha = 0.3)+
  geom_point(aes(par, rateo2, fill = factor(midge)), 
             shape = 21,
             size = 2)+
  geom_line(data = nep_sum %>% filter(model != "none.rds"), 
            aes(par,
                mean, 
                linetype = factor(midge),
                color = model),
            size = 0.5)+
  scale_fill_manual(values = c("gray90","black"))+
  scale_linetype_manual(values = c(2, 1))+
  scale_color_manual(values = c("goldenrod","dodgerblue","firebrick"),
                     guide = F)+
  theme(legend.position = "top")
