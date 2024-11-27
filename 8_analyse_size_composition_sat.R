# calc and plot size structure in 2018 plankton 
# moved to THESIS directory (in figures_conclusion)


library(tidyverse)
library(latex2exp)

source("weighted_moving_average.R")
wmv <- weighted_moving_average   # reassign function bc laziness



## read output from extract_warm_and_cold_clusters.R
ts_dat <- read_csv2("output/Output_3k_cluster_analysis_thesis.csv", na = c("", "Inf", "-Inf","NA"))


PFT_dat <- ts_dat %>%
  select(date, `mean PICO warm [mg / m^3]`:`mean MICRO front [mg / m^3]`) %>%
  pivot_longer(-date, names_to = c("type", "regime"), names_prefix = "mean ", names_sep = " ", values_to = "mean_vals")

# stacked area plot
ggplot(data = PFT_dat, aes(x = date, y = mean_vals, group = type)) +
  geom_area(aes(fill = type)) +
  facet_grid(.~regime) +
  scale_x_date(name = "",
               date_breaks = "1 month",
               date_labels = "%b",
               limits = as.Date(c("2018-04-01", "2018-09-15"))) +
  theme_bw()
# # This was needed for troubleshooting.... maybe useful later
# # add nano and pico
# size_dat <-  ts_dat %>%
#   select(date, `mean PICO warm [mg / m^3]`:`mean MICRO front [mg / m^3]`) 
# names(size_dat) = gsub(pattern = " \\[mg / m\\^3\\]", replacement = "", x = names(size_dat)) # get rif of unit, this messes up pivoting 
# 
# size_dat <- size_dat %>%
#   mutate(`mean PICO+NANO cold` = `mean PICO cold`+`mean NANO cold`,
#          `mean PICO+NANO warm` = `mean PICO warm`+`mean NANO warm`,
#          `mean PICO+NANO front` = `mean PICO front`+`mean NANO front`
#   ) %>%
#   select(-c(`mean PICO cold`,`mean NANO cold`,`mean PICO warm`,`mean NANO warm`,`mean PICO front`,`mean NANO front`))  %>%
#   pivot_longer(-date, names_to = c("type", "regime"), names_prefix = "mean ", names_sep = " ", values_to = "mean_vals")

# add nano and pico
size_dat <-  ts_dat %>%
  select(date, `mean PICO warm [mg / m^3]`:`mean MICRO front [mg / m^3]`)  %>%
  mutate(`mean PICO+NANO cold [mg / m^3]` = `mean PICO cold [mg / m^3]`+`mean NANO cold [mg / m^3]`,
         `mean PICO+NANO warm [mg / m^3]` = `mean PICO warm [mg / m^3]`+`mean NANO warm [mg / m^3]`,
         `mean PICO+NANO front [mg / m^3]` = `mean PICO front [mg / m^3]`+`mean NANO front [mg / m^3]`
         ) %>%
  select(-c(`mean PICO cold [mg / m^3]`,`mean NANO cold [mg / m^3]`,`mean PICO warm [mg / m^3]`,`mean NANO warm [mg / m^3]`,`mean PICO front [mg / m^3]`,`mean NANO front [mg / m^3]`))  %>%
  pivot_longer(-date, names_to = c("type", "regime"), names_prefix = "mean ", names_sep = " ", values_to = "mean_vals")

# plot stacked bar plot (micro, nano+pico)
ggplot(data = size_dat, aes(x = date, y = mean_vals, group = type)) + #filter(size_dat, type ==c("MICRO", "PICO+NANO"))
  geom_area(aes(fill = type)) +
  facet_grid(.~regime) +
  scale_x_date(name = "",
               date_breaks = "1 month",
               date_labels = "%b",
               limits = as.Date(c("2018-04-01", "2018-09-15"))) +
  coord_cartesian(expand = F, ylim = c(0,7)) +
  labs(title = "Phytoplankton size structure over time", fill = "Size class", 
       y = TeX("Chl \\textit{a} [mg m$^{-3}$]")) + 
  theme_bw()
ggsave(paste0("output/plots/PFT_timeseries_3k_absolute.pdf"))

# make plot relative

plot_dat <- size_dat %>%
  group_by(date, regime) %>%
  mutate(freq = mean_vals/sum(mean_vals))

# plot stacked bar plot (micro, nano+pico)
ggplot(data = plot_dat, aes(x = date, y = freq*100, group = type)) + #filter(size_dat, type ==c("MICRO", "PICO+NANO"))
  geom_area(aes(fill = type)) +
  facet_grid(.~regime) +
  scale_x_date(name = "",
               date_breaks = "1 month",
               date_labels = "%b",
               limits = as.Date(c("2018-04-01", "2018-09-15"))) +
  coord_cartesian(expand = F) +
  labs(title = "Relative shares of phytoplankton size classes", fill = "Size class",
       y = "Proportion [%]") +
  theme_bw()
ggsave(paste0("output/plots/PFT_timeseries_3k_relative.pdf"))

# micro / pico+nano

plot_dat <- ts_dat %>%
  select(date, `mean PICO warm [mg / m^3]`:`mean MICRO front [mg / m^3]`)  %>%
  mutate(`MICRO_PICO+NANO cold [mg / m^3]` = `mean MICRO cold [mg / m^3]` / (`mean PICO cold [mg / m^3]`+`mean NANO cold [mg / m^3]`),
         `MICRO_PICO+NANO warm [mg / m^3]` = `mean MICRO warm [mg / m^3]` / (`mean PICO warm [mg / m^3]`+`mean NANO warm [mg / m^3]`),
         `MICRO_PICO+NANO front [mg / m^3]` = `mean MICRO front [mg / m^3]` / (`mean PICO front [mg / m^3]`+`mean NANO front [mg / m^3]`)
  ) %>%
  select(c(date, `MICRO_PICO+NANO cold [mg / m^3]`, `MICRO_PICO+NANO front [mg / m^3]`, `MICRO_PICO+NANO warm [mg / m^3]`))  %>%
  pivot_longer(-date, names_to = c("type", "regime"), names_prefix = "mean ", names_sep = " ", values_to = "mean_vals") %>%
  group_by(regime) %>%
  mutate(run_mean = wmv(mean_vals, weights=rep(1, 1095), width = 7)) # apply moving average to smooth curves a bit





ggplot(data = plot_dat, aes(x = date, y = run_mean, color = regime)) + #filter(size_dat, type ==c("MICRO", "PICO+NANO"))
  geom_line() +
  geom_line(aes(y = mean_vals), linewidth = 0.1) +
  facet_grid(.~regime) +
  scale_x_date(name = "",
               date_breaks = "1 month",
               date_labels = "%b",
               limits = as.Date(c("2018-04-01", "2018-09-15"))) +
  scale_color_brewer(aesthetics = c("color", "fill"), type = "qual", palette = "Dark2", 
                     limits = c("front", "warm", "cold"))+
  coord_cartesian(expand = F, ylim = c(0,15)) +
  labs(title = "proportions of phytoplankton size classes", color = "cluster",
       y = "micro / (pico+nano)") +
  theme_bw()
ggsave(paste0("output/plots/PFT_timeseries_3k_relative2.pdf"))




## fehlt noch: vergleich zu paper1 (obs, lokal und nur summer und autumn) und paper2 (long ts fur beide trajectory origins)

