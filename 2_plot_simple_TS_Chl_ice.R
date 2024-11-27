## plot time series
## author: Vanessa Lampe
## date: 2023-10-23
## version: 0
## git: 

library(tidyverse)
library(readxl)
library(writexl)
library(ggpubr)
source("weighted_moving_average.R")
wmv <- weighted_moving_average   # reassign function bc laziness



## read output from extract_warm_and_cold_clusters.R
ts_dat <- read_csv2("output/Output_3k_cluster_analysis_thesis.csv", na = c("", "Inf", "-Inf","NA"))




#_______________________________________________________________________________
# Simple plots

# TS of mean total Chl in FS
ggplot(data = ts_dat, aes(x = date, y = `mean chl total [mg / m^3]`)) + 
  geom_line()

# compare SST in clusters 
ggplot(data = ts_dat, aes(x = date)) + 
  geom_line(aes(y = `mean temp warm [ºC]`, col = "warm")) +
  geom_ribbon(aes(ymin = `mean temp warm [ºC]` - `sd temp warm [ºC]`, 
                  ymax = `mean temp warm [ºC]` + `sd temp warm [ºC]`, 
                  fill = "warm"), alpha = .3) +
  geom_line(aes(y = `mean temp cold [ºC]`, col = "cold")) +
  geom_ribbon(aes(ymin = `mean temp cold [ºC]` - `sd temp cold [ºC]`, 
                  ymax = `mean temp cold [ºC]` + `sd temp cold [ºC]`, 
                  fill = "cold"), alpha = .3) +
  geom_line(aes(y = `mean temp front [ºC]`, col = "front")) +
  geom_ribbon(aes(ymin = `mean temp front [ºC]` - `sd temp front [ºC]`, 
                  ymax = `mean temp front [ºC]` + `sd temp front [ºC]`, 
                  fill = "front"), alpha = .3) +
  labs(y = "STT [ºC]", colour = "cluster") +
  scale_x_date(name = "",
               date_breaks = "1 month",
               date_labels = "%b") +
#  scale_y_continuous(limits = c(0, NA), minor_breaks = F, name = "Chl a conc. [µg / L]") +
  scale_color_manual(aesthetics =c("colour", "fill"), values = c("warm" = "red", "cold" = "blue","front" = "grey"), 
                     breaks = c("cold", "front", "warm")) +
  guides(fill = "none") +
  theme_bw()




# compare raw Chl in clusters 
ggplot(data = ts_dat, aes(x = date)) + 
  geom_line(aes(y = `mean chl warm [mg / m^3]`, col = "warm")) +
  geom_line(aes(y = `mean chl cold [mg / m^3]`, col = "cold")) +
  geom_line(aes(y = `mean chl front [mg / m^3]`, col = "front")) +
  labs(y = expression(paste("Chl ", italic("a"), " mg ", m^{-3})), colour = "cluster") +
 # labs(y = "Chl *a* [mg m^(-3)] ", colour = "cluster") +
  scale_x_date(name = "",
               date_breaks = "1 month",
               date_labels = "%b") +
  #  scale_y_continuous(limits = c(0, NA), minor_breaks = F, name = "Chl a conc. [µg / L]") +
  scale_color_manual(values = c("warm" = "red", "cold" = "blue","front" = "grey"), 
                     breaks = c("cold", "front", "warm")) 
  #theme(axis.title.y = ggtext::element_markdown())



# compare chl in clusters, weighed moving average
wmeans <- ts_dat
wmeans$w_w_chl <- wmv(values = wmeans$`mean chl warm [mg / m^3]`,
                      weights = wmeans$`n.obs. chl warm [pixels]`,
                      width = 15)
wmeans$w_c_chl <- wmv(values = wmeans$`mean chl cold [mg / m^3]`, 
                      weights = wmeans$`n.obs chl cold [pixels]`,
                      width = 15)
wmeans$w_f_chl <- wmv(values = wmeans$`mean chl front [mg / m^3]`, 
                      weights = wmeans$`n.obs. chl front [pixels]`,
                      width = 15)


ggplot(data = wmeans, aes(x = date)) + 
  geom_line(aes(y = (`w_w_chl`), col = "warm")) +
  geom_line(aes(y = `w_c_chl`, col = "cold")) +
  geom_line(aes(y = `w_f_chl`, col = "front")) +
  scale_color_manual(values = c("warm" = "red", "cold" = "blue","front" = "grey"), 
                     breaks = c("cold", "front", "warm")) +
  labs(y = expression(paste("Chl ", italic("a"), " mg ", m^{-3})), colour = "cluster") 
# total dazu in schwarz


# compare PFT compositon in clusters 

# rearrange data a bit

PFT_dat <- ts_dat %>%
  select(date, `mean PICO warm [mg / m^3]`:`mean MICRO front [mg / m^3]`) %>%
  pivot_longer(-date, names_to = c("type", "regime"), names_prefix = "mean ", names_sep = " ", values_to = "mean_vals")

chl_dat <- ts_dat %>%
  select(date, `mean chl warm [mg / m^3]`, `mean chl cold [mg / m^3]`, `mean chl front [mg / m^3]`) %>%
  pivot_longer(-date, names_to = c("type", "regime"), names_prefix = "mean ", names_sep = " ", values_to = "mean_vals")

# line plots
ggplot(data = PFT_dat, aes(x = date, y = mean_vals, col = type)) +
  geom_line() +
  facet_grid(.~regime) +
  scale_x_date(name = "",
               date_breaks = "1 month",
               date_labels = "%b",
               limits = as.Date(c("2018-04-01", "2018-09-15"))) 


# stacked area plot
ggplot(data = PFT_dat, aes(x = date, y = mean_vals, group = type)) +
  geom_area(aes(fill = type)) +
  facet_grid(.~regime) +
  scale_x_date(name = "",
               date_breaks = "1 month",
               date_labels = "%b",
               limits = as.Date(c("2018-04-01", "2018-09-15"))) +
geom_line(data = chl_dat)
# -> Markus sagt: pft algos sind noch nicht komplett ausgereift. Daher kein match total chl zu sum(pfts)
# ich konzentriere mich weiter auf total chl


# compare Ice in clusters 
ggplot(data = ts_dat, aes(x = date)) + 
  geom_line(aes(y = `ice cover warm region [%]`*100, col = "warm")) +
  geom_ribbon(aes(ymin = (`ice cover warm region [%]` - `ice cover warm region sd [%]`)*100, 
                  ymax = (`ice cover warm region [%]` + `ice cover warm region sd [%]`)*100, 
                  fill = "warm"), alpha = .3) +
  geom_line(aes(y = `ice cover cold region [%]`*100, col = "cold")) +
  geom_ribbon(aes(ymin = (`ice cover cold region [%]` - `ice cover cold region  sd [%]`)*100, 
                  ymax = (`ice cover cold region [%]` + `ice cover cold region  sd [%]`)*100, 
                  fill = "cold"), alpha = .3) +
  geom_line(aes(y = `ice cover front region [%]`*100, col = "front")) +
  geom_ribbon(aes(ymin = (`ice cover front region [%]` - `ice cover front region sd [%]`)*100, 
                  ymax = (`ice cover front region [%]` + `ice cover front region sd [%]`)*100, 
                  fill = "front"), alpha = .3) +
  geom_abline(intercept = 15, slope = 0, linetype = "dotted") + 
  labs(y = "sea ice cover [%]", colour = "cluster") +
  scale_x_date(name = "",
               date_breaks = "1 month",
               date_labels = "%b") +
  scale_y_continuous(breaks = c(0, 15, 25, 50, 75, 100)) +
  scale_color_manual(aesthetics =c("colour", "fill"), values = c("warm" = "red", "cold" = "blue","front" = "grey"), 
                     breaks = c("cold", "front", "warm")) +
  guides(fill = "none") +
  coord_cartesian(ylim=c(-0.7, 100), expand = F) +
  theme_bw()


