## plot time series
## author: Vanessa Lampe
## date: 2020-04-01
## version: 1.2
## git: https://git.geomar.de/vanessa-lampe/chlorophyll-time-series

library(tidyverse)
library(readxl)
library(writexl)
library(ggpubr)
source("weighted_moving_average.R")


## read output from extract_warm_and_cold_clusters.R
ts_dat <- read_excel("Output_warm_and_cold_cluster_analysis.xlsx")


# plot daily chl in whole (total) sampling area:
ggplot(data = ts_dat, aes(x = as.Date(date), y = `mean chl total [mg / m^3]`)) +
  geom_line(aes(color = "g")) +
  scale_x_date(name = "",
               date_breaks = "1 year",
               date_labels = "%Y",
               date_minor_breaks = "1 month") +
  scale_y_continuous(limits = c(0, NA), minor_breaks = F, name = "Chl a conc. [µg / L]") +
  scale_color_manual(values = c("g" = "darkgreen", "w" = "red","c" = "blue"), 
                     breaks = c("g", "c", "w"),labels = c("g" = "total", 
                                                          "w" = "warm",
                                                          "c" = "cold"), 
                     name = "regime", drop = F) +
  theme_bw() +
  theme(legend.position=c(1, 1), 
        legend.justification = c(1,1), 
        legend.box.background = element_rect(fill = NA, color = "black", 
                                             size = 0.7)
  ) ->chltotday
chltotday

#plot daily chl in warm and cold sections
ggplot(data = ts_dat, aes(x = as.Date(date))) +
  geom_line(aes(y = `mean chl warm [mg / m^3]`, color = "w")) +
  geom_line(aes(y = `mean chl cold [mg / m^3]`, color  = "c")) +
  scale_x_date(name = "",
               date_breaks = "1 year",
               date_labels = "%Y",
               date_minor_breaks = "1 month") +
  scale_y_continuous(limits = c(0, NA), name = "Chl a conc. [µg / L]", 
                     minor_breaks = F) +
  scale_color_manual(values = c("g" = "darkgreen", "w" = "red","c" = "blue"), 
                     breaks = c("g", "c", "w"),labels = c("g" = "total", 
                                                          "w" = "warm",
                                                          "c" = "cold"), 
                     name = "regime", drop = F) +
  theme_bw() +
  theme(legend.position=c(1, 1), 
        legend.justification = c(1,1), 
        legend.box.background = element_rect(fill = NA, color = "black", 
                                             size = 0.7)
  ) -> chlwcday
chlwcday

# plot daily SST in warm and cold sections
ggplot(data = ts_dat, aes(x = as.Date(date))) +
  geom_line(aes(y = `mean temp cold [ºC]`, color = "c")) +
  geom_line(aes(y = `mean temp warm [ºC]`, color = "w")) +
  scale_x_date(name = "",
               date_breaks = "1 year",
               date_labels = "%Y",
               date_minor_breaks = "1 month", 
               expand = c(0,0)) +
  scale_y_continuous(name = "mean temp [ºC]", minor_breaks = F) +
  scale_color_manual(values = c("g" = "darkgreen", "w" = "red","c" = "blue"), 
                     breaks = c("g", "c", "w"), labels = c("g" = "total", 
                                                           "w" = "warm",
                                                           "c" = "cold"), 
                     name = "regime", drop = F) +
  theme_bw() +
  theme(legend.position=c(1, 1), 
        legend.justification = c(1,1), 
        legend.box.background = element_rect(fill = NA, color = "black", 
                                             size = 0.7)
  ) -> tempwcday
tempwcday

ggarrange(chltotday +  scale_color_manual(limits = c("g", "w", "c"), values = c("g" = "darkgreen", "w" = "red","c" = "blue"), breaks = c("g", "c", "w"),
                                          labels = c("g" = "total", "w" = "warm","c" = "cold"), name = "regime"),  
          chlwcday, tempwcday, nrow = 3, common.legend = T,
          align = "v") %>%
  ggexport(filename = "figures/daily_chl_sst.pdf")
# vielleicht noch shades?

#_______________________________________________________________________________
# look at weekly means
weekly_means <- ts_dat %>%
  mutate(weekday = weekdays(date)) %>%
  group_by(Week = lubridate::floor_date(date, unit = "week")) %>%
  summarise(min_date  = min(date), max_date = max(date),
            `mean temp warm [ºC]` = mean(`mean temp warm [ºC]`, na.rm = T),
            `mean temp cold [ºC]` = mean(`mean temp cold [ºC]`, na.rm = T),
            `mean chl total [mg / m^3]` = mean(`mean chl total [mg / m^3]`, na.rm = T),
            `mean chl warm [mg / m^3]` = mean(`mean chl warm [mg / m^3]`, na.rm = T),
            `min chl warm [mg / m^3]` = min(`min chl warm [mg / m^3]`, na.rm = T), 
            `max chl warm [mg / m^3]` = max(`max chl warm [mg / m^3]`, na.rm = T),
            `mean chl cold [mg / m^3]`  = mean(`mean chl cold [mg / m^3]`, na.rm = T),
            `min chl cold [mg / m^3]` = min(`min chl cold [mg / m^3]`, na.rm = T),
            `max chl cold [mg / m^3]` = max(`max chl cold [mg / m^3]`, na.rm = T)
            )

# plot weekly means of SST in warm and cold sections
ggplot(data = weekly_means, aes(x = as.Date(Week))) +
  scale_x_date(name = "",
               date_breaks = "1 year",
               date_labels = "%Y",
               date_minor_breaks = "1 month", expand = c(0,0)) +
  scale_y_continuous(name = "mean temp [ºC]", minor_breaks = F) +
  geom_line(aes(y = `mean temp warm [ºC]`, color = "w")) +
  geom_line(aes(y = `mean temp cold [ºC]`, color = "c")) +
  scale_color_manual(values = c("w" = "red","c" = "blue"), breaks = c("c", "w"),
                     labels = c("w" = "warm","c" = "cold"), name = "regime") +
  theme_bw() +
  theme(legend.position=c(1, 1), 
        legend.justification = c(1,1), 
        legend.box.background = element_rect(fill = NA, color = "black", 
                                             size = 0.7)
  ) -> tempwcweek
tempwcweek

# plot weekly means of chl in warm and cold section
ggplot(data = weekly_means, aes(x = as.Date(Week))) +
  scale_x_date(name = "",
               date_breaks = "1 year",
               date_labels = "%Y",
               date_minor_breaks = "1 month") +
  scale_y_continuous(name = "Chl a conc. [µg / L]", minor_breaks = F) +
  geom_line(aes(y = `mean chl warm [mg / m^3]`, color = "w")) +
  geom_line(aes(y = `mean chl cold [mg / m^3]`, color = "c")) +
  scale_color_manual(values = c("w" = "red","c" = "blue"), breaks = c("c", "w"),
                     labels = c("w" = "warm","c" = "cold"), name = "regime") +
  theme_bw() +
  theme(legend.position=c(1, 1), 
        legend.justification = c(1,1), 
        legend.box.background = element_rect(fill = NA, color = "black", 
                                             size = 0.7)
  ) -> chlwcweek
chlwcweek

# plot weekly means of chl in total sampling area
ggplot(data = weekly_means, aes(x = as.Date(Week))) +
  scale_x_date(name = "",
               date_breaks = "1 year",
               date_labels = "%Y",
               date_minor_breaks = "1 month") +
  scale_y_continuous(name = "Chl a conc. [µg / L]", minor_breaks = F) +
  geom_line(aes(y = `mean chl total [mg / m^3]`, color = "g")) +
  scale_color_manual(values = c("g" = "darkgreen"), breaks = c("g"),
                     labels = c("g" = "total"), name = "regime") +
  theme_bw() +
  theme(legend.position=c(1, 1), 
        legend.justification = c(1,1), 
        legend.box.background = element_rect(fill = NA, color = "black", size = 0.7)
  ) -> chltotweek
chltotweek

ggarrange(chltotweek+  scale_color_manual(limits = c("g", "w", "c"), values = c("g" = "darkgreen", "w" = "red","c" = "blue"), breaks = c("g", "c", "w"),
                                          labels = c("g" = "total", "w" = "warm","c" = "cold"), name = "regime", drop = F),  
          chlwcweek, tempwcweek, nrow = 3, common.legend = T, align = "v") %>%
  ggexport(filename = "figures/weekly_chl_sst.pdf")

#_______________________________________________________________________________
# look at weighted moving averages
wmv <- weighted_moving_average   # reassign function bc lazyness

wmeans <- ts_dat
wmeans$w_w_chl <- wmv(values = wmeans$`mean chl warm [mg / m^3]`,
                      weights = wmeans$`n.obs. chl warm [pixels]`,
                      width = 15)
wmeans$w_c_chl <- wmv(values = wmeans$`mean chl cold [mg / m^3]`, 
                      weights = wmeans$`n.obs chl cold [pixels]`,
                      width = 15)
wmeans$tot_chl_nobs <- base::rowSums(cbind(wmeans$`n.obs. chl warm [pixels]`, 
                                           wmeans$`n.obs chl cold [pixels]`), na.rm = T)
wmeans$tot_chl_nobs[wmeans$tot_chl_nobs == 0] <- NA
wmeans$w_tot_chl <- wmv(values = wmeans$`mean chl total [mg / m^3]`, 
                        weights = wmeans$tot_chl_nobs,
                        width = 15)

# plot weighed (by n.obs) moving average of chl in warm and cold section
ggplot(wmeans, aes(x = as.Date(date))) +
  scale_x_date(name = "",
               date_breaks = "1 year",
               date_labels = "%Y",
               date_minor_breaks = "1 month", expand = c(0,0)) +
  scale_y_continuous(name = "Chl a conc. [µg / L]", minor_breaks = F) +
  geom_line(aes(y = w_c_chl, color = "c")) +
  geom_line(aes(y = w_w_chl, color = "w")) +
  scale_color_manual(values = c("w" = "red","c" = "blue"), breaks = c("c", "w"),
                     labels = c("w" = "warm","c" = "cold"), name = "regime") +
  theme_bw() +
  theme(legend.position=c(1, 1), 
       legend.justification = c(1,1), 
       legend.box.background = element_rect(fill = NA, color = "black", size = 0.7)
       )  -> wmeans_wc
wmeans_wc
ggsave(wmeans_wc, filename = "figures/weighed_moving_average_chl_wc.pdf")

# plot weighed (by n.obs) moving average of chl in total sampling area
ggplot(wmeans, aes(x = as.Date(date))) +
  scale_x_date(name = "",
               date_breaks = "1 year",
               date_labels = "%Y",
               date_minor_breaks = "1 month") +
  scale_y_continuous(name = "Chl a conc. [µg / L]", minor_breaks = F) +
  geom_line(aes(y = w_tot_chl, color = "g")) +
  scale_color_manual(values = c("g" = "darkgreen"), labels = c("g" = "total"),
                     name = "regime") +
  theme_bw() +
  theme(legend.position = c(1,1),
        legend.justification = c(1,1),
        legend.box.background = element_rect(fill = NA, color = "black", size = 0.7)) -> wmeans_chl_total
wmeans_chl_total

ggarrange(wmeans_chl_total +  scale_color_manual(limits = c("g", "w", "c"), values = c("g" = "darkgreen", "w" = "red","c" = "blue"), breaks = c("g", "c", "w"),
                                                 labels = c("g" = "total", "w" = "warm","c" = "cold"), name = "regime"),
          wmeans_wc, tempwcweek, nrow = 3, common.legend = T, align = "v") %>%
  ggexport(filename = "figures/wmv_chl_sst.pdf")
