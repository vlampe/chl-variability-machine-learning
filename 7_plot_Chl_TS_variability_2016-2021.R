# analogous to 3_plot_Chl_TS_variability.R, plot median and variability of the long 
# (2016-2021) Chl time series. Do we see the same patterns of high/low variability
# in other years, too?

# show time series of Chl in the clusters, with uncertainties. 
# here, uncertainties are shown as 1) IQR and 2) 95% Confidence Intervals. 

library(tidyverse)
require(reticulate)
library(latex2exp)
reticulate::py_version()
use_python("/opt/homebrew/Cellar/python@3.11/3.11.5/Frameworks/Python.framework/Versions/3.11/Resources/Python.app/Contents/MacOS/Python")
np      <- import("numpy",   convert=FALSE)
py_csv  <- import("csv",     convert=FALSE)
diffKDE <- import("diffKDE", convert=FALSE)


# normalise data?
normalise <- FALSE # one of c(TRUE or FALSE)

# choose transformation of chl data
chl_trans = "no_trans" # one of c("log_trans", "no_trans")

# load Rds with raw chl data (created in 1_extract_clusters.R)
chl_list <- readRDS(file = "output/3k_raw_chl_conc_2016-2021.rds")

results <- tibble(`date` = NA,
                  `cluster` = NA,
                  `min` = NA,
                  `max` = NA,
                  `IQR_l` = NA,
                  `IQR_u` = NA,
                  `CI_l` = NA,
                  `CI_u` = NA,
                  `SD` = NA,
                  `mean` = NA, 
                  `median` = NA,
                  `mode` = NA, 
                  `n.obs` = NA
)[-1,]

var <- "CHL" 
regimes <- c("warm", "front", "cold")

max_d <- lapply(chl_list$CHL, length)[[1]]

# enter daily loop
for(d in c(1:max_d)){
  doy <- d
  datum <- as.Date(d-1, origin = "2016-01-01") # sonst passt zählung nicht
  
  # enter cluster loop
  for(clus in regimes){
    
    # extract values from list 
    chl_vals <- chl_list[[var]][[clus]][[d]]
    
    # if chl_vals is empty, skip and attach row to results table
    if (length(chl_vals) == 0){
      out <- tibble(`date` = datum,
                    `cluster` = clus,
                    `min` = NA,
                    `max` = NA,
                    `IQR_l` = NA,
                    `IQR_u` = NA,
                    `CI_l` = NA,
                    `CI_u` = NA,
                    `SD` = NA,
                    `mean` = NA, 
                    `median` = NA,
                    `mode` = NA, 
                    `n.obs` = 0
      )
    } else { # if chl_vals is not empty: do this
      
      
      # normalise?
      if (normalise == TRUE){
        norm_min <- 0.15
        chl_vals = pmax(norm_min, chl_vals)/norm_min
      }
      
      # log-transform?
      if (chl_trans=="log_trans"){
        chl_vals = log10(chl_vals)
      } 
      
      
      
      # minimum and maximum
      min_val = min(chl_vals)
      max_val = max(chl_vals)
      
      # IQR 
      #Q1 <-  quantile(chl_vals)[[2]]
      #Q3 <-  quantile(chl_vals)[[4]]
      
      # mean, median, mode ( BESSER AUS cdfs ermitteln )
      mean_val <- mean(chl_vals)
      #mode_val <- DescTools::Mode(chl_vals)
      #median_val <- median(chl_vals)
      sd_val <- sd(chl_vals)
      
      # diffusion KDE of sample data (min is 0. max is 60, bc the max value for chl pixel was 57ish)
      # n set to 6000, so step size is ~0.01 (this is required I think for Chl concentrations) but leads to regular failing of KDE (not sure why exactly... )
      kde_min = 0.01
      kde_max = 60
      
      if (normalise == TRUE){
        kde_min = pmax(norm_min, kde_min)/norm_min
        kde_max = pmax(norm_min, kde_max)/norm_min
      }
      
      if (chl_trans == "log_trans"){
        kde_min = log10(kde_min)
        kde_max = log10(kde_max)
      }
      
      
      
      # estimate PDF
      pdf_list <-  py_to_r(diffKDE$KDE(np$array(chl_vals), kde_min, kde_max, as.integer(6000))) 
      # diffusion KDE of sample data (min is 0. max is 60, bc the max value for chl pixel was 57ish)
      # n set to 6000, so step size is ~0.01 (this is required I think for Chl concentrations) but leads to regular failing of KDE (not sure why exactly... )
      
      pdf_x <- pdf_list[[2]] # convert numpy array to R array
      pdf_y <- pdf_list[[1]]
      
      # # fix failing KDEs: set max based on maximum value, but keep stepsize at 0.01. since this also not always works, also try diffKDE with different parameters. 
      # if(any(is.nan(pdf_y))){ # try again with less steps, if KDE failed
      #   
      #   try(pdf_list <- py_to_r(diffKDE$KDE(np$array(chl_vals), kde_min, kde_max)))
      #   try(pdf_list <-  py_to_r(diffKDE$KDE(np$array(chl_vals))))
      #   try(pdf_list <-  py_to_r(diffKDE$KDE(np$array(chl_vals), xmin = min_val-(0.5*min_val)))) 
      #   try(pdf_list <-  py_to_r(diffKDE$KDE(np$array(chl_vals), kde_min, round(max_val+0.5*max_val,2))))
      #   
      #   
      #   pdf_x <- pdf_list[[2]] # convert numpy array to R array
      #   pdf_y <- pdf_list[[1]]
      # }
      
      # fix failing KDEs: set max based on maximum value, but keep stepsize at 0.01. since this also not always works, also try diffKDE with different parameters. 
      success <- FALSE # since any(is.nan(pdf_y)) == TRUE
      i <- 0
      kde_try_list <-  (list(try(pdf_list <- py_to_r(diffKDE$KDE(np$array(chl_vals), kde_min, kde_max))),
                             try(pdf_list <-  py_to_r(diffKDE$KDE(np$array(chl_vals), xmin = min_val-(0.5*min_val)))),
                             try(pdf_list <-  py_to_r(diffKDE$KDE(np$array(chl_vals), kde_min, round(max_val+0.5*max_val,2)))),
                             try(pdf_list <-  py_to_r(diffKDE$KDE(np$array(chl_vals))))
                             ))
      while (!success) {
        i <- i+1
        # do something: search through kde_try_list
        
        if( typeof(kde_try_list[[i]][[1]])!="double" ) { next } # skip if error prevented calculation
        
        pdf_x <- kde_try_list[[i]][[2]] # convert numpy array to R array
        pdf_y <- kde_try_list[[i]][[1]]
        
        # check for success
        success <- (any(is.nan(pdf_y)) == FALSE)
        
      }
      
      
      if(any(is.nan(pdf_y))){
        print(d) 
        error("oh no") # this is not a real function but causes the loop to break, if pdf_y is still NA... 
      }
      #plot(pdf_x, pdf_y, 'l', col="red", xlim = c(0,4) )
      
      
      # calculate CDF
      h <-  diff(pdf_x[1:2])
      cdf <- cumsum(pdf_y)*h
      #lines(pdf_x, cdf, 'l', col="blue")
      
      # calc statistics from CDF
      # in y vector of cdf, find position where y is closest to .025 and .975 for CI,
      # (0.5 for median), 
      # 0.25 and 0.75 for Q1 and Q3
      
      pos <- which.min(abs(cdf - 0.5))
      median_val <- pdf_x[pos]
      # abline(v = pdf_x[pos])
      
      pos <- which.min(abs(cdf - 0.025))
      CI_l <- pdf_x[pos]
      pos <- which.min(abs(cdf - 0.975))
      CI_u <- pdf_x[pos]
      
      pos <- which.min(abs(cdf - 0.25))
      Q1 <- pdf_x[pos]
      pos <- which.min(abs(cdf - 0.75))
      Q3 <- pdf_x[pos]
      
      # mode = where is pdf max?
      pos <- which.max(pdf_y)
      mode_val <- pdf_x[pos]
      
      # and save calculated values 
      out <- tibble(`date` = datum,
                    `cluster` = clus,
                    `min` = min_val,
                    `max` = max_val,
                    `IQR_l` = Q1,
                    `IQR_u` = Q3,
                    `CI_l` = CI_l,
                    `CI_u` = CI_u,
                    `SD` = sd_val,
                    `mean` = mean_val, 
                    `median` = median_val,
                    `mode` = mode_val, 
                    `n.obs` = length(chl_vals)
      )
    }
    
    # append out to results table 
    results <- rbind(results, out)
    
  }
}
# save output
write_csv2(x = results %>% arrange(`date`) , file = paste0("output/Chl_norm", normalise, "_", chl_trans, "_timeseries_variability_3k_2016-2021.csv"))

# ______________________________________________________________________________
## PLOTTING

# normalised data?
normalise <- FALSE # one of c(TRUE or FALSE)

# choose transformation of chl data
chl_trans = "no_trans" # one of c("log_trans", "no_trans")
#...............................................................................

yunit <- ifelse(normalise, "no dim", "mg m$^{-3}$")



results <- read_csv2(paste0("output/Chl_norm", normalise, "_", chl_trans, "_timeseries_variability_3k_2016-2021.csv"))

ggplot(data = results, aes(x = `date`, y= `median`, colour = `cluster`)) + 
  geom_line() +
  geom_line(aes(y = `mean`), linetype = "dotted") +
  geom_ribbon(aes(ymin = `IQR_l`, 
                  ymax = `IQR_u`, 
                  fill = `cluster`), alpha = .3, linetype = 0) +
  geom_ribbon(aes(ymin = `CI_l`, 
                  ymax = `CI_u`, 
                  fill = `cluster`), alpha = .1, linetype = 0) +
  scale_x_date(name = "",
               date_breaks = "3 month",
               date_minor_breaks = "1 month",
               date_labels = "%b %y") +
  #scale_y_continuous(breaks = c(0.1, 1, 2, 4, 8, 12, 16)) +
  # labs(y = expression(paste("Chl ", italic("a"), " mg ", m^{-3})), colour = "cluster") +
  labs(y = TeX(paste0("Chl \\textit{a}", " [", yunit, "]")), colour = "cluster",
       title = paste0("normalised: ", normalise, " trans: ", chl_trans)) +
  scale_color_manual(aesthetics =c("colour", "fill"), values = c("warm" = "red", "cold" = "blue","front" = "cornsilk4"), 
                     breaks = c("cold", "front", "warm")) +
  # guides(fill = "none") +
  # coord_cartesian(ylim=c(0, 4), xlim = as.Date(c("2018-03-01", "2018-10-15")), expand = F) +
  #coord_trans(y = "log10", ylim=c(0.1, NA), xlim = as.Date(c("2018-03-01", "2018-10-15")), expand = F) +
  #annotation_logticks(scaled = F) +
#  facet_grid(cluster~.) +
  theme_bw()
ggsave(paste0("output/plots/Chl_norm", normalise, "_", chl_trans, "_timeseries_variability_3k_2016-2021.pdf"))


ggplot(data = results, aes(x = `date`, y = `IQR_u` - `IQR_l`, colour = `cluster`)) +
  geom_line() +
  labs(x = "Date", y = TeX(paste0("IQR", " [", yunit, "]")), title = paste0("normalised: ", normalise, " trans: ", chl_trans)) +
  #facet_grid(cluster~.) +
  scale_x_date(name = "",
               date_breaks = "1 month",
               date_labels = "%b") +
  scale_color_manual(aesthetics =c("colour", "fill"), values = c("warm" = "red", "cold" = "blue","front" = "cornsilk4"), 
                     breaks = c("cold", "front", "warm")) +
  coord_cartesian(ylim=c(0, 7), xlim = as.Date(c("2016-03-01", "2021-10-15")), expand = F) +
  theme_bw()
ggsave(paste0("output/plots/Chl_norm", normalise, "_", chl_trans, "_timeseries_IQR_3k_2016-2021.pdf"))


# plot number of observations
ggplot(data = results, aes(x = `date`)) + 
  geom_line(aes(y= `n.obs`)) +
  scale_x_date(name = "",
               date_breaks = "1 month",
               date_labels = "%b") +
  facet_grid(cluster~.) +
  theme_bw()


# plot IQR and number of observations
ggplot(data = results, aes(x = `date`)) + 
  geom_line(aes(y= (`IQR_u`-`IQR_l`)*10000)) +
  geom_line(aes(y = `n.obs`), linetype = "dotted", color = "red") + 
  scale_x_date(name = "",
               date_breaks = "1 month",
               date_labels = "%b") +
  facet_grid(cluster~.) +
  theme_bw()
# in some years, n.obs is very high when iqr is small



## _____ new plots for long ts

results2 <- results %>%
  mutate(`year` = year(`date`),
         year_date = as.Date(paste(day(`date`),month(`date`),"2000", sep = "-"), format = "%d-%m-%Y"),
         IQR_l = pmax(IQR_l, 0.01))


p <- ggplot(data = results2, aes(x = `year_date`, y= `median`, colour = `cluster`)) + 
  geom_line(linewidth = 0.2) +
  geom_ribbon(aes(ymin = `IQR_l`, 
                  ymax = `IQR_u`, 
                  fill = `cluster`), alpha = .1, linetype = 0) +
  scale_x_date(name = "",
               date_breaks = "1 month",
               date_labels = "%b") +
  scale_y_continuous(breaks = c(0.1, 0.5, 1, 2, 4, 8)) +
  labs(y = TeX(paste0("Chl \\textit{a}", " [", yunit, "]")), colour = "cluster",
       title = paste0("normalised: ", normalise, " trans: ", chl_trans)) +
  scale_color_manual(aesthetics =c("colour", "fill"), values = c("warm" = "red", "cold" = "blue","front" = "cornsilk4"), 
                     breaks = c("cold", "front", "warm")) +
#  scale_color_brewer(aesthetics = c("color", "fill"), type = "qual", palette = "Dark2", 
#                     limits = c("front", "warm", "cold"))+
#  coord_cartesian(ylim=c(0, 4), xlim = as.Date(c("2000-03-01", "2000-10-15")), expand = F) +
  coord_trans(y = "log10", ylim=c(0.1, NA), xlim = as.Date(c("2000-03-01", "2000-10-15")), expand = F, clip = "on")+
  annotation_logticks(scaled = F) +
  facet_grid(year~.) +
  theme_bw() +
  theme(legend.position = "bottom")
p + labs(title = NULL)
ggsave(paste0("output/plots/Chl_norm", normalise, "_", chl_trans, "_stack_timeseries_IQR_3k_2016-2021.pdf"),
       width = 15, height = 16, units = "cm")

## revised plot, after comments of R1 and R2: 

## edit with reviewers comments 
tsdat <- read_csv2("output/Output_3k_cluster_analysis_2016-2021.csv") %>%
  dplyr::select(c("julian day", "date", "ice cover warm region [%]", 
           "ice cover front region [%]", "ice cover cold region [%]")) %>%
  mutate(`year` = year(`date`),
         `month` = month(`date`),
         year_date = as.Date(paste(day(`date`),month(`date`),"2000", sep = "-"), format = "%d-%m-%Y")) %>%
  pivot_longer(cols = c("ice cover warm region [%]", "ice cover front region [%]", 
                        "ice cover cold region [%]"),
               names_pattern = "\\b(warm|cold|front)\\b",
               names_to = "cluster", 
               values_to = "ice cover [%]"
  ) %>%
  group_by(month, year, cluster) %>%
  summarise(avg_ice = mean(`ice cover [%]`)*100, 
            year_date = first(year_date)) 


p_rev1 <- ggplot(data = results2, aes(x = `year_date`, y= `median`, colour = `cluster`)) + 
  geom_abline(intercept = 1, slope = 0, linetype = "solid", color = "gold", linewidth = 0.3) +
  geom_ribbon(aes(ymin = `IQR_l`, 
                  ymax = `IQR_u`, 
                  fill = `cluster`), alpha = .1, linetype = 0) +
  geom_line(linewidth = 0.2) +
  geom_line(aes(y=`mean`), linetype = "dotted") +
  geom_text(data = filter(tsdat, cluster == "cold"), aes(x = year_date+11, y = 4.5, 
                              label = signif(avg_ice, digits=2),
                              colour = cluster), show.legend = F, size = 3, vjust = 0)+ 
  geom_text(data = filter(tsdat, cluster == "front"), 
            aes(x = year_date+20, y = 4.5, label = signif(round(avg_ice,1), digits=2), colour = cluster),
            show.legend = F, size = 3, vjust = 0)+ 
  scale_x_date(name = "",
               date_breaks = "1 month",
               date_labels = "%b") +
  scale_y_continuous(breaks = c(0.1, 0.5, 1, 2, 4, 8)) +
  labs(y = TeX(paste0("Chl \\textit{a}", " [", yunit, "]")), colour = "cluster",
       title = paste0("normalised: ", normalise, " trans: ", chl_trans)) +
  scale_color_manual(aesthetics =c("colour", "fill"), values = c("warm" = "red", "cold" = "blue","front" = "cornsilk4"), 
                     breaks = c("cold", "front", "warm")) +
  #  scale_color_brewer(aesthetics = c("color", "fill"), type = "qual", palette = "Dark2", 
  #                     limits = c("front", "warm", "cold"))+
  #  coord_cartesian(ylim=c(0, 4), xlim = as.Date(c("2000-03-01", "2000-10-15")), expand = F) +
  coord_trans(y = "log10", ylim=c(0.1, 8), xlim = as.Date(c("2000-03-01", "2000-10-09")), expand = F, clip = "on")+
  annotation_logticks(scaled = F) +
  facet_grid(year~.) +
  theme_bw() +
  theme(legend.position = "bottom")
p_rev1 + labs(title = NULL)
ggsave(paste0("output/plots/Chl_norm", normalise, "_", chl_trans, "_stack_timeseries_IQR_3k_2016-2021_rev1.pdf"),
       width = 15, height = 16, units = "cm")


library(colorblindcheck)
palette_check(rainbow(n=3),
              , plot = TRUE)

RColorBrewer::brewer.pal(3,"Set1")
palette_check(RColorBrewer::brewer.pal(3,"Dark2"),
              , plot = TRUE)

palette_check(RColorBrewer::brewer.pal(3,"Set1"),
              , plot = TRUE)

## anomalies from yearly mean

out <- results2 %>%
  group_by(year, cluster) %>%
  summarise(avg_year = mean(`mean`, na.rm = T)) %>%
  spread(year, avg_year) %>%
  mutate_at(2:7, round, 2)

results3 <- results2 %>%
  group_by(year, cluster) %>%
  mutate(avg_year = mean(`mean`, na.rm = T),
         anomaly = mean - avg_year,
         year = factor(year, levels = c(2016:2021)))
# attention: these are not yearly averages, but averages over the period RS data is available! Mid-March -- October

p <- ggplot(results3, aes(year_date, anomaly, color = year)) +
  geom_line(linewidth = 0.4) +
  scale_color_manual(values = RColorBrewer::brewer.pal(6, "Dark2")) +
    coord_cartesian(xlim = as.Date(c("2000-03-01", "2000-10-15")), expand = F) +
  scale_x_date(name = "",
               date_breaks = "1 month",
               date_labels = c("M", "A", "M", "J", "J", "A", "S", "O")) + # %b for 3-letter-abbrv.
  labs(x = "Month", y = TeX(paste0("Chl \\textit{a} anomaly [", yunit, "]")),
       color = "year") +
  facet_wrap(~cluster) +
  theme_bw() +
  theme(legend.position = "bottom")
p
ggsave(paste0("output/plots/Chl_norm", normalise, "_", chl_trans, "_anomaly_yearly_mean_3k_2016-2021.pdf"),
       width = 15, height = 7.5, units = "cm")


p_rev1 <- ggplot(results3, aes(year_date, anomaly, color = year)) +
  geom_line(linewidth = 0.4) +
  scale_color_manual(values = RColorBrewer::brewer.pal(6, "Dark2")) +
  coord_cartesian(xlim = as.Date(c("2000-03-01", "2000-10-15")), expand = F) +
  scale_x_date(name = "",
               date_breaks = "1 month",
              # date_labels = c("M", "A", "M", "J", "J", "A", "S", "O")) + # %b for 3-letter-abbrv.
              date_labels = "%b") +
  labs(x = "Month", y = TeX(paste0("Chl \\textit{a} anomaly [", yunit, "]")),
       color = "year") +
  facet_grid(cluster~.) +
  theme_bw() +
  theme(legend.position = "bottom",
        legend.direction = "horizontal")
p_rev1 + 
  labs(title = NULL) +
  guides(color=guide_legend(nrow=1, title.position = "left", title = "Year"))
ggsave(paste0("output/plots/Chl_norm", normalise, "_", chl_trans, "_anomaly_yearly_mean_3k_2016-2021_rev1.pdf"),
       width = 15, height = 14, units = "cm")

# save vor Verteidigung
p + theme(text = element_text(color = "white"),
        panel.background = element_rect(fill = NA, color = "white"),
        panel.grid = element_line(color = "grey30"),
        panel.border = element_rect(color = "white"),
        plot.background = element_rect(fill = "black", color = "black"),
        axis.text = element_text(color = "white"),
        axis.ticks = element_line(color = "white"),
        legend.background = element_rect(fill = "black"),
        strip.background = element_rect(fill = NA, color = "white"),
        strip.text = element_text(color = "white"))
ggsave(paste0("~/Documents/Verteidigung/plots/figures_paper3/Chl_norm", normalise, "_", chl_trans, "_anomaly_yearly_mean_3k_2016-2021.pdf"),
       width = 15, height = 7.5, units = "cm")



# new plot: day of year of max chla throughout the years

test <- results3 %>%
  mutate(doy = julian(year_date, origin=as.Date("2000-01-01"))+1) %>%
  ungroup() %>%
  slice_max(mean, by=c(cluster, year))
# this results in 2020 finding the max in autumn (based on 4 pixels), which is probably just an artifact; therefore exclude sept 2020 first!

test <- results3 %>%
  mutate(doy = julian(year_date, origin=as.Date("2000-01-01"))+1) %>%
  filter(doy < 260) %>%
  ungroup() %>%
  slice_max(mean, by=c(cluster, year))

maxplot <- ggplot(data=test, aes(x=year_date, y=year, color=cluster, shape = cluster)) +
  geom_point(size = 2) +
  labs(y=NULL)+
  geom_text(aes(label=round(mean,2)), nudge_y = 0.25, show.legend=F, size = 3) +
  scale_color_manual(aesthetics =c("colour"), values = c("warm" = "red", "cold" = "blue","front" = "cornsilk4"), 
                     breaks = c("cold", "front", "warm")) +
  scale_x_date(name = "",
               date_breaks = "2 week", date_minor_breaks = "1 day",
               date_labels = "%b %d") +
  theme_bw()
maxplot
# cool


# when does the mean Chl-a first cross the bloom threshold?
test2 <- results3 %>%
  mutate(doy = julian(year_date, origin=as.Date("2000-01-01"))+1) %>%
  filter(mean >= 1) %>%
  ungroup() %>%
  slice_min(doy, by=c(cluster, year))

firstplot <- ggplot(data = test2, aes(x=year_date, y=year, color=cluster, shape = cluster)) +
  geom_point(size=2) +
  labs(y=NULL)+
  scale_color_manual(aesthetics =c("colour"), values = c("warm" = "red", "cold" = "blue","front" = "cornsilk4"), 
                     breaks = c("cold", "front", "warm")) +
  scale_x_date(name = "",
               date_breaks = "1 month", date_minor_breaks = "1 day",
               date_labels = "%b %d") +
  theme_bw()
  # geom_jitter(width=1,height = 0)
firstplot

library(ggpubr)
ggarrange(firstplot, maxplot, common.legend = T, legend = "bottom", align = "hv", labels = "auto")
ggsave(paste0("output/plots/First_and_max_bloom_per_cluster_2016-2021.pdf"),
       height = 8, width = 15, units = "cm",  device = pdf)

  