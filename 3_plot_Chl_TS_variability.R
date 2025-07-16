## plot more complex Chl time series
## author: Vanessa Lampe
## date: 2023-10-23
## version: 0
## git: 

# show time series of Chl in the clusters, with uncertainties. 
# here, uncertainties are shown as 1) IQR and 2) 95% Confidence Intervals. 
#setwd("Documents/microARC/Manuscripts/Last_chapter_clustering_manuskript_skripts/skripts/")

library(tidyverse)
require(reticulate)
library(latex2exp)

library(extrafont)
library(Cairo)

reticulate::py_version()
use_python("/opt/homebrew/Cellar/python@3.11/3.11.11/Frameworks/Python.framework/Versions/3.11/Resources/Python.app/Contents/MacOS/Python")
np      <- import("numpy",   convert=FALSE)
py_csv  <- import("csv",     convert=FALSE)
diffKDE <- import("diffKDE", convert=FALSE)


# normalise data?
normalise <- FALSE # one of c(TRUE or FALSE)

# choose transformation of chl data
chl_trans = "log_trans" # one of c("log_trans", "no_trans")

# load Rds with raw chl data (created in 1_extract_clusters.R)
chl_list <- readRDS(file = "output/3k_raw_chl_conc.rds")

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

# enter daily loop
for(d in c(1:365)){
  doy <- d
  datum <- as.Date(d-1, origin = "2018-01-01") # sonst passt zählung nicht
  
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
        chl_vals = log10(chl_vals/0.1) # normalise by down-rounded minimum Chla conc in data
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
        kde_min = log10(kde_min/0.1)
        kde_max = log10(kde_max/0.1)
      }
      
      
      
      # estimate PDF
      pdf_list <-  py_to_r(diffKDE$KDE(np$array(chl_vals), kde_min, kde_max, as.integer(6000))) 
      # diffusion KDE of sample data (min is 0. max is 60, bc the max value for chl pixel was 57ish)
      # n set to 6000, so step size is ~0.01 (this is required I think for Chl concentrations) but leads to regular failing of KDE (not sure why exactly... )

      pdf_x <- pdf_list[[2]] # convert numpy array to R array
      pdf_y <- pdf_list[[1]]
      
      # fix failing KDEs: set max based on maximum value, but keep stepsize at 0.01. since this also not always works, also try diffKDE with different parameters. 
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
      
      while(any(is.nan(pdf_y))){ # try again with less steps, if KDE failed
        
        try(pdf_list <- py_to_r(diffKDE$KDE(np$array(chl_vals), kde_min, kde_max)))
          if(!any(is.nan(pdf_list[[1]]))) break
        try(pdf_list <-  py_to_r(diffKDE$KDE(np$array(chl_vals))))
          if(!any(is.nan(pdf_list[[1]]))) break
        try(pdf_list <-  py_to_r(diffKDE$KDE(np$array(chl_vals), xmin = min_val-(0.5*min_val)))) 
          if(!any(is.nan(pdf_list[[1]]))) break
        try(pdf_list <-  py_to_r(diffKDE$KDE(np$array(chl_vals), kde_min, round(max_val+0.5*max_val,2))))
          if(!any(is.nan(pdf_list[[1]]))) break
        } 
      pdf_x <- pdf_list[[2]] # convert numpy array to R array
      pdf_y <- pdf_list[[1]]
     

      
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
write_csv2(x = results %>% arrange(`date`) , file = paste0("output/Chl_norm", normalise, "_", chl_trans, "_timeseries_variability_3k.csv"))

# ______________________________________________________________________________
## PLOTTING

# normalised data?
normalise <- FALSE # one of c(TRUE or FALSE)

# choose transformation of chl data
chl_trans = "no_trans" # one of c("log_trans", "no_trans")
#...............................................................................

yunit <- ifelse(normalise, "no dim", "mg m$^{-3}$")

# but log transformation requires the value to be unitless. therefore we assume:  
# If 𝑥=0.5 is measured in some units, say, seconds, then taking the log actually means ln(0.5𝑠/1𝑠)=ln(0.5)
if (chl_trans == "log_trans"){yunit = "no dim"}



results <- read_csv2(paste0("output/Chl_norm", normalise, "_", chl_trans, "_timeseries_variability_3k.csv"))

p <- ggplot(data = results, aes(x = `date`, y= `median`, colour = `cluster`)) + 
  geom_line() +
  geom_line(aes(y = `mean`), linetype = "dotted") +
  geom_ribbon(aes(ymin = `IQR_l`, 
                  ymax = `IQR_u`, 
                  fill = `cluster`), alpha = .6, linetype = 0) +
  geom_ribbon(aes(ymin = `CI_l`, 
                  ymax = `CI_u`, 
                  fill = `cluster`), alpha = .4, linetype = 0) +
  scale_x_date(name = "",
               date_breaks = "1 month",
               date_labels = "%b") +
  scale_y_continuous(breaks = c(0.1, 0.5, 1, 2, 4, 8)) +
  # labs(y = expression(paste("Chl ", italic("a"), " mg ", m^{-3})), colour = "cluster") +
  labs(y = TeX(paste0("Chl \\textit{a}", " [", yunit, "]")), colour = "cluster",
       title = paste0("normalised: ", normalise, " trans: ", chl_trans)) +
  scale_color_manual(aesthetics =c("colour", "fill"), values = c("warm" = "red", "cold" = "blue","front" = "cornsilk4"), 
                     breaks = c("cold", "front", "warm")) +
 # guides(fill = "none") +
  coord_cartesian(ylim=c(0, 4), xlim = as.Date(c("2018-03-01", "2018-10-15")), expand = F) +
  coord_trans(y = "log10", ylim=c(0.1, NA), xlim = as.Date(c("2018-03-01", "2018-10-15")), expand = F) +
  annotation_logticks(scaled = F) +
  facet_grid(cluster~.) +
  theme_bw() +
  theme(legend.position = "bottom")
p + labs(title = NULL)
ggsave(paste0("output/plots/Chl_norm", normalise, "_", chl_trans, "_timeseries_variability_3k.pdf"), 
       width = 15, height = 14, units = "cm")

# plot with 1 mg m-3 line as vidual help for bloom threshold, as requested by R1
p_rev1 <- ggplot(data = results, aes(x = `date`, y= `median`, colour = `cluster`)) + 
  geom_abline(intercept = 1, slope = 0, linetype = "solid", color = "gold", linewidth = 0.3) +
  geom_line() +
  geom_line(aes(y = `mean`), linetype = "dotted") +
  geom_ribbon(aes(ymin = `IQR_l`, 
                  ymax = `IQR_u`, 
                  fill = `cluster`), alpha = .3, linetype = 0) +
  geom_ribbon(aes(ymin = `CI_l`, 
                  ymax = `CI_u`, 
                  fill = `cluster`), alpha = .1, linetype = 0) +
  scale_x_date(name = "",
               date_breaks = "1 month",
               date_labels = "%b") +
  scale_y_continuous(breaks = c(0.1, 0.5, 1, 2, 4, 8)) +
  # labs(y = expression(paste("Chl ", italic("a"), " mg ", m^{-3})), colour = "cluster") +
  labs(y = TeX(paste0("Chl \\textit{a}", " [", yunit, "]")), colour = "cluster",
       title = paste0("normalised: ", normalise, " trans: ", chl_trans)) +
  scale_color_manual(aesthetics =c("colour", "fill"), values = c("warm" = "red", "cold" = "blue","front" = "cornsilk4"), 
                     breaks = c("cold", "front", "warm")) +
  # guides(fill = "none") +
  coord_cartesian(ylim=c(0, 4), xlim = as.Date(c("2018-03-01", "2018-10-15")), expand = F) +
  coord_trans(y = "log10", ylim=c(0.1, NA), xlim = as.Date(c("2018-03-01", "2018-10-15")), expand = F) +
  annotation_logticks(scaled = F) +
  facet_grid(cluster~.) +
  theme_bw() +
  theme(legend.position = "bottom")
p_rev1 + labs(title = NULL)
ggsave(paste0("output/plots/Chl_norm", normalise, "_", chl_trans, "_timeseries_variability_3k_rev1.pdf"), 
       width = 15, height = 14, units = "cm")


# save for Verteidigung (Also adjust alpha of CIs... before plotting, and change back to 0.3 and 0.1)
p + 
  annotation_logticks(scaled = F, color = "white") +
  scale_color_manual(aesthetics =c("colour", "fill"), 
                     values = c("warm" = "darkorange3", "cold" = "cyan","front" = "grey"), 
                     breaks = c("cold", "front", "warm")) +
  theme(text = element_text(color = "white"),
        panel.background = element_rect(fill = NA, color = "white"),
        panel.grid = element_line(color = "grey30"),
        panel.border = element_rect(color = "white"),
        plot.background = element_rect(fill = "black", color = "black"),
        axis.text = element_text(color = "white"),
        axis.ticks = element_line(color = "white"),
        legend.background = element_rect(fill = "black"),
        strip.background = element_rect(fill = NA, color = "white"),
        strip.text = element_text(color = "white")
        )
ggsave(paste0("~/Documents/Verteidigung/plots/figures_paper3/Chl_norm", normalise, "_", chl_trans, "_timeseries_variability_3k_verteidigung.pdf"), 
       width = 15, height = 14, units = "cm")

# plot means together
ggplot(data = results, aes(x = `date`, y= `mean`, colour = `cluster`)) + 
  geom_line() +
  scale_x_date(name = "",
               date_breaks = "1 month",
               date_labels = "%b") +
  scale_y_continuous(breaks = c(0.1, 0.5, 1, 2, 4, 8)) +
  # labs(y = expression(paste("Chl ", italic("a"), " mg ", m^{-3})), colour = "cluster") +
  labs(y = TeX(paste0("Chl \\textit{a}", " [", yunit, "]")), colour = "cluster",
       title = paste0("normalised: ", normalise, " trans: ", chl_trans)) +
  scale_color_manual(aesthetics =c("colour", "fill"), values = c("warm" = "red", "cold" = "blue","front" = "cornsilk4"), 
                     breaks = c("cold", "front", "warm")) +
  # guides(fill = "none") +
#  coord_cartesian(ylim=c(0, 4), xlim = as.Date(c("2018-03-01", "2018-10-15")), expand = F) +
  coord_trans(y = "log10", ylim=c(0.1, NA), xlim = as.Date(c("2018-03-01", "2018-10-15")), expand = F) +
  annotation_logticks(scaled = F) +
  theme_bw() +
  theme(legend.position = "bottom")


## plot IQR over time
ggplot(data = results, aes(x = `date`, y = `IQR_u` - `IQR_l`, colour = `cluster`)) +
  geom_line() +
  labs(x = "Date", y = TeX(paste0("IQR", " [", yunit, "]")), title = paste0("normalised: ", normalise, " trans: ", chl_trans)) +
  #facet_grid(cluster~.) +
  scale_x_date(name = "",
               date_breaks = "1 month",
               date_labels = "%b") +
  scale_color_manual(aesthetics =c("colour", "fill"), values = c("warm" = "red", "cold" = "blue","front" = "cornsilk4"), 
                                            breaks = c("cold", "front", "warm")) +
  coord_cartesian(ylim=c(0, 7), xlim = as.Date(c("2018-03-01", "2018-10-15")), expand = F) +
  theme_bw()
ggsave(paste0("output/plots/Chl_norm", normalise, "_", chl_trans, "_timeseries_IQR_3k.pdf"))


## get post-bloom chl-a concentrations 
tst <- results %>%
  filter((date > as.Date("2018-06-25")) & (date < as.Date("2018-09-15"))) %>%
  group_by(cluster) %>%
  summarise(mean = mean(mean, na.rm = T),
            median = mean(median, na.rm = T), 
            CI_l = mean(CI_l, na.rm = T), 
            n = n())


#...............................................................................
# create common figure of CIR, CIR of log10 Chl, and CIR of log10, normalised by concentration

## 1) plot CI range over time

# normalised data?
normalise <- FALSE # one of c(TRUE or FALSE)
# choose transformation of chl data
chl_trans = "no_trans" # one of c("log_trans", "no_trans")
yunit <- ifelse(normalise, "no dim", "mg m$^{-3}$")

# but log transformation requires the value to be unitless. therefore we assume:  
# If 𝑥=0.5 is measured in some units, say, seconds, then taking the log actually means ln(0.5𝑠/1𝑠)=ln(0.5)
if (chl_trans == "log_trans"){yunit = "no dim"}
results <- read_csv2(paste0("output/Chl_norm", normalise, "_", chl_trans, "_timeseries_variability_3k.csv"))


p1 <- ggplot(data = results, aes(x = `date`, y = `CI_u` - `CI_l`, colour = `cluster`)) +
  geom_line() +
  labs(x = "Date", y = TeX(paste0("$\\Delta$CI95 [", yunit, "]")), title = TeX("Width of CI95 for untransformed Chl-\\textit{a}")) + # $Q_{0.975} - Q_{0.025
  #facet_grid(cluster~.) +
  scale_x_date(name = "",
               date_breaks = "1 month",
               date_labels = "%b") +
  scale_color_manual(aesthetics =c("colour", "fill"), values = c("warm" = "red", "cold" = "blue","front" = "cornsilk4"), 
                     breaks = c("cold", "front", "warm")) +
  coord_cartesian(ylim=c(0, NA), xlim = as.Date(c("2018-03-01", "2018-10-15")), expand = F) +
  theme_bw() +
  theme(legend.position = "bottom", 
        title = element_text(size = 8))

p1
ggsave(paste0("output/plots/Chl_norm", normalise, "_", chl_trans, "_timeseries_CI95_3k.pdf"),
       height = 8, width = 15, units = "cm")

## 2) plot CI range over time for log10 transformed data

# normalised data?
normalise <- FALSE # one of c(TRUE or FALSE)
# choose transformation of chl data
chl_trans = "log_trans" # one of c("log_trans", "no_trans")
yunit <- ifelse(normalise, "no dim", "mg m$^{-3}$")

# but log transformation requires the value to be unitless. therefore we assume:  
# If 𝑥=0.5 is measured in some units, say, seconds, then taking the log actually means ln(0.5𝑠/1𝑠)=ln(0.5)
if (chl_trans == "log_trans"){yunit = "no dim"}
results <- read_csv2(paste0("output/Chl_norm", normalise, "_", chl_trans, "_timeseries_variability_3k.csv"))


p2 <- ggplot(data = results, aes(x = `date`, y = `CI_u` - `CI_l`, colour = `cluster`)) +
  geom_line() +
  labs(x = "Date", y = TeX(paste0("$\\Delta$CI95 [", yunit, "]")), title = TeX("Width of CI95 for log$_{10}$(Chl-\\textit{a}/0.1 mg m$^{-3}$)")) +
  #facet_grid(cluster~.) +
  scale_x_date(name = "",
               date_breaks = "1 month",
               date_labels = "%b") +
  scale_color_manual(aesthetics =c("colour", "fill"), values = c("warm" = "red", "cold" = "blue","front" = "cornsilk4"), 
                     breaks = c("cold", "front", "warm")) +
  coord_cartesian(ylim=c(0, NA), xlim = as.Date(c("2018-03-01", "2018-10-15")), expand = F) +
  theme_bw() +
  theme(legend.position = "bottom", 
        title = element_text(size = 8))

p2
ggsave(paste0("output/plots/Chl_norm", normalise, "_", chl_trans, "_timeseries_CI95_3k.pdf"),
       height = 8, width = 15, units = "cm")

## 3) plot CI range over time for log10 transformed data, divided by log10(conc)

# normalised data?
normalise <- FALSE # one of c(TRUE or FALSE)
# choose transformation of chl data
chl_trans = "log_trans" # one of c("log_trans", "no_trans")
yunit <- ifelse(normalise, "no dim", "mg m$^{-3}$")

# but log transformation requires the value to be unitless. therefore we assume:  
# If 𝑥=0.5 is measured in some units, say, seconds, then taking the log actually means ln(0.5𝑠/1𝑠)=ln(0.5)
if (chl_trans == "log_trans"){yunit = "no dim"}
results <- read_csv2(paste0("output/Chl_norm", normalise, "_", chl_trans, "_timeseries_variability_3k.csv"))


p3 <- 
  ggplot(data = results, aes(x = `date`, y = (`CI_u` - `CI_l`)/`mean`, colour = `cluster`)) +
  geom_line() +
  labs(x = "Date", y = TeX("normalised $\\Delta$CI95"), title = TeX("Normalised width of CI95 for log$_{10}$(Chl-\\textit{a}/0.1 mg m$^{-3}$)")) +
  #facet_grid(cluster~.) +
  scale_x_date(name = "",
               date_breaks = "1 month",
               date_labels = "%b") +
  scale_color_manual(aesthetics =c("colour", "fill"), values = c("warm" = "red", "cold" = "blue","front" = "cornsilk4"), 
                     breaks = c("cold", "front", "warm")) +
  coord_cartesian(ylim=c(NA, NA), xlim = as.Date(c("2018-03-01", "2018-10-15")), expand = F) +
  theme_bw() +
  theme(legend.position = "bottom", 
        title = element_text(size = 8))

p3

ggpubr::ggarrange(p1, p2, p3, ncol = 1, common.legend = T, legend = "bottom", 
                  align = "hv", labels = "auto")
ggsave(paste0("output/plots/Chl_3in1_timeseries_CI95_3k.pdf"),
        height = 15, width = 10, units = "cm", device = cairo_pdf)


# save for Verteidigung
  ggpubr::ggarrange(p1 + 
                      scale_color_manual(aesthetics =c("colour", "fill"), 
                                         values = c("warm" = "darkorange3", "cold" = "cyan","front" = "grey"), 
                                         breaks = c("cold", "front", "warm")) +
                      theme(text = element_text(color = "white"),
                            panel.background = element_rect(fill = NA, color = "white"),
                            panel.grid = element_line(color = "grey30"),
                            panel.border = element_rect(color = "white"),
                            plot.background = element_rect(fill = "black", color = "black"),
                            axis.text = element_text(color = "white"),
                            axis.ticks = element_line(color = "white"),
                            legend.background = element_rect(fill = "black"),
                            strip.background = element_rect(fill = NA, color = "white"),
                            strip.text = element_text(color = "white")
                      ),
                    p2+ 
                      scale_color_manual(aesthetics =c("colour", "fill"), 
                                         values = c("warm" = "darkorange3", "cold" = "cyan","front" = "grey"), 
                                         breaks = c("cold", "front", "warm")) +
                      theme(text = element_text(color = "white"),
                            panel.background = element_rect(fill = NA, color = "white"),
                            panel.grid = element_line(color = "grey30"),
                            panel.border = element_rect(color = "white"),
                            plot.background = element_rect(fill = "black", color = "black"),
                            axis.text = element_text(color = "white"),
                            axis.ticks = element_line(color = "white"),
                            legend.background = element_rect(fill = "black"),
                            strip.background = element_rect(fill = NA, color = "white"),
                            strip.text = element_text(color = "white")
                      ), 
                    p3+ 
                      scale_color_manual(aesthetics =c("colour", "fill"), 
                                         values = c("warm" = "darkorange3", "cold" = "cyan","front" = "grey"), 
                                         breaks = c("cold", "front", "warm")) +
                      theme(text = element_text(color = "white"),
                            panel.background = element_rect(fill = NA, color = "white"),
                            panel.grid = element_line(color = "grey30"),
                            panel.border = element_rect(color = "white"),
                            plot.background = element_rect(fill = "black", color = "black"),
                            axis.text = element_text(color = "white"),
                            axis.ticks = element_line(color = "white"),
                            legend.background = element_rect(fill = "black"),
                            strip.background = element_rect(fill = NA, color = "white"),
                            strip.text = element_text(color = "white")
                      ), 
                    ncol = 1, common.legend = T, legend = "bottom", 
                    align = "hv", labels = "auto")
  ggsave(paste0("~/Documents/Verteidigung/plots/figures_paper3/Chl_3in1_timeseries_CI95_3k.pdf"),
         height = 15, width = 10, units = "cm", device = cairo_pdf, bg = "black")


#...............................................................................

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
  geom_line(aes(y = `n.obs`), linetype = "dotted") + 
  scale_x_date(name = "",
               date_breaks = "1 month",
               date_labels = "%b") +
  facet_grid(cluster~.) +
  theme_bw()


# using the normalised log trans results:
# plot the number of observations per cluster
ggplot(data = results, aes(x = `date`, y = n.obs, colour = `cluster`)) +
  geom_line() +
  labs(x = "Date", y = TeX("Number of Chl-$\\textit{a}$ observations [pixels]")) +
  facet_grid(cluster~.) +
  scale_x_date(name = "",
               date_breaks = "1 month",
               date_labels = "%b") +
  scale_color_manual(aesthetics =c("colour", "fill"), values = c("warm" = "red", "cold" = "blue","front" = "cornsilk4"), 
                     breaks = c("cold", "front", "warm")) +
  coord_cartesian(ylim=c(NA, NA), xlim = as.Date(c("2018-03-01", "2018-10-15")), expand = F) +
  theme_bw() +
  theme(legend.position = "bottom", 
        title = element_text(size = 8))
ggsave(paste0("output/plots/appx_Chl_obs_2018_rev1.pdf"),
       width = 15, height = 14, units = "cm")
# get average n.obs per cluster from March to October
tst <- results %>%
  filter((date > as.Date("2018-03-01")) & (date < as.Date("2018-10-01"))) %>%
  group_by(cluster) %>%
  summarise(mean.nobs = mean(n.obs, na.rm = T), # there are no nans in n.obs
            mean.nobs2 = mean(n.obs), 
            n = n())
tst





## plot coefficient of variation:
# normalised, log-transformed IQR/median or IQR/mean  -> multiply this by 100, get % deviation 

normalise = TRUE
chl_trans = "log_trans"

results <- read_csv2(paste0("output/Chl_norm", normalise, "_", chl_trans, "_timeseries_variability_3k.csv"))


ggplot(data = results, aes(x = `date`, y = ((`IQR_u` - `IQR_l`)/`mean`)*100, colour = `cluster`)) +
  geom_line() +
  labs(x = "Date", y = "CV [%]", title = "(IQR/mean)*100") +
  coord_cartesian(ylim=c(0, 220), xlim = as.Date(c("2018-03-01", "2018-10-15")), expand = F) +
  #facet_grid(cluster~.) +
  scale_x_date(name = "",
               date_breaks = "1 month",
               date_labels = "%b") +
  scale_color_manual(aesthetics =c("colour", "fill"), values = c("warm" = "red", "cold" = "blue","front" = "grey"), 
                     breaks = c("cold", "front", "warm")) +
  theme_bw()
ggsave(paste0("output/plots/Chl_norm", normalise, "_", chl_trans, "_timeseries_CV_3k.pdf"))



## ___________________________________________________________________________
## plot scatterplot of CI-Range vs chl-conc
normalise = FALSE
chl_trans = "no_trans"

results <- read_csv2(paste0("output/Chl_norm", normalise, "_", chl_trans, "_timeseries_variability_3k.csv"))
p1 <- 
  ggplot(data = results, aes(x = mean, y = (`CI_u` - `CI_l`), color = cluster)) +
  geom_point() +
#  geom_smooth(method = "lm")+ # show linear regression lines
  scale_color_manual(aesthetics =c("colour", "fill"), values = c("warm" = "red", "cold" = "blue","front" = "cornsilk4"), 
                     breaks = c("cold", "front", "warm")) +
  labs(x = TeX("Mean Chl-\\textit{a} concentration [mg m$^{-3}$]"),
       y = TeX("$\\Delta$CI95 [mg m$^{-3}$]")) +
  theme_bw() +
  theme(legend.position = "bottom",
        text = element_text(family = "sans"))
p1
ggsave(p1, paste0("output/plots/Chl_norm", normalise, "_", chl_trans, "_Var_rel_Conc_3k.pdf"),
       height = 8, width = 8.5, units = "cm")

normalise = FALSE
chl_trans = "log_trans"

results <- read_csv2(paste0("output/Chl_norm", normalise, "_", chl_trans, "_timeseries_variability_3k.csv"))
p2 <- ggplot(data = results, aes(x = mean, y = (`CI_u` - `CI_l`), color = cluster)) +
  geom_point() +
#  geom_smooth(method = "lm")+
  scale_color_manual(aesthetics =c("colour", "fill"), values = c("warm" = "red", "cold" = "blue","front" = "cornsilk4"), 
                     breaks = c("cold", "front", "warm")) +
  labs(#x = TeX("Mean log$_{10} \\left(\\frac{Chl \\textit{a}}{0.1~mg~m^{-3}}\\right)$"),
       x = expression('Mean log'[10]*bgroup('(', frac(Chl*~italic(a), 0.1*~mg~m^{phantom() - 3}) * phantom(.), ')')),
       y = TeX("$\\Delta$CI95")) +
  theme_bw() +
  theme(legend.position = "bottom",
        text = element_text(family = "sans"))
p2
ggsave(p2, paste0("output/plots/Chl_norm", normalise, "_", chl_trans, "_Var_rel_Conc_3k.pdf"),
       height = 9, width = 8.5, units = "cm")

p_a <- ggpubr::ggarrange(p1, p2, common.legend = T, legend = "bottom", align = "hv", labels = "auto")
# ggsave(paste0("output/plots/Chl_norm", normalise, "_no_AND_log_trans_Var_rel_Conc_3k.pdf"),
#        height = 9, width = 15, units = "cm",  device = "pdf") # Delta will not show up in pdf, use cairo_pdf (but makes parentheses ugly)

ggsave(paste0("output/plots/Chl_norm", normalise, "_no_AND_log_trans_Var_rel_Conc_3k.pdf"),
       height = 9, width = 15, units = "cm",  device = cairo_pdf)

## ___________________________________________________________________________
## plot scatterplot of CI-Range vs chl-conc; with linear regression lines as suggested by Reviewer (R2, 19)


normalise = FALSE
chl_trans = "no_trans"

results <- read_csv2(paste0("output/Chl_norm", normalise, "_", chl_trans, "_timeseries_variability_3k.csv"))
lmod_cold <- lm((`CI_u` - `CI_l`) ~ mean - 1, data = results, subset = cluster=="cold")
lmod_front <- lm((`CI_u` - `CI_l`) ~ mean - 1, data = results, subset = cluster=="front")
lmod_warm <- lm((`CI_u` - `CI_l`) ~ mean - 1, data = results, subset = cluster=="warm")
# summary(lm)
# predict(lmod_cold, data.frame(mean=c(1:5)))
confint(lmod_cold, 'mean', level = 0.95)[[1]] # get 95% CI for the estimate of "mean" 


pred_lm <- data.frame(mean=c(0:5)) %>%
  mutate(predCI_cold = mean * lmod_cold$coefficients[['mean']],
         predCI_front = mean * lmod_front$coefficients[['mean']],
         predCI_warm = mean * lmod_warm$coefficients[['mean']]) %>%
  pivot_longer(cols = -mean, names_prefix = "predCI_") %>%
  mutate(lower_val = case_when(name == "cold" ~ mean * confint(lmod_cold, 'mean', level = 0.95)[[1]],
                               name == "warm" ~ mean * confint(lmod_warm, 'mean', level = 0.95)[[1]],
                               name == "front" ~ mean * confint(lmod_front, 'mean', level = 0.95)[[1]]),
         upper_val = case_when(name == "cold" ~ mean * confint(lmod_cold, 'mean', level = 0.95)[[2]],
                               name == "warm" ~ mean * confint(lmod_warm, 'mean', level = 0.95)[[2]],
                               name == "front" ~ mean * confint(lmod_front, 'mean', level = 0.95)[[2]])
         )

 p1 <- 
  ggplot() +
   geom_ribbon(data=pred_lm, aes(x = mean, ymin = lower_val, ymax = upper_val, 
                                 fill = name), alpha = 0.3) +
   geom_line(data=pred_lm, aes(x=mean, y=value, color = name)) +
  geom_point(data = results, aes(x = mean, y = (`CI_u` - `CI_l`), color = cluster)) +
  #  geom_smooth(method = "lm")+ # show linear regression lines
  scale_color_manual(aesthetics =c("colour", "fill"), values = c("warm" = "red", "cold" = "blue","front" = "cornsilk4"), 
                     breaks = c("cold", "front", "warm")) +
  labs(x = TeX("Mean Chl \\textit{a} concentration [mg m$^{-3}$]"),
       y = TeX("$\\Delta$CI95 [mg m$^{-3}$]")) +
  coord_cartesian(xlim=c(0,4.8), expand = F)+
  theme_bw() +
  theme(legend.position = "bottom",
        text = element_text(family = "sans"))
p1


normalise = FALSE
chl_trans = "log_trans"

results <- read_csv2(paste0("output/Chl_norm", normalise, "_", chl_trans, "_timeseries_variability_3k.csv"))
lmod_cold <- lm((`CI_u` - `CI_l`) ~ mean - 1, data = results, subset = cluster=="cold")
lmod_front <- lm((`CI_u` - `CI_l`) ~ mean - 1, data = results, subset = cluster=="front")
lmod_warm <- lm((`CI_u` - `CI_l`) ~ mean - 1, data = results, subset = cluster=="warm")

# summary(lmod_cold)
pred_lm <- data.frame(mean=seq(0,1.55, by=0.1)) %>%
  mutate(predCI_cold = mean * lmod_cold$coefficients[['mean']],
         predCI_front = mean * lmod_front$coefficients[['mean']],
         predCI_warm = mean * lmod_warm$coefficients[['mean']]) %>%
  pivot_longer(cols = -mean, names_prefix = "predCI_") %>%
  mutate(lower_val = case_when(name == "cold" ~ mean * confint(lmod_cold, 'mean', level = 0.95)[[1]],
                               name == "warm" ~ mean * confint(lmod_warm, 'mean', level = 0.95)[[1]],
                               name == "front" ~ mean * confint(lmod_front, 'mean', level = 0.95)[[1]]),
         upper_val = case_when(name == "cold" ~ mean * confint(lmod_cold, 'mean', level = 0.95)[[2]],
                               name == "warm" ~ mean * confint(lmod_warm, 'mean', level = 0.95)[[2]],
                               name == "front" ~ mean * confint(lmod_front, 'mean', level = 0.95)[[2]]))

p2 <- ggplot() +
  geom_ribbon(data=pred_lm, aes(x = mean, ymin = lower_val, ymax = upper_val, 
                                fill = name), alpha = 0.3) +
  geom_line(data=pred_lm, aes(x=mean, y=value, color = name)) +
  geom_point(data = results, aes(x = mean, y = (`CI_u` - `CI_l`), color = cluster)) +
  #  geom_smooth(method = "lm")+
  scale_color_manual(aesthetics =c("colour", "fill"), values = c("warm" = "red", "cold" = "blue","front" = "cornsilk4"), 
                     breaks = c("cold", "front", "warm")) +
  labs(#x = TeX("Mean log$_{10} \\left(\\frac{Chl \\textit{a}}{0.1~mg~m^{-3}}\\right)$"),
    x = expression('Mean log'[10]*bgroup('(', frac(Chl*~italic(a), 0.1*~mg~m^{phantom() - 3}) * phantom(.), ')')),
    y = TeX("$\\Delta$CI95")) +
  coord_cartesian(xlim=c(0,1.55), expand = F)+
  theme_bw() +
  theme(legend.position = "bottom",
        text = element_text(family = "sans"))
p2

p_a <- ggpubr::ggarrange(p1, p2, common.legend = T, legend = "bottom", align = "hv", labels = "auto")
p_a
ggsave(paste0("output/plots/Chl_norm", normalise, "_no_AND_log_trans_Var_rel_Conc_3k_update.pdf"),
       height = 9, width = 15, units = "cm",  device = pdf)



## ___________________________________________________________________________
# try: shape for cluters, color for date
ggplot(data = results, aes(x = mean, y = (`CI_u` - `CI_l`), shape = cluster, color = date)) +
  geom_point() +
 # geom_abline(intercept = 1, slope = 1) +
  scale_color_date(low = "yellow", high = "blue") +
  labs(x = TeX("Mean log$_{10}$(Chl \\textit{a} concentration / 0.1 mg m$^{-3}$)"),
       y = TeX("$Q_{0.975} - Q_{0.025}$ (of log$_{10}$(Chl \\textit{a} concentration / 0.1 mg m$^{-3}$))")) +
  facet_wrap(~cluster)+
  theme_bw()



## create table to compare variability in total FS and clustered subregions 
normalise = FALSE
chl_trans = "no_trans"

results_c <- read_csv2(paste0("output/Chl_norm", normalise, "_", chl_trans, "_timeseries_variability_3k.csv"))
results_uc <- read_csv2(paste0("output/Chl_norm", normalise, "_", chl_trans, "_timeseries_variability_0k.csv"))

tab <- bind_rows(results_c, results_uc) %>%
  mutate(cluster = factor(cluster, levels = c("total", "cold", "front", "warm")),
         CIR = CI_u - CI_l,
         IQR = IQR_u - IQR_l,
         CV = (CIR / `mean`)*100) %>%
  group_by(cluster) %>%
  summarise(av_sd = mean(SD, na.rm = T),
            av_IQR = mean(IQR, na.rm = T),
            av_CIR = mean(CIR, na.rm = T),
            av_CV = mean(CV, na.rm = T)) %>%
  mutate_at(2:5, round, 2) 
  
  

