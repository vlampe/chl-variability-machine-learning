## plot more complex Chl time series
## without cluster distinction
## because each data point must belong to one of the clusters, their union covers the whole region

## author: Vanessa Lampe
## date: 2023-10-23
## version: 0
## git: 

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
  
    
  # extract values from list 
  chl_vals <- c(chl_list[[var]][[regimes[1]]][[d]], 
                chl_list[[var]][[regimes[2]]][[d]],
                chl_list[[var]][[regimes[3]]][[d]])
  
  # if chl_vals is empty, skip and attach row to results table
  if (length(chl_vals) == 0){
    out <- tibble(`date` = datum,
                  `cluster` = "total",
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
    
    # fix failing KDEs: set max based on maximum value, but keep stepsize at 0.01. since this also not always works, also try diffKDE with different parameters. 
    if(any(is.nan(pdf_y))){ # try again with less steps, if KDE failed
      
      try(pdf_list <- py_to_r(diffKDE$KDE(np$array(chl_vals), kde_min, kde_max)))
      try(pdf_list <-  py_to_r(diffKDE$KDE(np$array(chl_vals))))
      try(pdf_list <-  py_to_r(diffKDE$KDE(np$array(chl_vals), xmin = min_val-(0.5*min_val)))) 
      try(pdf_list <-  py_to_r(diffKDE$KDE(np$array(chl_vals), kde_min, round(max_val+0.5*max_val,2))))
      
      
      pdf_x <- pdf_list[[2]] # convert numpy array to R array
      pdf_y <- pdf_list[[1]]
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
                  `cluster` = "total",
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
# save output
write_csv2(x = results %>% arrange(`date`) , file = paste0("output/Chl_norm", normalise, "_", chl_trans, "_timeseries_variability_0k.csv"))

# ______________________________________________________________________________
## PLOTTING

# normalised data?
normalise <- FALSE # one of c(TRUE or FALSE)

# choose transformation of chl data
chl_trans = "no_trans" # one of c("log_trans", "no_trans")
#...............................................................................

yunit <- ifelse(normalise, "no dim", "mg m$^{-3}$")



results <- read_csv2(paste0("output/Chl_norm", normalise, "_", chl_trans, "_timeseries_variability_0k.csv"))

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
               date_breaks = "1 month",
               date_labels = "%b") +
  #scale_y_continuous(breaks = c(0.1, 1, 2, 4, 8, 12, 16)) +
  # labs(y = expression(paste("Chl ", italic("a"), " mg ", m^{-3})), colour = "cluster") +
  labs(y = TeX(paste0("Chl \\textit{a}", " [", yunit, "]")), colour = "cluster",
#        title = paste0("normalised: ", normalise, " trans: ", chl_trans)) +
       title = NULL) +
  scale_color_manual(aesthetics =c("colour", "fill"), values = c("total" = "black", "warm" = "red", "cold" = "blue","front" = "grey"), 
                     breaks = c("total", "cold", "front", "warm")) +
   guides(fill = "none", color = "none") +
#  coord_cartesian(ylim=c(0, 6), xlim = as.Date(c("2018-03-01", "2018-10-15")), expand = F) +
  coord_trans(y = "log10", ylim=c(0.1, NA), xlim = as.Date(c("2018-03-01", "2018-10-15")), expand = F) +
  annotation_logticks(scaled = F) +
  scale_y_continuous(breaks = c(0.1, 0.5, 1, 5, 10, 15)) +
#  facet_grid(cluster~.) +
  theme_bw()
ggsave(paste0("output/plots/Chl_norm", normalise, "_", chl_trans, "_timeseries_variability_0k.pdf"),
              width = 15, height = 8, units = "cm")

