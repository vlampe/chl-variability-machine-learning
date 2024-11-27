# use SST clusters to analyse distribution and time series of phytoplankton
# functional types

## author: Vanessa Lampe
## date: 2021-03-12/ 2021-05-07
## version: 1 (updated pft data sets (May 2021 update))
## git: https://git.geomar.de/vanessa-lampe/chlorophyll-time-series

library(raster)
library(sf)
library(tidyverse) 
library(writexl)
library(readxl)
library(scales)
library(ggpubr)
library("rnaturalearth")
library("rnaturalearthdata")

## data sources
# https://resources.marine.copernicus.eu/?option=com_csw&view=details&product_id=OCEANCOLOUR_GLO_CHL_L3_REP_OBSERVATIONS_009_085
# data path
datapath <- "updated_pft_data/"
#read warm/cold cluter output
warmcold_out <- read_excel("Output_warm_and_cold_cluster_analysis.xlsx")%>%
  mutate(year = lubridate::year(date))


# first find out which data are available
filenames <- list.files(path = "pft_data/", pattern = "*.nc")
fnamelist <- str_split(filenames, pattern = "_")
filepfts <- sapply(fnamelist, "[", 4)
fileyears <- sapply(fnamelist, "[", 5)
fileyears <- str_remove(fileyears, pattern = ".nc")
 
files <- tibble(filename = filenames, pft = filepfts, year = fileyears)

# set PFTs and years to analyse
years <- unique(files$year)
PFTs <- unique(files$pft)

# load polys
polys <- readRDS("polygons_2009-2018.rds")

# prepare df for results
results <- tibble()

# should daily plots (maps) be produced?
producePlot = T
if(producePlot) {
  world <- ne_countries(scale = "medium", returnclass = "sf")
  stereproj <- "+proj=stere +lat_0=90 +lat_ts=75 +lon_0=0 +k=1 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0"
}

# create folders for daily maps
for(pft in PFTs){
  if(!dir.exists(paste0("figures/daily_maps/updated_", pft))&producePlot) {
    dir.create(paste0("figures/daily_maps/updated_", pft))
  }
}


# enter yearly loop

for(y in years){
  
  for(pft in PFTs) {
    
    dat <- brick(paste0(datapath, "dataset-oc-glo-bio-multi-l3-pft_4km_daily-rep_",pft,"_", y, ".nc"),
                                  values = T)
    nc.dates <- getZ(dat)
    
    #get min and max values for the color scale in the plots
    # minValue(dat)
    # test <- dat
    # test <- setMinMax(test)
    # min(minValue(test), na.rm = T)
    # max(maxValue(test), na.rm = T)
  # enter loop for each day of year in nc file
    for(d in 1:length(nc.dates)){
      
      # get daily data and crop to the study area 
      pft_d <- dat[[d]]
      pft_d <- crop(pft_d, extent(-10, 12, 77, 81))
      today <- getZ(pft_d)
      
      # get jd from warmcold_out
      jd <- warmcold_out$`julian day`[as.Date(warmcold_out$date) == today]
      
      # get polys
      wc_poly <- polys[[y]][[as.character(jd)]]
      warm_poly <- wc_poly$warm
      cold_poly <- wc_poly$cold
      
      # 1. get mean, min, max and n.obs for total area
      chl <- getValues(pft_d)
      chl_mean_total <- mean(chl, na.rm = T)
      chl_max_total <- max(chl, na.rm = T)
      chl_min_total <- min(chl, na.rm = T)
      nobs_total <- length(na.omit(chl))
      
      # 2. get mean, min, max and n.obs for warm area
      pft_d_warm <- mask(pft_d, warm_poly)
      chl <- getValues(pft_d_warm)
      chl_mean_warm <- mean(chl, na.rm = T)
      chl_max_warm <- max(chl, na.rm = T)
      chl_min_warm <- min(chl, na.rm = T)
      nobs_warm <- length(na.omit(chl))
      
      # 3. get mean, min, max and n.obs for cold area
      pft_d_cold <- mask(pft_d, cold_poly)
      chl <- getValues(pft_d_cold)
      chl_mean_cold <- mean(chl, na.rm = T)
      chl_max_cold <- max(chl, na.rm = T)
      chl_min_cold <- min(chl, na.rm = T)
      nobs_cold <- length(na.omit(chl))
      
      # put output to out
      out <- tibble("Julian day" = jd,
                    date = today, 
                    PFT = pft, 
                    chl_mean_total = chl_mean_total,
                    chl_max_total = chl_max_total,
                    chl_min_total = chl_min_total,
                    nobs_total = nobs_total, 
                    chl_mean_cold = chl_mean_cold, 
                    chl_max_cold = chl_max_cold, 
                    chl_min_cold = chl_min_cold, 
                    nobs_cold = nobs_cold,
                    chl_mean_warm = chl_mean_warm, 
                    chl_max_warm = chl_max_warm, 
                    chl_min_warm = chl_min_warm, 
                    nobs_warm = nobs_warm
                    )
      
      results <- bind_rows(results, out)
      
      # produce a daily plot
      if(producePlot){
        
        
        # if there was no chlorophyll today, use empty dummy ## Funktoiniert noch nicht!!
        if(is.nan(chl_mean_total)) {
          chl_today = F
          nrows <- 10
          #create empty dummy sf
          dummy_sf <- st_sf(id = 1:nrows, geometry = st_sfc(lapply(1:nrows, function(x) st_geometrycollection())),
                            crs = "+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0")
          names(dummy_sf) <- c("FILL", "geometry")
          dummy_sf$FILL <- NA
          
        } else {
          chl_today <- T
          names(pft_d) <- "FILL"
          chl_sf <- pft_d %>%
            as("SpatialPixelsDataFrame") %>%
            as("SpatialPolygonsDataFrame") %>%
            as("sf")
        }
        
        warmpoly <- spTransform(warm_poly, CRSobj = stereproj)
        coldpoly <- spTransform(cold_poly, CRSobj = stereproj)
        
        pft_plot <- 
          ggplot() +
          labs(#title = paste0("Sea surface temperature and cold-warm separation on ", thisdate),
            #subtitle = "NOAA High-resolution Blended Analysis",
            y = "", 
            x = "") +
          # add chl layer
            {if(chl_today) geom_sf(data = chl_sf, aes(fill = FILL, colour = FILL), lwd = 0)}+
            {if(!chl_today) geom_sf(data = dummy_sf, aes(fill = FILL, colour = FILL), lwd = 0)}+
          # geom_sf(data = chl_sf, aes(fill = FILL, colour = FILL), lwd = 0) + 
          # geom_raster(data = sst_df, aes(x = x, y = y, fill = sst), alpha = 1) +
          
          # add colour scale for sill and colour
          scale_color_gradient2(aesthetics = c("color", "fill"),
                                #  colors = pals::ocean.algae(4),
                                low = "blue", mid = "limegreen", high = "red", midpoint = 0.1,
                                limits = c(0.1, 10),
                                breaks = c(0.1, 1, 10),
                                labels = format(c(0.1, 1, 10)),
                                trans = "log", 
                                oob = scales::squish , na.value = NA,
                                name = paste0(pft, "\nchl-a \n[mg m^-3]")) +
          
          # add land masses
          geom_sf(data = world, colour = "darkgrey", size = 0.2, fill = "grey") +
          
          # add white underlines for warm/cold polys
          geom_polygon(data = coldpoly, aes(x = long, y = lat, group = group), color = "white", fill = NA, linetype = "solid", size = 1.5) +
          geom_polygon(data = warmpoly, aes(x = long, y = lat, group = group), color = "white", fill = NA, linetype = "solid", size = 1.5) +
          
          # add warm/cold polys
          geom_polygon(data = coldpoly, aes(x = long, y = lat, group = group), color = "blue", fill = NA, linetype = "dashed", size = 1) +
          geom_polygon(data = warmpoly, aes(x = long, y = lat, group = group), color = "red", fill = NA, linetype = "dashed", size = 1) +
          
          # set y breaks
          scale_y_continuous(breaks = 75:83) +
          
          # set coord system
          coord_sf(crs = stereproj, xlim = c(-4e5, 4e5), ylim = c(-14.7e5, -9.3e5), #xlim war vorher 6e5
                   expand = F) +
          #add date label
          geom_label(data = data.frame(x = -4e5, y = -9.3e5, label = paste0(today)), 
                     aes(x = x, y = y, label = label), hjust = "bottom", vjust = "right", label.r = unit(0, "lines"))+
          #set theme
          theme_bw() +
          theme(legend.position = "right",
                panel.ontop=TRUE, 
                panel.border = element_rect(size = .5, colour = "black"),
                panel.background = element_blank(), 
                panel.grid = element_line(color = "darkgrey", size = 0.1),
                legend.direction = "vertical",
                legend.box.background = element_rect(colour = "black", size = 0.5, fill = "white")
          )
        
        # save daily plot
        jd_save <- str_pad(jd, 4, pad = "0")  
        cowplot::ggsave2(plot = pft_plot, filename = paste0("figures/daily_maps/updated_", pft,"/", jd_save, "_", today, ".png"),
                         width = 5, height = 2.5) 
      }
      
    }
  }
}

# save all data into excel sheet
write_xlsx(x = results %>% arrange(`Julian day`) , path = "updated_Output_PFT_SSTcluster_analysis.xlsx")


#
results <- read_xlsx(path = "updated_Output_PFT_SSTcluster_analysis.xlsx") 

## create time series plots
source("weighted_moving_average.R")

dat <- results %>%
  group_by(PFT)


# plot daily chl in whole (total) sampling area and SST sections:

# total
ggplot(data = dat, aes(x = as.Date(date), y = chl_mean_total, color = PFT)) +
  geom_line() +
  scale_x_date(name = "",
               date_breaks = "1 year",
               date_labels = "%Y",
               date_minor_breaks = "1 month") +
  scale_y_continuous(limits = c(0, NA), minor_breaks = F, name = "Chl a conc. [µg / L]") +
  theme_bw() +
  theme(legend.position=c(1, 1), 
        legend.justification = c(1,1), 
        legend.box.background = element_rect(fill = NA, color = "black", 
                                             size = 0.7)) ->chltotday
chltotday

# SST sections
ggplot(data = dat, aes(x = as.Date(date))) +
  geom_line(aes(y = chl_mean_warm, color = "w")) +
  geom_line(aes(y = chl_mean_cold, color = "c")) +
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
  facet_wrap(~PFT, dir = "v", scales = "free") +
  theme_bw() +
  theme(legend.position=c(1, 1), 
        legend.justification = c(1,1), 
        legend.box.background = element_rect(fill = NA, color = "black", 
                                             size = 0.7)) -> chlwcday
chlwcday

ggarrange(chltotday,  
          chlwcday, nrow = 2, common.legend = F,
          align = "v") %>%
  ggexport(filename = "figures/updated_daily_chl_sst_PFT.pdf")

# look at weekly means
weekly_means <- results %>%
  mutate(weekday = weekdays(date)) %>%
  group_by(Week = lubridate::floor_date(date, unit = "week"), PFT) %>%
  summarise(min_date  = min(date), max_date = max(date),
            `chl_mean_total` = mean(`chl_mean_total`, na.rm = T),
            `chl_mean_warm` = mean(`chl_mean_warm`, na.rm = T),
            `chl_mean_cold` = mean(`chl_mean_cold`, na.rm = T))

# total
ggplot(data = weekly_means, aes(x = as.Date(Week), y = chl_mean_total, color = PFT)) +
  geom_line() +
  scale_x_date(name = "",
               date_breaks = "1 year",
               date_labels = "%Y",
               date_minor_breaks = "1 month") +
  scale_y_continuous(limits = c(0, NA), minor_breaks = F, name = "Chl a conc. [µg / L]") +
  theme_bw() +
  theme(legend.position=c(1, 1), 
        legend.justification = c(1,1), 
        legend.box.background = element_rect(fill = NA, color = "black", 
                                             size = 0.7)) -> chltotweek
chltotweek

# SST sections
ggplot(data = weekly_means, aes(x = as.Date(Week))) +
  geom_line(aes(y = chl_mean_warm, color = "w")) +
  geom_line(aes(y = chl_mean_cold, color = "c")) +
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
  facet_wrap(~PFT, dir = "v", scales = "free") +
  theme_bw() +
  theme(legend.position=c(1, 1), 
        legend.justification = c(1,1), 
        legend.box.background = element_rect(fill = NA, color = "black", 
                                             size = 0.7)) -> chlwcweek
chlwcweek

ggarrange(chltotweek,  
          chlwcweek, nrow = 2, common.legend = F, align = "v") %>%
  ggexport(filename = "figures/updated_weekly_chl_sst_PFT.pdf")

# look at weighted moving averages
wmv <- weighted_moving_average   # reassign function bc lazyness
wmeans <- results %>%
  arrange(PFT)

wmeans$w_t_chl <- wmv(values = wmeans$chl_mean_total,
                      weights = wmeans$nobs_total, 
                      width = 15)
wmeans$w_w_chl <- wmv(values = wmeans$chl_mean_warm,
                      weights = wmeans$nobs_warm, 
                      width = 15)
wmeans$w_c_chl <- wmv(values = wmeans$chl_mean_cold,
                      weights = wmeans$nobs_cold, 
                      width = 15)

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
  facet_wrap(~PFT, scales = "free", dir = "v") +
  theme_bw() +
  theme(legend.position=c(1, 1), 
        legend.justification = c(1,1), 
        legend.box.background = element_rect(fill = NA, color = "black", size = 0.7)
  ) -> wmeans_wc
wmeans_wc

# plot weighed (by n.obs) moving average of chl in total sampling area
ggplot(wmeans, aes(x = as.Date(date))) +
  scale_x_date(name = "",
               date_breaks = "1 year",
               date_labels = "%Y",
               date_minor_breaks = "1 month") +
  scale_y_continuous(name = "Chl a conc. [µg / L]", minor_breaks = F) +
  geom_line(aes(y = w_t_chl, color = PFT)) +
  theme_bw() +
  theme(legend.position = c(1,1),
        legend.justification = c(1,1),
        legend.box.background = element_rect(fill = NA, color = "black", size = 0.7))-> wmeans_chl_total
wmeans_chl_total

ggarrange(wmeans_chl_total ,
          wmeans_wc, nrow = 2, common.legend = F, align = "v") %>%
  ggexport(filename = "figures/updated_wmv_chl_sst_PFT.pdf")


# plot wmeans_wc with watermass as panel and PFT as color, easier to understand

wmeans <- results %>%
  arrange(PFT)

wmeans$w_t_chl <- wmv(values = wmeans$chl_mean_total,
                      weights = wmeans$nobs_total, 
                      width = 15)
wmeans$w_w_chl <- wmv(values = wmeans$chl_mean_warm,
                      weights = wmeans$nobs_warm, 
                      width = 15)
wmeans$w_c_chl <- wmv(values = wmeans$chl_mean_cold,
                      weights = wmeans$nobs_cold, 
                      width = 15)


wmeans2 <- wmeans %>%
  pivot_longer(cols = c(w_t_chl, w_w_chl, w_c_chl),
               names_to = "cluster",
               values_to = "chl conc")
wmeans2$cluster <- factor(wmeans2$cluster, levels = c('w_t_chl', 'w_c_chl', 'w_w_chl'),
                             labels = c("total", "Arctic", "Atlantic")) 

# wmeans3 <-wmeans %>%
#   select(`Julian day`, date, PFT, w_t_chl, w_w_chl, w_c_chl) %>%
#   group_by(`Julian day`, date, PFT) %>%
#   gather(key = "cluster", value = "chl conc", -c(`Julian day`, date, PFT))


ggplot(filter(wmeans2, cluster != 'total'), aes(x = as.Date(date) ,y =`chl conc`, color = PFT)) +
  scale_x_date(name = "",
               date_breaks = "1 year",
               date_labels = "%Y",
               date_minor_breaks = "1 month", expand = c(0,0)) +
  scale_y_continuous(name = "Chl a conc. [µg / L]", minor_breaks = F) +
  geom_line() +
  facet_wrap(~cluster, scales = "free", dir = "v") +
  theme_bw() +
  theme(legend.position=c(1, 1), 
        legend.justification = c(1,1), 
        legend.box.background = element_rect(fill = NA, color = "black", size = 0.7)
  ) -> p
p

ts_dat <- read_excel("Output_warm_and_cold_cluster_analysis.xlsx") 

ts_dat$wmv_chl_warm <- wmv(values = ts_dat$`mean chl warm [mg / m^3]`,
                           weights = ts_dat$`n.obs. chl warm [pixels]`,
                           width = 15)
ts_dat$wmv_chl_cold <- wmv(values = ts_dat$`mean chl cold [mg / m^3]`,
                           weights = ts_dat$`n.obs chl cold [pixels]`,
                           width = 15)

ts_dat <- ts_dat %>%
  select(`julian day`, date, wmv_chl_warm, wmv_chl_cold) %>%
  pivot_longer(cols = c(wmv_chl_warm, wmv_chl_cold),
               names_to = "cluster", values_to = "chl conc") %>%
  arrange(cluster, `julian day`)

ts_dat$cluster <- factor(ts_dat$cluster, 
                         levels = c("wmv_chl_cold", 
                                    "wmv_chl_warm"),
                         labels = c("Arctic", "Atlantic"))

p +
  geom_line(data = ts_dat, aes(y = `chl conc`), color = "darkgreen", alpha = 0.4) +
  coord_cartesian(ylim = c(0, 6))

## in terminal, cd to the pft subdirectories, then
# ffmpeg -framerate 4 -pattern_type glob -i '*.png' -vf pad="width=ceil(iw/2)*2:height=ceil(ih/2)*2" -c:v libx264 -pix_fmt yuv420p output_diatoms.mp4
# ffmpeg -framerate 4 -pattern_type glob -i '*.png' -vf pad="width=ceil(iw/2)*2:height=ceil(ih/2)*2" -c:v libx264 -pix_fmt yuv420p output_haptophytes.mp4