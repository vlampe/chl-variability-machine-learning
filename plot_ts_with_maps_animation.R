# make time series plots with maps, but as an animation
# body of plotting skript is from plot_ts_with_maps.R
setwd("~/Documents/microARC/Manuscripts/Last_chapter_clustering_manuskript_skripts/skripts/")



library(tidyverse)
library(readxl)
library(writexl)
library(ggpubr)
source("weighted_moving_average.R")
wmv <- weighted_moving_average   # reassign function bc laziness

# for maps
library(raster)
library(pals)
library(sf)
library("rnaturalearth")
library("rnaturalearthdata")
library(ggrepel)
library(scales)
library(ggpattern)
library(latex2exp)

library(Cairo)

ts_dat <- read_csv2("output/Output_3k_cluster_analysis_thesis_2nd_run_no_lines.csv", na = c("", "Inf", "-Inf","NA"))

# load ice and sst data
ice.nc <- brick("data/METOFFICE-GLO-SST-L4-REP-OBS-SST_1695911795678.nc", values = T,
                varname = 'sea_ice_fraction')
sst.nc <- brick("data/METOFFICE-GLO-SST-L4-REP-OBS-SST_1695911795678.nc", values = T,
                varname = 'analysed_sst')

polys <- readRDS(file = "output/cluster_3k_2018_2nd_run_no_lines.rds")
stereproj <- "+proj=stere +lat_0=90 +lat_ts=75 +lon_0=0 +k=1 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0"
# use polygons and nc files
world <- ne_countries(scale = "medium", returnclass = "sf")

dates <- as.Date(getZ(sst.nc))


# 
# doy <- 208
# dates <- as.Date(getZ(sst.nc))
# date_highlight <- dates[doy]


# for (date_highlight in dates[182:212]){
#   sst_ts_d <- sst_ts + 
#     geom_vline(xintercept = date_highlight, linetype = "dashed") 
#   sst_ts_d
# }


# set alpha for confidence intervals: 0.3 for light, 0.7 for dark
alpha_ci <-  0.5

# prepare ts plots
sst_ts <- 
  ggplot(data = ts_dat, aes(x = date)) + 
  geom_line(aes(y = `mean temp warm [ºC]`, col = "warm")) +
  geom_ribbon(aes(ymin = `mean temp warm [ºC]` - `sd temp warm [ºC]`, 
                  ymax = `mean temp warm [ºC]` + `sd temp warm [ºC]`, 
                  fill = "warm"), alpha = alpha_ci) +
  geom_line(aes(y = `mean temp cold [ºC]`, col = "cold")) +
  geom_ribbon(aes(ymin = `mean temp cold [ºC]` - `sd temp cold [ºC]`, 
                  ymax = `mean temp cold [ºC]` + `sd temp cold [ºC]`, 
                  fill = "cold"), alpha = alpha_ci) +
  geom_line(aes(y = `mean temp front [ºC]`, col = "front")) +
  geom_ribbon(aes(ymin = `mean temp front [ºC]` - `sd temp front [ºC]`, 
                  ymax = `mean temp front [ºC]` + `sd temp front [ºC]`, 
                  fill = "front"), alpha = alpha_ci) +
  # geom_vline(xintercept = date_highlight, linetype = "dashed") +
  labs(y = "SST [ºC]", colour = "cluster") +
  scale_x_date(name = "",
               date_breaks = "1 month",
               date_labels = c("J", "F", "M", "A", "M", "J", "J", "A", "S", "O", "N", "D")) + # %b for 3-letter-abbrv.
  #  scale_y_continuous(limits = c(0, NA), minor_breaks = F, name = "Chl a conc. [µg / L]") +
  scale_color_manual(aesthetics =c("colour", "fill"), values = c("warm" = "red", "cold" = "blue","front" = "cornsilk4"), 
                     breaks = c("cold", "front", "warm")) +
  guides(fill = "none", color = "none") +
  coord_cartesian(expand = F) +
  theme_bw() +
  theme(axis.title = element_text(size = 9))

ice_ts <- ggplot(data = ts_dat, aes(x = date)) + 

  geom_line(aes(y = `ice cover cold region [%]`*100, col = "cold")) +
  geom_ribbon(aes(ymin = (`ice cover cold region [%]` - `ice cover cold region  sd [%]`)*100, 
                  ymax = (`ice cover cold region [%]` + `ice cover cold region  sd [%]`)*100, 
                  fill = "cold"), alpha = alpha_ci) +
  geom_line(aes(y = `ice cover front region [%]`*100, col = "front")) +
  geom_ribbon(aes(ymin = (`ice cover front region [%]` - `ice cover front region sd [%]`)*100, 
                  ymax = (`ice cover front region [%]` + `ice cover front region sd [%]`)*100, 
                  fill = "front"), alpha = alpha_ci) +
  geom_line(aes(y = `ice cover warm region [%]`*100, col = "warm")) +
  geom_ribbon(aes(ymin = (`ice cover warm region [%]` - `ice cover warm region sd [%]`)*100, 
                  ymax = (`ice cover warm region [%]` + `ice cover warm region sd [%]`)*100, 
                  fill = "warm"), alpha = alpha_ci) +
  geom_abline(intercept = 15, slope = 0, linetype = "dotted") + 
  # geom_vline(xintercept = date_highlight, linetype = "dashed") +
  labs(y = "SIC [%]", colour = "cluster") +
  scale_x_date(name = "",
               date_breaks = "1 month",
               date_labels = c("J", "F", "M", "A", "M", "J", "J", "A", "S", "O", "N", "D")) +
  scale_y_continuous(breaks = c(0, 15, 25, 50, 75, 100)) +
  scale_color_manual(aesthetics =c("colour", "fill"), values = c("warm" = "red", "cold" = "blue","front" = "cornsilk4"), 
                     breaks = c("cold", "front", "warm")) +
  guides(fill = "none") +
  coord_cartesian(ylim=c(-0.7, 100), expand = F) +
  theme_bw() +
  theme(axis.title = element_text(size = 9))

counter <- 0

dates[182:212]

# all of july: 182:212

# enter loop over plotting days
for (doy in c(182:244)){ # all of july
  
  
  date_highlight <- dates[doy]
  counter <- counter+1
  
  # SST ts plot
  sst_ts_d <- sst_ts + 
    geom_vline(xintercept = date_highlight, linetype = "dashed", color = "white", linewidth = 1.1) 
  
  # ice ts plot
  ice_ts_d <- ice_ts +
    geom_vline(xintercept = date_highlight, linetype = "dashed", color = "white", linewidth = 1.1) 
  
  
  
  ## maps
  
  
  warmpoly <- spTransform(polys[[doy]][["warm"]], CRSobj = stereproj)
  coldpoly <- spTransform(polys[[doy]][["cold"]], CRSobj = stereproj)
  frontpoly <- spTransform(polys[[doy]][["front"]], CRSobj = stereproj)
  
  
  # SST 
  sst_d <- sst.nc[[doy]]
  names(sst_d) <- "FILL"
  # convert raster to specialfeatures (dieser umweg ist leider notwendig)
  sst_sf <- sst_d %>%
    as("SpatialPixelsDataFrame") %>%
    as("SpatialPolygonsDataFrame") %>%
    as("sf")
  
  
  
  sst_map <- 
    ggplot() +
    # SST layer
    geom_sf(data = sst_sf, aes(fill = FILL, color = FILL)) +# , linewidth = 0, alpha = 1
    # add land masses
    geom_sf(data = world, colour = "darkgrey", size = 0.2, fill = "oldlace") +
    # # add white underlines for warm/cold polys
    # geom_polygon(data = coldpoly, aes(x = long, y = lat, group = group), color = "white", fill = NA, linetype = "solid", size = 1.5) +
    # geom_polygon(data = warmpoly, aes(x = long, y = lat, group = group), color = "white", fill = NA, linetype = "solid", size = 1.5) +
    # 
    # # add warm/cold polys with patterns
    # geom_polygon_pattern(data = coldpoly, aes(x = long, y = lat, group = group), color = "cyan", fill = NA, pattern_fill = "cyan", pattern_density = 0.1, pattern_colour = NA, linetype = "dashed", linewidth = 0.7) +
    # geom_polygon_pattern(data = warmpoly, aes(x = long, y = lat, group = group), color = "orange3", fill = NA, pattern_fill = "orange3", pattern_density = 0.1, pattern_colour = NA, linetype = "dashed", linewidth = 0.7) +
    # geom_polygon_pattern(data = frontpoly, aes(x = long, y = lat, group = group), color = "grey", fill = NA, pattern_fill = "grey", pattern_density = 0.1, pattern_colour = NA, linetype = "dashed", linewidth = 0.7) +
    # add warm/cold polys without patterns
    geom_polygon(data = coldpoly, aes(x = long, y = lat, group = group), color = "cyan", fill = NA, linewidth = .8, linetype = 'dashed') +
    geom_polygon(data = warmpoly, aes(x = long, y = lat, group = group), color = "darkorange3", fill = NA, linewidth = .8, linetype = 'dashed') +
    geom_polygon(data = frontpoly, aes(x = long, y = lat, group = group), color = "grey", fill = NA, linewidth = .8, linetype = 'dashed') +
    
    # set CRS
    coord_sf(crs = stereproj, xlim = c(-5e5, 5e5), ylim = c(-16.7e5, -5.3e5), #xlim war vorher 6e5
             expand = F) +
    # add colour scale for sill and colour
    scale_color_gradientn(aesthetics = c("color", "fill"),
                          guide = "colourbar",
                          colors = rev(rainbow(7)),
                          limits = c(-2.5, 8.5)+273.13,
                          breaks = c(-2, 0, 2, 4, 6, 8, 10)+273.13,
                          labels = format(c(-2, 0, 2, 4, 6, 8, 10)), # +273.13
                          oob = scales::squish , na.value = NA,
                          name = "SST \n[ºC]") +
    scale_x_continuous(breaks = seq(-30, 30, by = 10)) +
    labs(title = "SST", x = NULL, y = NULL) +
    guides(fill = guide_colorbar(title = "SST [ºC]", title.position = "top", title.hjust = 0.5),
           color = guide_colorbar(title = "SST [ºC]", title.position = "top", title.hjust = 0.5)) +
    theme_bw() +
    theme(legend.position = "bottom",
          legend.box.margin = margin(t = 0),
          title = element_text(size = 8),
          panel.ontop = T,
          panel.grid = element_line(color = "snow", linewidth = 0.1)
    )
  
  # sst_map
  
  
  # ICE
  ice_d <- ice.nc[[doy]]
  names(ice_d) <- "FILL"
  # convert raster to specialfeatures (dieser umweg ist leider notwendig)
  ice_sf <- ice_d %>%
    as("SpatialPixelsDataFrame") %>%
    as("SpatialPolygonsDataFrame") %>%
    as("sf")
  
  ice_map <- 
    ggplot() +
    # ICE layer
    geom_sf(data = ice_sf, aes(fill = FILL, color = FILL)) +# , linewidth = 0, alpha = 1
    # add land masses
    geom_sf(data = world, colour = "darkgrey", size = 0.2, fill = "oldlace") +
    # # add white underlines for warm/cold polys
    # geom_polygon(data = coldpoly, aes(x = long, y = lat, group = group), color = "white", fill = NA, linetype = "solid", size = 1.5) +
    # geom_polygon(data = warmpoly, aes(x = long, y = lat, group = group), color = "white", fill = NA, linetype = "solid", size = 1.5) +
    # 
    # # add warm/cold polys
    # geom_polygon_pattern(data = coldpoly, aes(x = long, y = lat, group = group), color = "cyan", fill = NA, pattern_fill = "cyan", pattern_density = 0.1, pattern_colour = NA, linetype = "dashed", linewidth = 0.7) +
    # geom_polygon_pattern(data = warmpoly, aes(x = long, y = lat, group = group), color = "orange3", fill = NA, pattern_fill = "orange3", pattern_density = 0.1, pattern_colour = NA, linetype = "dashed", linewidth = 0.7) +
    # geom_polygon_pattern(data = frontpoly, aes(x = long, y = lat, group = group), color = "grey", fill = NA, pattern_fill ="grey", pattern_density = 0.1, pattern_colour = NA,linetype = "dashed", linewidth = 0.7) +
    geom_polygon(data = coldpoly, aes(x = long, y = lat, group = group), color = "cyan", fill = NA, linewidth = 1, linetype = 'dashed') +
    geom_polygon(data = warmpoly, aes(x = long, y = lat, group = group), color = "darkorange3", fill = NA, linewidth = 1, linetype = 'dashed') +
    geom_polygon(data = frontpoly, aes(x = long, y = lat, group = group), color = "grey", fill = NA, linewidth = 1, linetype = 'dashed') +
    
    # set CRS
    coord_sf(crs = stereproj, xlim = c(-5e5, 5e5), ylim = c(-16.7e5, -5.3e5), #xlim war vorher 6e5
             expand = F) +
    scale_color_gradient2(aesthetics = c("color", "fill"),
                          low = "blue", mid = "white", high = "red", midpoint = 0.15,
                          breaks = c(0.15, 0.3, 0.6, 0.9),
                          labels = c(0.15, 0.3, 0.6, 0.9)*100,
                          limits = c(0, 1),
                          name = "SIC \n[%]") +
    guides(fill = guide_colorbar(title = "SIC [%]", title.position = "top", title.hjust = 0.5),
           color = guide_colorbar(title = "SIC [%]", title.position = "top", title.hjust = 0.5)) +
    scale_x_continuous(breaks = seq(-30, 30, by = 10)) +
    labs(title =paste0("SIC - ", date_highlight), x = NULL, y = NULL) +
    theme_bw() +
    theme(legend.position = "bottom",
          legend.box.margin = margin(t = 0),
          title = element_text(size = 8),
          panel.ontop = T,
          panel.grid = element_line(color = "snow", linewidth = 0.1)
          )
  
  
  # CHL 
  #dates <- as.Date(getZ(sst.nc))
  date <- dates[doy]
  chl.nc <- brick(paste0('data/', stringr::str_pad(month(date), 2, pad = "0"), '_cmems_obs-oc_glo_bgc-plankton_my_l3-multi-4km_P1D.nc'), values=T, 
                  varname = "CHL")
  
  chl_index <- which(getZ(chl.nc) == date)
  
  
  chl_d <- chl.nc[[chl_index]]
  
  
  names(chl_d) <- "FILL"
  # convert raster to specialfeatures (dieser umweg ist leider notwendig)
  chl_sf <- chl_d %>%
    as("SpatialPixelsDataFrame") %>%
    as("SpatialPolygonsDataFrame") %>%
    as("sf")
  
  
  chl_map <- 
    ggplot() +
    # chl layer
    geom_sf(data = chl_sf, aes(fill = FILL, color = FILL)) +# , linewidth = 0, alpha = 1
    # add land masses
    geom_sf(data = world, colour = "darkgrey", size = 0.2, fill = "oldlace") +
    # add white underlines for warm/cold polys
    # geom_polygon(data = coldpoly, aes(x = long, y = lat, group = group), color = "white", fill = NA, linetype = "solid", size = 1.5) +
    # geom_polygon(data = warmpoly, aes(x = long, y = lat, group = group), color = "white", fill = NA, linetype = "solid", size = 1.5) +
    # 
    # add warm/cold polys
    # geom_polygon_pattern(data = coldpoly, aes(x = long, y = lat, group = group), color = "cyan", fill = NA, pattern_fill = "cyan", pattern_density = 0.1, pattern_colour = NA, linetype = "dashed", linewidth = 0.7) +
    # geom_polygon_pattern(data = warmpoly, aes(x = long, y = lat, group = group), color = "orange3", fill = NA, pattern_fill = "orange3", pattern_density = 0.1, pattern_colour = NA, linetype = "dashed", linewidth = 0.7) +
    # geom_polygon_pattern(data = frontpoly, aes(x = long, y = lat, group = group), color = "grey", fill = NA, pattern_fill = "grey", pattern_density = 0.1, pattern_colour = NA,linetype = "dashed", linewidth = 0.7) +
    geom_polygon(data = coldpoly, aes(x = long, y = lat, group = group), color = "cyan", fill = NA, linewidth = 1, linetype = 'dashed') +
    geom_polygon(data = warmpoly, aes(x = long, y = lat, group = group), color = "darkorange3", fill = NA, linewidth = 1, linetype = 'dashed') +
    geom_polygon(data = frontpoly, aes(x = long, y = lat, group = group), color = "grey", fill = NA, linewidth = 1, linetype = 'dashed') +
    
    # set CRS
    coord_sf(crs = stereproj, xlim = c(-5e5, 5e5), ylim = c(-16.7e5, -5.3e5), #xlim war vorher 6e5
             expand = F) +
    # scale_fill_gradientn(aesthetics = c("color", "fill"), colours = rev(pals::ocean.algae(10)), #, limits = c(0,6), pals::ocean.algae(10) oce::oceColorsPalette(120)
    #                      na.value = "", name = "Chl \n[mg /m^3]", trans = "log10", 
    #                      breaks = c(0.1, 0.3, 1, 3, 10), limits = c(0.1, 10)) +
    scale_fill_gradientn(aesthetics = c("color", "fill"), 
                         colours = rev(pals::ocean.algae(100)), #, limits = c(0,6), pals::ocean.algae(10) oce::oceColorsPalette(120)
                         na.value = "", 
                         name = "Chl \n[mg /m^3]", 
                         trans = "log10", 
                         breaks = c(0.1, 0.3, 1, 3, 10), 
                         limits = c(0.1, 4),
                         oob = scales::squish
    ) +
    scale_x_continuous(breaks = seq(-30, 30, by = 10)) +
    labs(title = expression("Chl-"~italic(a)), x = NULL, y = NULL) +
    guides(fill = guide_colorbar(title = TeX("Chl [mg m$^{-3}$]"), 
                                 title.position = "top", title.hjust = 0.5),
           color = guide_colorbar(title = TeX("Chl [mg m$^{-3}$]"), 
                                  title.position = "top", title.hjust = 0.5)) +
    theme_bw() +
    theme(legend.position = "bottom",
          legend.box.margin = margin(t = 0),
          title = element_text(size = 8),
          panel.ontop = T, 
          panel.grid = element_line(color = "snow", linewidth = 0.1)
          )# size = 9

  

  
  # combine plots 
  gg <-  ggarrange(ggarrange(sst_ts_d + 
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
                                     legend.background = element_rect(fill = "black")
                             ),
                      ice_ts_d + 
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
                              legend.background = element_rect(fill = "black")), 
                      ncol = 2, common.legend = TRUE, legend = "bottom", align = "hv"), 
            ggarrange(sst_map +
                        theme(text = element_text(color = "white"),
                              panel.background = element_rect(fill = NA, color = "white"),
                              #panel.grid = element_line(color = "grey30"),
                              panel.border = element_rect(color = "white"),
                              plot.background = element_rect(fill = "black", color = "black"),
                              axis.text = element_text(color = "white"),
                              axis.ticks = element_line(color = "white"),
                              legend.background = element_rect(fill = "black")), 
                      ice_map +
                        theme(text = element_text(color = "white"),
                              panel.background = element_rect(fill = NA, color = "white"),
                              # panel.grid = element_line(color = "grey30"),
                              panel.border = element_rect(color = "white"),
                              plot.background = element_rect(fill = "black", color = "black"),
                              axis.text = element_text(color = "white"),
                              axis.ticks = element_line(color = "white"),
                              legend.background = element_rect(fill = "black")), 
                      chl_map + 
                        theme(text = element_text(color = "white"),
                              panel.background = element_rect(fill = NA, color = "white"),
                              # panel.grid = element_line(color = "grey30"),
                              panel.border = element_rect(color = "white"),
                              plot.background = element_rect(fill = "black", color = "black"),
                              axis.text = element_text(color = "white"),
                              axis.ticks = element_line(color = "white"),
                              legend.background = element_rect(fill = "black")), 
                      ncol = 3, align = "hv"),
            nrow = 2, heights = c(1.5,2)) # dauert
  
    cairo_pdf(filename = paste0("~/Documents/Verteidigung/plots/figures_paper3/animation_map2/",str_pad(counter, 3, pad = "0"),"_SST_Ice_ts_with_daily_maps_",date_highlight,"_DARK_cairo.pdf"),
              width = 1.772*3, height = 1.654*3, # in px at 300 dpi (default dpi is 100, so mulitply)
              bg = "black")
    print(gg)
  
    dev.off() # close cairo device
      # looks very good, all characters are fine, but in middle plot a unwanted line appears ... 
    
  # save
  # ggsave(paste0("output/plots/animation/",str_pad(counter, 3, pad = "0"),"_SST_Ice_ts_with_daily_maps_",date_highlight,"_DARK_whitepanels.pdf"),
  #        height = 14, width = 15, units = "cm", bg = "black", device = cairo_pdf, antialias = "default")
  # # looks good in file, but when making gif, the special characters get replaced
  # 
  # ggsave(paste0("output/plots/animation/",str_pad(counter, 3, pad = "0"),"_SST_Ice_ts_with_daily_maps_",date_highlight,"_DARK_whitepanels.png"),
  #        height = 14, width = 15, units = "cm", bg = "black", device = "png", dpi = 300)
  # # looks good in file, but colors are faded 
  



}



# to make gif: 
# magick -verbose -delay 50 -density 300 *.pdf gifname.gif
# see if using other parameters helps?
