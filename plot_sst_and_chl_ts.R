## plot daily maps of SST and Chl satellite observations
## to construct animation


library(raster)
library(pals)
library(tidyverse)
library(readxl)
library(sf)
library("rnaturalearth")
library("rnaturalearthdata")
library(ggrepel)
library(scales)


# use polygons and nc files
world <- ne_countries(scale = "medium", returnclass = "sf")
allpolys <- readRDS("polygons_2009-2018.rds")

overview <- read_excel("Output_warm_and_cold_cluster_analysis.xlsx")

years <- names(allpolys)


## create a basemap to ease plotting
# load data for basic map
basic_map <- readRDS("basic_map.RDS")
# projection 
stereproj <- basic_map$sterecrs

# plot empty baseplot
baseplot <- ggplot() +
  labs(
    y = "", 
    x = "") +
  # add land masses
  geom_sf(data = world, colour = "darkgrey", size = 0.2, fill = "grey") +
  
  # set y breaks
  scale_y_continuous(breaks = 75:83) +
  
  # set coord system
  coord_sf(crs = stereproj, xlim = c(-4e5, 4e5), ylim = c(-14.7e5, -9.3e5), #xlim war vorher 6e5
           expand = F) +
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




for(y in years) {
  print(paste0("begin plotting for y = ", y))
  polys <- allpolys[[paste0(y)]]
  
  # import data
  sst.nc <- brick(x = paste0("sst.day.mean.", y, ".nc"), values = T, 
                  varname = "sst")
  sst_crs <- sst.nc@crs@projargs
  
  chl_f_name <- if_else(y %in% c(2009:2011),
                        "erdMH1chla1day_", 
                        "nesdisVHNSQchlaDaily_")
  chl_var_name <- if_else(y %in% c(2009:2011),
                          "chlorophyll",
                          "chlor_a")
  chla.nc <- brick(x = paste0(chl_f_name, y, ".nc"), values = T, 
                   varname = paste0(chl_var_name))
  chla.nc@z[[1]] <- as.Date(as.POSIXct(chla.nc@z[[1]], origin = "1970-01-01 00:00:00", tz = "UTC"))
  
  # make a Julian Day - Date list (for lookups of chl a indices)
  # in case sst and chl file do not line up
  dates.sst <- data.frame(sst_index = 1:length(getZ(sst.nc)),
                          Date = getZ(sst.nc), 
                          julianday = julian(x = getZ(sst.nc), 
                                             origin = as.Date("2009-01-01")))
  dates.chl <- data.frame(chl_index = 1:length(getZ(chla.nc)),
                          Date = getZ(chla.nc), 
                          julianday = julian(x = getZ(chla.nc), 
                                             origin = as.Date("2009-01-01")))
  
  dates <- left_join(dates.sst, dates.chl, by = c("Date", "julianday")) %>%
    # in 2012, some days have two chl measurements (12pm or 23pm), the night should be excluded, else the following loop will fail
    distinct(sst_index, .keep_all = T)
  
  # now enter loop for each day of year
  for(d in 1:nrow(dates)){ 
    
    ## get polys
    jd <- paste0(dates$julianday[dates$sst_index==d])
    warmpoly <- spTransform(polys[[jd]]$warm, CRSobj = stereproj)
    coldpoly <- spTransform(polys[[jd]]$cold, CRSobj = stereproj)
    
    
    ## extract sst values for d
    sst_d <- rotate(sst.nc[[d]])
    sst_d <- crop(sst_d, extent(-25, 25, 76, 82)) # extent(-25, 25, 76, 82) is extent is for whole map; extent(-10, 12, 77, 81) is for polys
    today <- getZ(sst_d)
    # plot(sst_d)
    names(sst_d) <- "FILL"
    # convert raster to specialfeatures (dieser umweg ist leider notwendig)
    sst_sf <- sst_d %>%
      as("SpatialPixelsDataFrame") %>%
      as("SpatialPolygonsDataFrame") %>%
      as("sf")
    rm(sst_d)
    
    
    
    
    ## extract chl values for d
    chl_indx <- dates$chl_index[dates$sst_index == d] #[1] in case there are two measurements per day, take the earlier one (only in 2012, days 82 and 83 i think)
    # is there any chl data available for today?
    chl_today <- !is.na(chl_indx)
    # if yes, go on
    if(chl_today){
      chl_d <- chla.nc[[chl_indx]]
      # plot(chl_d)
      values(chl_d)[values(chl_d) < 0] = NA  # for 2012-2015 and 2019, NAvalues are not recognised correctly. R thinks NAvalue is -999, but fills table with -32767
      # plot(chl_d)
      #if all values are NA, set chl_today to F and break
      if(all(is.na(values(chl_d)))) {
        chl_today <- F
        } else {
      
          names(chl_d) <- "FILL"
          chl_sf <- chl_d %>%
            as("SpatialPixelsDataFrame") %>%
            as("SpatialPolygonsDataFrame") %>%
            as("sf") 
        }
    }
    
   sst_plot <- 
    ggplot() +
      labs(#title = paste0("Sea surface temperature and cold-warm separation on ", thisdate),
        #subtitle = "NOAA High-resolution Blended Analysis",
        y = "", 
        x = "") +
      # add sst layer
      geom_sf(data = sst_sf, aes(fill = FILL, colour = FILL), lwd = 0) +
      # geom_raster(data = sst_df, aes(x = x, y = y, fill = sst), alpha = 1) +
      
      # add colour scale for sill and colour
      scale_color_gradientn(aesthetics = c("color", "fill"),
                            colors = rev(rainbow(7)),
                            limits = c(-3.5, 11.5),
                            breaks = c(-2, 0, 2, 4, 6, 8, 10),
                            labels = format(c(-2, 0, 2, 4, 6, 8, 10)),
                            oob = scales::squish , na.value = NA,
                            name = "SST \n[\u00B0C]") +
      
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
   
   # if there is no chl observation for today, use an empty sf to create a plot with legend
    dummy_sf <- sst_sf  
    dummy_sf$FILL <- NA
    
   chl_plot <- 
      ggplot() +
        labs(#title = paste0("Sea surface temperature and cold-warm separation on ", thisdate),
          #subtitle = "NOAA High-resolution Blended Analysis",
          y = "", 
          x = "") +
        # add chl layer
        {if(chl_today) geom_sf(data = chl_sf, aes(fill = FILL, colour = FILL), lwd = 0)}+ 
        {if(!chl_today) geom_sf(data = dummy_sf, aes(fill = FILL, colour = FILL), lwd = 0)}+
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
                              name = "Chl-a \n[mg m^-3]") +
        
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
      
     jd <- str_pad(jd, 4, pad = "0")  
      ## then arrange and save plots
     ggpubr::ggarrange(sst_plot, chl_plot, 
                       ncol = 1, nrow = 2, align = "hv")
      
     cowplot::ggsave2(filename = paste0("figures/daily_maps/both/", jd, "_", today, ".png"),
                      width = 5, height = 5)    
      
     cowplot::ggsave2(plot = sst_plot, filename = paste0("figures/daily_maps/sst/", jd, "_", today, ".png"),
                      width = 5, height = 2.5) 
     cowplot::ggsave2(plot = chl_plot, filename = paste0("figures/daily_maps/chl/", jd, "_", today, ".png"),
                      width = 5, height = 2.5) 
     
    rm(sst_plot, chl_plot,
       chl_sf, sst_sf, dummy_sf,
       chl_indx, chl_today)
    } # end daily loop
  
} # end yearly loop


## because I forgot including leading zeros to the filenames, ffmpeg fails to
## order the pngs properly. (has been resolved for future runs)
## solution: create input list for ffmpeg

# og_fnames <- list.files("figures/daily_maps/both/")
# og_fnum <- str_extract(og_fnames, "^[:digit:]*")
# appen_fname <- str_extract(og_fnames, "_[:digit:]*-[:digit:]*-[:digit:]*.png")
# 
# new_fnum <- str_pad(og_fnum, 4, pad = "0")
# new_fname <- paste0(new_fnum, appen_fname)

## ensure nothing got mixed up
# tst <- tibble(og_fnames, og_fnum, appen_fname, new_fnum, new_fname)

## rename the files
# setwd("figures/daily_maps/both/")
# file.rename(from = og_fnames, to = new_fname)


# then create video in terminal
#ffmpeg -framerate 4 -pattern_type glob -i '*.png' -vf scale=790:-1 -c:v libx264 -pix_fmt yuv420p output2.mp4