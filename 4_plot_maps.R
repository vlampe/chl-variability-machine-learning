## some plotting ... how does ice correspond to polygons?
## make some maps

library(raster)
library(pals)
library(tidyverse)
library(sf)
library("rnaturalearth")
library("rnaturalearthdata")
library(ggrepel)
library(scales)
library(ggpattern)

# load ice data
ice.nc <- brick("data/METOFFICE-GLO-SST-L4-REP-OBS-SST_1695911795678.nc", values = T,
                varname = 'sea_ice_fraction')
sst.nc <- brick("data/METOFFICE-GLO-SST-L4-REP-OBS-SST_1695911795678.nc", values = T,
                varname = 'analysed_sst')

polys <- readRDS(file = "output/cluster_3k_2018.rds")



# use polygons and nc files
world <- ne_countries(scale = "medium", returnclass = "sf")

# ## create a basemap to ease plotting
# # load data for basic map
# basic_map <- readRDS("basic_map.RDS")
# 
# # projection 
# stereproj <- basic_map$sterecrs
stereproj <- "+proj=stere +lat_0=90 +lat_ts=75 +lon_0=0 +k=1 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0"
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
baseplot



# to plot contents of ncdf files, I need to convert a layer to a sf object first

doy <- 200

warmpoly <- spTransform(polys[[doy]][["warm"]], CRSobj = stereproj)
coldpoly <- spTransform(polys[[doy]][["cold"]], CRSobj = stereproj)
frontpoly <- spTransform(polys[[doy]][["front"]], CRSobj = stereproj)
plot(warmpoly)

# SST 
sst_d <- sst.nc[[doy]]
names(sst_d) <- "FILL"
# convert raster to specialfeatures (dieser umweg ist leider notwendig)
sst_sf <- sst_d %>%
  as("SpatialPixelsDataFrame") %>%
  as("SpatialPolygonsDataFrame") %>%
  as("sf")



ggplot() +
  # SST layer
  geom_sf(data = sst_sf, aes(fill = FILL, color = FILL), lwd = 0) +
  # add land masses
  geom_sf(data = world, colour = "darkgrey", size = 0.2, fill = "grey") +
  # add white underlines for warm/cold polys
  geom_polygon(data = coldpoly, aes(x = long, y = lat, group = group), color = "white", fill = NA, linetype = "solid", size = 1.5) +
  geom_polygon(data = warmpoly, aes(x = long, y = lat, group = group), color = "white", fill = NA, linetype = "solid", size = 1.5) +
  
  # add warm/cold polys
  geom_polygon_pattern(data = coldpoly, aes(x = long, y = lat, group = group), color = "blue", fill = NA, pattern_fill = "blue", pattern = 'circle', linetype = "dashed", size = 1) +
  geom_polygon_pattern(data = warmpoly, aes(x = long, y = lat, group = group), color = "red", fill = NA, pattern_fill = "red", pattern = 'circle', linetype = "dashed", size = 1) +
  geom_polygon_pattern(data = frontpoly, aes(x = long, y = lat, group = group), color = "grey", fill = NA, pattern = 'circle',linetype = "dashed", size = 1) +
  
  # set CRS
  coord_sf(crs = stereproj, xlim = c(-4e5, 4e5), ylim = c(-14.7e5, -9.3e5), #xlim war vorher 6e5
           expand = F) +
  # add colour scale for sill and colour
  scale_color_gradientn(aesthetics = c("color", "fill"),
                        colors = rev(rainbow(7)),
                        limits = c(-3.5, 11.5)+273.13,
                        breaks = c(-2, 0, 2, 4, 6, 8, 10)+273.13,
                        labels = format(c(-2, 0, 2, 4, 6, 8, 10)), # +273.13
                        oob = scales::squish , na.value = NA,
                        name = "SST \n[ºC]") 



# ICE
ice_d <- ice.nc[[doy]]
names(ice_d) <- "FILL"
# convert raster to specialfeatures (dieser umweg ist leider notwendig)
ice_sf <- ice_d %>%
  as("SpatialPixelsDataFrame") %>%
  as("SpatialPolygonsDataFrame") %>%
  as("sf")

ggplot() +
  # ICE layer
  geom_sf(data = ice_sf, aes(fill = FILL, color = FILL), lwd = 0) +
  # add land masses
  geom_sf(data = world, colour = "darkgrey", size = 0.2, fill = "grey") +
  # add white underlines for warm/cold polys
  geom_polygon(data = coldpoly, aes(x = long, y = lat, group = group), color = "white", fill = NA, linetype = "solid", size = 1.5) +
  geom_polygon(data = warmpoly, aes(x = long, y = lat, group = group), color = "white", fill = NA, linetype = "solid", size = 1.5) +
  
  # add warm/cold polys
  geom_polygon_pattern(data = coldpoly, aes(x = long, y = lat, group = group), color = "blue", fill = NA, pattern_fill = "blue", pattern = 'circle', linetype = "dashed", size = 1) +
  geom_polygon_pattern(data = warmpoly, aes(x = long, y = lat, group = group), color = "red", fill = NA, pattern_fill = "red", pattern = 'circle', linetype = "dashed", size = 1) +
  geom_polygon_pattern(data = frontpoly, aes(x = long, y = lat, group = group), color = "grey", fill = NA, pattern = 'circle',linetype = "dashed", size = 1) +
  
  # set CRS
  coord_sf(crs = stereproj, xlim = c(-4e5, 4e5), ylim = c(-14.7e5, -9.3e5), #xlim war vorher 6e5
           expand = F) +
  scale_color_gradient2(aesthetics = c("color", "fill"),
                        low = "blue", mid = "white", high = "red", midpoint = 0.15,
                        breaks = c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9),
                        limits = c(0, 1),
                        name = "Ice\n (%)")
  
  # scale_fill_gradientn(aesthetics = c("color", "fill"), colours = oce::oceColorsPalette(120), limits = c(0,1),
  #                      na.value = "", name = "Ice\n (%)")
  # scale_fill_gradientn(aesthetics = c("color", "fill"), colours = pals::ocean.ice(100), limits = c(0,1),
  #                      breaks = c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9),
  #                      na.value = "", name = "Ice\n (%)")


# CHL 

dates <- as.Date(getZ(sst.nc))
date <- dates[doy]
chl_index <- which(getZ(chl.nc) == date)

chl.nc <- brick(paste0('data/', stringr::str_pad(month(date), 2, pad = "0"), '_cmems_obs-oc_glo_bgc-plankton_my_l3-multi-4km_P1D.nc'), values=T, 
                varname = "CHL")
chl_d <- chl.nc[[chl_index]]


names(chl_d) <- "FILL"
# convert raster to specialfeatures (dieser umweg ist leider notwendig)
chl_sf <- chl_d %>%
  as("SpatialPixelsDataFrame") %>%
  as("SpatialPolygonsDataFrame") %>%
  as("sf")


ggplot() +
  # chl layer
  geom_sf(data = chl_sf, aes(fill = FILL, colour = FILL), lwd = 0) +
  # add land masses
  geom_sf(data = world, colour = "darkgrey", size = 0.2, fill = "grey") +
  # add white underlines for warm/cold polys
  geom_polygon(data = coldpoly, aes(x = long, y = lat, group = group), color = "white", fill = NA, linetype = "solid", size = 1.5) +
  geom_polygon(data = warmpoly, aes(x = long, y = lat, group = group), color = "white", fill = NA, linetype = "solid", size = 1.5) +

  # add warm/cold polys
  geom_polygon_pattern(data = coldpoly, aes(x = long, y = lat, group = group), color = "blue", fill = NA, pattern_fill = "blue", pattern = 'circle', linetype = "dashed", size = 1) +
  geom_polygon_pattern(data = warmpoly, aes(x = long, y = lat, group = group), color = "red", fill = NA, pattern_fill = "red", pattern = 'circle', linetype = "dashed", size = 1) +
  geom_polygon_pattern(data = frontpoly, aes(x = long, y = lat, group = group), color = "grey", fill = NA, pattern = 'circle',linetype = "dashed", size = 1) +

  # set CRS
  coord_sf(crs = stereproj, xlim = c(-4e5, 4e5), ylim = c(-14.7e5, -9.3e5), #xlim war vorher 6e5
           expand = F) +
  scale_fill_gradientn(aesthetics = c("color", "fill"), colours = oce::oceColorsPalette(120), #, limits = c(0,6), pals::ocean.algae(10)
                      na.value = "", name = "Chl\n (mg /m^3)", trans = "log10")


