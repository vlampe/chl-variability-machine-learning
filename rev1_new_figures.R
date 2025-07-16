## create new figures based on the anonymuous reviewers suggestions.

setwd("~/Documents/microARC/Manuscripts/Last_chapter_clustering_manuskript_skripts/Revision 1/")

library(tidyverse)
#library(ggOceanMaps)
library(raster)

# 1. overview map of the study area, as suggest by R2, 9
## load data
# load fine land data
NEDPath <- "~/shapefiles/10m_physical" # Natural Earth Data location, downloaded from https://www.naturalearthdata.com/downloads/10m-physical-vectors/

world <- sf::read_sf(paste(NEDPath, "ne_10m_land.shp", sep = "/"))
islands <- sf::read_sf(paste(NEDPath, "ne_10m_minor_islands.shp", sep = "/"))
world <- rbind(world, islands)



#larger product for Fram Strait Map
ice.nc <- brick("~/Documents/THESIS_clean/Thesis_extra_Figures/maps of bathy and ice/METOFFICE-GLO-SST-L4-REP-OBS-SST_multi-vars_49.97W-49.97E_65.03N-84.97N_2018-01-01-2018-12-31.nc", values = T,
                varname = 'sea_ice_fraction')


## aggregate by month

#ice
groups <- month(as_datetime(ice.nc@z$`Date/time`))
month_mean_ice <- raster::stackApply(ice.nc, indices = groups, fun = raster::mean, na.rm = T)
names(month_mean_ice) <- month(unique(groups), label = T) 


## make high resolution bathy map of fram strait

# load fine bathymetry
bathy_raster_path <- "~/shapefiles/IBCAO_v4_2_13_200m.nc" # downloaded from https://www.gebco.net/data_and_products/gridded_bathymetry_data/arctic_ocean/
bathy_raster_path <- "~/shapefiles/IBCAO_v4_2_13_200m_ice.nc" # with greenland ice eleveation


bathy_raster <- raster::raster(bathy_raster_path)
crs(bathy_raster) <- "EPSG:3996" # the projection of IBCAO data is 3996, not 3995
b <- as(extent(-20, 20, 72, 86), 'SpatialPolygons') 
crs(b) <- "EPSG:4326"
b <- sp::spTransform(b, "EPSG:3996")

# crop bathy_raster to smaller size
bathy_raster_c <- raster::crop(bathy_raster, b)

g_raster <- cut(bathy_raster_c, breaks = c(-8000,0,50,100,200,300,500,1000,1500,2000,4000,6000,8000)*-1)
# test3 <- rasterToPolygons(g_raster, dissolve = T)

terra::contour(g_raster) 
#test7 <- terra::as.contour(as(g_raster, 'SpatRaster'))
g_spatRas <- terra::as.polygons(as(g_raster, 'SpatRaster'))
crs(g_spatRas) <- "EPSG:3996"
g_sf <- sf::st_as_sf(g_spatRas)
g_sf$layer <- factor(g_sf$layer)

extent(g_raster)

# test5 <- g_raster %>% # dauert lange oder geht vllt gar nicht
#   as("SpatialPixelsDataFrame") %>%
#   as("SpatialPolygonsDataFrame") %>%
#   as("sf")

# plot fine IBCAO data
b_plot <- ggplot() +
  geom_sf(data = g_sf, aes(fill = layer, color = layer)) +
  # geom_sf(data = world, color = "cornsilk3", fill="cornsilk3") +
  scale_fill_manual(aesthetics = c("fill", "color"),
                    name = "Depth (m)",
                    breaks = c(1:12),
                    values = c(oce::oceColorsGebco(n = 11, region = "water"), "cornsilk3"),
                    labels = rev(c("land", "0-50", "50-100","100-200","200-300","300-500","500-1000","1000-1500","1500-2000","2000-4000","4000-6000", "6000-8000"))) +
  # scale_fill_manual(aesthetics = c("fill", "color"),
  #                   name = "Depth (m)",
  #                   breaks = c(2:12),
  #                   values = oce::oceColorsGebco(n = 11, region = "water"),
  #                   labels = rev(c("0-50", "50-100","100-200","200-300","300-500","500-1000","1000-1500","1500-2000","2000-4000","4000-6000", "6000-8000"))) + 
  scale_x_continuous(breaks = seq(-60,60,by=15), limits = c(-681400, 681400)) +
  scale_y_continuous(breaks = seq(70,88,by=2), limits = c(-1771800, -519600))+
  coord_sf(crs = "EPSG:3996", expand = F) +
  ggOceanMaps::theme_map(grid.col = "#606060", grid.size = 0.1) +
  theme(panel.ontop = T) +
  ggspatial::annotation_scale(location = "br") #+ 
#  ggspatial::annotation_north_arrow(location = "tr", which_north = "true") 
b_plot

crs(world)



## get 15% ice contours for Apr, May and Sept


plot(month_mean_ice[[9]])
crs(month_mean_ice)
# apr
ice_apr_spatVector <- terra::as.contour(as(month_mean_ice[[4]], 'SpatRaster'), levels = c(0.15, 0.8))
ice_apr_sf <- sf::st_as_sf(ice_apr_spatVector)
ice_apr_sf$level <- factor(ice_apr_sf$level)
# may
ice_may_spatVector <- terra::as.contour(as(month_mean_ice[[5]], 'SpatRaster'), levels = c(0.15, 0.8))
ice_may_sf <- sf::st_as_sf(ice_may_spatVector)
ice_may_sf$level <- factor(ice_may_sf$level)
# sep
ice_sep_spatVector <- terra::as.contour(as(month_mean_ice[[9]], 'SpatRaster'), levels = c(0.15, 0.8))
ice_sep_sf <- sf::st_as_sf(ice_sep_spatVector)
ice_sep_sf$level <- factor(ice_sep_sf$level)
ice_sep_sf <- sf::st_transform(ice_sep_sf, "EPSG:3996")

b_plot + 
  geom_sf(data = ice_apr_sf, aes(linetype = level), color = "red") +
  geom_sf(data = ice_sep_sf, aes(linetype = level), color = "darkred") +
  geom_sf(data = world, color = "cornsilk3", fill="cornsilk3") +
  scale_linetype_manual(name = "SIC", 
                        breaks = c(0.15, 0.8), 
                        values = c("dotted", "longdash")
  ) +
  coord_sf(crs = "EPSG:3996", expand = F) 
# ggsave("fram_bathy.pdf", device = "pdf", )



# study area box 
library(sf)
polys <- readRDS(file = "../skripts/chl-variability-machine-learning/output/cluster_3k_2018.rds")

ix <- 1
# sa <- st_as_sf(polys[[ix]]$warm) %>%
#   st_union(st_as_sf(polys[[ix]]$cold)) %>%
#   st_union(st_as_sf(polys[[ix]]$front)) 
# produces random spots, probably where clusters do not precisely match up?

sf_use_s2(FALSE)
sa <- st_as_sf(polys[[ix]]$warm) %>%
  st_union(st_as_sf(polys[[ix]]$cold)) %>%
  st_union(st_as_sf(polys[[ix]]$front)) %>%
  summarise(geometry = st_union(geometry)) 


# sa <- as(extent(-15, 15, 75, 85), 'SpatialPolygons')
# crs(sa) <- "EPSG:4326"
# sa <- st_as_sfc(sa)
# sa



# sa <- polys[[ix]]$warm %>% 
#   left_join(y=polys[[ix]]$warm)
# does not work
# test <- as_tibble(polys[[ix]]$warm)

b_plot + 
  geom_sf(data = ice_apr_sf, aes(linetype = level), color = "red", linewidth = 0.7) +
  geom_sf(data = ice_sep_sf, aes(linetype = level), color = "darkred", linewidth = 0.7) +
  geom_sf(data = world, color = "pink", fill="cornsilk3", linewidth=0.75) +
  scale_linetype_manual(name = "SIC", 
                        breaks = c(0.15, 0.8), 
                        values = c("dotted", "longdash")
  ) +
  geom_sf(data = sa, color = "black" , alpha = 0, linewidth = 0.8 ) +
  coord_sf(crs = "EPSG:3996", expand = F) 
# Natural Earth data and IBCAO do not seem to overlap perfectly in NE Greenland...

## add selected HAUSGARTEN Stations
stations <- read_csv2("../skripts/chl-variability-machine-learning/stationlist_2018.csv") %>%
  st_as_sf(crs="EPSG:4326", coords=c("long", "lat"))

pa <- b_plot + 
  geom_sf(data = ice_may_sf, aes(linetype = level), color = "orangered", linewidth = 0.7) +
  geom_sf(data = ice_sep_sf, aes(linetype = level), color = "darkred", linewidth = 0.7) +
  geom_sf(data = world, color = NA, fill="cornsilk3") +
  geom_sf(data = stations, color = "gold", shape=17)+
  scale_linetype_manual(name = "SIC (%)", 
                        breaks = c(0.15, 0.8), labels = c(15, 80),
                        values = c("dotted", "longdash")
  ) +
  geom_sf(data = sa, color = "black" , alpha = 0, linewidth = 0.8 ) +
  guides(fill = guide_legend(reverse = T, label.position = "left"), 
         color = guide_legend(reverse = T, label.position = "left"),
         linetype = guide_legend(label.position = "left")) +
  coord_sf(crs = "EPSG:3996", expand = F) 
pa
# ggsave("map_fs_ice.pdf", device = "pdf")

cairo_pdf("map_fs_ice.pdf", width = 6, height = 4.2, bg = NA)
pa
dev.off()

# make panarctic map with similar colors
library(ggOceanMaps)
# dp <- as(extent(-20, 20, 72, 82), 'SpatialPolygons')
# crs(dp) <- "EPSG:4326"
# dp <- st_as_sf(dp)
# 
# 
# dp <- as(extent(g_sf), "SpatialPolygons")
# crs(g_sf)
# crs(dp) <- crs(g_sf)# "EPSG:3996"
# dp <- st_as_sf(dp)

# scale_x_continuous(breaks = seq(-60,60,by=15), limits = c(-681400, 681400)) +
#   scale_y_continuous(breaks = seq(70,88,by=2), limits = c(-1871800, -519600))+
  
coords <- matrix( data = 
  c(-681400, -1771800, # LL
  -681400, -519600, # UL
  681400, -519600, #UR
  681400, -1771800, # LR
  -681400, -1771800), # LL, close polygon
  ncol=2,
  byrow = TRUE
)
dp <- st_polygon(list(coords))  
dp <- st_sf(geometry = st_sfc(dp), crs = crs(g_sf))
crs(dp) 
  
  # sa <- as(extent(-15, 15, 75, 85), 'SpatialPolygons')
  # crs(sa) <- "EPSG:4326"
  # sa <- st_as_sfc(sa)
  
p <- basemap(60, bathymetry = T, glaciers = T, shapefiles = "Arctic", land.col = "cornsilk3") + 
  scale_fill_manual(aesthetics = c("fill"),
                    name = "Depth (m)",
                    values = rev(c(oce::oceColorsGebco(n =8, region = "water")))) +
  # geom_sf(data = sf::st_as_sf(dp), aes(),color = "darkred" , alpha = 0, linewidth = 1 ) + 
  geom_sf(data = dp, aes(),color = "gold" , alpha = 0, linewidth = 1 ) + 
    theme(legend.position = "right")
p  

p <- p +
  geom_text(
    data =
      attributes(p)$map.grid$lat.grid.lines %>% 
      sf::st_coordinates() %>% 
      as.data.frame() %>% 
      dplyr::filter(X == -180) %>%
      dplyr::add_row(X = -180, Y = 90, L1 = 3, L2 = 1) %>% 
      transform_coord(proj.out = attributes(p)$crs, bind = TRUE),
    aes(x = lon.proj, y = lat.proj, label = paste0(abs(Y), "\u00B0N")),
    color = "grey20"
  ) +
  geom_text(
    data =
      attributes(p)$map.grid$lon.grid.lines %>% 
      sf::st_coordinates() %>% 
      as.data.frame() %>% 
      dplyr::filter(Y == 60) %>%
      dplyr::mutate(Y = 64) %>% 
      transform_coord(proj.out = attributes(p)$crs, bind = TRUE),
    aes(x = lon.proj, y = lat.proj, 
        label = ifelse(sign(X) == 1, paste0(abs(X), "\u00B0E"), 
                       paste0(abs(X), "\u00B0W"))
    ),
    color = "grey20"
  ) 
p


cairo_pdf("panarctic_bathy.pdf", width=4, height = 3.5, bg = NA)
p
dev.off()

# ggsave("panarctic_bathy.pdf", device = "pdf")