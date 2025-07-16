## illustration of right-skewed Chl-a Values

# this is to adress R2 comment:
# Comment R2, 15:
#   L294 ++Right-skewed means and medians: would be easier to see if a figure was added to
# the appendix that illustrated this point specifically



# plot pdfs on example day of sst, ice and Chl
library(raster)
library(sf)
require(latex2exp)
require(reticulate)
library("rnaturalearth")
library("rnaturalearthdata")
use_python("/opt/homebrew/Cellar/python@3.11/3.11.11/Frameworks/Python.framework/Versions/3.11/Resources/Python.app/Contents/MacOS/Python")
np      <- import("numpy",   convert=FALSE)
py_csv  <- import("csv",     convert=FALSE)
diffKDE <- import("diffKDE", convert=FALSE)


# get data
## read output from extract_warm_and_cold_clusters.R
ts_dat <- read_csv2("output/Output_3k_cluster_analysis_thesis.csv", na = c("", "Inf", "-Inf","NA"))


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



# load ice and sst data
ice.nc <- brick("data/METOFFICE-GLO-SST-L4-REP-OBS-SST_1695911795678.nc", values = T,
                varname = 'sea_ice_fraction')
sst.nc <- brick("data/METOFFICE-GLO-SST-L4-REP-OBS-SST_1695911795678.nc", values = T,
                varname = 'analysed_sst')

polys <- readRDS(file = "output/cluster_3k_2018.rds")
stereproj <- "+proj=stere +lat_0=90 +lat_ts=75 +lon_0=0 +k=1 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0"
# use polygons and nc files
world <- ne_countries(scale = "medium", returnclass = "sf")

# set doy for focus: 
# doy <- 208
# doy <- 171
# doy <- 134
# doy <- 139
doy <- 144

dates <- as.Date(getZ(sst.nc))
date_highlight <- dates[doy]

# set alpha for confidence intervals: 0.3 for light, 0.7 for dark
alpha_ci <-  0.7


warmpoly <- spTransform(polys[[doy]][["warm"]], CRSobj = stereproj)
coldpoly <- spTransform(polys[[doy]][["cold"]], CRSobj = stereproj)
frontpoly <- spTransform(polys[[doy]][["front"]], CRSobj = stereproj)


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
  geom_sf(data = chl_sf, aes(fill = FILL, color = FILL), lwd = 0) +
  # add land masses
  geom_sf(data = world, colour = "grey", size = 0.2, fill = "lightgrey") +
  # add white underlines for warm/cold polys
  # geom_polygon(data = coldpoly, aes(x = long, y = lat, group = group), color = "white", fill = NA, linetype = "solid", size = 1.5) +
  # geom_polygon(data = warmpoly, aes(x = long, y = lat, group = group), color = "white", fill = NA, linetype = "solid", size = 1.5) +
  # 
  # add warm/cold polys
  geom_polygon(data = frontpoly, aes(x = long, y = lat, group = group), color = "cornsilk4", fill = NA, linetype = "solid", linewidth = 0.7) +
  geom_polygon(data = coldpoly, aes(x = long, y = lat, group = group), color = "blue", fill = NA, linetype = "solid", linewidth = 0.7) +
  geom_polygon(data = warmpoly, aes(x = long, y = lat, group = group), color = "red", fill = NA, linetype = "solid", linewidth = 0.7) +
  geom_polygon(data = frontpoly, aes(x = long, y = lat, group = group), color = "cornsilk4", fill = NA, linetype = "dashed", linewidth = 0.5) +
  # geom_polygon_pattern(data = coldpoly, aes(x = long, y = lat, group = group), color = "blue", fill = NA, pattern_fill = "blue", pattern_density = 0.1, pattern_colour = NA, linetype = "dashed", linewidth = 0.7) +
  # geom_polygon_pattern(data = warmpoly, aes(x = long, y = lat, group = group), color = "red", fill = NA, pattern_fill = "red", pattern_density = 0.1, pattern_colour = NA, linetype = "dashed", linewidth = 0.7) +
  # geom_polygon_pattern(data = frontpoly, aes(x = long, y = lat, group = group), color = "cornsilk4", fill = NA, pattern_fill = "cornsilk4", pattern_density = 0.1, pattern_colour = NA,linetype = "dashed", linewidth = 0.7) +
  
  # set CRS
  coord_sf(crs = stereproj, xlim = c(-5e5, 5e5), ylim = c(-16.7e5, -5.3e5), #xlim war vorher 6e5
           expand = F) +
  # scale_fill_gradientn(aesthetics = c("color", "fill"), colours = rev(pals::ocean.algae(10)), #, limits = c(0,6), pals::ocean.algae(10) oce::oceColorsPalette(120)
  #                      na.value = "", name = "Chl \n[mg /m^3]", trans = "log10", 
  #                      breaks = c(0.1, 0.3, 1, 3, 10), limits = c(0.1, 10)) +
  scale_fill_gradientn(aesthetics = c("color", "fill"), 
                       colours = rev(pals::ocean.algae(5)), #, limits = c(0,6), pals::ocean.algae(10) oce::oceColorsPalette(120)
                       na.value = "", 
                       name = "Chl \n[mg /m^3]", 
                       trans = "log10", 
                       breaks = c(0.1, 0.3, 1, 3, 10), 
                       limits = c(0.1, 4),
                       oob = scales::squish
  ) +
  scale_x_continuous(breaks = seq(-30, 30, by = 10)) +
  labs(title = expression("Chl-"*italic(a)), x = NULL, y = NULL) +
  guides(fill = guide_colorbar(title = TeX("Chl-$a$ [mg m$^{-3}$]"), 
                               title.position = "top", title.hjust = 0.5),
         color = guide_colorbar(title = TeX("Chl-$a$ [mg m$^{-3}$]"), 
                                title.position = "top", title.hjust = 0.5)) +
  theme_bw() +
  theme(legend.position = "bottom",
        legend.margin=margin(0,0,0,0),
        legend.box.margin=margin(-10,0,0,0),
        title = element_text(size = 8))# size = 9
chl_map


# extract sst, sic, chl per cluster 
raw_vals <- list()

raw_vals[["Chl"]][["cold"]] <- getValues(mask(chl_d, polys[[doy]][["cold"]]))[!is.na(getValues(mask(chl_d, polys[[doy]][["cold"]])))]
raw_vals[["Chl"]][["front"]] <- getValues(mask(chl_d, polys[[doy]][["front"]]))[!is.na(getValues(mask(chl_d, polys[[doy]][["front"]])))]
raw_vals[["Chl"]][["warm"]] <- getValues(mask(chl_d, polys[[doy]][["warm"]]))[!is.na(getValues(mask(chl_d, polys[[doy]][["warm"]])))]




vars <- names(raw_vals)
clusters <- names(raw_vals[[1]])

KDEs_list <- list()

for (v in vars){
  
  ## min max per var
  switch (v,
          "Chl" = {kde_min = 0
          kde_max = 20},
          "SST" = {kde_min = -5
          kde_max = 15},
          "SIC" = {kde_min = 0
          kde_max = 100}
  )
  
  for (c in clusters){
    
    # get pdf using diffKDE
    
    pdf_list <- py_to_r(diffKDE$KDE(np$array(raw_vals[[v]][[c]]), kde_min, kde_max))
    
    if (v == "SIC" & (c == "warm" | c == "front")) { 
      # repeat KDE with larger T final iteration time to avoid ragged look
      pdf_list <- py_to_r(diffKDE$KDE(np$array(raw_vals[[v]][[c]]), kde_min, kde_max, 1004L, 20L, 0.05))
    }
    
    #plot(pdf_list[[2]],pdf_list[[1]], xlim = c(0,25), ylim = c(0, 5), type = "l")
    
    pdf_x <- pdf_list[[2]] # convert numpy array to R array
    pdf_y <- pdf_list[[1]]
    
    # KDEs_list[[v]][[c]][["x"]] <- pdf_x
    # KDEs_list[[v]][[c]][["y"]] <- pdf_y
    
    KDEs_list[[v]][[c]] <- tibble("x"= pdf_x, "y" = pdf_y)
    
  }
}


df1 <- bind_rows(KDEs_list[[1]], .id = c("cluster")) %>%
  mutate(variable = names(KDEs_list)[1] )
# df2 <- bind_rows(KDEs_list[[2]], .id = c("cluster")) %>%
#   mutate(variable = names(KDEs_list)[2] )
# df3 <- bind_rows(KDEs_list[[3]], .id = c("cluster")) %>%
#   mutate(variable = names(KDEs_list)[3] )

# KDE_table <-rbind(df1, df2, df3)

var_dat <- results %>%
  filter(date == date_highlight)

raw_vals_df <- data.frame(vals = raw_vals$Chl$cold, cluster = "cold") %>%
  bind_rows(data.frame(vals = raw_vals$Chl$front, cluster = "front"))%>%
  bind_rows(data.frame(vals = raw_vals$Chl$warm, cluster = "warm"))

ggplot(data = df1, aes(x = x, y = y, fill = `cluster`, color = `cluster`)) + # fill = `cluster` color = `cluster`
  geom_rug(data = raw_vals_df, aes(x = vals, y = NULL), linewidth = 0.1) +
  geom_rect(data = var_dat, aes(xmin = `IQR_l`, 
                  xmax = `IQR_u`, ymin = 0, ymax = Inf, x = NULL, y = NULL), alpha = .4, linetype = 0) +
  geom_rect(data = var_dat, aes(xmin = `CI_l`,
                  xmax = `CI_u`, ymin = 0, ymax = Inf, x = NULL, y = NULL), alpha = .2, linetype = 0) +
  geom_line() +
  geom_vline(aes(xintercept = mean), var_dat, linetype = "dotted") +
  geom_vline(aes(xintercept = median), var_dat, linetype = "solid") +
  labs(y = "prob. density") +
  scale_color_manual(aesthetics =c("colour", "fill"), values = c("warm" = "red", "cold" = "blue","front" = "darkgrey"), 
                     breaks = c("cold", "front", "warm")) +
  facet_wrap(~cluster, scales = "free") +
  scale_x_continuous(limits=c(0,7), name = TeX("Chl \\textit{a} (mg m$^{-3}$)")) +
  # scale_x_log10(limits=c(0.1, 10)) +
  theme_bw() +
  theme(legend.position = "none")
# ggsave(paste0("output/plots/appendix_rightskewedChla.pdf"),
#        height = 6, width = 15, units = "cm")


# alternative, less shaded area
ggplot(data = df1, aes(x = x, y = y, fill = `cluster`, color = `cluster`)) + # fill = `cluster` color = `cluster`
  geom_rug(data = raw_vals_df, aes(x = vals, y = NULL), linewidth = 0.1, length = unit(0.03, "npc"),
           outside = F, sides = "b") +
  geom_rect(data = var_dat, aes(xmin = `IQR_l`, 
                                xmax = `IQR_u`, ymin = -Inf, ymax = 0, x = NULL, y = NULL), alpha = .4, linetype = 0) +
  geom_rect(data = var_dat, aes(xmin = `CI_l`,
                                xmax = `CI_u`, ymin = -Inf, ymax = 0, x = NULL, y = NULL), alpha = .2, linetype = 0) +
  geom_line() +
  geom_vline(aes(xintercept = mean), var_dat, linetype = "dotted") +
  geom_vline(aes(xintercept = median), var_dat, linetype = "solid") +
  labs(y = "prob. density") +
  scale_color_manual(aesthetics =c("colour", "fill"), values = c("warm" = "red", "cold" = "blue","front" = "darkgrey"), 
                     breaks = c("cold", "front", "warm")) +
  facet_wrap(~cluster, scales = "free") +
  scale_x_continuous(limits=c(0,7), name = TeX("Chl \\textit{a} (mg m$^{-3}$)")) +
  scale_y_continuous(expand = expansion(mult=0.12)) +
  # coord_cartesian(expand = F, ylim = c(-1, NA)) +
  # scale_x_log10(limits=c(0.1, 10)) +
  theme_bw() +
  theme(legend.position = "none")
ggsave(paste0("output/plots/appendix_rightskewedChla.pdf"),
       height = 6, width = 15, units = "cm")



