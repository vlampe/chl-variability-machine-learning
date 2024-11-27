# read and analyse ncdf files of SST and chl a daily values
# extract clusters of warm and cold sections of sampling site
# calculate variables of interest

# author: Vanessa Lampe
# based on https://git.geomar.de/vanessa-lampe/chlorophyll-time-series

#setwd("Users/vlampe/Documents/Last_chapter_clustering_manuskript_skripts/skripts/")
#setwd("Documents/Last_chapter_clustering_manuskript_skripts/skripts/")
library(raster)
library(tidyverse) 
# 
# # python -m motuclient --motu https://my.cmems-du.eu/motu-web/Motu --service-id ARCTIC_MULTIYEAR_PHY_002_003-TDS --product-id cmems_mod_arc_phy_my_topaz4_P1D-m --longitude-min -24.189033492189047 --longitude-max 16.911226280619612 --latitude-min 74.88458885266205 --latitude-max 81.6996554900983 --date-min "2018-01-01 00:00:00" --date-max "2018-12-31 23:59:59" --depth-min 0 --depth-max 0 --variable siconc --variable sithick --variable so --variable thetao --out-dir <OUTPUT_DIRECTORY> --out-name <OUTPUT_FILENAME> --user <USERNAME> --pwd <PASSWORD>
#   
# #python3 -m motuclient --motu https://my.cmems-du.eu/motu-web/Motu --service-id OCEANCOLOUR_GLO_BGC_L3_MY_009_103-TDS --product-id cmems_obs-oc_glo_bgc-plankton_my_l3-multi-4km_P1D --longitude-min -15 --longitude-max 15 --latitude-min 75 --latitude-max 85 --date-min '2018-01-01 00:00:00' --date-max '2018-12-31 23:59:59' --variable CHL --variable CHL_uncertainty --variable DIATO --variable DIATO_uncertainty --variable DINO --variable DINO_uncertainty --variable GREEN --variable GREEN_uncertainty --variable HAPTO --variable HAPTO_uncertainty --variable MICRO --variable MICRO_uncertainty --variable NANO --variable NANO_uncertainty --variable PICO --variable PICO_uncertainty --variable PROCHLO --variable PROCHLO_uncertainty --variable PROKAR --variable PROKAR_uncertainty --variable flags --out-dir '~' --out-name 'oceancolor.ncdf' --user <user> --pwd <password>
# 
# setwd('Documents/Last_chapter_clustering_manuskript_skripts /skripts/')
# 
# # load sst data
# sst.nc <- brick("data/METOFFICE-GLO-SST-L4-REP-OBS-SST_1695911795678.nc", values = T,
#                 varname = 'analysed_sst')
# 
# chl.nc <- brick('data/06_cmems_obs-oc_glo_bgc-plankton_my_l3-multi-4km_P1D_1695911037858.nc', values=T, 
#                 varname = "CHL")
# 
# sst_crs <- sst.nc@crs@projargs
# sst_crs == chl.nc@crs@projargs
# plot(sst.nc[[100]])
# 
# 
# 
# plot(chl.nc[[8]])
# 
# vals =getValues(chl.nc[[8]]) # need to remove nans first
# vals_nona <- na.omit(vals)
# max(vals, na.rm =T)
# 
# hist(vals, nclass=100)
# hist(vals_nona, nclass=100)
# 
# # lists
# 
# chl_test_list <- list(c(1:3), c(4:7))
# 
# # ______________________________________________________________________________


# load sst data
sst.nc <- brick("data/METOFFICE-GLO-SST-L4-REP-OBS-SST_1695911795678.nc", values = T,
                varname = 'analysed_sst')
# load ice data
ice.nc <- brick("data/METOFFICE-GLO-SST-L4-REP-OBS-SST_1695911795678.nc", values = T,
                varname = 'sea_ice_fraction')
# load chl data in loop, because each month is in own nc file 

sst_crs <- sst.nc@crs@projargs


dates <- as.Date(getZ(sst.nc))
mons <- format(dates,"%m")
doys <- c(1:length(dates))
# round(julian(dates, origin = as.Date("2018-01-01")))

# prepare data output
results <- tibble(`julian day` = NA,
                  `date` = NA, 
                  `mean temp warm [ºC]` = NA, 
                  `sd temp warm [ºC]` = NA, 
                  `area warm [km2]` = NA,
                  `mean temp front [ºC]` = NA, 
                  `sd temp front [ºC]` = NA,
                  `area front [km2]` = NA,
                  `mean temp cold [ºC]` = NA, 
                  `sd temp cold [ºC]` = NA, 
                  `area cold [km2]` = NA, 
                  `mean chl total [mg / m^3]` = NA,
                  `mean chl warm [mg / m^3]` = NA,
                  `min chl warm [mg / m^3]` = NA,
                  `max chl warm [mg / m^3]` = NA,
                  `n.obs. chl warm [pixels]` = NA,
                  `mean chl cold [mg / m^3]` = NA,
                  `min chl cold [mg / m^3]` = NA,
                  `max chl cold [mg / m^3]` = NA,
                  `n.obs chl cold [pixels]` = NA,
                  `mean chl front [mg / m^3]` = NA,
                  `min chl front [mg / m^3]` = NA,
                  `max chl front [mg / m^3]` = NA,
                  `n.obs. chl front [pixels]` = NA,
                  `mean PICO warm [mg / m^3]` = NA,
                  `mean NANO warm [mg / m^3]` = NA,
                  `mean MICRO warm [mg / m^3]` = NA,
                  `mean PICO cold [mg / m^3]` = NA,
                  `mean NANO cold [mg / m^3]` = NA,
                  `mean MICRO cold [mg / m^3]` = NA,
                  `mean PICO front [mg / m^3]` = NA,
                  `mean NANO front [mg / m^3]` = NA,
                  `mean MICRO front [mg / m^3]` = NA,
                  `ice cover warm region [%]` = NA, # könnte man wieder aktivieren, muss aber erst die Daten bekommen
                  `ice cover cold region [%]` = NA ,
                  `ice cover front region [%]` = NA,
                  `ice cover warm region sd [%]` = NA, 
                  `ice cover cold region sd [%]` = NA ,
                  `ice cover front region sd [%]` = NA
)[-1,]

# empty list to fill with spatial polygons (warm and cold clusters)
polys <- list()
chl_list <- list()

# loop over days
old_mon <- 0 # counter for months; needed for chl data loading

for(doy in doys){
  
  # get date
  date <- dates[doy]
  mon <- mons[doy]
  
  
  # extent of chl.nc and sst.nc are equal (beacuse I downloaded them like this)
  #extent(chl.nc[[doy]])
  
  sst_d <- sst.nc[[doy]]
  ice_d <- ice.nc[[doy]] 
  
  # convert Kelvin to Celsius
  values(sst_d) <- values(sst_d) - 273.13
  # plot(sst_d)
  
  
  # 1.: CLUSTERING of SST
  km <- kmeans(na.omit(as.vector(sst_d)), centers = 3)  # get vector of sst values, without NAs
  df <- data.frame(val = as.vector(sst_d),              # make data.frame to
                   clus = NA)
  df$clus[which(!is.na(df$val))] <- km$cluster          # align clusters to values that are not NA
  
  # identity of clusters is numeric, the larger group is 1 and the smaller 2. 
  # now, is cluster 1 = cold or 1 = warm?
  df$group <- NA
  

  
  # # this is for 2 groups:
  # warm1 <- mean(df$val[df$clus==1], na.rm = T) > mean(df$val[df$clus==2], na.rm = T) # if true, 1 is warmer than 2
  # df$group[df$clus==1] <- ifelse(warm1, "warm", "cold")
  # df$group[df$clus==2] <- ifelse(warm1, "cold", "warm")
  
  # df$group <- as.factor(df$group)
  # # put values into original Raster
  # clusters <- setValues(sst_d, as.matrix(df$clus, ncol = ncol(sst_d)))   
  
  
  # for 3 groups:
  df$group <- "NA"
  means <- c(mean(df$val[df$clus==1], na.rm = T), mean(df$val[df$clus==2], na.rm = T), mean(df$val[df$clus==3], na.rm = T))
  df$group[df$clus==which.max(means)] <- "warm"
  df$group[df$clus==which.min(means)] <- "cold"
  # which(means != min(means) & means != max(means)) which index of means is neither max nor min?
  df$group[df$clus==which(means != min(means) & means != max(means))] <- "front"
  
  df$group <- as.factor(df$group)
  # put values into original Raster
  clusters <- setValues(sst_d, as.matrix(df$clus, ncol = ncol(sst_d)))  
  #  plot(clusters) # look at magnificent plots
  
  res(clusters) <- 0.05 # https://stackoverflow.com/questions/56115830/why-rastertopolygons-is-creating-horizontal-lines
  
  # 2.: create polygon masks for both clusters
  
  # # create polygon masks for warm and cold sections, for 2 clusters
  # warmpoly <- rasterToPolygons(x= clusters, fun = function(x){x==ifelse(warm1, 1, 2)}, dissolve = T)
  # coldpoly <- rasterToPolygons(x= clusters, fun = function(x){x==ifelse(warm1, 2, 1)}, dissolve = T)
  
  # create polygon masks for warm and cold sections, for 3 clusters
  warmpoly <- rasterToPolygons(x = clusters, fun = function(x){x==which.max(means)}, dissolve = T)
  coldpoly <- rasterToPolygons(x = clusters, fun = function(x){x==which.min(means)}, dissolve = T)
  frontpoly <- rasterToPolygons(x = clusters, fun = function(x){x==which(means != min(means) & means != max(means))}, dissolve = T)
  
  # 3.: use these masks to assign chl a data to clusters
  
  # load correct chl nc, if not done already
  if(doy==1 | mon != old_mon){
    chl.nc <- brick(paste0('data/', mon, '_cmems_obs-oc_glo_bgc-plankton_my_l3-multi-4km_P1D.nc'), values=T, 
                              varname = "CHL")
    pico.nc <- brick(paste0('data/', mon, '_cmems_obs-oc_glo_bgc-plankton_my_l3-multi-4km_P1D.nc'), values=T, 
                    varname = "PICO")
    nano.nc <- brick(paste0('data/', mon, '_cmems_obs-oc_glo_bgc-plankton_my_l3-multi-4km_P1D.nc'), values=T, 
                     varname = "NANO")
    micro.nc <- brick(paste0('data/', mon, '_cmems_obs-oc_glo_bgc-plankton_my_l3-multi-4km_P1D.nc'), values=T, 
                     varname = "MICRO")
    
  }
  
  
  # find index of date in chl.nc
  chl_index <- which(getZ(chl.nc) == date)
  
  chl_d <- chl.nc[[chl_index]]
  pico_d <- pico.nc[[chl_index]]
  nano_d <- nano.nc[[chl_index]]
  micro_d <- micro.nc[[chl_index]]
  
  # plot(chl_d)
 
  # apply sst masks on chl data
  cold_chl <- mask(chl_d, coldpoly)
  cold_pico <- mask(pico_d, coldpoly)
  cold_nano <- mask(nano_d, coldpoly)
  cold_micro <- mask(micro_d, coldpoly)
  
  warm_chl <- mask(chl_d, warmpoly)
  warm_pico <- mask(pico_d, warmpoly)
  warm_nano <- mask(nano_d, warmpoly)
  warm_micro <- mask(micro_d, warmpoly)
  
  front_chl <- mask(chl_d, frontpoly)
  front_pico <- mask(pico_d, frontpoly)
  front_nano <- mask(nano_d, frontpoly)
  front_micro <- mask(micro_d, frontpoly)
  
  
  # 4. calculate statistics for overview table 
  # chl a mean concentration in whole area
  chl_total <- mean(getValues(chl_d), na.rm = T)
  # mean(na.omit(as.vector(chl_d))) 
  
  # mean chl warm
  chl_ex_warm <- getValues(warm_chl)
  mean_chl_warm <- mean(chl_ex_warm, na.rm = T)
  min_chl_warm <- min(chl_ex_warm, na.rm = T)
  max_chl_warm <- max(chl_ex_warm, na.rm = T)
  
  mean_pico_warm <- mean(getValues(warm_pico), na.rm = T)
  mean_nano_warm <- mean(getValues(warm_nano), na.rm = T)
  mean_micro_warm <- mean(getValues(warm_micro), na.rm = T)
  
  # mean chl cold
  chl_ex_cold <- getValues(cold_chl)
  mean_chl_cold <- mean(chl_ex_cold, na.rm = T)
  min_chl_cold <- min(chl_ex_cold, na.rm = T)
  max_chl_cold <- max(chl_ex_cold, na.rm = T)
  
  mean_pico_cold <- mean(getValues(cold_pico), na.rm = T)
  mean_nano_cold <- mean(getValues(cold_nano), na.rm = T)
  mean_micro_cold <- mean(getValues(cold_micro), na.rm = T)
  
  # mean chl front
  chl_ex_front <- getValues(front_chl)
  mean_chl_front <- mean(chl_ex_front, na.rm = T)
  min_chl_front <- min(chl_ex_front, na.rm = T)
  max_chl_front <- max(chl_ex_front, na.rm = T)
  
  mean_pico_front <- mean(getValues(front_pico), na.rm = T)
  mean_nano_front <- mean(getValues(front_nano), na.rm = T)
  mean_micro_front <- mean(getValues(front_micro), na.rm = T)
  
  # nobs warm
  nobs_chl_warm <- length(na.omit(chl_ex_warm))
  # nobs cold
  nobs_chl_cold <- length(na.omit(chl_ex_cold))
  # nobs front
  nobs_chl_front <- length(na.omit(chl_ex_front))
  
  
  
  
  # other, chl independent, variables of interest
  # area warm poly
  area_temp_warm <- area(warmpoly) / 1000000 #km2
  # mean temp warm
  mean_warm_temp <- mean(na.omit(as.vector(mask(sst_d, warmpoly)))) 
  # sd temp warm
  sd_warm_temp <- sd(na.omit(as.vector(mask(sst_d, warmpoly)))) 
  # area cold poly
  area_temp_cold <- area(coldpoly) / 1000000 #km2
  # mean temp cold
  mean_cold_temp <- mean(na.omit(as.vector(mask(sst_d, coldpoly)))) 
  # sd temp cold
  sd_cold_temp <- sd(na.omit(as.vector(mask(sst_d, coldpoly)))) 
  # area front poly
  area_temp_front <- area(frontpoly) / 1000000 #km2
  # mean temp front
  mean_front_temp <- mean(na.omit(as.vector(mask(sst_d, frontpoly)))) 
  # sd temp front
  sd_front_temp <- sd(na.omit(as.vector(mask(sst_d, frontpoly)))) 
  
  
  # e.g. ice
  mean_ice_warm <- mean(as.vector(mask(ice_d, warmpoly)), na.rm = T)
  mean_ice_cold <- mean(as.vector(mask(ice_d, coldpoly)), na.rm = T)
  mean_ice_front <- mean(as.vector(mask(ice_d, frontpoly)), na.rm = T)
  
  sd_ice_warm <- sd(as.vector(mask(ice_d, warmpoly)), na.rm = T)
  sd_ice_cold <- sd(as.vector(mask(ice_d, coldpoly)), na.rm = T)
  sd_ice_front <- sd(as.vector(mask(ice_d, frontpoly)), na.rm = T)
  # stations in each of the clusters
  
  
  
  # 5. save results
  ## save everything to dataframe results
  
  out <- tibble(`julian day` = doy,
                `date` = date, 
                `mean temp warm [ºC]` = mean_warm_temp,
                `sd temp warm [ºC]` = sd_warm_temp,
                `area warm [km2]` = area_temp_warm,
                `mean temp front [ºC]` = mean_front_temp, 
                `sd temp front [ºC]` = sd_front_temp, 
                `area front [km2]` = area_temp_front,
                `mean temp cold [ºC]` = mean_cold_temp, 
                `sd temp cold [ºC]` = sd_cold_temp, 
                `area cold [km2]` = area_temp_cold, 
                `mean chl total [mg / m^3]` = chl_total,
                `mean chl warm [mg / m^3]` = mean_chl_warm,
                `min chl warm [mg / m^3]` = min_chl_warm,
                `max chl warm [mg / m^3]` = max_chl_warm,
                `n.obs. chl warm [pixels]` = nobs_chl_warm,
                `mean chl cold [mg / m^3]` = mean_chl_cold,
                `min chl cold [mg / m^3]` = min_chl_cold,
                `max chl cold [mg / m^3]` = max_chl_cold,
                `n.obs chl cold [pixels]` = nobs_chl_cold,
                `mean chl front [mg / m^3]` = mean_chl_front,
                `min chl front [mg / m^3]` = min_chl_front,
                `max chl front [mg / m^3]` = max_chl_front,
                `n.obs. chl front [pixels]` = nobs_chl_front,
                `mean PICO warm [mg / m^3]` = mean_pico_warm,
                `mean NANO warm [mg / m^3]` = mean_nano_warm,
                `mean MICRO warm [mg / m^3]` = mean_micro_warm,
                `mean PICO cold [mg / m^3]` = mean_pico_cold,
                `mean NANO cold [mg / m^3]` = mean_nano_cold,
                `mean MICRO cold [mg / m^3]` = mean_micro_cold,
                `mean PICO front [mg / m^3]` = mean_pico_front,
                `mean NANO front [mg / m^3]` = mean_nano_front,
                `mean MICRO front [mg / m^3]` = mean_micro_front,
                `ice cover warm region [%]` = mean_ice_warm, # könnte man wieder aktivieren, muss aber erst die Daten bekommen
                `ice cover cold region [%]` = mean_ice_cold ,
                `ice cover front region [%]` = mean_ice_front,
                `ice cover warm region sd [%]` = sd_ice_warm, # könnte man wieder aktivieren, muss aber erst die Daten bekommen
                `ice cover cold region  sd [%]` = sd_ice_cold ,
                `ice cover front region sd [%]` = sd_ice_front
               # `ice cover front region [%]` = NA
  )
  
  # write to output
  # overview table
  results <- rbind(results, out)
  # cluster polygons
  polys[[paste0(doy)]] <- list(warm = warmpoly,
                               front = frontpoly,
                               cold = coldpoly)
  # Chl values by clusters
  chl_list[["CHL"]][["warm"]][[paste0(doy)]] <- getValues(warm_chl)[!is.na(getValues(warm_chl))]
  chl_list[["CHL"]][["front"]][[paste0(doy)]] <- getValues(front_chl)[!is.na(getValues(front_chl))]
  chl_list[["CHL"]][["cold"]][[paste0(doy)]] <- getValues(cold_chl)[!is.na(getValues(cold_chl))]
  chl_list[["PICO"]][["warm"]][[paste0(doy)]] <- getValues(warm_pico)[!is.na(getValues(warm_pico))]
  chl_list[["PICO"]][["front"]][[paste0(doy)]] <- getValues(front_pico)[!is.na(getValues(front_pico))]
  chl_list[["PICO"]][["cold"]][[paste0(doy)]] <- getValues(cold_pico)[!is.na(getValues(cold_pico))]
  chl_list[["NANO"]][["warm"]][[paste0(doy)]] <- getValues(warm_nano)[!is.na(getValues(warm_nano))]
  chl_list[["NANO"]][["front"]][[paste0(doy)]] <- getValues(front_nano)[!is.na(getValues(front_nano))]
  chl_list[["NANO"]][["cold"]][[paste0(doy)]] <- getValues(cold_nano)[!is.na(getValues(cold_nano))]
  chl_list[["MICRO"]][["warm"]][[paste0(doy)]] <- getValues(warm_micro)[!is.na(getValues(warm_micro))]
  chl_list[["MICRO"]][["front"]][[paste0(doy)]] <- getValues(front_micro)[!is.na(getValues(front_micro))]
  chl_list[["MICRO"]][["cold"]][[paste0(doy)]] <- getValues(cold_micro)[!is.na(getValues(cold_micro))]
  
  print(doy)
  mon_old = mon
}

file_appx = "_2nd_run_no_lines"
# save output
write_csv2(x = results %>% arrange(`date`) , file = paste0("output/Output_3k_cluster_analysis_thesis", file_appx, ".csv"))
saveRDS(object = polys, file = paste0("output/cluster_3k_2018", file_appx, ".rds"))
saveRDS(object = chl_list, file = paste0("output/3k_raw_chl_conc", file_appx, ".rds"))



# expand to pft groups

# expand to 3 clusters



