# get time series of chl over longer time:

# get data 
# Chl: OC CCI, level4 v5, 8-day mean
# https://rsg.pml.ac.uk/thredds/ncss/CCI_ALL-v5.0-8DAY?var=chlor_a&north=82&west=-25&east=17&south=75&disableProjSubset=on&horizStride=1&time_start=2016-01-01T00%3A00%3A00Z&time_end=2021-12-31T00%3A00%3A00Z&timeStride=1&accept=netcdf

# SST and Ice
# python -m motuclient --motu https://nrt.cmems-du.eu/motu-web/Motu --service-id SST_GLO_SST_L4_NRT_OBSERVATIONS_010_001-TDS --product-id METOFFICE-GLO-SST-L4-NRT-OBS-SST-V2 --longitude-min -25 --longitude-max 17 --latitude-min 75 --latitude-max 82 --date-min "2016-01-01 00:00:00" --date-max "2021-12-31 23:59:59" --variable analysed_sst --variable sea_ice_fraction --out-dir <OUTPUT_DIRECTORY> --out-name <OUTPUT_FILENAME> --user <USERNAME> --pwd <PASSWORD>



#setwd("Users/vlampe/Documents/Last_chapter_clustering_manuskript_skripts/skripts/")
#setwd("Documents/Last_chapter_clustering_manuskript_skripts/skripts/")
library(raster)
library(tidyverse) 




# load sst data
sst.nc <- brick("data/METOFFICE-GLO-SST-L4-NRT-OBS-SST-V2_1698858096055_2016-2021.nc", values = T,
                varname = 'analysed_sst')
# load ice data
ice.nc <- brick("data/METOFFICE-GLO-SST-L4-NRT-OBS-SST-V2_1698858096055_2016-2021.nc", values = T,
                varname = 'sea_ice_fraction')
# load chl data
chl.nc <- brick("data/CCI_ALL-v5.0-8DAY_2016-2021.nc", values = T,
                varname = "chlor_a")




# make a Julian Day - Date list (for look-ups of chl a indices)
# in case sst and chl file do not line up
dates.sst <- data.frame(sst_index = 1:length(getZ(sst.nc)),
                        Date = as.Date(getZ(sst.nc)), 
                        julianday = julian(x = as.Date(getZ(sst.nc)), 
                                           origin = as.Date("2016-01-01")-1))
dates.chl <- data.frame(chl_index = 1:length(getZ(chl.nc)),
                        Date = as.Date(getZ(chl.nc)), 
                        julianday = julian(x = as.Date(getZ(chl.nc)), 
                                           origin = as.Date("2016-01-01")-1))

dates <- left_join(dates.sst, dates.chl, by = c("Date", "julianday")) %>%
  fill(chl_index, .direction = "down")
# %>%  distinct(sst_index, .keep_all = T)

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
                  `ice cover warm region [%]` = NA, 
                  `ice cover cold region [%]` = NA ,
                  `ice cover front region [%]` = NA,
                  `ice cover warm region sd [%]` = NA, 
                  `ice cover cold region sd [%]` = NA ,
                  `ice cover front region sd [%]` = NA
)[-1,]

# empty list to fill with spatial polygons (warm and cold clusters)
polys <- list()
chl_list <- list()



for(d in 1:nrow(dates)){ 
 
  # get date
  date <- dates$Date[d]
  
  # get data for today
  sst_d <- sst.nc[[d]]
  # convert Kelvin to Celsius
  values(sst_d) <- values(sst_d) - 273.13
  
  ice_d <- ice.nc[[d]] 
  
  chl_d <- chl.nc[[dates$chl_index[d]]]
  
  # 1.: CLUSTERING of SST
  km <- kmeans(na.omit(as.vector(sst_d)), centers = 3)  # get vector of sst values, without NAs
  df <- data.frame(val = as.vector(sst_d),              # make data.frame to
                   clus = NA)
  df$clus[which(!is.na(df$val))] <- km$cluster          # align clusters to values that are not NA
  
  # identity of clusters is numeric, the larger group is 1 and the smaller 2. 
  # now, is cluster 1 = cold or 1 = warm?
  df$group <- NA
  
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
  
  
  # 2.: create polygon masks for both clusters
  
  # create polygon masks for warm and cold sections, for 3 clusters
  warmpoly <- rasterToPolygons(x = clusters, fun = function(x){x==which.max(means)}, dissolve = T)
  coldpoly <- rasterToPolygons(x = clusters, fun = function(x){x==which.min(means)}, dissolve = T)
  frontpoly <- rasterToPolygons(x = clusters, fun = function(x){x==which(means != min(means) & means != max(means))}, dissolve = T)
  
  # 3.: use these masks to assign chl a data to clusters
  
  # apply sst masks on chl data
  cold_chl <- mask(chl_d, coldpoly)
  warm_chl <- mask(chl_d, warmpoly)
  front_chl <- mask(chl_d, frontpoly)
  
  # 4. calculate statistics for overview table 
  # chl a mean concentration in whole area
  chl_total <- mean(getValues(chl_d), na.rm = T)
  # mean(na.omit(as.vector(chl_d))) 
  
  # mean chl warm
  chl_ex_warm <- getValues(warm_chl)
  mean_chl_warm <- mean(chl_ex_warm, na.rm = T)
  min_chl_warm <- min(chl_ex_warm, na.rm = T)
  max_chl_warm <- max(chl_ex_warm, na.rm = T)
  

  # mean chl cold
  chl_ex_cold <- getValues(cold_chl)
  mean_chl_cold <- mean(chl_ex_cold, na.rm = T)
  min_chl_cold <- min(chl_ex_cold, na.rm = T)
  max_chl_cold <- max(chl_ex_cold, na.rm = T)
  
  # mean chl front
  chl_ex_front <- getValues(front_chl)
  mean_chl_front <- mean(chl_ex_front, na.rm = T)
  min_chl_front <- min(chl_ex_front, na.rm = T)
  max_chl_front <- max(chl_ex_front, na.rm = T)
  
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
  
  out <- tibble(`julian day` = d,
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
  polys[[paste0(d)]] <- list(warm = warmpoly,
                               front = frontpoly,
                               cold = coldpoly)
  # Chl values by clusters
  chl_list[["CHL"]][["warm"]][[paste0(d)]] <- getValues(warm_chl)[!is.na(getValues(warm_chl))]
  chl_list[["CHL"]][["front"]][[paste0(d)]] <- getValues(front_chl)[!is.na(getValues(front_chl))]
  chl_list[["CHL"]][["cold"]][[paste0(d)]] <- getValues(cold_chl)[!is.na(getValues(cold_chl))]

  
  print(d)

}

# save output
write_csv2(x = results %>% arrange(`date`) , file = "output/Output_3k_cluster_analysis_2016-2021.csv")
saveRDS(object = polys, file = "output/cluster_3k_2016-2021.rds")
saveRDS(object = chl_list, file = "output/3k_raw_chl_conc_2016-2021.rds")







