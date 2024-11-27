## what is causing variability in Chl?
# try some correlation analysis

library(tidyverse)

# results_chl_var <- read_csv2("output/Chl_timeseries_variability_3k.csv") # cant find file
results_chl_var <- read_csv2("output/Chl_no_trans_timeseries_variability_3k.csv")

ts_dat <- read_csv2("output/Output_3k_cluster_analysis_thesis.csv", na = c("", "Inf", "-Inf","NA"))



# gather ts dat 
temp_dat <- ts_dat %>%
  dplyr::select(date, `mean temp warm [ºC]`, `mean temp cold [ºC]`, `mean temp front [ºC]`) %>%
  pivot_longer(-date, names_prefix = "mean temp ", names_to = "cluster", values_to = "temp [ºC]")%>%
  mutate(cluster = str_extract(cluster, pattern = "^(warm|cold|front)")) # remove awkward name appendix

temp_sd_dat <- ts_dat %>%
  dplyr::select(date, `sd temp warm [ºC]`, `sd temp cold [ºC]`, `sd temp front [ºC]`) %>%
  pivot_longer(-date, names_prefix = "sd temp ", names_to = "cluster", values_to = "sd temp [ºC]")%>%
  mutate(cluster = str_extract(cluster, pattern = "^(warm|cold|front)")) # remove awkward name appendix

ice_dat <- ts_dat %>%
  dplyr::select(date, `ice cover warm region [%]`, `ice cover cold region [%]`, `ice cover front region [%]`) %>%
  pivot_longer(-date, names_prefix = "ice cover ", names_to = "cluster", values_to = "ice cover [%]")  %>%
  mutate(cluster = str_extract(cluster, pattern = "^(warm|cold|front)")) # remove awkward name appendix

ice_sd_dat <- ts_dat %>%
  dplyr::select(date, `ice cover warm region sd [%]`, `ice cover cold region  sd [%]`, `ice cover front region sd [%]`) %>%
  pivot_longer(-date, names_prefix = "ice cover ", names_to = "cluster", values_to = "ice cover sd [%]")  %>%
  mutate(cluster = str_extract(cluster, pattern = "^(warm|cold|front)")) # remove awkward name appendix



# quantify variability in Chl
chl_var <- results_chl_var %>%
  mutate(IQR = `IQR_u`-`IQR_l`,
         CI_r = CI_u - CI_l) %>%
  dplyr::select(date:max, SD:CI_r) %>%
  left_join(temp_dat, by = c('date', 'cluster')) %>%
  left_join(temp_sd_dat, by = c('date', 'cluster')) %>%
  left_join(ice_dat, by = c('date', 'cluster')) %>%
  left_join(ice_sd_dat, by = c('date', 'cluster')) %>%
  drop_na() #%>%
#  group_by(cluster)


test <- chl_var %>%
  dplyr::select(-cluster, -date) 
# rename columns
nums <- c(1:length(test))
cnames <-  paste0("X", nums, "_",colnames(test))
colnames(test) <- cnames


test %>%
  as.matrix %>%
  cor %>% 
  ## set redundant to `NA`
  `[<-`(lower.tri(., TRUE), NA) %>%
  ## back to tibble
  as_tibble(rownames="var1") %>%
  ## long format, dropping redundant  
  pivot_longer(cols=-1, names_to="var2", values_to="rho", values_drop_na=TRUE) %>%
  mutate(var1 = factor(var1, levels = cnames),
         var2 = factor(var2, levels = cnames)) %>%
  
  ## descending sort most correlated pairs
  #arrange(rho) %>%
  ggplot(aes(x=var1, y=var2, fill=rho))+
    geom_tile() + 
    geom_text(aes(label=round(rho, 2))) +
    scale_fill_distiller(palette = "RdBu") 
  

