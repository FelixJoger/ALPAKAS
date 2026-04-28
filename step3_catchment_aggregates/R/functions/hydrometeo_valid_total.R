### determination of discharge and meteorological observations within valid hydrological years
# defined by availability of both discharge and meteorological data for at least 80% of days within that hydrological year

#-------------------------------------------------------------------------------

# function for determination of valid hydrometeorological observations
compute_hydrometeo_valid_total <- function(meteo_prod_class){
  
  print(meteo_prod_class)
  
  # read valid discharge data
  discharge_valid_total <- readRDS(file = "step2_buffer_approximations/temp/discharge_valid_total.RDS")
  dir.create(paste0(path_catch_agg_temp, catch_del_name_r, "/"))
  
  meteo_list <- list()
  for (i in 1:nrow(catch_repr)) {
    
    ### read daily meteorological aggregates and bind to data frame
    meteo_list[[i]] <- read.csv(file = paste0(path_catch_agg_out, catch_del_name_r,
                                              "/meteorological_time_series/daily/AlpAKaS_meteo_daily_", catch_repr$AlpAKaS_ID[i], ".csv")) %>%
      mutate(AlpAKaS_ID = catch_repr$AlpAKaS_ID[i],
             # for national data products, Tmean is approximated by (Tmin + Tmax) / 2, wherever missing
             temperature_mean_nat = coalesce(temperature_mean_nat, (temperature_min_nat + temperature_max_nat) / 2),
             # select respective product class
             precipitation = .data[[paste0("precipitation_", meteo_prod_class)]],
             temperature_min = .data[[paste0("temperature_min_", meteo_prod_class)]],
             temperature_mean = .data[[paste0("temperature_mean_", meteo_prod_class)]],
             temperature_max = .data[[paste0("temperature_max_", meteo_prod_class)]]) %>%
      relocate(AlpAKaS_ID, .before = 1) %>%
      dplyr::select(AlpAKaS_ID, date, precipitation, temperature_min, temperature_mean, temperature_max)
  }
  meteo_total <- bind_rows(meteo_list) %>%
    mutate(date = as.Date(date),
           # add hydrological years
           month = month(date),
           hydro_year = if_else(month >= 10, year(date) + 1, year(date)))
  
  # join observations of valid hydrological years of discharge and meteorological data and determine all observations laying within valid years
  hydrometeo_valid_hydro_year <- readRDS(file = paste0("step2_buffer_approximations/temp/hydrometeo_valid_hydro_years_",
                                                       meteo_prod_class, ".RDS")) %>%
    filter(AlpAKaS_ID %in% unique(meteo_total$AlpAKaS_ID))
  hydrometeo_valid_total <- discharge_valid_total %>%
    full_join(meteo_total, by = c("AlpAKaS_ID", "date", "month", "hydro_year")) %>%
    inner_join(hydrometeo_valid_hydro_year %>% dplyr::select(AlpAKaS_ID, hydro_year), by = c("AlpAKaS_ID", "hydro_year"))
  
  saveRDS(hydrometeo_valid_total, file = paste0(path_catch_agg_temp, catch_del_name_r, "/hydrometeo_valid_total_", meteo_prod_class, ".RDS"))
}


### main function to determine valid observation of hydrometeorological data for each meteorological product class
hydrometeo_valid <- function(catch_del_name_r, catch_repr){
  
  # apply function to all meteorological product classes
  compute_hydrometeo_valid_total(meteo_prod_class = "eobs")
  compute_hydrometeo_valid_total(meteo_prod_class = "era5land")
  compute_hydrometeo_valid_total(meteo_prod_class = "nat")

}