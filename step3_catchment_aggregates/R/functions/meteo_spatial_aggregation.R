### aggregation of gridded meteorological data products within catchment representations to derive meteorological time series

#-------------------------------------------------------------------------------

# general function to compute catchment aggregates from gridded data products
meteo_catch_mean <- function(catch_repr_i, meteo_prod, var, years) {
  
  # define file path to read processed geospatial data products
  path_proc_data_prod <- paste0(path_catch_agg_temp, "proc_data_prod/")

  #-----------------------------------------------------------------------------
  
  catch_mean_long_list <- list() # empty list for joining all variables
  for (v in 1:length(var)){
    
    var_v <- var[v]
    print(var_v)
    var_v_catch_mean <- data.frame()
    
    # conditional time period, depending on MeteoSwiss variables
    if (var_v %in% c("RhiresD", "TabsD")){
      years <- 1961:2024
    } else if (var_v %in% c("TminD", "TmaxD")){
      years <- 1971:2024
    }
    
    for (y in 1:length(years)){
      
      year_y <- years[y]
      print(year_y)
      var_v_y_catch_mean <- data.frame()
      
      # read yearly raster files
      var_v_y_bbox <- rast(paste0(path_proc_data_prod, meteo_prod, "/", var_v, "/", meteo_prod, "_",
                                  var_v, "_raster_bbox_", year_y, ".nc"))
      # compute mean over exact overlapping grid cell portions and write into long format
      var_v_y_catch_mean <- exactextractr::exact_extract(var_v_y_bbox, catch_repr_i, "mean") %>%
        pivot_longer(cols = everything(), names_to = "date", values_to = var_v)
      
      # for ERA5-Land data assign hourly dates, otherwise daily dates
      if (meteo_prod == "era5land") {
        
        # generate time vector for different raster layers and overwrite in data frame, time will get shifted by one hour when assigned
        if (year_y == 1950){
          start_date <- as.POSIXct(paste0(year_y, "-01-01 01:00:00", tz = "UTC")) # for 1950-01-01 no data for 00:00 exist
        } else{
          start_date <- as.POSIXct(paste0(year_y, "-01-01 00:00:00", tz = "UTC"))
        }
        var_v_y_catch_mean$date <- seq(start_date, by = "1 hour", length.out = nrow(var_v_y_catch_mean))
        
      } else {
        
        # extract dates and overwrite in data frame
        var_v_y_catch_mean$date <- as.Date(time(var_v_y_bbox))
        
      }
      
      # bind results of all years for variable
      var_v_catch_mean <- rbind(var_v_catch_mean, var_v_y_catch_mean)
    }
    
    catch_mean_long_list[[v]] <- var_v_catch_mean
    
  }
  
  # bind data of all variables from list
  catch_mean_full <- reduce(catch_mean_long_list, full_join, by = c("date"))
  
  return(catch_mean_full)
  
}


### main function to perform spatial aggregation of meteorological products
meteo_catch_agg <- function(AlpAKaS_ID_i, catch_repr_i) {
  
  country_code_i <- catch_repr_i$country_code
  # exclude certain Austrian springs which catchments are better covered by MeteoSwiss data
  if (catch_repr_i$spring_name %in% c("Fidelisquelle", "Obwaldquelle")) {
    country_code_i <- "CH"
  }

  ### apply function for data aggregation and select and rename columns
  
  
  ### national datasets ###
  
  ### SPARTACUS ###
  if (country_code_i == "AT") {
    
    print("SPARTACUS")
    catch_mean_nat <- meteo_catch_mean(catch_repr_i = catch_repr_i, meteo_prod = "spartacus",
                                       var = c("RR", "TN", "TX"), years = 1961:2024)
    catch_mean_nat <- catch_mean_nat %>%
      mutate(temperature_mean_nat = NA) %>%
      dplyr::select(date, "RR", "TN", "temperature_mean_nat",  "TX") %>%
      rename(
        precipitation_nat = RR,
        temperature_min_nat = TN,
        temperature_max_nat = TX
      )
  
  ### MeteoSwiss ###
  } else if (country_code_i %in% c("CH", "LI")) {
    
    print("MeteoSwiss")
    catch_mean_nat <- meteo_catch_mean(catch_repr_i = catch_repr_i, meteo_prod = "meteoswiss",
                                       var = c("RhiresD", "TabsD", "TminD", "TmaxD"), years = NULL)
    catch_mean_nat <- catch_mean_nat %>%
      dplyr::select(date, "RhiresD", "TminD", "TabsD", "TmaxD") %>%
      rename(
        precipitation_nat = RhiresD,
        temperature_min_nat = TminD,
        temperature_mean_nat = TabsD,
        temperature_max_nat = TmaxD
      )

  ### HYRAS ###
  } else if (country_code_i == "DE") {
    
    print("HYRAS")
    catch_mean_nat <- meteo_catch_mean(catch_repr_i = catch_repr_i, meteo_prod = "hyras",
                                       var = c("pr", "tas", "tasmin", "tasmax"), years = 1951:2024)
    catch_mean_nat <- catch_mean_nat %>%
      dplyr::select(date, pr, tasmin, tas, tasmax) %>%
      rename(
        precipitation_nat = pr,
        temperature_min_nat = tasmin,
        temperature_mean_nat = tas,
        temperature_max_nat = tasmax
      )
    
  ### SAFRAN ###
  } else if (country_code_i == "FR") {
    
    print("SAFRAN")
    catch_mean_nat <- meteo_catch_mean(catch_repr_i = catch_repr_i, meteo_prod = "safran",
                                       var = c("PRENEI", "PRELIQ", "T", "TINF_H", "TSUP_H"), years = 1958:2024)
    # sum liquid and solid precipitation
    catch_mean_nat <- catch_mean_nat %>%
      mutate(precipitation_nat = PRENEI + PRELIQ) %>%
      dplyr::select(date, precipitation_nat, TINF_H, T, TSUP_H) %>%
      rename(
        temperature_min_nat = TINF_H,
        temperature_mean_nat = T,
        temperature_max_nat = TSUP_H
      )
    
  ### SLOCLIM ###
  } else if (country_code_i == "SI") {
    
    print("SLOCLIM")
    catch_mean_nat <- meteo_catch_mean(catch_repr_i = catch_repr_i, meteo_prod = "sloclim",
                                       var = c("pcp", "tmax_h", "tmin_h"), years = 1950:2018)
    catch_mean_nat <- catch_mean_nat %>%
      mutate(temperature_mean_nat = NA) %>%
      dplyr::select(date, "pcp", "tmin_h", "temperature_mean_nat", "tmax_h") %>%
      rename(
        precipitation_nat = pcp,
        temperature_min_nat = tmin_h,
        temperature_max_nat = tmax_h
      )
  
  } else {
    
    # empty columns
    catch_mean_nat <- data.frame(date = as.Date(character()),
                                 precipitation_nat = as.numeric(),
                                 temperature_min_nat = as.numeric(),
                                 temperature_mean_nat = as.numeric(),
                                 temperature_max_nat = as.numeric())
    
  }
  
  
  ### EOBS ###
  print("E-OBS")
  catch_mean_eobs <- meteo_catch_mean(catch_repr_i = catch_repr_i, meteo_prod = "eobs",
                                      var = c("rr", "tg", "tn", "tx"), years = 1950:2024)
  catch_mean_eobs <- catch_mean_eobs %>%
    dplyr::select(date, "rr", "tn", "tg", "tx") %>%
    rename(
      precipitation_eobs = rr,
      temperature_min_eobs = tn,
      temperature_mean_eobs = tg,
      temperature_max_eobs = tx
    )
  
  
  ### ERA5-Land ###
  print("ERA5-Land")
  catch_mean_era5land <- meteo_catch_mean(catch_repr_i = catch_repr_i, meteo_prod = "era5land",
                                          var = c("2m_temperature","total_precipitation"), years = 1950:2025)
  catch_mean_era5land <- catch_mean_era5land %>%
    dplyr::select(date, "total_precipitation", "2m_temperature") %>%
    rename(
      precipitation_era5land = total_precipitation,
      temperature_mean_era5land = `2m_temperature`
    )
  
  ### calculate daily aggregates for ERA5-Land
  # precipitation data is accumulated within each day in ERA5land, thus, the precipitation at 00UTC represents the summed precipitation fallen on the previous day
  # to obtain the precipitation fallen within the previous hour, the precipitation data of the previous hour has to be subtracted for all time steps != 01UTC
  # https://confluence.ecmwf.int/pages/viewpage.action?pageId=197702790
  
  # daily precipitation
  catch_mean_era5land_daily_prec <- catch_mean_era5land %>%
    dplyr::select(date, precipitation_era5land) %>%
    filter(hour(date) == 0) %>%
    mutate(date = as.Date(date))
  
  # calculate daily min, mean and max temperature
  catch_mean_era5land_daily_temp <- catch_mean_era5land %>%
    mutate(date = as.Date(date)) %>%
    group_by(date) %>%
    summarise(
      temperature_min_era5land = min(temperature_mean_era5land),
      temperature_max_era5land = max(temperature_mean_era5land),
      temperature_mean_era5land = mean(temperature_mean_era5land),
      .groups = "drop"
    ) %>% 
    dplyr::select(date, temperature_min_era5land, temperature_mean_era5land, temperature_max_era5land)
  
  catch_mean_era5land_daily <- catch_mean_era5land_daily_prec %>%
    full_join(catch_mean_era5land_daily_temp, by = c("date"))
  
  # join daily global and national meteorological data
  catch_mean_daily <- catch_mean_eobs %>%
    full_join(catch_mean_era5land_daily, by = c("date")) %>%
    full_join(catch_mean_nat, by = c("date")) %>%
    # cap time series at end of 2024
    filter(date <= as.Date("2024-12-31")) %>%
    mutate(across(where(is.numeric), ~ round(.x, 3)))
  
  # convert cumulative precipitation to hourly increments
  # generally, ERA5land also uses the end-of period convention for the hourly values
  # (precipitation at UTC+01 represents the precipitation fallen between UTC00 and UTC+01 after de-accumulating)
  catch_mean_hourly <- catch_mean_era5land %>%
    arrange(date) %>%
    mutate(
      tp_diff = precipitation_era5land - lag(precipitation_era5land),
      precipitation_era5land = ifelse(hour(date) == 1, precipitation_era5land, tp_diff),
      precipitation_era5land = ifelse(precipitation_era5land < 0, 0, precipitation_era5land) # very small negative values can occur which are set to 0
    ) %>%
    dplyr::select(-tp_diff) %>%
    # cap time series at end of 2024
    filter(as.Date(date) <= as.Date("2024-12-31")) %>%
    # make sure all dates are in correct format
    mutate(date = format(date, "%Y-%m-%d %H:%M:%S"),
           across(where(is.numeric), ~ round(.x, 3)))
  
  return(list(catch_mean_daily = catch_mean_daily, catch_mean_hourly = catch_mean_hourly))
  
}