### determination of valid hydrological years
# defined by availability of both discharge and meteorological data for at least 80% of days within that hydrological year

#-------------------------------------------------------------------------------

# function for computation of valid hydrological years
compute_hydrometeo_valid <- function(meteo_prod_class, 
                                     temp_avail_meteo_prod,
                                     discharge_valid_hydro_year,
                                     discharge_valid_total){
  
  print(meteo_prod_class)
  
  ### compute valid hydrological years for meteorological data: dependent on data products
  if (meteo_prod_class %in% c("eobs", "era5land")){
    
    temp_avail_meteo_prod_p <- temp_avail_meteo_prod %>%
      filter(prod == meteo_prod_class)
    # create date sequence per site
    meteo_df <- station_meta %>%
      mutate(date = map2(as.Date(temp_avail_meteo_prod_p$start_date), as.Date(temp_avail_meteo_prod_p$end_date), seq, by = "day")) %>%
      unnest(date) %>%
      transmute(AlpAKaS_ID, date, available = TRUE)
    
  } else if (meteo_prod_class == "nat"){
    
    meteo_df <- station_meta %>%
      # join start and end dates based on country_code
      left_join(temp_avail_meteo_prod %>% dplyr::select(country_code, start_date, end_date), by = "country_code") %>%
      filter(!is.na(start_date)) %>%
      # create date sequence per site
      mutate(date = map2(as.Date(start_date), as.Date(end_date), seq, by = "day")) %>%
      unnest(date) %>%
      transmute(AlpAKaS_ID, date, available = TRUE)
  }
  
  ### fill missing values in first and last hydrological years of meteorological time series with NA values for completeness analysis
  meteo_filled_list <- list()
  for (i in 1:n_distinct(meteo_df$AlpAKaS_ID)) {
    
    AlpAKaS_ID_i <- unique(meteo_df$AlpAKaS_ID)[i]
    print(paste(i, AlpAKaS_ID_i, sep = ": "))
    
    meteo_df_i <- meteo_df %>%
      filter(AlpAKaS_ID == AlpAKaS_ID_i)
    
    # extract start and end dates of first and last hydrological years
    first_year <- year(min(meteo_df_i$date))
    last_year  <- year(max(meteo_df_i$date))
    first_month <- month(min(meteo_df_i$date))
    last_month  <- month(max(meteo_df_i$date))
    if (first_month < 10){
      date_first <- as.Date(paste0(first_year - 1, "-10-01"))
    } else {
      date_first <- as.Date(paste0(first_year, "-10-01"))
    }
    if (last_month < 10){
      date_last <- as.Date(paste0(last_year, "-09-30"))
    } else {
      date_last <- as.Date(paste0(last_year + 1, "-09-30"))
    }
    
    # extract missing dates within the dataset
    full_dates <- tibble(date = seq(date_first, date_last, by = "day"))
    fill_dates <- full_dates %>% anti_join(meteo_df_i, by = "date") %>%
      mutate(AlpAKaS_ID = AlpAKaS_ID_i)
    
    # fill missing dates and save filled time series in list
    meteo_df_i_filled <- bind_rows(meteo_df_i, fill_dates) %>%
      arrange(date)
    meteo_filled_list[[i]] <- meteo_df_i_filled
  }
  
  # bind files and add information about hydrological year and seasons
  meteo_filled_total <- do.call(rbind, meteo_filled_list) %>%
    mutate(
      month = month(date),
      hydro_year = if_else(month >= 10, year(date) + 1, year(date)),
    )
  
  # determine ratio of data completeness per hydrological year and only keep years with ≥80% observations
  meteo_filled_hydro_year <- meteo_filled_total %>%
    group_by(AlpAKaS_ID, hydro_year) %>%
    summarise(
      n_days = n(),
      n_valid = sum(!is.na(available)),
      .groups = "drop") %>%
    mutate(valid_fraction = n_valid / n_days) 
  
  meteo_valid_hydro_year <- meteo_filled_hydro_year %>%
    filter(valid_fraction >= 0.8)
  
  # determine period and number of valid years for each site 
  meteo_valid_summary <- meteo_valid_hydro_year %>%
    group_by(AlpAKaS_ID) %>%
    summarise(
      start_year_valid = as.Date(ifelse(n() > 0, paste0(min(hydro_year) - 1, "-10-01"), NA)),
      end_year_valid = as.Date(ifelse(n() > 0, paste0(max(hydro_year), "-09-30"), NA)),
      n_years_valid = n(),
      .groups = "drop")
  
  saveRDS(file = paste0(path_buff_approx_temp, "/meteo_valid_hydro_year_", meteo_prod_class, ".RDS"), meteo_valid_hydro_year)
  saveRDS(file = paste0(path_buff_approx_temp, "/meteo_valid_summary_", meteo_prod_class, ".RDS"), meteo_valid_summary)
  
  
  #-----------------------------------------------------------------------------
  ### determine valid hydrological years based the data availability of both discharge and meteorological time series
  
  # join valid hydrological years of discharge and meteorological data
  hydrometeo_valid_hydro_year <- discharge_valid_hydro_year %>%
    inner_join(meteo_valid_hydro_year, by = c("AlpAKaS_ID", "hydro_year"))
  
  # determine period and number of valid years for each site used for calculation of hydrometeorological indices
  hydrometeo_valid_summary <- hydrometeo_valid_hydro_year %>%
    group_by(AlpAKaS_ID) %>%
    summarise(
      ind_start_date = as.Date(paste0(min(hydro_year) - 1, "-10-01")),
      ind_end_date = as.Date(paste0(max(hydro_year), "-09-30")),
      ind_years_valid = n(),
      .groups = "drop")
  
  # determine all discharge observations laying within valid years and the mean discharge for each site (required for buffer method)
  hydrometeo_valid_total_discharge <- discharge_valid_total %>%
    inner_join(hydrometeo_valid_hydro_year %>% dplyr::select(AlpAKaS_ID, hydro_year), by = c("AlpAKaS_ID", "hydro_year"))
  hydrometeo_valid_mean_discharge <- hydrometeo_valid_total_discharge %>%
    group_by(AlpAKaS_ID) %>%
    summarise(q_mean = mean(discharge, na.rm = TRUE),
              .groups = "drop")
  
  saveRDS(file = paste0(path_buff_approx_temp, "/hydrometeo_valid_hydro_years_", meteo_prod_class, ".RDS"), hydrometeo_valid_hydro_year)
  saveRDS(file = paste0(path_buff_approx_temp, "/hydrometeo_valid_summary_", meteo_prod_class, ".RDS"), hydrometeo_valid_summary)
  saveRDS(file = paste0(path_buff_approx_temp, "/hydrometeo_valid_total_discharge_",
                        meteo_prod_class, ".RDS"), hydrometeo_valid_total_discharge)
  saveRDS(file = paste0(path_buff_approx_temp, "/hydrometeo_valid_mean_discharge_",
                        meteo_prod_class, ".RDS"), hydrometeo_valid_mean_discharge)
  
}


### main function for determination of valid hydrological years
hydrometeo_valid_hy <- function(station_meta){
  
  ### compute valid hydrological years for discharge data: only dependent on hydragraphs
  
  # fill missing values in first and last hydrological years of discharge time series with NA values for completeness analysis
  discharge_filled_list <- list()
  for (i in 1:nrow(station_meta)) {
  
    AlpAKaS_ID_i <- station_meta$AlpAKaS_ID[i]
    print(paste(i, AlpAKaS_ID_i, sep = ": "))
    
    file_res_i <- read.csv(paste0("step1_discharge_time_series/output_data/daily/AlpAKaS_discharge_daily_", AlpAKaS_ID_i,".csv")) %>%
      mutate(date = ymd(date)) %>%
      # remove artefacts, outliers, segments not identified as main segments in changepoint detection (keep NA values)
      filter(qc_flag != TRUE, is.na(cpd_segment_main) | cpd_segment_main != FALSE) %>%
      dplyr::select(date, discharge)
    
    # extract start and end dates of first and last hydrological years
    first_year <- year(min(file_res_i$date))
    last_year  <- year(max(file_res_i$date))
    first_month <- month(min(file_res_i$date))
    last_month  <- month(max(file_res_i$date))
    if (first_month < 10){
      date_first <- as.Date(paste0(first_year - 1, "-10-01"))
    } else {
      date_first <- as.Date(paste0(first_year, "-10-01"))
    }
    if (last_month < 10){
      date_last <- as.Date(paste0(last_year, "-09-30"))
    } else {
      date_last <- as.Date(paste0(last_year + 1, "-09-30"))
    }
    
    # extract missing dates within the dataset
    full_dates <- tibble(date = seq(date_first, date_last, by = "day"))
    fill_dates <- full_dates %>% anti_join(file_res_i, by = "date")
    
    # fill missing dates and save filled time series in list
    file_res_i_filled <- bind_rows(file_res_i, fill_dates) %>%
      arrange(date) %>%
      mutate(AlpAKaS_ID = AlpAKaS_ID_i)
    discharge_filled_list[[i]] <- file_res_i_filled
  }
    
  # bind files and add information about hydrological year and seasons
  discharge_filled_total <- do.call(rbind, discharge_filled_list) %>%
    mutate(
      month = month(date),
      hydro_year = if_else(month >= 10, year(date) + 1, year(date)),
      season = case_when(
        month %in% c(12, 1, 2) ~ "djf",
        month %in% c(3, 4, 5) ~ "mam",
        month %in% c(6, 7, 8) ~ "jja",
        month %in% c(9, 10, 11) ~ "son"
      ))
  
  # determine ratio of data completeness per hydrological year and only keep years with ≥80 % observations
  discharge_filled_hydro_year <- discharge_filled_total %>%
    group_by(AlpAKaS_ID, hydro_year) %>%
    summarise(
      n_days = n(),
      n_valid = sum(!is.na(discharge)),
      mean_discharge = mean(discharge, na.rm = TRUE),
      AlpAKaS_ID = unique(AlpAKaS_ID),
      .groups = "drop") %>%
    mutate(valid_fraction = n_valid / n_days)
  
  discharge_valid_hydro_year <- discharge_filled_hydro_year %>%
    filter(valid_fraction >= 0.8)
  
  # determine period and number of valid years for each site 
  discharge_valid_summary <- discharge_valid_hydro_year %>%
    group_by(AlpAKaS_ID) %>%
    summarise(
      start_year_valid = as.Date(ifelse(n() > 0, paste0(min(hydro_year) - 1, "-10-01"), NA)),
      end_year_valid = as.Date(ifelse(n() > 0, paste0(max(hydro_year), "-09-30"), NA)),
      n_years_valid = n(),
      .groups = "drop")
  
  # determine all observations laying within valid hydrological years
  discharge_valid_total <- discharge_filled_total %>%
    inner_join(discharge_valid_hydro_year %>% dplyr::select(AlpAKaS_ID, hydro_year), by = c("AlpAKaS_ID", "hydro_year"))
  
  saveRDS(file = paste0(path_buff_approx_temp, "/discharge_valid_total.RDS"), discharge_valid_total)
  saveRDS(file = paste0(path_buff_approx_temp, "/discharge_valid_hydro_year.RDS"), discharge_valid_hydro_year)
  saveRDS(file = paste0(path_buff_approx_temp, "/discharge_valid_summary.RDS"), discharge_valid_summary)
  
  
  #-----------------------------------------------------------------------------
  # read summary of data availability of different meteorological products
  temp_avail_meteo_prod <- read.csv(paste0(path_buff_in, "/temp_avail_meteo_prod.csv"))

  # apply function to all meteorological product classes
  compute_hydrometeo_valid(meteo_prod_class = "eobs",
                           temp_avail_meteo_prod = temp_avail_meteo_prod,
                           discharge_valid_hydro_year = discharge_valid_hydro_year,
                           discharge_valid_total = discharge_valid_total)
  compute_hydrometeo_valid(meteo_prod_class = "era5land",
                           temp_avail_meteo_prod = temp_avail_meteo_prod,
                           discharge_valid_hydro_year = discharge_valid_hydro_year,
                           discharge_valid_total = discharge_valid_total)
  compute_hydrometeo_valid(meteo_prod_class = "nat",
                           temp_avail_meteo_prod = temp_avail_meteo_prod,
                           discharge_valid_hydro_year = discharge_valid_hydro_year,
                           discharge_valid_total = discharge_valid_total)

}
