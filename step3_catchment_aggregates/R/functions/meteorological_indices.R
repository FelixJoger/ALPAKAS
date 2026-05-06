### computation of climatic indices based on observations within valid hydrological years
# defined by availability of both discharge and meteorological data for at least 80% of days within that hydrological year
### this script largely builds on code from the original CAMELS paper (Addor et al., 2017) and has been adapted
# to fit the data structure used in this project (see sections seasonality and extreme precipitation events)

#-------------------------------------------------------------------------------

# function for computation of climatic indices
compute_meteo_indices <- function(meteo_prod_class){
  
  print(meteo_prod_class)
  
  # read overlapping valid hydrometeorological observations for respective catchment type and meteorological product
  hydrometeo_valid_total <- readRDS(file = paste0(path_catch_agg_temp, catch_del_name_r, "/hydrometeo_valid_total_", meteo_prod_class, ".RDS"))
  
  meteo_summary <- hydrometeo_valid_total %>%
    group_by(AlpAKaS_ID) %>%
    arrange(date) %>%
    summarise({
      
      ### mean precipitation and mean temperature
      p_mean = mean(precipitation, na.rm = TRUE)
      t_mean = mean(temperature_mean, na.rm = TRUE)
      
      ### mean hydrological year minimum/maximum of daily mean temperature and
      # frequency of days with mean temperature below <= 0 °C
      hydroannual_extremes <- tibble(
        hydro_year = hydro_year,
        temperature_mean = temperature_mean) %>%
        group_by(hydro_year) %>%
        summarise(
          tgmax = max(temperature_mean, na.rm = TRUE),
          tgmin = min(temperature_mean, na.rm = TRUE),
          cold_days = sum(temperature_mean <= 0, na.rm = TRUE),
          .groups = "drop")
      tg_min_hy_mean <- mean(hydroannual_extremes$tgmin, na.rm = TRUE)
      tg_max_hy_mean <- mean(hydroannual_extremes$tgmax, na.rm = TRUE)
      cold_days_freq <- mean(hydroannual_extremes$cold_days, na.rm = TRUE)
      
      ### seasonality
      # nls() does not fit each year separately but the entire data at once, thus, entirely missing hydrological years are not problematic here
      s_p_first_guess <- 90 - which.max(rapply(split(precipitation, format(date, '%m')), mean, na.rm = TRUE)) * 30
      s_p_first_guess <- s_p_first_guess %% 360 # convert to a value between 0 and 360
      
      fit_temp = nls(temperature_mean ~ mean(temperature_mean, na.rm = TRUE) + delta_t * sin(2 * pi * (yday(date) - s_t) / 365.25),
                     data = cur_data(), start = list(delta_t = 5, s_t = -90))
      fit_prec = nls(precipitation ~ mean(precipitation, na.rm = TRUE) * (1 + delta_p * sin(2 * pi * (yday(date) - s_p) / 365.25)),
                     data = cur_data(), start = list(delta_p = 0.4, s_p = s_p_first_guess))
      
      s_p <- summary(fit_prec)$par['s_p', 'Estimate']
      delta_p <- summary(fit_prec)$par['delta_p', 'Estimate']
      s_t <- summary(fit_temp)$par['s_t','Estimate']
      delta_t <- summary(fit_temp)$par['delta_t', 'Estimate']
      
      # seasonality and timing of precipitation
      delta_p_star <- delta_p * sign(delta_t) * cos(2 * pi * (s_p - s_t) / 365.25)
      
      ### fraction of precipitation falling as snow
      if(any(temperature_mean <= 0 & precipitation > 0, na.rm = TRUE)){
        
        f_s_daily <- sum(precipitation[temperature_mean <= 0], na.rm = TRUE) / sum(precipitation, na.rm = TRUE)
        
      } else {
        
        f_s_daily <- 0
        
      }
      
      ### extreme precipitation indices
      # if there are several segments within the time series (one or several hydrological years missing), compute extreme events for each segment
      # and join results the underlying mean, however, was computed from the full available record rather than separately for each segment
      group_id <- cumsum(c(1, diff(date) > 1))
      segments_day <- split(date, group_id)
      segments_prec <- split(precipitation, group_id)
      
      # make sure previous data point is not considered as continuous event after gaps of one or several hydrological years
      hp_list <- list()
      for (i in seq_along(segments_prec)) {
        
        prec_i <- segments_prec[[i]]
        
        # frequency and duration of high intensity precipitation events
        hp_i <- prec_i >= 5 * mean(precipitation, na.rm = TRUE)
        hp_i[is.na(hp_i)] <- F # if no precipitation data available, consider it is not an event
        hp_list[[i]] <- hp_i
      }
      hp <- unlist(hp_list)
      
      hp_length <- nchar(strsplit(paste(ifelse(hp, 'H', '-'), collapse = ''), '-')[[1]]) # compute number of consecutive high precipitation days
      hp_length <- hp_length[hp_length > 0]
      if(sum(hp_length) != sum(hp)){stop('Unexpected total number of high precip days')}
      
      if(length(hp_length) > 0){ # at least one high precipitation event in the provided time series
        
        hp_freq <- sum(hp) / length(hp) * 365.25
        hp_dur <- mean(hp_length)
        hp_sea <- rapply(split(hp[hp], season[hp], drop = TRUE), length)
        
        if(max(rank(hp_sea) %% 1 != 0)){ # if tie between seasons with the most days with high precipitation, set timing to NA
          
          hp_timing <- NA
          
        } else{
          
          hp_timing <- names(hp_sea)[which.max(hp_sea)]
          
        }
        
      } else { # not a single high precipitation event in the provided time series
        
        hp_freq <- 0
        hp_dur <- 0
        hp_timing <- NA
        
      }
      
      # frequency and duration of low intensity precipitation events
      lp <- precipitation < 1
      lp[is.na(lp)] <- F # if no precip data available, consider it is not an event
      lp_length <- nchar(strsplit(paste(ifelse(lp, 'L', '-'), collapse = ''), '-')[[1]]) # compute number of consecutive low precip days
      lp_length <- lp_length[lp_length > 0]
      
      lp_freq <- sum(lp) / length(lp) * 365.25
      lp_dur <- mean(lp_length)
      lp_sea <- rapply(split(lp[lp], season[lp], drop = TRUE), length)
      
      if(max(rank(lp_sea) %% 1 != 0)){ # if tie between seasons with the most days with low precipitation, set timing to NA
        
        lp_timing <- NA
        
      } else{
        
        lp_timing <- names(lp_sea)[which.max(lp_sea)]
        
      }
      
      
      ### collect all computed indices in tibble
      tibble(
        p_mean,
        p_seasonality = delta_p_star,
        frac_snow = f_s_daily,
        high_prec_freq = hp_freq,
        high_prec_dur = hp_dur,
        high_prec_timing = hp_timing,
        low_prec_freq = lp_freq,
        low_prec_dur = lp_dur,
        low_prec_timing = lp_timing,
        t_mean,
        tg_min_hy_mean,
        tg_max_hy_mean,
        cold_days_freq
      )
      
    }, .groups = "drop")
  
  
  return(meteo_summary)
  
}


### main function to compute meteorological indices
meteo_catch_attr <- function(catch_del_name_r){

  # apply function to all meteorological product classes
  meteo_eobs_catch_attr <- compute_meteo_indices(meteo_prod_class = "eobs")
  meteo_era5land_catch_attr <- compute_meteo_indices(meteo_prod_class = "era5land")
  meteo_nat_catch_attr <- compute_meteo_indices(meteo_prod_class = "nat")
  
  #-----------------------------------------------------------------------------
  ### join and save meteorological catchment attributes
  meteo_catch_attr <- meteo_eobs_catch_attr %>%
    left_join(meteo_era5land_catch_attr, by = "AlpAKaS_ID", suffix = c("_eobs", "_era5land")) %>%
    left_join(meteo_nat_catch_attr %>% rename_with(~ paste0(.x, "_nat"), -AlpAKaS_ID), by = "AlpAKaS_ID") %>%
    mutate(across(where(is.numeric), ~ round(.x, 3)))
  write.csv(meteo_catch_attr, paste0(path_catch_agg_out, catch_del_name_r, "/static_attributes/ALPAKAS_meteorology_attributes.csv"), row.names = FALSE)
  
}