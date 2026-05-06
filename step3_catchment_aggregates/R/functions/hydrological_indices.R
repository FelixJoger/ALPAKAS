### computation of hydrological indices based on observations within valid hydrological years
# defined by availability of both discharge and meteorological data for at least 80% of days within that hydrological year

#-------------------------------------------------------------------------------

# function for computation of cross-correlational indices
compute_cc_indices <- function(meteo_prod_class){
  
  print(meteo_prod_class)
  
  # read overlapping valid hydrometeorological observations for respective catchment type and meteorological product
  hydrometeo_valid_total <- readRDS(file = paste0(path_catch_agg_temp, catch_del_name_r, "/hydrometeo_valid_total_", meteo_prod_class, ".RDS"))
  
  # missing observations in time series are filled with NA values, as this is required for the cross-correlational analysis
  hydrometeo_valid_total_filled <- hydrometeo_valid_total %>%
    group_by(AlpAKaS_ID) %>%
    complete(date = seq(min(date), max(date), by = "1 day")) %>%
    arrange(AlpAKaS_ID, date) %>%
    mutate(hydro_year = if_else(month(date) >= 10, year(date) + 1, year(date)))
  
  # use the entire precipitation time series and set discharge data with <80% completeness per hydrological year to NA
  cc_summary <- hydrometeo_valid_total_filled %>%
    group_by(AlpAKaS_ID) %>%
    arrange(date) %>%
    # interpolation is performed for a maximum gap of five days as the standard setting in KarstID
    mutate(discharge = na.spline(discharge, method = "monoH.FC", maxgap = 5, na.rm = FALSE)) %>%
    summarise({
      
      # na.pass computes correlations only for valid pairs per lag
      cc = ccf(precipitation, discharge, lag.max = 100, plot = FALSE, na.action = na.pass)
      
      # keep only lags <= 0 and set r to zero, if negative
      lags_all <- drop(cc$lag)
      acf_all  <- drop(cc$acf)
      lags_ok  <- lags_all[lags_all <= 0]
      acf_ok   <- acf_all[lags_all <= 0]
      
      tibble(
        cc_peak_lag = lags_ok[which.max(acf_ok)],
        cc_peak_r   = max(max(acf_ok, na.rm = TRUE), 0)
      )
      
    }, .groups = "drop")
  
  return(cc_summary)
  
}


# function to consider a minimum fraction of available discharge data
rollmean_min <- function(x, k, min_frac = 0.8) {
  rollapply(x, width = k, align = "center", partial = FALSE,
            FUN = function(w) {
              ok <- sum(!is.na(w))
              if (ok >= ceiling(min_frac * k)) mean(w, na.rm = TRUE) else NA_real_
            }
  )
}


### main function to compute hydrological indices
hydro_catch_attr <- function(catch_del_name_r, catch_repr){
  
  # indicators based solely on the hydrographs are computed for the valid hydrological years within the time period of E-OBS/Era5land data
  hydrometeo_valid_total_eobs <- readRDS(file = paste0(path_catch_agg_temp, "catchment_approx/hydrometeo_valid_total_eobs.RDS")) %>%
    filter(AlpAKaS_ID %in% unique(catch_repr$AlpAKaS_ID))
  
  # missing years in time series are filled with NA values, as this is required for the autocorrelational and cross-correlational analysis
  hydrometeo_valid_total_eobs_filled <- hydrometeo_valid_total_eobs %>%
    group_by(AlpAKaS_ID) %>%
    complete(date = seq(min(date), max(date), by = "1 day")) %>%
    arrange(AlpAKaS_ID, date) %>%
    mutate(hydro_year = if_else(month(date) >= 10, year(date) + 1, year(date)))
  
  ### statistical indicators based on hydrographs
  discharge_summary <- hydrometeo_valid_total_eobs_filled %>%
    group_by(AlpAKaS_ID) %>%
    summarise(q_mean = mean(discharge, na.rm = TRUE),
              q_min = min(discharge, na.rm = TRUE),
              q_10 = quantile(discharge, probs = 0.10, na.rm = TRUE),
              q_90 = quantile(discharge, probs = 0.90, na.rm = TRUE),
              q_max = max(discharge, na.rm = TRUE),
              q_sd = sd(discharge, na.rm = TRUE),
              CV = q_sd / q_mean,
              SVC = q_90 / q_10,
              .groups = "drop")
  
  # use the entire precipitation time series and set discharge data with <80% completeness per hydrological year to NA
  ac_summary <- hydrometeo_valid_total_eobs_filled %>%
    group_by(AlpAKaS_ID) %>%
    arrange(date) %>%
    # interpolation is performed for a maximum gap of five days as the standard setting in KarstID
    mutate(discharge = na.spline(discharge, method = "monoH.FC", maxgap = 5, na.rm = FALSE)) %>%
    summarise({
      
      
      ### autocorrelational analysis
      # na.pass keeps the series and computes each lag using the available non-missing pairs
      ac <- acf(discharge,
                lag.max = 125, # see Mangin (1984)
                type = "correlation",
                plot = FALSE,
                na.action = na.pass)
      # drop lag 0
      lags <- as.numeric(ac$lag[-1])
      acs <- drop(ac$acf[-1])
  
      memory_idx <- which(acs <= 0.2)[1]
  
      
      ### regulation time (proxy)
      # both KarstID (Larocque) and XLKarst (Tukey) require continuous discharge time series as only then regulation time
      # from spectral analysis is physically meaningful; thus, a proxy formula by XLKarst can be used to assess regulation time
      # σ₂₅₀ = standard deviation of the 250-day moving-average filtered series
      
      sd_raw <- sd(discharge, na.rm = TRUE)
      
      # 250-day centered moving mean with requirement of sufficient non-NA values per window (>=80%)
      Q_ma <- rollmean_min(discharge, k = 250, min_frac = 0.8)
      sd_ma <- sd(Q_ma, na.rm = TRUE)
      
      # proxy for regulation time
      reg_time <- 125 * (sd_ma / sd_raw)^2
      
      
      ### collect all computed indices in tibble
      tibble(
        ac_memory_lag = memory_idx,
        regulation_time = reg_time
      )
      
    }, .groups = "drop")
  
  hydrograph_site_attr <- discharge_summary %>%
    left_join(ac_summary, by = "AlpAKaS_ID")

  
  ### cross-correlational analysis
  # apply function to all meteorological product classes
  cc_eobs_catch_attr <- compute_cc_indices(meteo_prod_class = "eobs")
  cc_era5land_catch_attr <- compute_cc_indices(meteo_prod_class = "era5land")
  cc_nat_catch_attr <- compute_cc_indices(meteo_prod_class = "nat")

  #-----------------------------------------------------------------------------
  ### join and save hydrological catchment attributes
  cc_catch_attr <- cc_eobs_catch_attr %>%
    left_join(cc_era5land_catch_attr, by = "AlpAKaS_ID", suffix = c("_eobs", "_era5land")) %>%
    left_join(cc_nat_catch_attr %>% rename_with(~ paste0(.x, "_nat"), -AlpAKaS_ID), by = "AlpAKaS_ID")
  hydrology_attr_list <- list(hydrograph_site_attr, cc_catch_attr)
  hydrology_attr <- reduce(hydrology_attr_list, left_join, by = "AlpAKaS_ID") %>%
    mutate(across(where(is.numeric), ~ round(.x, 3)))
  write.csv(hydrology_attr, paste0(path_catch_agg_out, catch_del_name_r, "/static_attributes/ALPAKAS_hydrology_attributes.csv"), row.names = FALSE)
  
}