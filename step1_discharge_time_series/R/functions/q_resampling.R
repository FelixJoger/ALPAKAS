### resampling of discharge time series to hourly and daily temporal resolutions

#-------------------------------------------------------------------------------

preproc_q_resample <- function(agg_class, input) {
  
  ### perform hourly and daily data aggregation of instantaneous data
  if (agg_class == "hourly") {
    
    date_unit <- "hour"
    date_seq <- "1 min"
    thresh_interp <- 1
    
  ### perform daily data aggregation of hourly (aggregated) data
  } else if (agg_class == "daily") {
    
    date_unit <- "day"
    date_seq <- "1 hour"
    thresh_interp <- 24
    
  }
  
  file_res_i <- input %>%
    # remove all previously identified outliers and artefacts
    filter(!qc_flag) %>%
    dplyr::select(date, discharge) %>%
    # ensure that no timestamps include seconds
    mutate(date = floor_date(date, unit = "minute"))
  
  # sequence of dates from first to last measurement
  time_grid_i <- data.frame(date = seq(min(file_res_i$date), max(file_res_i$date), by = date_seq))
  # join observed values and interpolate missing values for trapezoidal integration
  file_res_i_interp <- time_grid_i %>%
    left_join(file_res_i, by = "date") %>%
    mutate(value_interp  = na.approx(discharge, x = date, na.rm = FALSE))
  
  # extract actual observation times
  real_obs_times_i <- file_res_i_interp %>%
    filter(!is.na(discharge)) %>%
    pull(date)
  
  ### only keep interpolated values meeting the predefined conditions
  for (n in 1:(length(real_obs_times_i) - 1)) {
    
    current_time_i <- real_obs_times_i[n+1]
    previous_time_i <- real_obs_times_i[n]
    
    # check if the previous actual observation was within a certain period before the measurement
    has_obs_prev <- isTRUE(previous_time_i >= current_time_i - hours(thresh_interp))
    
    # if the condition is TRUE, retain interpolated values; otherwise, remove all interpolated values up to the
    # previous observation value to ensure interpolation is only applied when preceding observations
    # are within 1 hour (for hourly resolution) or 24 hours (for daily resolution)
    if (!has_obs_prev) {
      file_res_i_interp <- file_res_i_interp %>%
        mutate(value_interp = if_else(
          date < current_time_i & date > previous_time_i & is.na(discharge),
          NA_real_,
          value_interp))
    }
  }
  
  # aggregate values within each hour or day, depending on the aggregation class (trapezoidal integration)
  file_res_i_final <- file_res_i_interp %>%
    mutate(date = floor_date(date, unit = date_unit)) %>%
    group_by(date) %>%
    summarise(discharge = mean(value_interp, na.rm = TRUE), .groups = "drop") %>%
    mutate(discharge = ifelse(is.nan(discharge), NA, discharge))
  
  
  return(file_res_i_final)
    
}