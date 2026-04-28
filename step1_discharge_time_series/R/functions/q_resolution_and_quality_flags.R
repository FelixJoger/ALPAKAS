### resolution and quality flags for the final hourly and daily data to the discharge time series

#-------------------------------------------------------------------------------

preproc_q_flags <- function(final_class, station_meta_i) {

  ### read summary files indicating periods for quality and temporal resolution flags
  
  # periods of lower dominant resolution (over longer time series segments)
  periods_res_dominant_i <- read.csv(paste0(path_q_in, "ts_segments_temp_res_dominant.csv")) %>%
    mutate(
      period_start_date = parse_date_time(period_start_date, orders = c("ymd_HMS", "ymd")),
      period_end_date = parse_date_time(period_end_date, orders = c("ymd_HMS", "ymd"))
    ) %>%
    filter(AlpAKaS_ID == AlpAKaS_ID_i)
  
  # periods of lower quality
  periods_poor_quality_i <- read.csv(paste0(path_q_in, "ts_segments_poor_quality.csv")) %>%
    mutate(
      period_start_date = parse_date_time(period_start_date, orders = c("ymd_HMS", "ymd")),
      period_end_date = parse_date_time(period_end_date, orders = c("ymd_HMS", "ymd"))
    ) %>%
    filter(AlpAKaS_ID == AlpAKaS_ID_i)
  
  # changepoints previously detected via script "q_automatic_changepoint_detection.R" and manually verified
  periods_valid_changepoints_i <- read.csv(paste0(path_q_in, "ts_segments_valid_changepoints.csv")) %>%
    mutate(
      period_start_date = parse_date_time(period_start_date, orders = c("ymd_HMS", "ymd")),
      period_end_date = parse_date_time(period_end_date, orders = c("ymd_HMS", "ymd"))
    ) %>%
    filter(AlpAKaS_ID == AlpAKaS_ID_i)

#-------------------------------------------------------------------------------

  ### final hourly time series
  if (final_class == "hourly") {
    
    # read time series after AOC1 and AOC2 from temp and initialize columns
    file_res_i <- readRDS(paste0(path_q_temp, final_class, "/", AlpAKaS_ID_i, ".rds")) %>%
      filter(!is.na(discharge)) %>%
      mutate(
        temp_res_dom = "hourly",
        dq_issue = NA_character_,
        cpd_segment_id = NA_integer_,
        cpd_segment_main = NA_character_
      )
    
    # filter quality flags for hourly data
    periods_poor_quality_i <- periods_poor_quality_i %>%
      filter(hourly == "TRUE")
  
  ### final daily time series
  } else if (final_class == "daily") {
    
    ### original time series at daily resolution
    if (station_meta_i$q_orig_res_daily) {
      
      # read original daily time series and initialize columns
      file_res_i <- read.csv(paste0(path_q_in, final_class, "/", AlpAKaS_ID_i, ".csv")) %>%
        filter(!is.na(discharge)) %>%
        mutate(date = ymd(date),
               temp_res = "<= day",
               temp_res_dom = "daily",
               qc_detect_method = if_else(qc_flag == TRUE, "manual", NA_character_),
               dq_issue = NA_character_,
               cpd_segment_id = NA_integer_,
               cpd_segment_main = NA_character_
        )
    
    ### daily aggregated time series
    } else {
      
      # read time series after daily aggregation based on instantaneous and hourly data and initialize columns
      file_res_i <- readRDS(paste0(path_q_temp, final_class, "/", AlpAKaS_ID_i, ".rds")) %>%
        filter(!is.na(discharge)) %>%
          mutate(
            qc_flag = if (!"qc_flag" %in% names(.)) FALSE else qc_flag,
            qc_type = if (!"qc_type" %in% names(.)) NA_character_ else qc_type,
            qc_detect_method = if (!"qc_detect_method" %in% names(.)) NA_character_ else qc_detect_method,
            temp_res = "<= day",
            temp_res_dom = "daily",
            dq_issue = NA_character_,
            cpd_segment_id = NA_integer_,
            cpd_segment_main = NA_character_
          )
      
    }
    
    # assign temporal resolution classes allowing for data gaps within a consistent resolution class
    file_res_i <- file_res_i %>%
      mutate(
        time_prev = lag(date),
        diff_prev_hr = as.numeric(difftime(date, time_prev, units = "hours")),
        temp_res = case_when(
          diff_prev_hr > 1 & diff_prev_hr <= 24 & lead(diff_prev_hr) > 1 & lead(diff_prev_hr) <= 24 ~ "<= day",
          diff_prev_hr > 24 & diff_prev_hr <= 168 & lead(diff_prev_hr) > 24 & lead(diff_prev_hr) <= 168 ~ "<= week",
          diff_prev_hr > 168 & lead(diff_prev_hr) > 168 ~ "> week",
          TRUE ~ temp_res
        )
      )
    # assign temporal resolution class to the first row specifically
    if(file_res_i$diff_prev_hr[2] == 24) {
      file_res_i$temp_res[1] <- "<= day"
    } else if (file_res_i$diff_prev_hr[2] > 24 & file_res_i$diff_prev_hr[2] <= 168) {
      file_res_i$temp_res[1] <- "<= week"
    } else {
      file_res_i$temp_res[1] <- "> week"
    }
    
    # convert dates to daily resolution to exclude days flagged days entirely
    periods_res_dominant_i <- periods_res_dominant_i %>%
      mutate(period_start_date = as.Date(period_start_date),
             period_end_date = as.Date(period_end_date))
    periods_poor_quality_i <- periods_poor_quality_i %>%
      # filter quality flags for daily data
      filter(daily == "TRUE") %>%
      mutate(period_start_date = as.Date(period_start_date),
             period_end_date = as.Date(period_end_date))
 
    }
  
  ### apply resolution and quality flags, if available
  for (p in seq_len(nrow(periods_res_dominant_i))) {
    file_res_i <- file_res_i %>%
      mutate(
        temp_res_dom = if_else(
        date >= periods_res_dominant_i$period_start_date[p] & date <= periods_res_dominant_i$period_end_date[p],
        periods_res_dominant_i$temp_res_dominant[p],
        temp_res_dom
      ))
  }
  
  for (p in seq_len(nrow(periods_poor_quality_i))) {
    file_res_i <- file_res_i %>%
      mutate(
        dq_issue = if_else(
        date >= periods_poor_quality_i$period_start_date[p] & date <= periods_poor_quality_i$period_end_date[p],
        periods_poor_quality_i$poor_data_quality[p],
        dq_issue
      ))
  }
  
  for (p in seq_len(nrow(periods_valid_changepoints_i))) {
    file_res_i <- file_res_i %>%
      mutate(
        cpd_segment_id = if_else(as.Date(date) >= periods_valid_changepoints_i$period_start_date[p] &
                                   as.Date(date) <= periods_valid_changepoints_i$period_end_date[p],
                                 periods_valid_changepoints_i$cpd_segment_id[p], cpd_segment_id),
        cpd_segment_main = if_else(as.Date(date) >= periods_valid_changepoints_i$period_start_date[p] &
                                     as.Date(date) <= periods_valid_changepoints_i$period_end_date[p],
                                   as.character(periods_valid_changepoints_i$cpd_segment_main[p]), cpd_segment_main))
  }
  
  # select and sort required columns
  file_res_i_final <- file_res_i %>%
    dplyr::select(date, discharge, qc_flag, qc_type, qc_detect_method, temp_res, temp_res_dom,
           dq_issue, cpd_segment_id, cpd_segment_main) %>%
    mutate(across(where(is.numeric), ~ round(.x, 3)))
  
  # make sure all dates are in correct format for hourly data
  if (final_class == "hourly") {
    file_res_i_final <- file_res_i_final %>%
      mutate(date = format(date, "%Y-%m-%d %H:%M:%S"))
  }
  
  return(file_res_i_final)
  
}