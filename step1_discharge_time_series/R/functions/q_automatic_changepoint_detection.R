### changepoint identification in time series based on daily (aggregated) data

#-------------------------------------------------------------------------------

preproc_q_cpd <- function(path_q_in, path_q_temp, station_meta) {
  
  cpd_segment_summary <- data.frame() # collect changepoint segments per station
  for (i in 1:nrow(station_meta)){
  
    AlpAKaS_ID_i <- station_meta$AlpAKaS_ID[i]
    print(AlpAKaS_ID_i)

    if (station_meta$q_orig_res_daily[i]) {
      
      ### read original time series at daily resolution and remove flagged observations
      file_res_i <- read.csv(paste0(path_q_in, "daily/", AlpAKaS_ID_i, ".csv")) %>%
        filter(!is.na(discharge),
               !qc_flag == TRUE) %>%
        mutate(date = ymd(date),
               temp_res_dom = "daily") %>%
        dplyr::select(date, discharge, temp_res_dom)
      
    } else {
      
      ### read daily aggregated time series
      file_res_i <- readRDS(paste0(path_q_temp, "daily/", AlpAKaS_ID_i, ".rds")) %>%
        mutate(temp_res_dom = "daily")
      
    }
    
    # read identified periods of lower dominant resolution (over longer time series segments)
    periods_res_dominant_i <- read.csv(paste0(path_q_in, "ts_segments_temp_res_dominant.csv")) %>%
      mutate(
        period_start_date =  as.Date(parse_date_time(period_start_date, orders = c("ymd_HMS", "ymd"))),
        period_end_date = as.Date(parse_date_time(period_end_date, orders = c("ymd_HMS", "ymd")))
      ) %>%
      filter(AlpAKaS_ID == AlpAKaS_ID_i)
    # apply resolution flags, if available
    for (p in seq_len(nrow(periods_res_dominant_i))) {
      file_res_i <- file_res_i %>%
        mutate(
          temp_res_dom = if_else(
            date >= periods_res_dominant_i$period_start_date[p] & date <= periods_res_dominant_i$period_end_date[p],
            periods_res_dominant_i$temp_res_dominant[p],
            temp_res_dom
          ))
    }
    
    # filter time series, interpolate data up to 4 weeks and delete larger gaps
    file_res_i <- file_res_i %>%
      filter(temp_res_dom == "daily") %>%
      complete(date = seq(min(date), max(date), by = "day")) %>%
      mutate(discharge = na.approx(discharge, maxgap = 28)) %>%
      filter(!is.na(discharge))
    
    # detect changes in mean and variance using PELT within the "changepoint" package
    # a pre-tuned penalty factor of 250 multiplied by the logarithm of the time series length is set
    cpt_result <- cpt.meanvar(file_res_i$discharge, method = "PELT", penalty = "Manual",
                              pen.value = 250 * log(length(file_res_i$discharge)))

    # extract dates of changepoints by indices
    cpt_dates <- file_res_i$date[cpts(cpt_result)]
    if (length(cpt_dates > 0)){
      df_CPD_dates_i <- data.frame(date = cpt_dates) %>%
        mutate(AlpAKaS_ID = AlpAKaS_ID_i)
    }
    
    # perform segmentation based on minimum, maximum, and changepoint dates
    segment_breaks <- c(min(file_res_i$date), cpt_dates + days(1), max(file_res_i$date))

    # assign each observation to a segment and compute a mean per segment
    file_res_i <- file_res_i %>%
      mutate(cpd_segment_id = cut(date, breaks = segment_breaks, include.lowest = TRUE)) %>%
      mutate(cpd_segment_id = as.integer(factor(cpd_segment_id)))
    
    # summarize identified segments
    segment_summary_i <- file_res_i %>%
      group_by(cpd_segment_id) %>%
      summarise(
        AlpAKaS_ID = AlpAKaS_ID_i,
        period_start_date = min(date),
        period_end_date = max(date),
        mean_discharge = mean(discharge, na.rm = TRUE),
        .groups = "drop") %>%
      dplyr::select(AlpAKaS_ID, period_start_date, period_end_date, cpd_segment_id, mean_discharge)
    if (nrow(segment_summary_i) > 1) {
      cpd_segment_summary <- bind_rows(cpd_segment_summary, segment_summary_i)
    }
    
  }
  

  return(cpd_segment_summary)

}