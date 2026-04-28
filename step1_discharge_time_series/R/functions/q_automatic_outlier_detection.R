### automatic outlier detection 1 and 2 for discharge time series

#-------------------------------------------------------------------------------

preproc_q_aoc <- function(res_class, AlpAKaS_ID_i, aoc_no, input) {
  
  ### aoc1: performed for instantaneous and hourly original time series
  if (aoc_no == "aoc1") {
    
    # read original files after manual data cleaning (flags)
    file_res_i_init <- read.csv(paste0(path_q_in, res_class, "/", AlpAKaS_ID_i, ".csv")) %>%
      mutate(date = parse_date_time(date, orders = c("ymd HMS", "ymd HM", "ymd")))
    # remove manually identified outliers and artefacts
    file_res_i_automatic <- file_res_i_init %>%
      filter(qc_flag != TRUE)
    # extract manually identified outliers and artefacts
    file_res_i_dq_manual <- file_res_i_init %>%
      filter(qc_flag == TRUE) %>%
      mutate(qc_detect_method = "manual",
             temp_res_prev = NA_character_,
             temp_res = NA_character_)


  ### aoc2: performed for hourly aggregated time series
  } else if (aoc_no == "aoc2") {
    
    # read hourly aggregated time series
    file_res_i_init <- input %>%
      na.omit(discharge) %>%
      mutate(qc_type = NA_character_,
             qc_flag = FALSE)
    # prepare file for autmatic outlier removal
    file_res_i_automatic <- file_res_i_init
    # extract empty data frame as no outliers or artefacts exist after aggregation
    file_res_i_dq_manual <- file_res_i_init[0, ] %>%
      mutate(qc_detect_method = "manual", 
             temp_res_prev = NA_character_,
             temp_res = NA_character_,
             qc_type = NA_character_,
             qc_flag = FALSE)
    
  }

  #-----------------------------------------------------------------------------
  
  ### start automatic outlier detection
  
  # initialize iteration loop that runs until no further outliers are identified
  df_remaining <- file_res_i_automatic
  all_outliers <- data.frame()  # to collect all detected outliers
  iteration <- 1
  results_list <- list()

  repeat {

    cat("Iteration:", iteration, "\n")

    # save copy before detection
    file_res_i <- df_remaining %>%
      dplyr::select(all_of(names(file_res_i_init)))

    # initialize column for temporal resolution classes
    file_res_i <- file_res_i %>%
      mutate(temp_res = "<= hour")

    # assign temporal resolution classes for each data points based on preceding value
    file_res_i <- file_res_i %>%
      mutate(
        # assign temporal resolution classes focusing on preceding value
        time_prev = lag(date),
        diff_prev_hr = as.numeric(difftime(date, time_prev, units = "hours")),
        temp_res_prev = if_else(!is.na(diff_prev_hr) & diff_prev_hr <= 1, "<= hour",
                                if_else(!is.na(diff_prev_hr) & diff_prev_hr > 1 & diff_prev_hr <= 24, "<= day",
                                        if_else(!is.na(diff_prev_hr) & diff_prev_hr > 24 & diff_prev_hr <= 168, "<= week",
                                                if_else(!is.na(diff_prev_hr) & diff_prev_hr > 168, "> week", NA_character_)))),

        # check whether temporal resolution towards preceding and following value fall in the same class
        has_unique_temp_res = temp_res_prev == lead(temp_res_prev),

        # assign temporal resolution classes allowing for data gaps within the same resolution class
        temp_res = case_when(
          diff_prev_hr > 1 & diff_prev_hr <= 24 & lead(diff_prev_hr) > 1 & lead(diff_prev_hr) <= 24 ~ "<= day",
          diff_prev_hr > 24 & diff_prev_hr <= 168 & lead(diff_prev_hr) > 24 & lead(diff_prev_hr) <= 168 ~ "<= week",
          diff_prev_hr > 168 & lead(diff_prev_hr) > 168 ~ "> week",
          TRUE ~ temp_res)
      )

    # assign the correct class to the first row specifically
    if(file_res_i$diff_prev_hr[2] == 1) {
      file_res_i$temp_res_prev[1] <- "<= hour"
      file_res_i$temp_res[1] <- "<= hour"
    } else if(file_res_i$diff_prev_hr[2] > 1 & file_res_i$diff_prev_hr[2] <= 24) {
      file_res_i$temp_res_prev[1] <- "<= day"
      file_res_i$temp_res[1] <- "<= day"
    } else if (file_res_i$diff_prev_hr[2] > 24 & file_res_i$diff_prev_hr[2] <= 168) {
      file_res_i$temp_res_prev[1] <- "<= week"
      file_res_i$temp_res[1] <- "<= week"
    } else {
      file_res_i$temp_res_prev[1] <- "> week"
      file_res_i$temp_res[1] <- "> week"
    }

    # calculate change rates to preceding observations
    file_res_i <- file_res_i %>%
      mutate(change_rate_prev = (discharge - lag(discharge)) / diff_prev_hr)

    # summarize number of observations with positive and negative change rates in resolution class "<= hour"
    temp_res_sum <- file_res_i %>%
      filter(temp_res_prev == "<= hour") %>%
      summarise(
        n_positive = sum(change_rate_prev > 0, na.rm = TRUE),
        n_negative = sum(change_rate_prev <= 0, na.rm = TRUE)
        )
    
    # if either number of positive values or the number of negative values is below 1000,
    # skip the rest of this loop iteration and continue with the next one (no outlier detection applied)
    if (temp_res_sum$n_positive < 1000 | temp_res_sum$n_negative < 1000){
      next
    }

    # compute percentiles of each data point in resolution class "<= hour" and assign bins based on deciles
    file_res_i <- file_res_i %>%
      mutate(
        percentile = case_when(temp_res_prev == "<= hour" ~ ecdf(discharge)(discharge) * 100, TRUE ~ NA_real_),
        percentile_bin = cut(percentile, breaks = seq(0, 100, by = 10), include.lowest = TRUE, right = FALSE,
                             labels = paste0(seq(0, 90, 10), "-", seq(10, 100, 10), "%")))

    # compute mean negative and mean positive change rates per bin
    avg_by_bin <- file_res_i %>%
      group_by(percentile_bin) %>%
      summarise(
        change_rate_neg_mean = mean(change_rate_prev[change_rate_prev <= 0], na.rm = TRUE),
        change_rate_pos_mean = mean(change_rate_prev[change_rate_prev > 0], na.rm = TRUE),
        .groups = "drop") %>%
      na.omit()
    
    file_res_i <- file_res_i %>%
      # join mean negative and mean positive change rates per bin
      left_join(avg_by_bin, by = "percentile_bin") %>%
      # detect outliers based on pre-tuned change rate factors
      mutate(sharp_dip = change_rate_prev < 10 * change_rate_neg_mean &
               lead(change_rate_prev) > 10 * change_rate_pos_mean,
             sharp_spike = change_rate_prev > 20 * change_rate_pos_mean &
               lead(change_rate_prev) < 20 * change_rate_neg_mean,
             # add categorical column for different types of outliers and a quality control flag for outliers
             qc_type = case_when(sharp_dip ~ "sharp_dip",
                                 sharp_spike ~ "sharp_spike"),
             qc_flag = (sharp_dip | sharp_spike) & has_unique_temp_res
      )
    
    # save result of iteration
    results_list[[paste0("iteration_", iteration)]] <- file_res_i

    # extract outliers detected in iteration
    new_outliers <- file_res_i %>% filter(qc_flag == TRUE)

    # add outliers detected in iteration to the full list
    all_outliers <- bind_rows(all_outliers, new_outliers)

    # drop outliers for next iteration
    df_remaining <- file_res_i %>%
      filter(qc_flag != TRUE | is.na(qc_flag))

    # stop if no additional outliers found
    if (nrow(new_outliers) == 0) break

    iteration <- iteration + 1

  } # closing repeat loop

  # prepare cleaned data and automatic outlier detection file
  file_res_i_cleaned <- df_remaining %>%
    dplyr::select(all_of(names(file_res_i_init)), temp_res_prev, temp_res) %>%
    mutate(qc_detect_method = NA_character_,
           qc_type = NA_character_,
           qc_flag = FALSE)
  
  # if no outliers were detected, initialize empty data frame with correct structure
  if (nrow(all_outliers) == 0){
    file_res_i_dq_automatic <- file_res_i_init[0, ] %>%
      mutate(qc_detect_method = NA_character_,
             temp_res_prev = NA_character_,
             temp_res = NA_character_,
             qc_type = NA_character_,
             qc_flag = FALSE)
  } else {
    file_res_i_dq_automatic <- all_outliers %>%
      dplyr::select(all_of(names(file_res_i_init)), temp_res_prev, temp_res, qc_type, qc_flag) %>%
      mutate(qc_detect_method = "automatic")
  }

  # bind outliers and artefacts to cleaned file and select required columns
  file_res_i_flagged <- rbind(file_res_i_cleaned, file_res_i_dq_manual, file_res_i_dq_automatic) %>%
    arrange(date) %>%
    dplyr::select(date, discharge, qc_flag, qc_type, qc_detect_method, temp_res)
  

  return(file_res_i_flagged)

}