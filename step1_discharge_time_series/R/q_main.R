### main pipeline that sources all other required scrips for discharge time series preprocessing

# libraries
# core / general utilities
library(here)
library(lubridate)
library(stringr)

# data manipulation
library(dplyr)
library(tidyr)
library(zoo)

# statistics / time series analysis
library(changepoint)

#-------------------------------------------------------------------------------
# check if working directory is set correctly (e.g. start R from within the project directory)
here::i_am("step1_discharge_time_series/R/q_main.R")

path_q_in <- "step1_discharge_time_series/input_data/"
path_q_func <- "step1_discharge_time_series/R/functions/"
path_q_temp <- "step1_discharge_time_series/temp/"
path_q_out <- "step1_discharge_time_series/output_data/"

source(paste0(path_q_func, "q_automatic_outlier_detection.R"))
source(paste0(path_q_func, "q_resampling.R"))
source(paste0(path_q_func, "q_automatic_changepoint_detection.R"))
source(paste0(path_q_func, "q_resolution_and_quality_flags.R"))

dir.create(paste0(path_q_temp, "daily"))
dir.create(paste0(path_q_temp, "hourly"))
#-------------------------------------------------------------------------------

# read meta data file
station_meta <- read.csv(file = "station_meta_input.csv")


### inst ###

station_meta_inst <- station_meta %>%
  filter(as.logical(q_orig_res_inst))

for (i in 1:nrow(station_meta_inst)){

  AlpAKaS_ID_i <- station_meta_inst$AlpAKaS_ID[i]
  print(AlpAKaS_ID_i)
  
  ### automatic outlier detection 1
  q_aoc1 <- preproc_q_aoc(res_class = "inst", AlpAKaS_ID_i = AlpAKaS_ID_i, aoc_no = "aoc1", input = NULL)
  
  ### hourly aggregation
  q_agg_hourly <- preproc_q_resample(agg_class = "hourly", input = q_aoc1)
  
  ### AOC2
  q_aoc2 <- preproc_q_aoc(res_class = "inst", AlpAKaS_ID_i = NULL, aoc_no = "aoc2", input = q_agg_hourly)
  
  ### daily aggregation
  q_agg_daily <- preproc_q_resample(agg_class = "daily", input = q_aoc2)
  
  # save intermediate results
  saveRDS(q_aoc2, paste0(path_q_temp, "hourly/", AlpAKaS_ID_i, ".RDS"))
  saveRDS(q_agg_daily, paste0(path_q_temp, "daily/", AlpAKaS_ID_i, ".RDS"))
  
}


### hourly ###

station_meta_hourly <- station_meta %>%
  filter(as.logical(q_orig_res_hourly))

for (i in 1:nrow(station_meta_hourly)){

  AlpAKaS_ID_i <- station_meta_hourly$AlpAKaS_ID[i]
  print(AlpAKaS_ID_i)
  
  ### automatic outlier detection 1
  q_aoc1 <- preproc_q_aoc(res_class = "hourly", AlpAKaS_ID_i = AlpAKaS_ID_i, aoc_no = "aoc1", input = NULL)
  
  # save intermediate results
  saveRDS(q_aoc1, paste0(path_q_temp, "hourly", "/", AlpAKaS_ID_i, ".RDS"))
  
  # only compute daily aggregation for stations for which not also daily time series exist
  if (station_meta_hourly$q_orig_res_daily[i] == "FALSE") {
    
    ### daily aggregation
    q_agg_daily <- preproc_q_resample(agg_class = "daily", input = q_aoc1)
    
    # save intermediate results
    saveRDS(q_agg_daily, paste0(path_q_temp, "daily/", AlpAKaS_ID_i, ".RDS"))
    
  }

}

#-------------------------------------------------------------------------------
### manually detect quality and resolution flags in time series and adapt summary files for detected periods accordingly
### perform automatic changepoint detection
cpd_segment_summary <- preproc_q_cpd(path_q_in = path_q_in, path_q_temp = path_q_temp, station_meta = station_meta)
### review identified changepoints manually and adapt summary file for detected periods accordingly
#-------------------------------------------------------------------------------

### apply quality and resolution flags
for (i in 1:nrow(station_meta)){
  
  station_meta_i = station_meta[i, ]
  AlpAKaS_ID_i <- station_meta_i$AlpAKaS_ID
  print(AlpAKaS_ID_i)
  
  ### assign quality and resolution flags and save final daily time series
  q_final_daily_flags <- preproc_q_flags(final_class = "daily", station_meta_i = station_meta_i)
  write.csv(file = paste0(path_q_out, "daily/AlpAKaS_discharge_daily_", AlpAKaS_ID_i, ".csv"), q_final_daily_flags, row.names = FALSE)
  
  if (station_meta$q_hourly_available[i]) {
    
    ### assign quality and resolution flags and save final hourly time series
    q_final_hourly_flags <- preproc_q_flags(final_class = "hourly", station_meta_i = station_meta_i)
    write.csv(file = paste0(path_q_out, "hourly/AlpAKaS_discharge_hourly_", AlpAKaS_ID_i, ".csv"), q_final_hourly_flags, row.names = FALSE)
    
  }
  
}

#-------------------------------------------------------------------------------