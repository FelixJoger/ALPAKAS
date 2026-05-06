### integration of metadata attributes to obtain final metadata file

#-------------------------------------------------------------------------------

metadata_integr <- function() {
  
  # read input metadata file 
  input_meta <- read.csv("station_meta_input.csv")
  
  #-----------------------------------------------------------------------------
  ### time series attributes
  # initialize columns for daily and hourly time series attributes
  input_meta_add <- input_meta %>%
    mutate(
      q_daily_start = NA_character_,
      q_daily_end = NA_character_,
      q_daily_obs = NA_integer_,
      q_daily_compl = NA_integer_,
      q_hourly_start = NA_character_,
      q_hourly_end = NA_character_,
      q_hourly_obs = NA_integer_,
      q_hourly_compl = NA_integer_,
    )
  
  ### daily ###
  for(i in 1:nrow(input_meta_add)){
    
    AlpAKaS_ID_i <- input_meta_add$AlpAKaS_ID[i]
    print(AlpAKaS_ID_i)
  
    # load processed data file and extract time series attributes
    file_res_i <- read.csv(paste0("step1_discharge_time_series/output_data/daily/ALPAKAS_discharge_daily_", AlpAKaS_ID_i,".csv"))
    input_meta_add <- input_meta_add %>%
      mutate(
        q_daily_start  = if_else(AlpAKaS_ID == AlpAKaS_ID_i, as.character(file_res_i$date[1]), q_daily_start),
        q_daily_end  = if_else(AlpAKaS_ID == AlpAKaS_ID_i, as.character(file_res_i$date[nrow(file_res_i)]), q_daily_end),
        q_daily_obs = if_else(AlpAKaS_ID == AlpAKaS_ID_i, nrow(file_res_i), q_daily_obs),
        q_daily_compl = if_else(AlpAKaS_ID == AlpAKaS_ID_i, round(nrow(file_res_i) / as.numeric(difftime(q_daily_end, q_daily_start, units = "days") + 1), 3), q_daily_compl)
        )
      
  }
  
  ### hourly ###
  input_meta_hourly <- input_meta_add %>%
    filter(q_hourly_available)
  
  for(i in 1:nrow(input_meta_hourly)){
    
    AlpAKaS_ID_i <- input_meta_hourly$AlpAKaS_ID[i]
    print(AlpAKaS_ID_i)
    
    # load processed data file and extract time series attributes
    file_res_i <- read.csv(paste0("step1_discharge_time_series/output_data/hourly/ALPAKAS_discharge_hourly_", AlpAKaS_ID_i,".csv"))
    input_meta_add <- input_meta_add %>%
      mutate(
        q_hourly_start = if_else(AlpAKaS_ID == AlpAKaS_ID_i, as.character(file_res_i$date[1]), q_hourly_start),
        q_hourly_end = if_else(AlpAKaS_ID == AlpAKaS_ID_i, as.character(file_res_i$date[nrow(file_res_i)]), q_hourly_end),
        q_hourly_obs = if_else(AlpAKaS_ID == AlpAKaS_ID_i, nrow(file_res_i), q_hourly_obs),
        q_hourly_compl = if_else(AlpAKaS_ID == AlpAKaS_ID_i, round(nrow(file_res_i) / as.numeric(difftime(q_hourly_end, q_hourly_start, units = "hours") + 1), 3), q_hourly_compl)
      )
  
  }
  
  
  # add time series attributes for valid hydrological years defined by the pre-computed availability of both discharge and meteorological data
  ### eobs ###
  hydrometeo_valid_summary_eobs <- readRDS(file = "step2_buffer_approximations/temp/hydrometeo_valid_summary_eobs.RDS")
  ### era5land ###
  hydrometeo_valid_summary_era5land <- readRDS(file = "step2_buffer_approximations/temp/hydrometeo_valid_summary_era5land.RDS")
  ### national ###
  hydrometeo_valid_summary_nat <- readRDS(file = "step2_buffer_approximations/temp/hydrometeo_valid_summary_nat.RDS")
  # join  attributes
  hydrometeo_valid_summary <- hydrometeo_valid_summary_eobs %>%
    left_join(hydrometeo_valid_summary_era5land, by = "AlpAKaS_ID", suffix = c("_eobs", "_era5land")) %>%
    left_join(hydrometeo_valid_summary_nat %>% rename_with(~ paste0(.x, "_nat"), -AlpAKaS_ID), by = "AlpAKaS_ID")
  
  # join attributes and make sure at least two valid hydrological years are available
  input_meta_add <- input_meta_add %>%
    left_join(hydrometeo_valid_summary, by = "AlpAKaS_ID") %>%
    filter(ind_years_valid_eobs >= 2)
  
  
  #-----------------------------------------------------------------------------
  ### add catchment areas of delineated catchments and catchment buffers
  # read catchment files and compute catchment areas
  catch_expert <- st_read(paste0("step3_catchment_aggregates/input_data/catchment_delineations/catchment_expert.geojson"), quiet = TRUE) %>%
    mutate(area_catch_expert = round(set_units(st_area(.), km^2) %>% as.numeric(), 3))
  catch_approx <- st_read(paste0(path_buff_approx_out, "/catchment_delineations/catchment_approx.geojson"), quiet = TRUE) %>%
    mutate(area_catch_approx = round(set_units(st_area(.), km^2) %>% as.numeric(), 3))
  
  # join attributes
  input_meta_add <- input_meta_add %>%
    mutate(catch_expert_available = AlpAKaS_ID %in% catch_expert$AlpAKaS_ID,
           catch_expert_shp = AlpAKaS_ID %in% catch_expert$AlpAKaS_ID[catch_expert$publish_shp]) %>%
    left_join(catch_expert %>% st_drop_geometry() %>%
                dplyr::select(AlpAKaS_ID, area_catch_expert), by = "AlpAKaS_ID") %>%
    left_join(catch_approx %>% st_drop_geometry() %>%
                dplyr::select(AlpAKaS_ID, area_catch_approx), by = "AlpAKaS_ID") %>%
    relocate(catch_expert_available, catch_expert_shp, area_catch_expert, area_catch_approx, .after = comment_catchment)
  
  
  # save final metadata file
  write.csv(input_meta_add, file = paste0(path_buff_approx_out, "/ALPAKAS_station_meta.csv"), row.names = FALSE)

}