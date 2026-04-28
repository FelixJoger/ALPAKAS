### main pipeline that sources all other required scrips for catchment-scale aggregations

# libraries
# core / general utilities
library(here)
library(lubridate)

# data manipulation
library(dplyr)
library(tidyr)
library(purrr)
library(readr)
library(tibble)
library(zoo)

# spatial vector data
library(sf)

# raster / spatial gridded data
library(terra)
library(exactextractr)

# statistics
library(stats)

#-------------------------------------------------------------------------------
# check if working directory is set correctly (e.g. start R from within the project directory)
here::i_am("step3_catchment_aggregates/R/catch_agg_main.R")

path_catch_agg_in <- "step3_catchment_aggregates/input_data/"
path_catch_agg_func <- "step3_catchment_aggregates/R/functions/"
path_catch_agg_temp <- "step3_catchment_aggregates/temp/"
path_catch_agg_out <- "step3_catchment_aggregates/output_data/"

source(paste0(path_catch_agg_func, "meteo_geometry_conversion.R"))
source(paste0(path_catch_agg_func, "static_attr_geometry_conversion.R"))
source(paste0(path_catch_agg_func, "meteo_spatial_aggregation.R"))
source(paste0(path_catch_agg_func, "static_attr_spatial_aggregation.R"))
source(paste0(path_catch_agg_func, "hydrometeo_valid_total.R"))
source(paste0(path_catch_agg_func, "meteorological_indices.R"))
source(paste0(path_catch_agg_func, "hydrological_indices.R"))

#-------------------------------------------------------------------------------
### perform geometry conversion of geodata products (only required once)
meteo_geom_convers()
static_attr_geom_convers()
#-------------------------------------------------------------------------------

# define name of catchment representations
catch_del_name <- c("catchment_approx", "catchment_expert")

for (r in 1:length(catch_del_name)) {

  catch_del_name_r <- catch_del_name[r]
  
  # read catchment file
  if (catch_del_name_r == "catchment_expert") {
    catch_repr <- st_read(paste0(path_catch_agg_in, "catchment_delineations/", catch_del_name_r, ".geojson"), quiet = TRUE)
  } else {
    catch_repr <- st_read(paste0("step2_buffer_approximations/output_data/catchment_delineations/", catch_del_name_r, ".geojson"), quiet = TRUE)
  }

  for (i in 1:nrow(catch_repr)){

    catch_repr_i <- catch_repr[i, ]

    AlpAKaS_ID_i <- catch_repr_i$AlpAKaS_ID
    print(AlpAKaS_ID_i)

    ### spatial aggregation of meteorological time series
    meteo_catch_res <- meteo_catch_agg(AlpAKaS_ID_i = AlpAKaS_ID_i, catch_repr_i = catch_repr_i)
    write.csv(meteo_catch_res[[1]], paste0(path_catch_agg_out, catch_del_name_r, "/meteorological_time_series/daily/AlpAKaS_meteo_daily_", AlpAKaS_ID_i, ".csv"), row.names = FALSE)
    write.csv(meteo_catch_res[[2]], paste0(path_catch_agg_out, catch_del_name_r, "/meteorological_time_series/hourly/AlpAKaS_meteo_hourly_", AlpAKaS_ID_i, ".csv"), row.names = FALSE)

  }

  ### spatial aggregation of static attributes
  static_attr_catch_agg(catch_del_name_r = catch_del_name_r, catch_repr = catch_repr)

  ### hydrometeorological indices directly derived from time series data
  hydrometeo_valid(catch_del_name_r = catch_del_name_r, catch_repr = catch_repr)
  meteo_catch_attr(catch_del_name_r = catch_del_name_r)
  hydro_catch_attr(catch_del_name_r = catch_del_name_r, catch_repr = catch_repr)
  
}

#-------------------------------------------------------------------------------