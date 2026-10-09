### main pipeline that sources all other required scrips for catchment approximations based on the buffer approach

# libraries
# core / general utilities
library(rlang)
library(glue)
library(here)
library(lubridate)
library(units)

# data manipulation
library(dplyr)
library(tidyr)
library(purrr)
library(readr)
library(tidyverse)

# spatial vector data
library(sf)

# raster / spatial gridded data
library(terra)
library(raster)
library(stars)
library(exactextractr)

# NetCDF handling
library(ncdf4)

# GIS / hydrology / terrain tools
library(whitebox)

# visualization
library(ggplot2)
library(tmap)
library(cols4all)


rm(list = ls())
cat("\014")

#-------------------------------------------------------------------------------
# check if working directory is set correctly (e.g. start R from within the project directory)
here::i_am("step2_buffer_approximations/R/buffer_approx_main.R")

path_buff_in <- here::here("step2_buffer_approximations", "input_data")
path_buff_func <- here::here("step2_buffer_approximations", "R", "functions")
path_buff_approx_temp <- here::here("step2_buffer_approximations", "temp")
path_buff_approx_out <- here::here("step2_buffer_approximations", "output_data")


# load functions
source(file.path(path_buff_func, "hydrometeo_valid_hydro_years.R"))
source(file.path(path_buff_func, "compute_recharge_era5.R"))
source(file.path(path_buff_func, "iterative_buffers.R"))
source(file.path(path_buff_func, "compute_topographic_catchment.R"))
source(file.path(path_buff_func, "metadata_integration.R"))

#-------------------------------------------------------------------------------

# read meta data file
station_meta <- read.csv(file.path("station_meta_input.csv"), fileEncoding = "Windows-1252")

### determine hydrological years with valid basis of data for both discharge and meteorological time series
hydrometeo_valid_hy(station_meta = station_meta)


### compute catchment approximations based on buffer approach
# compute recharge per hydrological year
#
# The following settings correspond to those used to derive the
# catchment approximations included in the ALPAKAS dataset:
# - soil moisture storage changes: included
# - snow storage changes: included
# - surface runoff: excluded
#
# Thus, recharge is calculated as:
# R = P - ET - ΔSM - ΔSWE
#
# where P is precipitation, ET is actual evapotranspiration,
# ΔSM is the change in soil moisture storage, and ΔSWE is the
# change in snow-water-equivalent storage.
if (file.exists(file.path(path_buff_approx_temp, "R_annual_1951_2024.nc"))){
  recharge_nc<- nc_open(file.path(path_buff_approx_temp, "R_annual_1951_2024.nc"))
} else{
  recharge_nc = compute_recharge_nc(
    path_buff_in = path_buff_in, 
    out_dir = path_buff_approx_temp,
    include_soil = TRUE, 
    include_snow = TRUE, 
    include_runoff = FALSE)
}

# compute buffers iteratively using topography (and tracer tests)
catchment_approx_gdf = compute_iterative_buffers(
  recharge_nc  = recharge_nc,
  station_meta = station_meta,
  path_buff_approx_temp = path_buff_approx_temp, 
  path_buff_in = path_buff_in,
  include_tracer = TRUE # use tracer tests if available to get the orientation of the catchment
  )

# save catchment approximations
st_write(catchment_approx_gdf, file.path(path_buff_approx_out, "./catchment_delineations/catchment_approx.geojson"), driver = "GeoJSON", append=FALSE)


### integrate additional attributes into metadata file
metadata_integr()

#-------------------------------------------------------------------------------