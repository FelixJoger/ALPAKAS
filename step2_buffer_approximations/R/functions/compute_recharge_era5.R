### Compute recharge from monthly ERA5-Land data
# main function returns ncdf4 object with variable "R" and time dimension hydrological_year


#-------------------------------------------------------------------------------

# --- Main Function ---

compute_recharge_nc <- function(
    path_buff_in,
    include_soil = TRUE,
    include_snow = TRUE,
    include_runoff = TRUE,
    out_dir = "./temp/"
) {
  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
  
  message("Loading NetCDF layers...")
  # Load ERA5 variables (NetCDF-backed SpatRaster)
  # Required variables: "tp", "e", "sro","sd", "swvl1","swvl2","swvl3","swvl4"
  era5_path <- file.path(path_buff_in, "meteo_data_prod/ERA5land/era5land_monthly_1950_2024.nc")
  P   <- rast(era5_path, subds = "tp") # total precipitation (m)
  AET <- rast(era5_path, subds = "e") # total evaporation (m)
  Qs  <- rast(era5_path, subds = "sro") # surface runoff (m)
  SD  <- rast(era5_path, subds = "sd") # snow depth water equivalent (m of water equivalent)
  
  # Load soil moisture layers
  # volumetric soil water layer 1-4 (m**3/m**3)
  SM_layers <- lapply(c("swvl1","swvl2","swvl3","swvl4"), function(v) rast(era5_path, subds = v))
  layer_depth <- c(0.07, 0.28-0.07, 1-0.28, 2.89-1)
  
  # Read bounding box
  bbox_rds_path <- file.path(path_buff_in, "bounding_box.RDS")
  ext_geom <- terra::ext(terra::vect(sf::st_transform(sf::st_as_sf(readRDS(bbox_rds_path)), 4326)))
  
  message("Cropping all layers to bounding box ...")
  # Crop all rasters to extent first
  P   <- crop(P, ext_geom)
  AET <- crop(AET, ext_geom)
  Qs  <- crop(Qs, ext_geom)
  SD  <- crop(SD, ext_geom)
  SM_layers <- lapply(SM_layers, function(x) crop(x, ext_geom))
  
  message("Reading time information from NetCDF...")
  nc <- nc_open(era5_path)
  time_vals  <- ncvar_get(nc, "valid_time")
  time_dates <- as.Date(as.POSIXct(time_vals, origin="1970-01-01", tz="UTC"))
  days_in_month <- days_in_month(time_dates)
  nc_close(nc)
  
  scale_factor <- 1000 * days_in_month  # m → mm and multiply by days in month

  
  message("Computing soil moisture difference...")
  if (include_soil) {

    SM_total <- init(P, NA)  # one layer per time step
    
    for (t in 1:nlyr(P)) {
      # Initialize weighted sum raster
      weighted_sum <- init(P[[1]], 0)
      
      # Loop over soil layers to get total soil moisture
      for (i in 1:length(SM_layers)) {
        layer <- SM_layers[[i]][[t]]          # single layer at timestep t
        weighted_sum <- weighted_sum + layer * layer_depth[i]
      }
      
      # Convert to mm
      SM_total[[t]] <- weighted_sum * 1000
    }

    # Difference in total soil moisture between consecutive time steps
    SM_diff <- SM_total[[2:nlyr(SM_total)]] - SM_total[[1:(nlyr(SM_total)-1)]]
    
    # Prepend NA for first timestep
    first_layer <- init(SM_total[[1]], NA)
    SM_diff <- c(first_layer, SM_diff)
    
  } else {
    # create zero raster with same dimensions and number of layers
    SM_diff <- init(P, 0)
  }
  
  message("Computing snow difference...")
  if (include_snow) {
    SD_scaled <- SD * scale_factor
    
    diff_layers <- SD_scaled[[2:nlyr(SD_scaled)]] - SD_scaled[[1:(nlyr(SD_scaled)-1)]]  
    
    # Prepend NA for first timestep
    first_layer <- init(SD_scaled[[1]], NA)
    SD_diff <- c(first_layer, diff_layers)
    
  } else {
    SD_diff <- init(P, 0)
  }
  
  if (!include_runoff) Qs <- init(P, 0)
  
  message("Computing monthly recharge...")
  # scale variables to create monthly sums in mm
  P   <- P   * scale_factor
  AET <- AET * scale_factor
  Qs  <- Qs  * scale_factor
  
  
  R <- init(P, NA)  # empty SpatRaster with same dimensions/layers
  for (i in 1:nlyr(P)) {
    R[[i]] <- P[[i]] + AET[[i]] - Qs[[i]] - SM_diff[[i]] - SD_diff[[i]]
  }

  message("Aggregating to hydrological year (Oct-Sep)...")
  time_idx <- 1:nlyr(R)
  time_dates_sub <- time_dates[time_idx]
  
  year  <- year(time_dates_sub)
  month <- month(time_dates_sub)
  hyd_year <- year + ifelse(month >= 10, 1, 0)
  hyd_year_counts <- table(hyd_year)
  
  # Keep only years with >=10 months
  valid_years <- as.integer(names(hyd_year_counts[hyd_year_counts >= 10]))

  # Subset raster to only valid hydrological years
  valid_idx <- which(hyd_year %in% valid_years)
  R_valid <- R[[valid_idx]]
  P_valid    <- P[[valid_idx]]
  AET_valid  <- AET[[valid_idx]]
  Qs_valid   <- Qs[[valid_idx]]
  SD_diff_valid   <- SD_diff[[valid_idx]]
  SM_diff_valid   <- SM_diff[[valid_idx]]
  hyd_year_valid <- hyd_year[valid_idx]
  
  # Get first and last valid hydrological years
  first_year <- min(hyd_year_valid)
  last_year  <- max(hyd_year_valid)
  
  R_annual <- tapp(R_valid, index = hyd_year_valid, fun = sum, na.rm = TRUE)
  P_annual    <- tapp(P_valid,   hyd_year_valid, sum, na.rm = TRUE)
  AET_annual  <- tapp(AET_valid, hyd_year_valid, sum, na.rm = TRUE)
  Qs_annual   <- tapp(Qs_valid,  hyd_year_valid, sum, na.rm = TRUE)
  SD_diff_annual   <- tapp(SD_diff_valid,  hyd_year_valid, sum, na.rm = TRUE)
  SM_diff_annual   <- tapp(SM_diff_valid,  hyd_year_valid, sum, na.rm = TRUE)
  

  crs(R_annual) <- "EPSG:4326"
  
  time(R_annual)   <- valid_years
  time(P_annual)   <- valid_years
  time(AET_annual) <- valid_years
  time(Qs_annual)  <- valid_years
  time(SD_diff_annual) <- valid_years
  time(SM_diff_annual) <- valid_years
  
  out_file <- sprintf("R_annual_%d_%d.nc", first_year, last_year)
  message("Saving annual NetCDF to: ", out_file)
  
  s_list <- list(
    R_annual,
    P_annual,
    AET_annual
  )
  var_names <- c("R", "P", "AET")
  
  longnames <- c(
    R = "Annual groundwater recharge (Oct–Sep)",
    P = "Precipitation sum (Oct–Sep)",
    AET = "Actual evapotranspiration sum (Oct–Sep)",
    Qs = "Surface runoff sum (Oct–Sep)",
    SD_diff = "Snow depth change (Oct–Sep)",
    SM_diff = "Soil moisture change (Oct–Sep)"
  )
  
  if (include_runoff) {
    s_list[[length(s_list) + 1]] <- Qs_annual
    var_names <- c(var_names, "Qs")
  }
  
  if (include_snow) {
    s_list[[length(s_list) + 1]] <- SD_diff_annual
    var_names <- c(var_names, "SD_diff")
  }
  
  if (include_soil) {
    s_list[[length(s_list) + 1]] <- SM_diff_annual
    var_names <- c(var_names, "SM_diff")
  }
  

  s <- do.call(sds, s_list)
  names(s) <- var_names
  units(s) <- rep("mm", length(s_list))
  longnames(s) <- longnames[var_names]
  
  writeCDF(s, filename = file.path(out_dir, out_file), 
           timename = "hydrological_year", overwrite=TRUE)
  
  
  message("Recharge NetCDF saved to: ", out_file)
  R <- nc_open(file.path(out_dir, out_file))
  
  return(R)
}
