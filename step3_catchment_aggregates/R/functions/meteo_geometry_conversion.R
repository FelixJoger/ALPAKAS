### conversion of gridded meteorological data products to yearly files of uniform coordinate reference system, units, and cropped to study region
# data products have to be downloaded to the respective folders prior running the script
# download links are provided along with the code sections for the respective data products

#-------------------------------------------------------------------------------

# function to rasterize, for each date, the point data for the respective variable (SAFRAN)
make_layer <- function(d) {
  sub <- safran_df %>%
    filter(.data$date == d, !is.na(.data[[var_v]])) %>%
    transmute(x = .data$X, y = .data$Y, val = .data[[var_v]])
  
  r <- r_template
  if (nrow(sub) == 0) {
    names(r) <- format(d, "%Y%m%d")
    return(r)
  }
  
  # terra::cellFromXY expects cbind(x, y)
  cells <- cellFromXY(r, as.matrix(sub[, c("x", "y")]))
  vals <- rep(NA_real_, ncell(r))
  vals[cells] <- sub$val
  values(r) <- vals
  names(r) <- format(d, "%Y%m%d")
  r
}


### main function to perform geometry conversion of meteorological products
meteo_geom_convers <- function() {
  
  # define file path for original data products used as input
  # data products have to be downloaded to the respective folders prior running the script
  path_input_data <- paste0(path_catch_agg_in, "meteo_data_prod/")
  
  # # define file path for processed geospatial data products
  path_proc_data_prod <- paste0(path_catch_agg_temp, "proc_data_prod/")
  
  # define target EPSG (ETRS89 / LAEA Europe)
  crs_proj <- 3035
  
  # read bounding box to crop data products
  bounding_box <- readRDS(file = paste0(path_catch_agg_in, "bounding_box.RDS"))
  
  #-----------------------------------------------------------------------------
  ### load meteorological raster files, transform them, crop them to the bounding box, and save them as yearly NetCDF files
  dir.create(path_proc_data_prod)
  
  
  ### national datasets ###
  
  ### AT: SPARTACUS 2.1 ###
  # https://data.hub.geosphere.at/dataset/spartacus-v2-1d-1km
  dir.create(paste0(path_proc_data_prod, "spartacus/"))
  
  # for each variable and each month of each year, read data files
  var <- c("RR", "TN", "TX")
  years <- 1961:2024
  
  for (v in 1:length(var)){
    var_v <- var[v]
    print(var_v)
    dir.create(paste0(path_proc_data_prod, "spartacus/", var_v))
    
    for (y in 1:length(years)){
      year_y <- years[y]
      print(year_y)
      
      file_names <- list.files(path = paste0(path_input_data, "SPARTACUS/", var_v, "/", year_y, "/"), pattern = "\\.nc$", full.names = TRUE)  
      var_v_y_rast <- rast(file_names)
      var_v_y_rast_proj <- project(var_v_y_rast, paste0("EPSG:", crs_proj))
      var_v_y_bbox <- crop(var_v_y_rast_proj, vect(bounding_box))
      # divide raster values by a scale factor of 10
      var_v_y_bbox <- var_v_y_bbox / 10
      
      writeCDF(var_v_y_bbox, paste0(path_proc_data_prod, "spartacus/", var_v, "/spartacus_", var_v, "_raster_bbox_", year_y, ".nc"))
    }
  }
  
  ### CH: MeteoSwiss ###
  # provided upon request
  dir.create(paste0(path_proc_data_prod, "meteoswiss/"))
  
  # for each variable and each month of each year, read data files
  var <- c("RhiresD", "TabsD", "TminD", "TmaxD")
  
  for (v in 1:length(var)){
    var_v <- var[v]
    print(var_v)
    dir.create(paste0(path_proc_data_prod, "meteoswiss/", var_v))
    
    # conditional, depending on variable
    if (var_v %in% c("RhiresD", "TabsD")){
      years <- 1961:2024
    } else if (var_v %in% c("TminD", "TmaxD")){
      years <- 1971:2024
    }
    if (var_v %in% c("RhiresD")){
      key_name <- "h"
    } else if (var_v %in% c("TminD", "TmaxD", "TabsD")){
      key_name <- "r"
    }
    
    for (y in 1:length(years)){
      year_y <- years[y]
      print(year_y)
      
      var_v_y_rast <- rast(paste0(path_input_data, "MeteoSwiss/", var_v, "/",
                                  var_v, "_ch01", key_name, ".swiss.lv95_", year_y, "01010000_", year_y, "12310000.nc"))
      var_v_y_rast_proj <- project(var_v_y_rast, paste0("EPSG:", crs_proj))
      var_v_y_bbox <- crop(var_v_y_rast_proj, vect(bounding_box))
      
      writeCDF(var_v_y_bbox, paste0(path_proc_data_prod, "meteoswiss/", var_v, "/meteoswiss_", var_v, "_raster_bbox_", year_y, ".nc"))
    }
  }
  
  ### DE: HYRAS 6.0 ###
  # https://opendata.dwd.de/climate_environment/CDC/grids_germany/daily/hyras_de/
  dir.create(paste0(path_proc_data_prod, "hyras/"))
  
  # for each variable, read yearly data files
  var <- c("pr", "tas", "tasmin", "tasmax")
  
  for (v in 1:length(var)){
    var_v <- var[v]
    print(var_v)
    dir.create(paste0(path_proc_data_prod, "hyras/", var_v))
    
    # conditional, depending on variable
    if (var_v %in% c("pr", "tas")){
      years <- 1951:2025
    } else if (var_v %in% c("tasmin", "tasmax")){
      years <- 1951:2024
    }
    
    for (y in 1:length(years)){
      year_y <- years[y]
      print(year_y)
  
      var_v_y_rast <- rast(paste0(path_input_data, "HYRASv6.0/", var_v, "/", var_v, "_hyras_1_", year_y, "_v6-0_de.nc"))
      var_v_y_bbox <- crop(var_v_y_rast, vect(bounding_box))
      # projecting is required as otherwise problems with illegal characters in file
      var_v_y_bbox <- project(var_v_y_bbox, paste0("EPSG:", crs_proj))
      
      writeCDF(var_v_y_bbox, paste0(path_proc_data_prod, "hyras/", var_v, "/hyras_", var_v, "_raster_bbox_", year_y, ".nc"))
    }
  }
  
  ### FR: SAFRAN Sim2 ###
  # https://meteo.data.gouv.fr/datasets/donnees-changement-climatique-sim-quotidienne/
  dir.create(paste0(path_proc_data_prod, "safran/"))
  
  # list and read all .csv.gz files
  file_names <- list.files(path = paste0(path_input_data, "SAFRAN"),
                           pattern = "\\.csv\\.gz$", full.names = TRUE)
  
  for (l in 1:length(file_names)){
    
    # read data files separately
    safran_df <- read_delim(file_names[l], delim = ";")
  
    # relevant variables: PRENEI (solid precipitation), PRELIQ (liquid precipitation), T (mean temperature),
    # TINF_H (minimum temperature of hourly values), TSUP_H (maximum temperature of hourly values)
    var = c("PRENEI", "PRELIQ", "T", "TINF_H", "TSUP_H")
    
    # only select relevant columns
    safran_df <- safran_df %>%
      dplyr::select(LAMBX, LAMBY, DATE, all_of(var)) %>%
      rename(X = LAMBX,
             Y = LAMBY,
             date = DATE)
    
    # ymd() and as.Date() are slow due to the large amount of data
    # avoid slow coercion by first formatting once per unique value and then match values
    u   <- unique(safran_df$date)
    uc  <- if (is.numeric(u)) sprintf("%08d", u) else as.character(u)
    ud  <- as.Date(uc, format = "%Y%m%d")
    safran_df$date <- ud[match(safran_df$date, u)]
    
    # convert hm to m
    safran_df <- safran_df %>%
      mutate(
        X = X * 100,
        Y = Y * 100)
    
    # determine unique x & y values
    xs <- sort(unique(safran_df$X))
    ys <- sort(unique(safran_df$Y))
    # extract grid spacing
    dx <- median(diff(xs))
    dy <- median(diff(ys))
    # define grid extension (grid points are always centered within grid cells)
    ext <- terra::ext(min(xs) - dx/2, max(xs) + dx/2,
                      min(ys) - dy/2, max(ys) + dy/2)
    r_template <- rast(ncols = (max(xs) - min(xs)) / dx + 1, nrows = (max(ys) - min(ys)) / dy + 1, extent = ext, crs = "EPSG:27572")
    dates <- sort(unique(safran_df$date))
    
    for (v in 1:length(var)){
      var_v <- var[v]
      print(var_v)
      dir.create(paste0(path_proc_data_prod, "safran/", var_v))
      
      years <- unique(year(safran_df$date))
    
      for (y in 1:length(years)){
    
        year_y <- years[y]
        print(year_y)
    
        dates_year <- dates[year(dates) == year_y]
        # apply function to rasterize point data and write raster layers into list
        var_v_y_rast <- lapply(dates_year, make_layer)
        # stack data
        var_v_y_rast <- rast(var_v_y_rast)
        # attach time coordinate
        time(var_v_y_rast) <- dates_year
        var_v_y_rast_proj <- project(var_v_y_rast, paste0("EPSG:", crs_proj))
        var_v_y_bbox <- crop(var_v_y_rast_proj, vect(bounding_box))
        
        writeCDF(var_v_y_bbox, paste0(path_proc_data_prod, "safran/", var_v, "/safran_", var_v, "_raster_bbox_", year_y, ".nc"))
      }
    }
  }
  
  ### SI: SLOCLIM ###
  # https://zenodo.org/records/4108543
  dir.create(paste0(path_proc_data_prod, "sloclim/"))
  
  # for each variable, read data files
  var <- c("pcp", "tmax_h", "tmin_h")
  
  for (v in 1:length(var)){
    var_v <- var[v]
    print(var_v)
    dir.create(paste0(path_proc_data_prod, "sloclim/", var_v))
    
    var_v_rast <- rast(paste0(path_input_data, "SLOCLIM/sloclim_", var_v, ".nc"))
    var_v_rast_proj <- project(var_v_rast, paste0("EPSG:", crs_proj))
    var_v_bbox <- crop(var_v_rast_proj, vect(bounding_box))
    
    # assign dates to raster file
    time(var_v_bbox) <- seq(from = as.Date("1950-01-01"), to = as.Date("2018-12-31"), by = "day")
    years <- unique(year(time(var_v_bbox)))
    
    for (y in 1:length(years)) {
      
      year_y <- years[y]
      print(year_y)
      
      # split raster file to yearly files
      var_v_y_bbox <- var_v_bbox[[year(time(var_v_bbox)) == year_y]]
      
      writeCDF(var_v_y_bbox, paste0(path_proc_data_prod, "sloclim/", var_v, "/sloclim_", var_v, "_raster_bbox_", year_y, ".nc"))
      
    }
  }
  

  ### E-OBS ###
  # https://surfobs.climate.copernicus.eu/dataaccess/access_eobs.php
  dir.create(paste0(path_proc_data_prod, "eobs/"))
  
  # for each variable, read data files
  var <- c("rr", "tg", "tn", "tx")
  for (v in 1:length(var)){
    var_v <- var[v]
    print(var_v)
    dir.create(paste0(path_proc_data_prod, "eobs/", var_v))
    
    var_v_raster <- rast(paste0(path_input_data, "EOBS/", var_v, "_ens_mean_0.1deg_reg_v31.0e.nc"))
    
    # perform processing on yearly basis
    time <- time(var_v_raster)
    years_floor <- unique(floor_date(time, "year"))
    
    for (y in 1:length(years_floor)) {
      
      years_floor_y <- years_floor[y]
      year_y <- year(years_floor_y)
      print(year_y)
      
      idx_year <- which(floor_date(time, "year") == years_floor_y)
      
      var_v_y <- var_v_raster[[idx_year]]
      var_v_y_proj <- project(var_v_y, paste0("EPSG:", crs_proj))
      # divide raster values by a scale factor depending on variable
      if (var_v == "rr"){
        var_v_y_proj <- var_v_y_proj / 10
      } else {
        var_v_y_proj <- var_v_y_proj / 100
      }
      var_v_y_proj_bbox <- crop(var_v_y_proj, bounding_box)
      
      writeCDF(var_v_y_proj_bbox, paste0(path_proc_data_prod, "eobs/", var_v, "/eobs_", var_v, "_raster_bbox_", year_y, ".nc"))
    }
  }
  
  
  ### ERA5-Land ###
  # https://cds.climate.copernicus.eu/datasets/reanalysis-era5-land?tab=overview
  ### hourly data
  dir.create(paste0(path_proc_data_prod, "era5land/"))
  
  # for each variable, read data files
  var <- c("2m_temperature","total_precipitation")
  years <- 1950:2025
  
  for (v in 1:length(var)){
    var_v <- var[v]
    print(var_v)
    dir.create(paste0(path_proc_data_prod, "era5land/", var_v))
    
    # perform processing on yearly basis
    for (y in 1:length(years)){
      year_y <- years[y]
      print(year_y)
      
      # list all files for the respective year
      file_names <- list.files(path = paste0(path_input_data, "ERA5land/", var_v, "/", year_y, "/"),
                               pattern = "\\.nc$", full.names = TRUE) 
      var_v_y_rast <- rast(file_names)
      # filter files for variable as some files contain data for several variables
      if (var_v == "2m_temperature"){
        var_v_y_rast <- var_v_y_rast[[grep("t2m", names(var_v_y_rast))]]
      } else if (var_v == "total_precipitation"){
        var_v_y_rast <- var_v_y_rast[[grep("tp", names(var_v_y_rast))]]
      }
      var_v_y_rast_proj <- project(var_v_y_rast, paste0("EPSG:", crs_proj))
      var_v_y_bbox <- crop(var_v_y_rast_proj, vect(bounding_box))
      # convert units
      if (var_v == "2m_temperature"){
        # convert Kelvin to C°
        var_v_y_bbox <- var_v_y_bbox - 273.15  
      } else if (var_v == "total_precipitation"){
        # convert depth in m to depth in mm
        var_v_y_bbox <- var_v_y_bbox * 1000
      }
  
      writeCDF(var_v_y_bbox, paste0(path_proc_data_prod, "era5land/", var_v, "/era5land_", var_v, "_raster_bbox_", year_y, ".nc"))
    }
  }
      
}