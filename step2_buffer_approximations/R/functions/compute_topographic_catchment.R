### Compute topographic catchments using whitebox
# main function returns sf object with one topographic geometry per AlpAKaS_ID

#-------------------------------------------------------------------------------   

# --- Functions ---

# Function to reproject the DEM to EPSG:3035 for further usage
process_copernicus_dem <- function(path_buff_in,
                                  path_dem_3035,
                                  crs_proj = 3035,
                                  bbox_file = "bounding_box.RDS") {
  
  
  # define file paths
  path_input_data <- file.path(path_buff_in, "static_data_prod/")
  path_proc_data_prod <- dirname(path_dem_3035)
  
  # ensure output directory exists
  if (!dir.exists(path_proc_data_prod)) {
    dir.create(path_proc_data_prod, recursive = TRUE)
  }
  
  # read bounding box
  bounding_box <- readRDS(file = file.path(path_buff_in, bbox_file))
  
  # locate DEM files
  files_elev <- list.files(
    file.path(path_input_data, "Copernicus_GLO30_DEM/COP30_elevation/"),
    pattern = "\\.tif$",
    full.names = TRUE
  )
  
  if (length(files_elev) == 0) {
    stop("No DEM files found in expected directory.")
  }
  
  # read rasters
  cop_dem_elev <- lapply(files_elev, terra::rast)
  
  # mosaic (mean aggregation)
  cop_dem_elev_merged <- do.call(terra::mosaic, c(cop_dem_elev, fun = mean))
  
  # project
  cop_dem_elev_merged_proj <- terra::project(
    cop_dem_elev_merged,
    paste0("EPSG:", crs_proj)
  )
  
  # crop to bounding box
  cop_dem_elev_bbox <- terra::crop(
    cop_dem_elev_merged_proj,
    terra::vect(bounding_box)
  )
  
  
  # write result
  terra::writeRaster(cop_dem_elev_bbox, filename = path_dem_3035, overwrite = TRUE)
  
}

 
  
# Function to compute the topographic catchment iteratively for one station using whitebox
single_topographic_catchment <- function(row, path_dem, visualize, coastline, out_dir) {  
  # row: A single row from a dataframe or spatial dataframe containing the gauging station coordinates and ids.
  # path_dem: The file path to the Digital Elevation Model (DEM) data, which is used to calculate the topographic catchment.
  # visualize: If TRUE, visualizations of the catchment are created; if FALSE, they are not
  
  alpakas_id <- row$AlpAKaS_ID
  cat("\n", strrep("=", 20), " ", row$spring_name, " ",alpakas_id, " ", strrep("=", 20), "\n")
  
  
  ## Check if the topographic catchment has already been delineated
  target_dir <- file.path(out_dir, "topographic_catchments")
  
  if (!dir.exists(target_dir)) {
    dir.create(target_dir, recursive = TRUE)
  }
  
  pattern <- paste0("^topographic_catchment_", alpakas_id, "_snap.*_acc.*\\.shp$")
  
  shp_files <- list.files(target_dir, pattern = pattern, full.names = TRUE)

  # Early exit if file exists
  if (length(shp_files) > 0 ) { 
    
    # Get modification times
    files_info <- file.info(shp_files)
    
    # Pick the most recently modified file
    recent_file <- shp_files[which.max(files_info$mtime)]
    
    message("File exists, loading cached watershed...")
    
    # Read shapefile
    catchment <- sf::st_read(recent_file, quiet = TRUE)
    
    crs_info <- sf::st_crs(catchment)
    
    if (is.na(crs_info)) {
      catchment <- sf::st_set_crs(catchment, 3035)
    } else if (is.na(crs_info$epsg) || crs_info$epsg != 3035) {
      catchment <- sf::st_transform(catchment, 3035)
    }
    
    # Rename abbreviated columns back to original names
    if (!"sprng_n" %in% names(catchment)) {
      catchment$sprng_n <- row$spring_name
    }

    catchment <- catchment %>%
      dplyr::rename(
        local_station_ID  = lcl__ID,
        spring_name       = sprng_n,
        accumulation_thr  = accmlt_,
        snap_distance     = snp_dst,
        snapped_x         = snppd_x,
        snapped_y         = snppd_y,
        snapped_accum     = snppd_c
      )
    
    # Return the catchment
    return(catchment)
  }

  
  ## Create paths
  make_path <- function(target_dir, alpakas_id, name, extension) {
    file.path(target_dir, paste0(name, "_", alpakas_id, ".", extension))
  }
  
  pourpoint_path <- make_path(target_dir, alpakas_id, "pourpoints", "shp")
  dem_crop_path  <- make_path(target_dir, alpakas_id, "temp_dem_cropped", "tif")
  breach_path    <- make_path(target_dir, alpakas_id, "dem_breached", "tif")
  filled_path    <- make_path(target_dir, alpakas_id, "dem_filled", "tif")
  flowacc_path   <- make_path(target_dir, alpakas_id, "D8FA", "tif")
  pointer_path   <- make_path(target_dir, alpakas_id, "D8pointer", "tif")
  hillshade_path <- make_path(target_dir, alpakas_id, "brush_hillshade","tif")
  stream_path    <- make_path(target_dir, alpakas_id, "raster_streams","tif")
  snapped_pp_path <- make_path(target_dir, alpakas_id, "snappedpp", "shp")
  
  # save pourpoints (gauging station coordinates) as shapefile
  pourpoints <- row %>%
    dplyr::select(local_station_ID, geometry) 
  names(pourpoints)[names(pourpoints) == "local_station_ID"] <- "lcl__ID"
  st_write(pourpoints,pourpoint_path,append=FALSE)
  
  
  ## Preprocess DEM
  # read DEM
  dem <- rast(path_dem)
  
  if (!file.exists(dem_crop_path)) {
    # Create a buffer around the point (size depends on the buffer radius) and crop the DEM
    if (is.na(row$buffer_radius_km)) {
      buffer_dist <- 10000 # Default value if buffer_radius is NA
    } else {
      buffer_dist <- min(row$buffer_radius_km *1000*10)
    }

    buffer <- st_buffer(row$geometry, dist = buffer_dist) # Distance in meters
    bbox <- st_bbox(buffer)
    buffer_extent <- extent(bbox$xmin, bbox$xmax, bbox$ymin, bbox$ymax)
    dem_masked <- crop(dem, buffer_extent)
    dem_masked <- mask(dem_masked, coastline)
    
    writeRaster(dem_masked, dem_crop_path, overwrite = TRUE)
  }

  if (!file.exists(breach_path)) {
    # Breach depressions
    wbt_breach_depressions_least_cost(
      dem = dem_crop_path,
      output = breach_path,
      dist = 5, # default 5 for a resolution of 30 m
      fill = TRUE
    )
  }

  if (!file.exists(filled_path)) {
    # Fill depressions and sinks
    wbt_fill_depressions_wang_and_liu(
      dem = breach_path,
      output = filled_path
    )
  }
  
  if (!file.exists(flowacc_path)) {
    # Generate D8 flow accumulation 
    #(each cell is routed to one of the 8 neighboring cells based on direction of steepest descent)
    wbt_d8_flow_accumulation(
      input = filled_path,
      output = flowacc_path
    )
  }
  
  if (!file.exists(pointer_path)) {
    # Generate D8 pointer file (direction) using clockwise, base-2 numeric naming convention: (1,2,4,8,18,32,62,128)
    wbt_d8_pointer(
      dem = filled_path,
      output = pointer_path
    )
  }
  
  if (visualize){
    # Create hillshade (only for visualization)
    if (!file.exists(hillshade_path)) {
      wbt_hillshade(dem = dem_crop_path,
                    output = hillshade_path,
                    azimuth = 115)
    }
    hillshade <- raster("topographic_catchments/brush_hillshade.tif")
  }
  

  ## Iteration
  # Set iteration parameters
  n <- 0 # Initialize iteration counter
  snap_distance <- 180 # Maximum search radius (meters) around station for locating a stream pixel
  accumulation <- 1000 # Flow accumulation threshold used to define stream cells (in number of contributing cells)
  river_catch <- FALSE # Flag indicating whether a too large (probably river) catchment was incorrectly delineated instead of a spring catchment
  failed_snap <- FALSE # # Flag indicating whether no stream was found within snap_distance of the spring location
  
  # Cell resolution (length of one side, in metres)
  res_x <- res(dem)[1]  # x-direction
  res_y <- res(dem)[2]  # y-direction
  
  expected_area <- row$estimated_recharge_area_km2 # estimated area of the spring catchment based on water balance
  expected_min_accumulation <- (expected_area * 1e6) / (res_x * res_y) # Minimum expected flow accumulation (in cells), derived from estimated catchment area and cell size
  
  # Continue Iteration until a topographic catchment has been delineated which is not larger than the
  
  while(n==0 || river_catch && n<=10 || failed_snap && n<=10 ){
    print(glue("\nIteration : {n}",
               "\ncurrent snap_distance {snap_distance}, current accumulation {accumulation}"))
    
    if (n>0){
      accumulation <- accumulation/2
    }
     
    file_name <- file.path(target_dir,glue("topographic_catchment_{alpakas_id}_snap{snap_distance}_acc{accumulation}.shp"))
    

    # extract streams from flow accumulation grid
    wbt_extract_streams(flow_accum = flowacc_path,
                        output = stream_path,
                        threshold = accumulation) #6000
    
    
    # Snap pour points to closest stream
    wbt_jenson_snap_pour_points(pour_pts = pourpoint_path,
                                streams = stream_path,
                                output = snapped_pp_path,
                                snap_dist = snap_distance) #careful with this! Know the units of your data
    
    
    ## Verify snap point
    # Get accumulation at pour point to assure it has not snapped not a river nearby
    snapped_pp <- st_read(snapped_pp_path, quiet = TRUE)
    flow_accum <- rast(flowacc_path)
    snapped_vect <- vect(snapped_pp)
    accum_value <- extract(flow_accum, snapped_vect)[,2]
    
    print(glue("Expected minimum accumualtion: {expected_min_accumulation}"))
    print(glue("Current accumulation: {accum_value}"))
    
    # if the current accumulation is lower than the threshold, no river has been found: risk of too small catchment
    if (accum_value < accumulation && n<10){
      failed_snap <- TRUE
      river_catch <- FALSE
      n <- n + 1
      
      print("New iteration with lower accumulation, no snapping river found")
      next
    }
    # if the accumulation is >3* expected accumulation (based on water balance), the algorithm has probably snapped to a river
    if (accum_value > expected_min_accumulation*3 && n<10 ){ 
      river_catch <- TRUE
      failed_snap <- FALSE
      n <- n+1
      
      print("New iteration with lower accumulation, potential river catchment")
      next
    }
      
    
    ## Create watershed for the snapped pour point
    # Create watershed delineation as raster
    wbt_watershed(d8_pntr = pointer_path,
                  pour_pts = snapped_pp_path,
                  output = file.path(target_dir,"brush_watersheds.tif"))
    
    ws <- raster(file.path(target_dir,"brush_watersheds.tif"))
    
    #visualize catchments
    if (visualize==TRUE){
      pp_sf <- st_as_sf(pp)
      
      p <- tm_shape(hillshade)+
        tm_raster(col.scale = tm_scale_continuous(values = "brewer.greys"), col.legend = tm_legend_hide()) +
        tm_shape(ws) +
        tm_raster(col.scale = tm_scale_categorical(), legend.show = TRUE, col_alpha = 0.5) +
        tm_shape(pp_sf) +
        tm_dots(col = "red")
      print(p)
    }
    
    # convert to polygon shapefile
    wsshape <- st_as_stars(ws) %>% st_as_sf(merge = T)
    
    # merge catchment polygons if multiple polygons have been created for one snappoint
    if (nrow(wsshape) > 1) {
      buffered_shape <- st_buffer(wsshape, dist = 0.01) # distance for merging
      combined_shape <- st_union(buffered_shape)
      combined_shape <- st_sf(geometry = combined_shape)
    } else {
      combined_shape <- wsshape
    }

    n <- n + 1
    river_catch <- FALSE 
    failed_snap <- FALSE

  }
  
  ## save shapefile
  combined_shape_sf <- combined_shape %>%
    mutate(
      local_station_ID = row$local_station_ID,
      spring_name = row$spring_name,
      accumulation_thr = accumulation,
      snap_distance = snap_distance,
      snapped_x = st_coordinates(snapped_pp)[1, "X"],
      snapped_y = st_coordinates(snapped_pp)[1, "Y"],
      snapped_accum = accum_value
    )
  
  combined_shape_sf <- sf::st_set_crs(combined_shape_sf, 3035)
  st_write(combined_shape_sf, file_name, append = FALSE)
  
  return(combined_shape_sf)
}

#-------------------------------------------------------------------------------   

# --- Main Function ---

compute_topographic_catchments <- function(
    all_data_sf,
    path_buff_in,
    path_buff_approx_temp,
    visualize = FALSE
) {
  
  # Reproject DEM to EPSG:3035
  dem_3035_path <- file.path(path_buff_approx_temp, "proc_data_prod/cop_dem_elev_raster_bbox.tif")
  if (!file.exists(dem_3035_path)){
    process_copernicus_dem(path_buff_in, dem_3035_path)
  }
  
  # Read coastline
  path_coastline <- file.path(path_buff_in,"static_data_prod/EEA_Coastline/EEA_Coastline_20170228.shp")
  coastline <- sf::st_read(path_coastline, quiet = TRUE)

  # Process row-wise catchments
  rows_spring_station <- all_data_sf %>%
    rowwise() %>%
    mutate(
      result = list({
        catchment <- single_topographic_catchment(
          cur_data(),
          dem_3035_path,
          visualize = visualize,
          coastline = coastline, out_dir = path_buff_approx_temp
        )
        
        # Add identifiers
        catchment$country_code <- country_code
        catchment$AlpAKaS_ID <- AlpAKaS_ID
        
        # Set geometry column
        sf::st_geometry(catchment) <- "topographic_geometry"
        
        catchment
      })
    ) %>%
    ungroup()
  
  # Flatten list-column into sf
  rows_spring_station_flat <- dplyr::bind_rows(rows_spring_station$result)
  
  # Select and format output
  output_gdf <- rows_spring_station_flat %>%
    dplyr::select(
      country_code,
      AlpAKaS_ID,
      local_station_ID,
      spring_name,
      accumulation_thr,
      snap_distance,
      snapped_x,
      snapped_y,
      snapped_accum,
      topographic_geometry
    )
  
  # Set geometry and CRS
  output_gdf <- sf::st_set_geometry(output_gdf, "topographic_geometry")
  output_gdf <- sf::st_set_crs(output_gdf, sf::st_crs(all_data_sf))
  filename <- file.path(path_buff_approx_temp, "all_topographic_catchments.geojson")
  st_write(output_gdf, filename,driver = "GeoJSON", append=FALSE, quiet=TRUE)
  
  return(output_gdf)
}
