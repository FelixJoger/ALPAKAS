### conversion of geospatial data products to uniform coordinate reference system, units, and cropped to study region
# data products have to be downloaded to the respective folders prior running the script
# download links are provided along with the code sections for the respective data products

#-------------------------------------------------------------------------------

### main function to perform geometry conversion of geospatial products
static_attr_geom_convers <- function() {

  # define file path for original data products used as input
  # data products have to be downloaded to the respective folders prior running the script
  path_input_data <- paste0(path_catch_agg_in, "static_data_prod/")
  
  # define file path for processed geospatial data products
  path_proc_data_prod <- paste0(path_catch_agg_temp, "proc_data_prod/")
  
  # define target EPSG (ETRS89 / LAEA Europe)
  crs_proj <- 3035
  
  # read bounding box to clip data products
  bounding_box <- readRDS(file = paste0(path_catch_agg_in, "bounding_box.RDS"))
  
  #-----------------------------------------------------------------------------
  ### load spatial geodata products, transform them, clip them to the bounding box, and save them as RDS spatial files or TIF-files
  
  
  ### topography and positioning ###
  
  ### Copernicus DEM ###
  # https://portal.opentopography.org/raster?opentopoID=OTSDEM.032021.4326.3; https://doi.org/10.5270/esa-c5d3d65
  ### elevation
  files_elev <- list.files(paste0(path_input_data, "Copernicus_GLO30_DEM/COP30_elevation/"), pattern = "\\.tif$", full.names = TRUE)
  cop_dem_elev <- lapply(files_elev, rast)
  cop_dem_elev_merged <- do.call(mosaic, c(cop_dem_elev, fun = mean))
  cop_dem_elev_merged_proj <- project(cop_dem_elev_merged, paste0("EPSG:", crs_proj))
  cop_dem_elev_bbox = crop(cop_dem_elev_merged_proj, vect(bounding_box))
  writeRaster(cop_dem_elev_bbox, filename = paste0(path_proc_data_prod, "cop_dem_elev_raster_bbox.tif"), overwrite = TRUE)
  
  ### slope
  files_slope <- list.files(paste0(path_input_data, "Copernicus_GLO30_DEM/COP30_slope/"), pattern = "\\.tif$", full.names = TRUE)
  cop_dem_slope <- lapply(files_slope, rast)
  cop_dem_slope_merged <- do.call(mosaic, c(cop_dem_slope, fun = mean))
  cop_dem_slope_merged_proj <- project(cop_dem_slope_merged, paste0("EPSG:", crs_proj))
  cop_dem_slope_bbox = crop(cop_dem_slope_merged_proj, vect(bounding_box))
  writeRaster(cop_dem_slope_bbox, filename = paste0(path_proc_data_prod, "cop_dem_slope_raster_bbox.tif"), overwrite = TRUE)
  
  ### aspect (northness & eastness)
  files_aspect <- list.files(paste0(path_input_data, "Copernicus_GLO30_DEM/COP30_aspect/"), pattern = "\\.tif$", full.names = TRUE)
  cop_dem_aspect <- lapply(files_aspect, rast)
  cop_dem_aspect_merged <- do.call(mosaic, c(cop_dem_aspect, fun = mean))
  cop_dem_aspect_merged_proj <- project(cop_dem_aspect_merged, paste0("EPSG:", crs_proj))
  cop_dem_aspect_raster_bbox = crop(cop_dem_aspect_merged_proj, vect(bounding_box))
  # convert degrees to radians in order to be able to calculate averages using northness and eastness
  cop_dem_aspect_raster_bbox_rad <- cop_dem_aspect_raster_bbox * pi / 180
  cop_dem_slope_raster_bbox_rad <- cop_dem_slope_raster_bbox * pi / 180
  # calculate northness and eastness weighted by slope to account for high uncertainty for areas with low slope
  # -> see: https://doi.org/10.1038/sdata.2018.40
  cop_dem_northn_raster_bbox <- cos(cop_dem_aspect_raster_bbox_rad) * sin(cop_dem_slope_raster_bbox_rad)
  cop_dem_eastn_raster_bbox  <- sin(cop_dem_aspect_raster_bbox_rad) * sin(cop_dem_slope_raster_bbox_rad)
  writeRaster(cop_dem_northn_raster_bbox, filename = paste0(path_proc_data_prod, "cop_dem_northn_raster_bbox.tif"), overwrite = TRUE)
  writeRaster(cop_dem_eastn_raster_bbox, filename = paste0(path_proc_data_prod, "cop_dem_eastn_raster_bbox.tif"), overwrite = TRUE)
  

  ### EU-Hydro ###
  # https://land.copernicus.eu/en/products/eu-hydro/eu-hydro-river-network-database
  # load specific geopackage layers of all relevant major river catchments and merge them
  gpkg_names <- c("danube", "elbe", "loire", "po", "rhine", "rhone")
  gpkg_paths <- paste0(path_input_data, "EU_hydro_gpkg_eu/euhydro_", gpkg_names, "_v013_GPKG/euhydro_",
                       gpkg_names, "_v013_GPKG/euhydro_", gpkg_names, "_v013.gpkg")
  layers_list <- map(gpkg_paths, ~ st_read(.x, layer = "River_Net_l", quiet = TRUE))
  eu_hydro_shape <- bind_rows(layers_list)
  # drop z dimension
  eu_hydro_shape <- st_zm(eu_hydro_shape, drop = TRUE, what = "ZM")
  eu_hydro_shape_bbox <- st_intersection(eu_hydro_shape, bounding_box)
  saveRDS(eu_hydro_shape_bbox, file = paste0(path_proc_data_prod, "eu_hydro_shape_bbox.RDS"))
  
  
  ### MOHP ###
  # https://www.hydroshare.org/resource/0d6999591fb048cab5ab71fcb690eadb/
  dir.create(paste0(path_catch_agg_temp, "proc_data_prod/mohp/"))
  mohp_prod <- c("divide_stream_distance", "lateral_position", "stream_distance")
  mohp_var <- c("dsd", "lp", "sd")
  
  for (m in 1:length(mohp_prod)){
    mohp_prod_m <- mohp_prod[m]
    mohp_var_m <- mohp_var[m]
    
    for (n in 1:9){
      mohp_raster_m_n <- rast(paste0(path_input_data, "macro_mohp_feature/", mohp_prod_m, "/mohp_europe_europemainland_", mohp_var_m, "_streamorder", n, "_30m.tif"))
      crs(mohp_raster_m_n) <- "EPSG:3035"
      mohp_raster_m_n_bbox = crop(mohp_raster_m_n, vect(bounding_box))
      writeRaster(mohp_raster_m_n_bbox, filename = paste0(path_proc_data_prod, "mohp/mohp_", mohp_var_m, n, "_raster_bbox.tif"))
    }
  }
  
  
  ### Koeppen-Geiger climate zones ###
  # https://www.gloh2o.org/koppen/
  kg_raster <- rast(paste0(path_input_data, "koppen_geiger_tif/1991_2020/koppen_geiger_0p00833333.tif"))
  # define larger extend of bounding box for initial cropping and transform to box to product crs
  bbox <- st_bbox(bounding_box)
  bounding_box_large <- st_polygon(list(matrix(c(
    bbox[1] - 100000, bbox[2] - 100000,
    bbox[1] - 100000, bbox[4] + 100000,
    bbox[3] + 100000, bbox[4] + 100000,
    bbox[3] + 100000, bbox[2] - 100000,
    bbox[1] - 100000, bbox[2] - 100000), ncol = 2, byrow = TRUE)))
  bounding_box_large <- st_sfc(bounding_box_large, crs = crs_proj) %>%
    st_transform(bounding_box_large, crs = 4326)
  kg_raster_bboxl <- crop(kg_raster, vect(bounding_box_large))
  kg_raster_proj_bboxl <- project(kg_raster_bboxl, paste0("EPSG:", crs_proj), method = "near") # make sure categorical values don´t get averaged
  kg_raster_proj_bbox <- crop(kg_raster_proj_bboxl, vect(bounding_box))
  writeRaster(kg_raster_proj_bbox, filename = paste0(path_proc_data_prod, "kg_raster_bbox.tif"))
  
  
  ### hydrogeology ###
  
  ### MEDKAM ###
  # https://www.whymap.org/whymap/EN/Maps_Data/Medkam/medkam_node_en.html
  medkam_shape <- st_read(paste0(path_input_data, "WHYMAP_MEDKAM/shp/MEDKAM_Karst__v1_poly.shp")) %>% 
    st_transform(medkam_shape, crs = crs_proj)
  medkam_shape_bbox <- st_intersection(medkam_shape, bounding_box)
  saveRDS(medkam_shape_bbox, file = paste0(path_proc_data_prod, "medkam_shape_bbox.RDS"))
  
  
  ### IHME1500 ###
  # https://www.bgr.bund.de/DE/Themen/Grundwasser/Projekte/Flaechen-Rauminformationen/Ihme1500/ihme1500.html
  ihme_shape <- st_read(paste0(path_input_data, "IHME1500_v12/shp/ihme1500__ec4060_v12_poly.shp")) %>%
    st_transform(ihme_shape, crs = crs_proj)
  ihme_shape_bbox <- st_intersection(ihme_shape, bounding_box)
  saveRDS(ihme_shape_bbox, file = paste0(path_proc_data_prod, "ihme_shape_bbox.RDS"))
  
  
  ### soil ###
  
  ### ESDD ###
  # https://esdac.jrc.ec.europa.eu/content/european-soil-database-derived-data
  # all variables of interest except root depth are present for both topsoil (T) and subsoil (S))
  # the topsoil and subsoil have to be weigthed by its depth; topsoil: first 30 cm; subsoil: root depth - 30 cm
  # available soil water content (TAWC) has to be summed up
  
  # read, crop and save root depth raster file
  esdd_depth_roots <- rast(paste0(path_input_data, "STU_EU_Layers/STU_EU_DEPTH_ROOTS.rst"))
  crs(esdd_depth_roots) <- "EPSG:3035"
  esdd_depth_roots_bbox = crop(esdd_depth_roots, vect(bounding_box))  # `vect()` converts sf to terra's vector format
  writeRaster(esdd_depth_roots_bbox, filename = paste0(path_proc_data_prod, "esdd_depth_roots_raster_bbox.tif"))
  
  # perform cell-wise operations to split the raster into raster layers for topsoil and subsoil depths
  # topsoil depth is 30 cm, unless the full depth is shallower
  esdd_depth_T_bbox <- clamp(esdd_depth_roots_bbox, upper = 30, values = TRUE)
  # subsoil depth is full depth - 30 cm, unless the full depth if shallower than 30 cm
  esdd_depth_S_bbox <- ifel(esdd_depth_roots_bbox >= 30, esdd_depth_roots_bbox - 30, 0)
  
  # read, crop and compute depth weighted variables within bbox
  vars_weight <- c("SAND", "SILT", "CLAY", "OC", "BD", "GRAVEL")
  var_names_weight <- c("sand", "silt", "clay", "oc", "bd", "gravel") 
  for(n in 1:length(vars_weight)){
    
    var <- vars_weight[n]
    # topsoil
    esdd_var_T <- rast(paste0(path_input_data, "STU_EU_Layers/STU_EU_T_", var, ".rst"))
    crs(esdd_var_T) <- "EPSG:3035"
    esdd_var_T_bbox = crop(esdd_var_T, vect(bounding_box))
    # subsoil
    esdd_var_S <- rast(paste0(path_input_data, "STU_EU_Layers/STU_EU_S_", var, ".rst"))
    crs(esdd_var_S) <- "EPSG:3035"
    esdd_var_S_bbox = crop(esdd_var_S, vect(bounding_box))
    # weighted mean
    esdd_var_weighted <- (esdd_var_T_bbox * esdd_depth_T_bbox + esdd_var_S_bbox * esdd_depth_S_bbox) / (esdd_depth_T_bbox + esdd_depth_S_bbox)
    writeRaster(esdd_var_weighted, filename = paste0(path_proc_data_prod, "esdd_", var_names_weight[n],"_raster_bbox.tif"))
  }
  
  # read, crop and compute summed water contents within bbox
  # topsoil
  esdd_tawc_T <- rast(paste0(path_input_data, "STU_EU_Layers/STU_EU_T_TAWC", ".rst"))
  crs(esdd_tawc_T) <- "EPSG:3035"
  esdd_tawc_T_bbox = crop(esdd_tawc_T, vect(bounding_box))
  # suboil
  esdd_tawc_S <- rast(paste0(path_input_data, "STU_EU_Layers/STU_EU_S_TAWC", ".rst"))
  crs(esdd_tawc_S) <- "EPSG:3035"
  esdd_tawc_S_bbox = crop(esdd_tawc_S, vect(bounding_box))
  # sum
  esdd_tawc_sum <- esdd_tawc_T_bbox + esdd_tawc_S_bbox
  writeRaster(esdd_tawc_sum, filename = paste0(path_proc_data_prod, "esdd_tawc_raster_bbox.tif"))
  
    
  ### EU-SHD (3D Soil Hydraulic Database of Europe at 250 m resolution) ###
  # https://esdac.jrc.ec.europa.eu/content/3d-soil-hydraulic-database-europe-1-km-and-250-m-resolution
  # load shapefile to identify required tiles for project area
  eu_shd_shape <- st_read(paste0(path_input_data, "EU_SoilHydroGrids_250m/grid_cells_250m_etrs89/grid_cells_250m_etrs89.shp"))
  eu_shd_tiles_required <- st_intersection(eu_shd_shape, bounding_box)
  
  # define required variables and layers
  vars <- c("KS", "THS")
  lays <- c("sl1", "sl2", "sl3", "sl4", "sl5", "sl6", "sl7")
  
  # merge raster tiles for entire study area for each variable and soil layer, separately
  for (v in 1:length(vars)) {
  
    for (l in 1:length(lays)) {
    
    # create an empty vector to collect the raster files
    eu_shd_rasters <- c()
    for (i in 1:nrow(eu_shd_tiles_required)){
    
      eu_shd_raster_i <- rast(paste0(path_input_data, "EU_SoilHydroGrids_250m/EU_SoilHydroGrids_250m/", eu_shd_tiles_required$grid_id[i],
                                     "/", vars[v], "_M_", lays[l], "_", eu_shd_tiles_required$grid_id[i], ".tif"))
      eu_shd_rasters <- c(eu_shd_rasters, eu_shd_raster_i)
    }
    
    eu_shd_rasters_merge <- do.call(merge, eu_shd_rasters)
    eu_shd_rasters_merge_proj <- project(eu_shd_rasters_merge, paste0("EPSG:", crs_proj))
    eu_shd_rasters_merge_proj_bbox = crop(eu_shd_rasters_merge_proj, vect(bounding_box))
    
    writeRaster(eu_shd_rasters_merge_proj_bbox, filename = paste0(path_input_data, "EU_SoilHydroGrids_250m/layers_bbox/eu_shd_",
                                                                  vars[v], "_", lays[l],".tif"), overwrite = TRUE)
    }
  }
  
  # perform cell-wise operations to compute depth weighted averages for each variable, respectively
  for (v in 1:length(vars)){
    
    # read all raster layers for each variable
    file_paths <- list.files(paste0(path_input_data, "EU_SoilHydroGrids_250m/layers_bbox/"),
                             pattern = paste0("^eu_shd_", vars[v],".*\\.tif$"), full.names = TRUE)
    eu_shd_rasters_list <- lapply(file_paths, rast)
    eu_shd_raster_stack <- rast(eu_shd_rasters_list)
  
    # convert units
    # KS = saturated hydraulic conductivity (cm * day-1 * 100) -> convert to cm * h-1
    # THS = saturated volumetric water content -> porosity (cm3 * cm-3 * 100) -> convert to cm3 * cm-3
    if(vars[v] == "KS"){
      eu_shd_raster_stack = eu_shd_raster_stack / 100 / 24
    } else if(vars[v] == "THS") {
      eu_shd_raster_stack = eu_shd_raster_stack / 100
    }
    
    # number of layers in stack
    n_lays <- nlyr(eu_shd_raster_stack)
    # layer depth boundaries
    depth_lays <- c(0, 5, 15, 30, 60, 100, 200)
    # thickness of six layers
    thick_lays <- diff(depth_lays)
    
    # perform trapezoidal integration for layer weighting within modeled depths:
    # invalid values should be ignored (no soil data at this depth as soil shallower)
    # list for numerators: trapezoid contribution for each depth interval
    num_list <- list()
    # list for denominators: thicknesses of valid intervals
    den_list <- list()
    # loop over depth intervals
    for (i in 1:(n_lays - 1)) {
      
      # raster with values at the shallower boundary
      left  <- eu_shd_raster_stack[[i]]
      # raster with values at the deeper boundary
      right <- eu_shd_raster_stack[[i+1]]
      # check if segments are valid (no NA values)
      valid <- !is.na(left) & !is.na(right)
      # compute trapezoid contribution for depth interval, if valid
      num_list[[i]] <- ifel(valid, (left + right) / 2 * thick_lays[i], NA)
      # compute trapezoid contribution for depth interval, if valid
      den_list[[i]] <- ifel(valid, thick_lays[i], NA)
    }
    
    # sum numerators and denominators and calculate weighted mean across all grid cells
    num   <- sum(rast(num_list), na.rm = TRUE)
    den <- sum(rast(den_list), na.rm = TRUE)
    eu_shd_raster <- num / den
  
    writeRaster(eu_shd_raster, filename = paste0(path_proc_data_prod, "eu_shd_", vars[v],"_raster_bbox.tif"))
  }
  

  ### GGT (Global 1-km Gridded Thickness of Soil, Regolith, and Sedimentary Deposit Layers) ###
  # https://www.earthdata.nasa.gov/data/catalog/ornl-cloud-global-soil-regolith-sediment-1304-1
  ggt_raster <- rast(paste0(path_input_data, "Global_Soil_Regolith_Sediment/data/average_soil_and_sedimentary-deposit_thickness.tif"))
  ggt_raster_proj <- project(ggt_raster, paste0("EPSG:", crs_proj))
  ggt_raster_proj_bbox = crop(ggt_raster_proj, vect(bounding_box))
  writeRaster(ggt_raster_proj_bbox, filename = paste0(path_proc_data_prod, "ggt_raster_bbox.tif"))
  

  ### Corine land cover ###
  # https://land.copernicus.eu/en/products/corine-land-cover/clc2018
  st_layers(paste0(path_input_data, "CLC/DATA/U2018_CLC2018_V2020_20u1.gpkg"))
  # load specific geopackage layer (excluding European islands)
  corine_landcover_shape <- st_read(paste0(path_input_data, "CLC/DATA/U2018_CLC2018_V2020_20u1.gpkg"), layer = "U2018_CLC2018_V2020_20u1")
  corine_landcover_shape_bbox = st_intersection(corine_landcover_shape, bounding_box)
  saveRDS(corine_landcover_shape_bbox, file = paste0(path_proc_data_prod, "clc_shape_bbox.RDS"))
  
  
  ### Glacier outlines ###
  ### RGI v7: provides glacier outlines for 2003
  # https://nsidc.org/data/nsidc-0770/versions/7
  rgi_shape <- st_read(paste0(path_input_data,
                              "Randolph_Glacier_Inventory_V7/RGI2000-v7.0-G-11_central_europe/RGI2000-v7.0-G-11_central_europe.shp")) %>%
    st_transform(rgi_shape, crs = crs_proj)
  rgi_shape_bbox <- st_intersection(rgi_shape, bounding_box)
  saveRDS(rgi_shape_bbox, file = paste0(path_proc_data_prod, "rgi_shape_bbox.RDS"))
  
  ### LIA: Reinthaler & Paul (2025) provides reconstructed glacier outlines for 1850
  # https://zenodo.org/records/14336827
  lia_shape <- st_read(paste0(path_input_data, "LIA_glacier_outlines_Alps/LIA_Alps_Reinthaler_Paul_2024.shp")) %>%
    st_transform(lia_shape, crs = crs_proj)
  lia_shape_bbox <- st_intersection(lia_shape, bounding_box)
  saveRDS(lia_shape_bbox, file = paste0(path_proc_data_prod, "lia_shape_bbox.RDS"))
  
  ### Sentinel: Paul et al. (2019) provides glacier outlines for 2015/16
  # https://doi.pangaea.de/10.1594/PANGAEA.909133
  sent_shape <- st_read(paste0(path_input_data, "c3s_gi_rgi11_s2_2015_v2/c3s_gi_rgi11_s2_2015_v2.shp")) %>%
    st_transform(sent_shape, crs = crs_proj)
  sent_shape_bbox <- st_intersection(sent_shape, bounding_box)
  saveRDS(sent_shape_bbox, file = paste0(path_proc_data_prod, "sent_shape_bbox.RDS"))

}