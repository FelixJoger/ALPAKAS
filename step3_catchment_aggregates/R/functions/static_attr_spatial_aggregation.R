### aggregation of geospatial data products within catchment representations to derive static catchment attributes
# use exact_extract() for aggregation to weight raster cells according to their overlap (coverage fraction) with the polygon

#-------------------------------------------------------------------------------

# function to clip polygons to catchments
clip_to_catchments <- function(data_sf, catch_repr) {
  clipped_list <- lapply(seq_len(nrow(catch_repr)), function(i) {
    st_intersection(data_sf, catch_repr[i, ])
  })
  bind_rows(clipped_list)
}

# function to compute catchment proportions for polygons by aggregation
aggregate_polygon_classes <- function(data_sf, catch_repr, class_col = "attribute_name") {
  # clip data product to catchments
  clipped <- clip_to_catchments(data_sf, catch_repr)
  # dissolve polygons with same category within buffers
  clipped %>%
    group_by(AlpAKaS_ID, .data[[class_col]], area) %>%
    summarise(do_union = TRUE, .groups = "drop") %>%
    # compute area proportions of each class
    mutate(
      area_polygon = as.numeric(st_area(.)),
      proportion_polygon = as.numeric(area_polygon / area) * 100 # in percent
    ) %>%
    st_drop_geometry() %>%
    rename(attribute_name = all_of(class_col)) %>%
    dplyr::select(AlpAKaS_ID, attribute_name, proportion_polygon)
}

# function to bind rows together and add name of corresponding sampling site for every raster cell
bind_extract_list <- function(extract_list, catch_repr, id_col = "AlpAKaS_ID") {
  map2_dfr(extract_list, catch_repr[[id_col]], function(df, id) {
    df[[id_col]] <- id
    df
  })
}

# function to extract raster values
extract_catchment_values <- function(r, catch_repr, id_col = "AlpAKaS_ID") {
  # create list of data frames for each catchment, respectively
  extract_list <- exactextractr::exact_extract(r, catch_repr, include_xy = TRUE)
  # bind rows and add name of corresponding sampling site for every raster cell
  bind_extract_list(extract_list, catch_repr, id_col = id_col)
}

# function to compute weighted mean of variable for each site from raster values
weighted_mean_by_catchment <- function(df, var_name, id_col = "AlpAKaS_ID") {
  df %>%
    group_by(.data[[id_col]]) %>%
    summarise(
      !!var_name := sum(value * coverage_fraction, na.rm = TRUE) /
        sum(coverage_fraction, na.rm = TRUE),
      .groups = "drop"
    )
}

# define function to calculate weighted percentiles
weighted_perc <- function(val, weigth, prob) {
  ord <- order(val)
  val <- val[ord]
  # normalize weights to sum to 1
  weigth <- weigth[ord] / sum(weigth)
  # compute cumulative weight
  cw <- cumsum(weigth)
  val[which(cw >= prob)[1]]
}


### main function to perform spatial aggregation of geospatial products
static_attr_catch_agg <- function(catch_del_name_r, catch_repr) {
  
  # define file path for processed geospatial data products
  path_proc_data_prod <- paste0(path_catch_agg_temp, "proc_data_prod/")
  
  # define file paths for saving static attributes
  path_catch_attr_out <- paste0(path_catch_agg_out, catch_del_name_r, "/static_attributes/")
  
  # compute catchment area in m2
  catch_repr <- catch_repr %>%
    mutate(area = st_area(geometry))
  
  # read metadata for attribute names and lookup table for CAMELS-CH CLC classes
  spatial_metadata <- read.csv(paste0(path_catch_agg_in, "spatial_metadata.csv"))
  # load lookup table for CAMELS-CH CLC classification
  clc_camels_ch_lookup <- read.csv(paste0(path_catch_agg_in, "clc_camels_ch_lookup.csv"))
  # load station meta location to extract site attributes
  station_meta_loc <- read.csv(file = "station_meta_input.csv") %>%
    filter(AlpAKaS_ID %in% catch_repr$AlpAKaS_ID) %>%
    st_as_sf(coords = c("station_easting", "station_northing"), crs = 3035)
  #-----------------------------------------------------------------------------
  ### load reprojected geospatial data products and perform data aggregation as well as variable computations
  
  ### topography and positioning ###
  
  ### COP-DEM ###
  print("COP-DEM")
  
  ### elevation
  cop_dem_elev_raster_bbox <- rast(paste0(path_proc_data_prod, "cop_dem_elev_raster_bbox.tif"))
  cop_dem_elev_raster_df <- extract_catchment_values(cop_dem_elev_raster_bbox, catch_repr)
  ### compute elevation statistics for each site, respectively
  # apply function weighted_perc for weighted percentiles
  cop_dem_elev_catch_attr <- cop_dem_elev_raster_df %>%
    group_by(AlpAKaS_ID) %>%
    summarise(
      elev_mean = sum(value * coverage_fraction, na.rm = TRUE) / sum(coverage_fraction, na.rm = TRUE),
      elev_min = min(value),
      elev_10 = weighted_perc(value, coverage_fraction, 0.1),
      elev_25 = weighted_perc(value, coverage_fraction, 0.25),
      elev_50 = weighted_perc(value, coverage_fraction, 0.5),
      elev_75 = weighted_perc(value, coverage_fraction, 0.75),
      elev_90 = weighted_perc(value, coverage_fraction, 0.9),
      elev_max = max(value),
      .groups = "drop")
  
  ### slope
  cop_dem_slope_raster_bbox <- rast(paste0(path_proc_data_prod, "cop_dem_slope_raster_bbox.tif"))
  cop_dem_slope_raster_df <- extract_catchment_values(cop_dem_slope_raster_bbox, catch_repr)
  ### compute slope statistics for each site, respectively
  cop_dem_slope_catch_attr <- cop_dem_slope_raster_df %>%
    group_by(AlpAKaS_ID) %>%
    summarise(
      slope_mean = sum(value * coverage_fraction, na.rm = TRUE) / sum(coverage_fraction, na.rm = TRUE),
      flat_area_perc = sum(coverage_fraction[value < 3], na.rm = TRUE) / sum(coverage_fraction, na.rm = TRUE) * 100,
      steep_area_perc = sum(coverage_fraction[value > 15], na.rm = TRUE) / sum(coverage_fraction, na.rm = TRUE) * 100,
      .groups = "drop")
  
  ### aspect (northness & eastness)
  ### northness
  cop_dem_northn_raster_bbox <- rast(paste0(path_proc_data_prod, "cop_dem_northn_raster_bbox.tif"))
  cop_dem_northn_raster_df <- extract_catchment_values(cop_dem_northn_raster_bbox, catch_repr)
  cop_dem_northn_agg <- weighted_mean_by_catchment(cop_dem_northn_raster_df, "northn_mean")
  
  ### eastness
  cop_dem_eastn_raster_bbox <- rast(paste0(path_proc_data_prod, "cop_dem_eastn_raster_bbox.tif"))
  cop_dem_eastn_raster_df <- extract_catchment_values(cop_dem_eastn_raster_bbox, catch_repr)
  cop_dem_eastn_agg <- weighted_mean_by_catchment(cop_dem_eastn_raster_df, "eastn_mean")
  
  cop_dem_aspect_catch_attr <- cop_dem_northn_agg %>%
    full_join(cop_dem_eastn_agg, by = "AlpAKaS_ID")
  
  
  ### EU-Hydro ###
  print("EU-Hydro")
  
  eu_hydro_shape_bbox <- readRDS(file = paste0(path_proc_data_prod, "eu_hydro_shape_bbox.RDS")) %>%
    dplyr::select(STRAHLER)
  
  # clip data product to catchments
  eu_hydro_catch_repr_df <- clip_to_catchments(eu_hydro_shape_bbox, catch_repr)
  # dissolve linestrings of all orders within each buffer to sum length of streams within catchments
  eu_hydro_catch_repr_stats <- eu_hydro_catch_repr_df %>%
    group_by(AlpAKaS_ID, area) %>%
    summarise(geometry = st_union(Shape), .groups = "drop")
  # calculate total length within catchment per catchment area
  eu_hydro_catch_attr <- eu_hydro_catch_repr_stats %>%
    mutate(total_length_km = st_length(geometry) / 10^3) %>%
    full_join(catch_repr %>% st_drop_geometry() %>% dplyr::select(AlpAKaS_ID), by = "AlpAKaS_ID") %>%
    mutate(drainage_density = as.numeric(total_length_km / (area / 10^6))) %>% 
    dplyr::select(AlpAKaS_ID, drainage_density) %>%
    mutate(across(drainage_density, ~ replace_na(.x, 0))) %>% # add 0 if not present
    st_drop_geometry()

  
  ### MOHP ###
  print("MOHP")
  
  mohp_prod <- c("divide_stream_distance", "lateral_position", "stream_distance")
  mohp_var <- c("dsd", "lp", "sd")
  
  # data frame to bind aggregated attributes computed for variables individually in loop
  mohp_site_attr <- data.frame(AlpAKaS_ID = station_meta_loc$AlpAKaS_ID)
  for (m in 1:length(mohp_prod)){
    
    mohp_prod_m <- mohp_prod[m]
    mohp_var_m <- mohp_var[m]
    
    for (n in 1:9){
      
      mohp_raster_m_n_bbox <- rast(paste0(path_proc_data_prod, "mohp/mohp_", mohp_var_m, n, "_raster_bbox.tif"))
      # extract raster values for each station position
      mohp_m_n_df_stats <- terra::extract(mohp_raster_m_n_bbox, terra::vect(station_meta_loc)) %>%
        setNames(c("AlpAKaS_ID", paste0(mohp_var_m, n))) %>%
        mutate(AlpAKaS_ID = station_meta_loc$AlpAKaS_ID)
      
      # fill NA values by nearest raster cell with non-NA value, if present
      na_idx <- which(is.na(mohp_m_n_df_stats[[paste0(mohp_var_m, n)]]))
      if (length(na_idx) > 0) {
        
        pts_na <- vect(station_meta_loc[na_idx, ])
        
        # extract all raster cells within 1 km around NA stations
        ext_buf <- terra::extract(
          mohp_raster_m_n_bbox, buffer(pts_na, width = 100),
          cells  = TRUE, xy = TRUE, df = TRUE) %>%
          rename(x_cell = x, y_cell = y)
        names(ext_buf)[2] <- "val"
        
        # extract station coordinates
        pxy <- as.data.frame(crds(pts_na))
        pxy$ID <- seq_len(nrow(pxy))
        
        # for each NA station, only keep non-NA cells, compute distance, and take nearest
        nearest_df <- ext_buf %>%
          filter(!is.na(val)) %>%
          left_join(pxy, by = "ID") %>%
          mutate(dist_m = sqrt((x_cell - x)^2 + (y_cell - y)^2)) %>%
          group_by(ID) %>%
          slice_min(dist_m, with_ties = FALSE) %>%
          ungroup()
        
        # fill in nearest grid cell value
        fill_idx <- na_idx[nearest_df$ID]
        mohp_m_n_df_stats[[paste0(mohp_var_m, n)]][fill_idx] <- nearest_df$val
        
      }
      mohp_site_attr <- left_join(mohp_site_attr, mohp_m_n_df_stats)
    }
  }
  
  
  ### Koeppen-Geiger ###
  print("Koeppen-Geiger")
  
  kg_raster_proj_bbox <- rast(paste0(path_proc_data_prod, "kg_raster_bbox.tif"))
  # create lookup table based on present IDs and join to extracted IDs
  kg_lookup <- data.frame(kg_ID = c(8, 9, 14, 15, 18, 19, 26, 27, 29, 30),
                          climate_zone  = c("Csa", "Csb", "Cfa", "Cfb", "Dsb", "Dsc", "Dfb", "Dfc", "ET", "EF"))
  # extract Koeppen-Geiger zone at stations and join climate zones
  kg_site_attr <- terra::extract(kg_raster_proj_bbox, terra::vect(station_meta_loc)) %>%
    setNames(c("AlpAKaS_ID", "kg_ID")) %>%
    mutate(AlpAKaS_ID = station_meta_loc$AlpAKaS_ID,
           kg_ID = as.numeric(as.character(kg_ID))) %>%
    left_join(kg_lookup, by = "kg_ID") %>%
    dplyr::select(-kg_ID)
  
  
  ### land cover ###
  
  ### Corine land cover ###
  print("CLC")
  
  clc_shape_bbox <- readRDS(file = paste0(path_proc_data_prod, "clc_shape_bbox.RDS")) %>% 
    mutate(Code_18 = as.numeric(Code_18)) %>%
    # join lookup table to group CLC classes according to CAMELS-CH
    left_join(clc_camels_ch_lookup, by = c("Code_18" = "clc_code"))
  
  clc_catch_repr_df_stats <- aggregate_polygon_classes(clc_shape_bbox, catch_repr)
  expected_attr <- c("crop_perc", "grass_perc", "shrub_perc", "dwood_perc", "mix_wood_perc", "ewood_perc",
                     "wetland_perc", "inwater_perc", "ice_perc", "loose_rock_perc", "rock_perc", "urban_perc")
  # write proportions into data frame of wide format including all classes and 0 if not present and reorder
  clc_catch_attr <- clc_catch_repr_df_stats %>%
    complete(AlpAKaS_ID, attribute_name = expected_attr, fill = list(proportion_polygon = 0)) %>%
    pivot_wider(id_cols = AlpAKaS_ID, names_from = attribute_name, values_from = proportion_polygon, values_fill = 0) %>%
    dplyr::select(AlpAKaS_ID, all_of(expected_attr))
  
  
  ### glacier ###
  print("glacier")
  
  ### LIA
  lia_shape_bbox <- readRDS(file = paste0(path_proc_data_prod, "lia_shape_bbox.RDS"))
  lia_catch_repr_df <- clip_to_catchments(lia_shape_bbox, catch_repr) %>%
    dplyr::select(AlpAKaS_ID, area) %>%
    mutate(extent_year = "glacier_1850_perc")
  
  ### RGI
  rgi_shape_bbox <- readRDS(file = paste0(path_proc_data_prod, "rgi_shape_bbox.RDS"))
  rgi_catch_repr_df <- clip_to_catchments(rgi_shape_bbox, catch_repr) %>%
    dplyr::select(AlpAKaS_ID, area) %>%
    mutate(extent_year = "glacier_2003_perc")
  
  ### Sentinel
  sent_shape_bbox <- readRDS(file = paste0(path_proc_data_prod, "sent_shape_bbox.RDS"))
  sent_catch_repr_df <- clip_to_catchments(sent_shape_bbox, catch_repr) %>%
    dplyr::select(AlpAKaS_ID, area) %>%
    mutate(extent_year = "glacier_2015_perc")
  
  glacier_catch_repr_df <- bind_rows(lia_catch_repr_df, rgi_catch_repr_df, sent_catch_repr_df) %>%
    mutate(extent_year = factor(extent_year, levels = c("glacier_1850_perc", "glacier_2003_perc", "glacier_2015_perc")))
  glacier_catch_repr_df_stats <- aggregate_polygon_classes(glacier_catch_repr_df, catch_repr, class_col = "extent_year") %>%
    rename(glacier_prop = proportion_polygon)
  glacier_catch_attr <- glacier_catch_repr_df_stats %>%
    pivot_wider(id_cols = AlpAKaS_ID, names_from = attribute_name, values_from = c(glacier_prop), values_fill = 0)
  # add information of sites without glacier proportion
  glacier_catch_attr <- tibble(AlpAKaS_ID = catch_repr$AlpAKaS_ID) %>%
    left_join(glacier_catch_attr, by = "AlpAKaS_ID") %>%
    mutate(across(where(is.numeric), ~ replace_na(., 0)))

  
  ### soil ###
  
  ### ESDD ###
  print("ESDD")
  
  spatial_metadata_esdd <- spatial_metadata %>%
    filter(product_name == "esdd")
  
  # data frame to bind aggregated attributes computed for variables individually in loop
  esdd_catch_attr <- catch_repr[, c("AlpAKaS_ID")] %>%
    st_drop_geometry()
  for(n in 1:nrow(spatial_metadata_esdd)){
    
    var_name <- spatial_metadata_esdd$attribute_name[n]
    var_ident <- spatial_metadata_esdd$attr_ident[n]
    
    esdd_var <- rast(paste0(path_proc_data_prod, "esdd_", var_ident, "_raster_bbox.tif"))
    esdd_var_df <- extract_catchment_values(esdd_var, catch_repr)
    esdd_var_df_mean <- weighted_mean_by_catchment(esdd_var_df, var_name)
    
    esdd_catch_attr <- esdd_catch_attr %>%
      left_join(esdd_var_df_mean, by = "AlpAKaS_ID")
  }
  
  
  ### EU-SHD ###
  print("EU-SHD")
  
  spatial_metadata_eu_shd <- spatial_metadata %>%
    filter(product_name == "eu_shd")
  
  # data frame to bind aggregated attributes computed for variables individually in loop
  eu_shd_catch_attr <- catch_repr[, c("AlpAKaS_ID")] %>%
    st_drop_geometry()
  for(n in 1:nrow(spatial_metadata_eu_shd)){
    
    var_name <- spatial_metadata_eu_shd$attribute_name[n]
    var_ident <- spatial_metadata_eu_shd$attr_ident[n]
    
    eu_shd_var <- rast(paste0(path_proc_data_prod, "eu_shd_", var_ident, "_raster_bbox.tif"))
    eu_shd_var_df <- extract_catchment_values(eu_shd_var, catch_repr)
    eu_shd_var_df_mean <- weighted_mean_by_catchment(eu_shd_var_df, var_name)
    
    eu_shd_catch_attr <- eu_shd_catch_attr %>%
      left_join(eu_shd_var_df_mean, by = "AlpAKaS_ID")
  }
  
  
  ### GGT ###
  print("GGT")
  
  spatial_metadata_ggt <- spatial_metadata %>%
    filter(product_name == "ggt")
  
  var_name <- spatial_metadata_ggt$attribute_name[1]
  var_ident <- spatial_metadata_ggt$attr_ident[1]
  
  ggt_raster_bbox <- rast(paste0(path_proc_data_prod, "ggt_raster_bbox.tif"))
  ggt_raster_df <- extract_catchment_values(ggt_raster_bbox, catch_repr)
  ggt_catch_attr <- weighted_mean_by_catchment(ggt_raster_df, var_name)
  
  
  ### hydrogeology ###
  
  ### MEDKAM ###
  print("MEDKAM")
  medkam_shape_bbox <- readRDS(file = paste0(path_proc_data_prod, "medkam_shape_bbox.RDS"))
  
  # join metadata
  spatial_metadata_medkam <- spatial_metadata %>%
    filter(product_name == "medkam")
  medkam_shape_bbox <- medkam_shape_bbox %>%
    left_join(spatial_metadata_medkam, by = c("MEDKAM_Cla" = "attr_ident"))
  medkam_catch_repr_df_stats <- aggregate_polygon_classes(medkam_shape_bbox, catch_repr)
  # write proportions into data frame of wide format including all classes and 0 if not present
  medkam_catch_attr <- medkam_catch_repr_df_stats %>%
    pivot_wider(id_cols = AlpAKaS_ID, names_from = attribute_name, values_from = proportion_polygon, values_fill = 0)
 
  
  ### IHME ###
  print("IHME")
  ihme_shape_bbox <- readRDS(file = paste0(path_proc_data_prod, "ihme_shape_bbox.RDS"))
  
  spatial_metadata_ihme <- spatial_metadata %>%
    filter(product_name == "ihme")
  
  ### lithology-level 4
  ihme_shape_bbox_litho4 <- ihme_shape_bbox %>%
    left_join(spatial_metadata_ihme, by = c("LEVEL4" = "attr_ident"))
  ihme_litho4_catch_repr_df_stats <- aggregate_polygon_classes(ihme_shape_bbox_litho4, catch_repr) %>%
    filter(!is.na(attribute_name))
  ihme_litho4_catch_repr_df_stats_wide <- ihme_litho4_catch_repr_df_stats %>%
    pivot_wider(id_cols = AlpAKaS_ID, names_from = attribute_name, values_from = c(proportion_polygon), values_fill = 0)
  # add categories not present in the data to data frame and reorder
  ihme_level4_cat <- c("AlpAKaS_ID", "calc_perc", "silic_perc", "magm_perc", "metam_perc", "calc_crse_sed_perc",
                       "calc_fine_sed_perc", "silic_crse_sed_perc", "silic_fine_sed_perc", "crse_sed_perc", "fine_sed_perc")
  missing_cols <- setdiff(ihme_level4_cat, names(ihme_litho4_catch_repr_df_stats_wide))
  ihme_litho4_catch_repr_df_stats_wide[missing_cols] <- 0
  ihme_litho4_catch_repr_df_stats_wide <- ihme_litho4_catch_repr_df_stats_wide %>%
    dplyr::select(all_of(ihme_level4_cat))
  
  ### aquifer type
  ihme_shape_bbox_aquif <- ihme_shape_bbox %>%
    mutate(AQUIF_CODE = as.character(AQUIF_CODE)) %>%
    left_join(spatial_metadata_ihme, by = c("AQUIF_CODE" = "attr_ident"))
  ihme_aquif_catch_repr_df_stats <- aggregate_polygon_classes(ihme_shape_bbox_aquif, catch_repr) %>%
    filter(!is.na(attribute_name))
  ihme_aquif_catch_repr_df_stats_wide <- ihme_aquif_catch_repr_df_stats %>%
    pivot_wider(id_cols = AlpAKaS_ID, names_from = attribute_name, values_from = c(proportion_polygon), values_fill = 0) %>%
    dplyr::select(AlpAKaS_ID, porous_high_prod_perc, porous_low_mod_prod_perc, fissured_high_prod_perc, fissured_low_mod_prod_perc,
           local_aquif_perc, non_aquif_perc)
  
  ### lithology-level 1: subset (karstified)
  ihme_shape_bbox_litho1 <- ihme_shape_bbox %>%
    left_join(spatial_metadata_ihme, by = c("LEVEL1" = "attr_ident")) %>%
    filter(!is.na(attribute_name)) %>%
    mutate(across(where(is.numeric), ~ replace_na(., 0)))
  ihme_litho1_catch_repr_df_stats <- aggregate_polygon_classes(ihme_shape_bbox_litho1, catch_repr) %>%
    filter(!is.na(attribute_name))
  ihme_litho1_catch_repr_df_stats_wide <- ihme_litho1_catch_repr_df_stats %>%
    pivot_wider(id_cols = AlpAKaS_ID, names_from = attribute_name, values_from = c(proportion_polygon), values_fill = 0) %>%
    dplyr::select(AlpAKaS_ID, lime_karst_perc, dol_lime_karst_perc, chalk_lime_karst_perc)
  
  ### join and save derived attributes
  ihme_catch_attr <- ihme_litho4_catch_repr_df_stats_wide %>%
    left_join(ihme_litho1_catch_repr_df_stats_wide, by = "AlpAKaS_ID") %>%
    left_join(ihme_aquif_catch_repr_df_stats_wide, by = "AlpAKaS_ID") %>%
    mutate(across(where(is.numeric), ~ replace_na(., 0)))
  
  
  #-----------------------------------------------------------------------------
  ### join, round, and save attributes by categories
  
  ### topography and positioning ###
  # only use order 1 and 2 of MOHP
  mohp_site_attr_sub <- mohp_site_attr %>%
    dplyr::select(AlpAKaS_ID, matches("1|2"))
  topography_attr_list <- list(cop_dem_elev_catch_attr, cop_dem_slope_catch_attr, cop_dem_aspect_catch_attr,
                               eu_hydro_catch_attr, mohp_site_attr_sub, kg_site_attr)
  topography_attr <- reduce(topography_attr_list, left_join, by = "AlpAKaS_ID") %>%
    mutate(across(where(is.numeric), ~ round(.x, 3))) %>%
    rename(ALPAKAS_ID = AlpAKaS_ID)
  write.csv(topography_attr, paste0(path_catch_attr_out, "AlpAKaS_topography_attributes.csv"), row.names = FALSE)
  
  ### land cover ###
  land_cover_attr <- clc_catch_attr %>%
    mutate(across(where(is.numeric), ~ round(.x, 3))) %>%
    rename(ALPAKAS_ID = AlpAKaS_ID)
  write.csv(land_cover_attr, paste0(path_catch_attr_out, "AlpAKaS_land_cover_attributes.csv"), row.names = FALSE)
  
  ### glacier ###
  glacier_attr <- glacier_catch_attr %>%
    mutate(across(where(is.numeric), ~ round(.x, 3))) %>%
    rename(ALPAKAS_ID = AlpAKaS_ID)
  write.csv(glacier_attr, paste0(path_catch_attr_out, "AlpAKaS_glacier_attributes.csv"), row.names = FALSE)
  
  ### soil ###
  soil_attr_list <- list(esdd_catch_attr, eu_shd_catch_attr, ggt_catch_attr)
  soil_attr <- reduce(soil_attr_list, left_join, by = "AlpAKaS_ID") %>%
    mutate(across(where(is.numeric), ~ round(.x, 3))) %>%
    rename(ALPAKAS_ID = AlpAKaS_ID)
  write.csv(soil_attr, paste0(path_catch_attr_out, "AlpAKaS_soil_attributes.csv"), row.names = FALSE)
  
  ### hydrogeology ###
  hydrogeology_attr_list <- list(medkam_catch_attr, ihme_catch_attr)
  hydrogeology_attr <- reduce(hydrogeology_attr_list, left_join, by = "AlpAKaS_ID") %>%
    mutate(across(where(is.numeric), ~ round(.x, 3))) %>%
    rename(ALPAKAS_ID = AlpAKaS_ID)
  write.csv(hydrogeology_attr, paste0(path_catch_attr_out, "AlpAKaS_hydrogeology_attributes.csv"), row.names = FALSE)
  

}