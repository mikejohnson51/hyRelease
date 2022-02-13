
#' Aggregated LSTM variables
#' @param gpkg a geopackage with aggregation units
#' @param catchment_name the layer name of the aggregation units
#' @param geo_dir the geogroids geo directory
#' @param years number of years of PPT and Temperture to process
#' @param precision the precision of the computations
#' @param out_file the CSV file path to write the results to
#' @return the out_file path
#' @export
#' @importFrom geogrids geo_cache_list geogrid_warp make_grid
#' @importFrom sf read_sf st_area
#' @importFrom dplyr  `%>%` filter arrange desc slice arrange pull left_join rename mutate group_by summarise ungroup
#' @importFrom zonal weighting_grid execute_zonal
#' @importFrom foreach `%dopar%` foreach
#' @importFrom doMC registerDoMC
#' @importFrom data.table setnames data.table as.data.table fwrite `:=`
#' @importFrom terra rast terrain classify window ext

aggregate_basin_attributes = function(gpkg,
                                 catchment_name,
                                 geo_dir,
                                 years = 3,
                                 precision = 9,
                                 out_file = NULL){

  .SD <- NULL

  files = geo_cache_list(geo_dir)

  if(is.null(out_file)){ out_file = file.path(dirname(gpkg),  "lstm_attributes.csv") }

  message("Loading aggregation units ...")
  cats = read_sf(gpkg, catchment_name)

  # Load PPT, TAVG, and Snow Thresholds
  message("Loading PPT, AVG, and SNOW Thresholds ...")
  message("Using ", years, " years of Gridmet Data")

  pr = filter(files, grepl('pr', fullname)) %>%
    filter(grepl(".nc", fullname)) %>%
    arrange(desc(fullname)) %>%
    slice(1:years) %>%
    arrange(fullname)

  tm = filter(files, grepl('tavg', fullname))  %>%
    filter(grepl(".nc", fullname)) %>%
    arrange(desc(fullname)) %>%
    slice(1:years) %>%
    arrange(fullname)

  gridmet_w = zonal::weighting_grid(pr$fullname[1], cats, "ID")

  suppressWarnings({

    doMC::registerDoMC()

    pr_list = foreach::foreach(i = 1:nrow(pr)) %dopar% {
      zonal::execute_zonal(file = pr$fullname[i], w = gridmet_w)
    }

    tm_list = foreach::foreach(i = 1:nrow(tm)) %dopar% {
      zonal::execute_zonal(file = tm$fullname[i], w = gridmet_w)
    }
  })

  ppt  = Reduce(function(dtf1, dtf2) merge(dtf1, dtf2, by = "ID", all.x = TRUE), pr_list)
  tavg = Reduce(function(dtf1, dtf2) merge(dtf1, dtf2, by = "ID", all.x = TRUE), tm_list)

  ppt  = setnames(ppt, c("ID", paste0('day', 1:(ncol(ppt)-1))))
  tavg = setnames(tavg, c("ID", paste0('day', 1:(ncol(tavg)-1))))

  stopifnot(identical(ppt$ID, tavg$ID))

  gridmet_out_id = ppt$ID
  ppt  = ppt[,-"ID"]
  tavg = tavg[,-"ID"]

  # Snow Fractions
  jennings = filter(files, grepl("jennings", fullname))
  snow_tif = geogrid_warp(file = jennings$fullname, grid = make_grid(tm$fullname[1]), r = "bilinear")

  jen = zonal::execute_zonal(file = snow_tif + 273.15, w = gridmet_w, join = FALSE)
  snow_day = tavg[, .SD  < jen$V1]
  snow = ppt * snow_day

  ff = function(x){
    cs = cumsum(x)
    cs = cs - cummax((x == 0) * cs)
    c(ifelse(diff(cs) < 0, cs, NA), cs[length(cs)])
  }

  highPPT = ppt[, .SD > (5*rowMeans(ppt))]
  lowPPT = ppt[, .SD < 1]

  high     <- data.table(t(highPPT))
  low      <- data.table(t(lowPPT))
  mod_cols <- names(high)

  high[ , (mod_cols) := lapply(.SD, ff), .SDcols = mod_cols]
  low[ , (mod_cols) := lapply(.SD, ff), .SDcols = mod_cols]

  traits = data.frame(ID = gridmet_out_id,
                       meanPPT = rowMeans(ppt),
                       snowFrac = rowSums(snow) / rowSums(ppt),
                       high_ppt_freq = rowSums(highPPT) / years,
                       low_ppt_freq =  rowSums(lowPPT) / years,
                       high_ppt_dur = colMeans(high, na.rm = TRUE),
                       low_ppt_dur  = colMeans(low,  na.rm = TRUE))


  message("Processing Terrain (slope, mean elevation) ...")

  DEM = terra::rast(filter(files, grepl("elevation", fullname))$fullname)

  t = c(DEM, 1000*tan(terrain(DEM, v = "slope", unit = "radians")))

  terrain = zonal::execute_zonal(t, cats, 'ID', join = FALSE) %>%
    setnames(c("ID", "elevation", "slope"))

  terrain = left_join(terrain, select(st_drop_geometry(cats), ID, areasqkm = area_sqkm))

  traits = left_join(traits, terrain, by = "ID")

  message("Processing soil information (sand, silt, clay, depth, permiability) ...")

  s = filter(files, grepl("sand-1m|silt-1m|clay-1m|rockdepm|GLIM_xx|permeability_permafrost", fullname)) %>%
    filter(grepl(".tif", fullname)) %>%
    dplyr::pull(fullname) %>%
    terra::rast()

  soils_w = zonal::weighting_grid(s, cats, "ID")

  soils = zonal::execute_zonal(s, w = soils_w, join = FALSE) %>%
    setnames(c("ID", names(s))) %>%
    mutate(`clay-1m-percent` = .data$`clay-1m-percent` * 100,
           `sand-1m-percent` = .data$`sand-1m-percent` * 100,
           `silt-1m-percent` = .data$`silt-1m-percent` * 100,
            rockdepm = rockdepm / 100,
           carbonate_rocks_frac = .data$GLIM_xx/ 100,
           GLIM_xx = NULL,
           soil_conductivity = -0.60 +  (0.0126*.data$`sand-1m-percent`) - (0.0064*.data$`clay-1m-percent`),
           soil_conductivity = exp(.data$soil_conductivity),
           geol_porostiy = 50.5 - (0.142*.data$`sand-1m-percent`)  - (0.037*.data$`clay-1m-percent`),
           max_water_content = (.data$geol_porostiy/100) * .data$rockdepm)


  soils = filter(files, grepl('average_soil_and_sedimentary-deposit_thickness.tif', fullname))$fullname %>%
    zonal::execute_zonal(cats, "ID", join = FALSE) %>%
    setnames(c("ID", "soil_depth_pelletier")) %>%
    left_join(soils, by = "ID")

  traits = left_join(traits, soils, by = "ID")

  # Energy
  message("Processing Energy Information (PET, AI) ...")

  energy = filter(files, grepl("ai_normal|pet_normal", .data$fullname)) %>%
    dplyr::pull(.data$fullname) %>%
    terra::rast()

  energy = zonal::execute_zonal(energy, w = gridmet_w, join = FALSE) %>%
    setnames(c("ID", "aridity", "pet_mean")) %>%
    mutate(pet_mean = .data$pet_mean  *.1)

  traits = left_join(traits, energy, by = "ID")

  ## Landuse

  message("Processing MODIS Landcover ...")

  modis_lc = filter(files, grepl('MCD12Q1.006/mosaics_cog/2019-01-01.tif$', .data$fullname))

  lc_tiff = terra::rast(modis_lc$fullname)

  mod_500 = zonal::weighting_grid(lc_tiff, cats, "ID")

  lu = zonal::execute_zonal(lc_tiff, w = mod_500, FUN = "freq", join = FALSE)

  lu_dom = zonal::execute_zonal(lc_tiff, w = mod_500, FUN = "mode", join = FALSE) %>%
    setnames(c("ID", "dom_lc"))

  forest = lu %>%
    filter(value %in% c(1:5)) %>%
    group_by(ID) %>%
    summarise(forest = sum(percentage, na.rm = TRUE)) %>%
    ungroup()

  traits = left_join(left_join(traits, forest, by = "ID"), lu_dom, by = "ID") %>%
    mutate(forest = ifelse(is.na(forest), 0, forest))

  ## GVF & LAI

  message("Generating and processing LAI and GVF data ...")

  lai = filter(files, grepl("LAI/summary/cogs/", fullname)) %>%
    filter(grepl("tif$", fullname)) %>%
    filter(grepl("diff|max", fullname))

  lai = zonal::execute_zonal(rast(lai$fullname), w = mod_500, FUN = "mean", join = FALSE) %>%
    setnames(c("ID", "lai_diff", "lai_max"))

  gvf = filter(files, grepl("GVF/summary/cogs/", fullname)) %>%
    filter(grepl("tif$", fullname)) %>%
    filter(grepl("diff|max", fullname))

  gvf = zonal::execute_zonal(rast(gvf$fullname), cats, "ID", FUN = "mean", join = FALSE) %>%
    setnames(c("ID", "gvf_diff", "gvf_max"))

  traits = left_join(left_join(traits, gvf, by = "ID"), lai, by = "ID")

  message("Writing output CSV to\n ", out_file)

  cols = names(traits)[2:ncol(traits)]
  traits = data.table::as.data.table(traits)
  traits[,(cols) := round(.SD, precision), .SDcols = cols]

  data.table::fwrite(traits, out_file, row.names = FALSE)

  return(out_file)
}
