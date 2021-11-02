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
#' @importFrom zonal weighting_grid execute_zonal execute_zonal_cat
#' @importFrom foreach `%dopar%` foreach
#' @importFrom parallel detectCores
#' @importFrom doParallel registerDoParallel
#' @importFrom data.table setnames data.table as.data.table fwrite
#' @importFrom terra rast terrain classify window ext

aggregate_lstm_params = function(gpkg,
                                 catchment_name,
                                 geo_dir,
                                 years = 3,
                                 precision = 9,
                                 out_file = NULL){

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

  `%dopar%` <- foreach::`%dopar%`
  no_cores  <- parallel::detectCores() - 1
  suppressWarnings({
    doParallel::registerDoParallel(cores = no_cores)

  pr_list = foreach::foreach(i = 1:nrow(pr), .combine = cbind) %dopar% {
    zonal::execute_zonal(file = pr$fullname[i], w = gridmet_w)[,-c("ID")]
  }

  tm_list = foreach::foreach(i = 1:nrow(tm), .combine = cbind) %dopar% {
    zonal::execute_zonal(file = tm$fullname[i], w = gridmet_w)[,-c("ID")]
  }

  })

  ppt  = setnames(pr_list, paste0('day', 1:ncol(pr_list)))
  tavg = setnames(tm_list, paste0('day', 1:ncol(tm_list)))

  traits = data.frame(ID = sort(cats$ID))

  ## Mean PPT
  traits$meanPPT = rowMeans(ppt)

  # Snow Fractions
  jennings = filter(files, grepl("jennings", fullname))

  snow_tif = geogrid_warp(jennings$fullname, make_grid(tm$fullname[1]))

  jen = zonal::execute_zonal(file = snow_tif, w = gridmet_w)
  jen$kelvins = jen$V1 + 273.15

  snow_day = tavg[, .SD  < jen$kelvins]
  snow = ppt * snow_day

  traits$snowFrac  = rowSums(snow) / rowSums(ppt)

  ## High PPT Freq
  message("Computing High PPT freq ...")
  ppt5x   = 5 * traits$meanPPT
  highPPT = ppt[, .SD > ppt5x]
  traits$high_ppt_freq = rowSums(highPPT) / years

  message("Computing Low PPT freq ...")
  lowPPT = ppt[, .SD < 1]
  traits$low_ppt_freq = rowSums(lowPPT) / years

  high     <- data.table(t(highPPT))
  mod_cols <- names(high)

  ff = function(x){
    cs = cumsum(x)
    cs = cs - cummax((x == 0) * cs)
    c(ifelse(diff(cs) < 0, cs, NA), cs[length(cs)])
  }

  high[ , (mod_cols) := lapply(.SD, ff), .SDcols = mod_cols]

  low <- data.table(t(lowPPT))
  low[ , (mod_cols) := lapply(.SD, ff), .SDcols = mod_cols]

  traits$high_ppt_dur = colMeans(high, na.rm = TRUE)
  traits$low_ppt_dur  = colMeans(low,  na.rm = TRUE)

  message("Processing Terrain (slope, mean elevation) ...")

  DEM = filter(files, grepl("elevation", fullname))$fullname
  DEM = terra::rast(DEM)
  t   = c(terra::terrain(DEM), DEM)
  terrain = zonal::execute_zonal(t, cats, 'ID') %>%
    setNames(c("ID", "elevation", "slope"))

  terrain$areasqkm = as.numeric(st_area(cats)/1e6)

  traits = left_join(traits, terrain, by = "ID")

  message("Processing soil information (sand, silt, clay, depth, permiability) ...")

  s = filter(files, grepl("sand-1m|silt-1m|clay-1m|rockdepm|GLIM_xx|permeability_permafrost", fullname)) %>%
    filter(grepl(".tif", fullname)) %>%
    dplyr::pull(fullname) %>%
    terra::rast()

  soils_w = zonal::weighting_grid(s$`clay-1m-percent`, cats, "ID")

  soils = zonal::execute_zonal(s, w = soils_w) %>%
    setNames(c("ID", names(s)))

  soils = soils %>%
    mutate(k = -0.60 +  (0.0126*`sand-1m-percent`) - (0.0064*`clay-1m-percent`),
           k = exp(k),
           volumetric_porosity = 50.5 - (0.142*`sand-1m-percent`)  - (0.037*`clay-1m-percent`),
           max_water_content = (volumetric_porosity/100) * rockdepm) %>%
  dplyr::rename(carbonates = GLIM_xx)

  traits = left_join(traits, soils, by = "ID")

  # Energy
  message("Processing Energy Information (PET, AI) ...")

  energy = filter(files, grepl("ai_normal|pet_normal", fullname)) %>%
    dplyr::pull(fullname) %>%
    terra::rast()

  energy = zonal::execute_zonal(energy, w = gridmet_w) %>%
    setnames(c("ID", names(energy)))

  traits = left_join(traits, energy, by = "ID")

  ## Landuse

  message("Processing MODIS Landcover ...")

  modis_lc = filter(files, grepl('MCD12Q1.006/1km/2019-01-01.tif', fullname))

  lc_tiff = terra::rast(modis_lc$fullname) %>%
    terra::crop(cats, snap = "out")

  lu = zonal::execute_zonal_cat(lc_tiff, w = soils_w)

  forest = lu %>%
    filter(value %in% c(1:5)) %>%
    group_by(ID) %>%
    summarise(forest = sum(percentage, na.rm = TRUE)) %>%
    ungroup()

  traits = left_join(traits, forest, by = "ID") %>%
    mutate(forest = ifelse(is.na(forest), 0, forest))

  ## GVF & LAI

  message("Generating and processing LAI and GVF data ...")

  modis_mapping = read.csv('/Users/mjohnson/github/geogrids/inst/modis_lc.csv')

  modis_mapping$ndvi_inf = c(0.81,0.86,0.88,
                             0.90,0.87,0.86,
                             0.86,0.75,0.76,
                             0.74,0.86,0.84,
                             0.86,0.82,NA,
                             0.86,NA)

  ndvi_denom = terra::classify(lc_tiff, select(modis_mapping, is = Class, becomes = ndvi_inf)) - 0.05

  ########
  ndvi = filter(files, grepl('MOD13A3.006/mosaics/landcover', fullname)) %>%
    arrange(desc(fullname)) %>%
    slice(1:(12*years)) %>%
    arrange(fullname)

  # apply scaling factor
  ndvi_rast = terra::rast(ndvi$fullname)
  terra::window(ndvi_rast) <- terra::ext(ndvi_denom)
  x = terra::clamp(ndvi_rast, lower = -2000, upper = 10000, values = FALSE)

  ndvi_rast_processd = (x * 0.0001) - 0.5

  indices<-rep(1:12, times=years)

  gvf = ndvi_rast_processd / ndvi_denom
  gvf = terra::clamp(gvf, lower = 0, upper = 1, values = FALSE)
  gvf.mean <- terra::tapp(gvf, indices, fun = mean)
  gvf.max = max(gvf.mean)
  gvf.min = min(gvf.mean)

  #########
  lai_files = filter(files, grepl('MOD15A2H.006/monthly_means', fullname)) %>%
    arrange(desc(fullname)) %>%
    slice(1:(12*years)) %>%
    arrange(fullname)

  # apply scaling factor

  lai_rast = terra::rast(lai_files$fullname) #* 0.1
  terra::window(lai_rast) <- terra::ext(ndvi_denom)
  lai_rast = lai_rast * 0.1
  lai_rast = terra::clamp(lai_rast, lower = 0, upper = 100, values = FALSE)

  lai.mean <- terra::tapp(lai_rast, indices, fun = mean)
  lai.max = max(lai.mean)
  lai.min = min(lai.mean)

  gvf_lai_rast = rast(c(gvf_max = gvf.max,
                        gvf_diff = gvf.max-gvf.min,
                        lai_max = lai.max,
                        lai_diff = lai.max-lai.min))

  plot(gvf_lai_rast)
  gvf_lai_stat = zonal::execute_zonal(gvf_lai_rast, cats, "ID") %>%
    setnames(c("ID", names(gvf_lai_rast)))

  traits = left_join(traits, gvf_lai_stat, by = "ID")

  message("Writing output CSV to\n ", out_file)

  cols = names(traits)[2:ncol(traits)]
  traits = data.table::as.data.table(traits)
  traits[,(cols) := round(.SD, precision), .SDcols = cols]

  data.table::fwrite(traits, out_file, row.names = FALSE)

  return(out_file)
}
