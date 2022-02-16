#' Aggregated CFE variables frome NWM
#' @param gpkg a geopackage with aggregation units
#' @param catchment_name the layer name of the aggregation units
#' @param flowline_name the layer name of the flowpath units
#' @param single_layer should only the top layer of a multilayer parameter be processed?
#' @param nwm_dir the NWM directory
#' @param precision the precision of the computations
#' @param out_file the location to write output
#' @return NULL
#' @export
#' @importFrom sf read_sf st_drop_geometry
#' @importFrom terra rast ext crs
#' @importFrom zonal weighting_grid
#' @importFrom dplyr left_join rename mutate select summarise across everything bind_cols right_join group_by inner_join
#' @importFrom data.table fwrite setnames
#' @importFrom tidyr unnest_longer
#' @importFrom RNetCDF open.nc var.get.nc
#' @importFrom stats weighted.mean complete.cases

# aggregate_nwm_params_reduce(
#   gpkg           =  '/Volumes/Transcend/ngen/CONUS-hydrofabric/ngen-release/01a/2021-11-01/spatial/hydrofabric.gpkg',
#   catchment_name = "catchments",
#   flowline_name  = "flowpaths",
#   single_layer   = TRUE,
#   nwm_dir        =  '/Volumes/Transcend/nwmCONUS-v216',
#   precision      =  9,
#   out_dir        =  "/Volumes/Transcend/ngen/CONUS-hydrofabric/ngen-release/01a//2021-11-01/parameters"
# )

aggregate_cfe = function(gpkg,
                         nwm_dir,
                         catchment_name,
                         flowline_name,
                         single_layer = FALSE,
                         precision = 9,
                         out_file = NULL) {
  .SD <- . <- .data <-  NULL

  template = list(
    ext = ext(
      -2303999.62876143,
      2304000.37123857,
      -1920000.70008381,
      1919999.29991619
    ),
    crs = 'PROJCS["Sphere_Lambert_Conformal_Conic",GEOGCS["GCS_Sphere",DATUM["D_Sphere",SPHEROID["Sphere",6370000.0,0.0]],PRIMEM["Greenwich",0.0],UNIT["Degree",0.0174532925199433]],PROJECTION["Lambert_Conformal_Conic"],PARAMETER["false_easting",0.0],PARAMETER["false_northing",0.0],PARAMETER["central_meridian",-97.0],PARAMETER["standard_parallel_1",30.0],PARAMETER["standard_parallel_2",60.0],PARAMETER["latitude_of_origin",40.000008],UNIT["Meter",1.0]];-35691800 -29075200 126180232.640845;-100000 10000;-100000 10000;0.001;0.001;0.001;IsHighPrecision'
  )

  cats = sf::read_sf(gpkg, catchment_name)

  soils = rast(file.path(nwm_dir, 'soilproperties_CONUS_FullRouting.nc'))
  terra::ext(soils) = template$ext
  terra::crs(soils) = template$crs

  nwm_w_1000m = zonal::weighting_grid(soils, geom = cats,  ID = "ID")

  ### soil_properties
  soil_mode_var = "bexp"
  soil_gm_var   = c("dksat", "psisat")
  soil_mean_var = c("slope", "smcmax", "smcwlt", "refkdt")

  x1 = map_vars(
    r =  floor(soils[[grepl(soil_mode_var, names(soils))]]),
    w = nwm_w_1000m,
    file = NULL,
    precision,
    single_layer = single_layer,
    FUN = "mode",
    prefix = "sp_"
  )

  x2 = map_vars(
    r =  soils[[grepl(paste0(soil_gm_var, collapse = "|"), names(soils))]],
    w = nwm_w_1000m,
    file = NULL,
    precision,
    single_layer = single_layer,
    FUN = "gm_mean",
    prefix = "sp_"
  )

  x3 = map_vars(
    r =  soils[[grepl(paste0(soil_mean_var, collapse = "|"), names(soils))]],
    w = nwm_w_1000m,
    file = NULL,
    precision,
    single_layer = single_layer,
    FUN = "mean",
    prefix = "sp_"
  )


  traits = left_join(x1, left_join(x2, x3, by = 'ID'), by = 'ID')

  ### WRF
  wrf = rast(file.path(nwm_dir, 'wrfinput_CONUS.nc'),
             subds = c("IVGTYP", "ISLTYP"))
  terra::ext(wrf) = template$ext
  terra::crs(wrf) = template$crs

  x4 = map_vars(
    r =  wrf,
    w = nwm_w_1000m,
    file = NULL,
    precision,
    single_layer = single_layer,
    FUN = "mode",
    prefix = "wf_"
  )

  traits = left_join(traits, x4, by = 'ID')

  ### FULL DOM

  # fulldom = rast(file.path(nwm_dir, 'Fulldom_CONUS_FullRouting.nc'), subds = c('LKSATFAC'))
  # cols = paste0("fd_", names(fulldom))
  # fulldom = execute_zonal(fulldom, geom = cats, "ID", FUN  = "mean")
  # names(fulldom) = c("ID", cols)
  # fulldom[,(cols) := round(.SD, precision), .SDcols=cols]
  # traits = left_join(traits, fulldom, by = 'ID')

  ####

  if (!is.null(flowline_name)) {
    crosswalk <- st_drop_geometry(read_sf(gpkg, flowline_name))

    crosswalk = select(crosswalk, .data$ID, .data$member_COMID) %>%
      mutate(comid = strsplit(.data$member_COMID, ",")) %>%
      tidyr::unnest_longer(col = c("comid")) %>%
      mutate(ID = .data$ID, comid = as.integer(.data$comid)) %>%
      filter(!duplicated(.))

    gwnc = RNetCDF::open.nc(file.path(nwm_dir, 'GWBUCKPARM_CONUS_FullRouting.nc'))

    vars      = c("Area_sqkm", "ComID", "Coeff",  "Zmax")
    vars_mode = c("Area_sqkm", "ComID", "Expon")

    gwparams_means = suppressMessages({
      lapply(vars, function(x)
        RNetCDF::var.get.nc(gwnc, x)) %>%
        bind_cols() %>%
        setnames(vars) %>%
        rename(comid = .data$ComID) %>%
        mutate(comid = as.integer(.data$comid)) %>%
        inner_join(select(crosswalk, .data$ID, .data$comid), by = 'comid') %>%
        filter(complete.cases(.)) %>%
        filter(!duplicated(.)) %>%
        group_by(.data$ID) %>%
        summarise(across(everything(), ~ round(
          weighted.mean(.x, w = .data$Area_sqkm, na.rm = TRUE),
          precision
        ))) %>%
        select(-.data$comid, -.data$Area_sqkm)
    })


    gwparams_mode = suppressMessages({
      lapply(vars_mode, function(x)
        RNetCDF::var.get.nc(gwnc, x)) %>%
        bind_cols() %>%
        setnames(vars_mode) %>%
        rename(comid = .data$ComID) %>%
        inner_join(select(crosswalk, .data$ID, .data$comid), by = 'comid') %>%
        filter(complete.cases(.)) %>%
        filter(!duplicated(.)) %>%
        group_by(.data$ID) %>%
        summarise(Expon = zonal::getmode(floor(.data$Expon)))
    })

    traits = left_join(gwparams_means, gwparams_mode, by = 'ID') %>%
      mutate(ID = gsub("wb-", "cat-", .data$ID)) %>%
      setnames(c('ID', paste0('gw_', names(.)[-1]))) %>%
      right_join(traits, by = 'ID')
  }

  data.table::fwrite(traits, out_file, row.names = FALSE)

}

aggregate_noahowp = function(gpkg,
                             nwm_dir,
                             catchment_name,
                             single_layer = FALSE,
                             precision = 9,
                             out_file = NULL) {
  .SD <- NULL

  template = list(
    ext = ext(
      -2303999.62876143,
      2304000.37123857,
      -1920000.70008381,
      1919999.29991619
    ),
    crs = 'PROJCS["Sphere_Lambert_Conformal_Conic",GEOGCS["GCS_Sphere",DATUM["D_Sphere",SPHEROID["Sphere",6370000.0,0.0]],PRIMEM["Greenwich",0.0],UNIT["Degree",0.0174532925199433]],PROJECTION["Lambert_Conformal_Conic"],PARAMETER["false_easting",0.0],PARAMETER["false_northing",0.0],PARAMETER["central_meridian",-97.0],PARAMETER["standard_parallel_1",30.0],PARAMETER["standard_parallel_2",60.0],PARAMETER["latitude_of_origin",40.000008],UNIT["Meter",1.0]];-35691800 -29075200 126180232.640845;-100000 10000;-100000 10000;0.001;0.001;0.001;IsHighPrecision'
  )

  cats = sf::read_sf(gpkg, catchment_name)

  soils = suppressWarnings(rast(file.path(
    nwm_dir, 'soilproperties_CONUS_FullRouting.nc'
  )))
  terra::ext(soils) = template$ext
  terra::crs(soils) = template$crs

  nwm_w_1000m = zonal::weighting_grid(soils, geom = cats,  ID = "ID")

  ### soil_properties
  soil_mean_var = c('cwpvt', 'vcmx25', 'mp', 'mfsno')
  x1 = map_vars(
    r =  soils[[grepl(paste(soil_mean_var, collapse = "|"), names(soils))]],
    w = nwm_w_1000m,
    file = NULL,
    precision,
    single_layer = single_layer,
    FUN = "mean",
    prefix = "sp_"
  )

  data.table::fwrite(x1, out_file, row.names = FALSE)

}
