#' Map Variables
#' @description  Wrapper to extract, rename and return or write a csv of zonal variables
#' @param r a file path, terra SpatRaster or raster raster object
#' @param w a weighting file computed with zonal::weighting_grid
#' @param file a CSV path to write data to, if NULL data.table object is returned
#' @param precision computational precision of the zonal calculations
#' @param FUN a summary function (see zonal::execute_zonal)
#' @param verbose should messages be emitted?
#' @return file path or data.table
#' @export
#' @importFrom zonal execute_zonal
#' @importFrom terra window
#' @importFrom data.table fwrite

map_vars = function(r, w, file = NULL, precision = 9,
                    FUN = "mean", single_layer = FALSE,
                    verbose = TRUE, prefix = NULL){


  cols = names(r)

  if(single_layer){
    cols1 = grep('stag=1', cols, value = TRUE)
    cols2 = cols[!grepl('stag=', cols)]
    cols = c(cols1, cols2)
    r = r[[cols]]
  }

  cols = gsub("_Time=1", "", names(r))

  out = execute_zonal(r, w = w, FUN = FUN)

  if(is.null(prefix)){
    names(out) = c("ID",  cols)
  } else {
    cols = paste0(prefix, cols)
    names(out) = c("ID",  cols)
  }

  out[,(cols) := round(.SD, precision), .SDcols = cols]
  terra::window(r) = NULL

  if(is.null(file)){
    return(out)
  } else {
    data.table::fwrite(out, file, row.names = FALSE)
    message(file)
    return(file)
  }
}

#' Aggregated NWM variables
#' @param gpkg a geopackage with aggregation units
#' @param catchment_name the layer name of the aggregation units
#' @param flowline_name the layer name of the flowpath units
#' @param nwm_dir the geogroids geo directory
#' @param precision the precision of the computations
#' @param out_dir the directory to write files to
#' @return NULL
#' @export
#' @importFrom sf read_sf st_drop_geometry
#' @importFrom terra rast ext crs
#' @importFrom zonal weighting_grid
#' @importFrom dplyr left_join rename mutate select summarise across everything bind_cols right_join group_by
#' @importFrom data.table fwrite setnames
#' @importFrom tidyr unnest_longer
#' @importFrom RNetCDF open.nc var.get.nc
#' @importFrom stats weighted.mean

# aggregate_nwm_params(
#   gpkg           =  '/Volumes/Transcend/ngen/CONUS-hydrofabric/ngen-release/01a/2021-11-01/spatial/hydrofabric.gpkg',
#   catchment_name = "catchments",
#   flowline_name  = "flowpaths",
#   nwm_dir        =  '/Volumes/Transcend/nwmCONUS-v216',
#   precision      =  9,
#   out_dir        =  NULL
# )

aggregate_nwm_params = function(gpkg,
                                nwm_dir,
                                catchment_name,
                                flowline_name,
                                precision = 9,
                                out_dir = NULL){

  if(is.null(out_dir)){ out_dir = dirname(gpkg)}

  cats = sf::read_sf(gpkg, catchment_name)

  geogrid_vars = c('ALBEDO12M', 'GREENFRAC', 'HGT_M',
                   'LAI12M', 'SNOALB',  'SOILTEMP')

  geogrid_vars_cat = c('LU_INDEX',
                       'SOILCBOT', 'SOILCTOP',
                       'SCB_DOM', 'SCT_DOM',
                       'SLOPECAT')

  template = list(ext = ext(-2303999.62876143, 2304000.37123857, -1920000.70008381, 1919999.29991619),
                  crs = '+proj=lcc +lat_0=40.0000076293945 +lon_0=-97 +lat_1=30 +lat_2=60 +x_0=0 +y_0=0 +R=6370000 +units=m +no_defs')


  geo      = rast(file.path(nwm_dir, 'geo_em_CONUS.nc'), subds = geogrid_vars)
  ext(geo) =   template$ext
  crs(geo) =   template$crs
  geo_cat      = rast(file.path(nwm_dir, 'geo_em_CONUS.nc'), subds = geogrid_vars_cat)
  ext(geo_cat) =   template$ext
  crs(geo)     =   template$crs

  nwm_w_1000m = zonal::weighting_grid(geo, cats, "ID")
  ###
  x1 = map_vars(r = geo,     w = nwm_w_1000m, file = NULL, precision, FUN = "mean")
  x2 = map_vars(r = geo_cat, w = nwm_w_1000m, file = NULL, precision, FUN = "mode")
  out = left_join(x1, x2)
  cols = names(out)[-1]
  out[,(cols) := round(.SD, precision), .SDcols = cols]

  file = file.path(out_dir, "geo_em_CONUS.nc")
  data.table::fwrite(out, file, row.names = FALSE)
  message(file)

  ###
  soils = rast(file.path(nwm_dir, 'soilproperties_CONUS_FullRouting.nc'))
  ext(soils) = template$ext
  crs(soils) = template$crs
  map_vars(r = soils, w = nwm_w_1000m,file = file.path(out_dir, "soilproperties_FullRouting.csv"), precision)

  ###
  wrf = rast(file.path(nwm_dir, 'wrfinput_CONUS.nc'), subds = c('ISLTYP','IVGTYP', 'HGT'))
  ext(wrf) = template$ext
  crs(wrf) = template$crs
  map_vars(r = wrf, w = nwm_w_1000m, file = file.path(out_dir, "wrfinput.csv"), precision)

  ###
  fulldom = rast(file.path(nwm_dir, 'Fulldom_CONUS_FullRouting.nc'), subds = c('landuse', "LKSATFAC"))
  cols = names(fulldom)
  fulldom = execute_zonal(fulldom, geom = cats, "ID")
  names(fulldom) = c("ID", cols)
  fulldom[,(cols) := round(.SD, precision), .SDcols=cols]
  data.table::fwrite(fulldom, file.path(out_dir, "Fulldom_FullRouting.csv"), row.names = FALSE)

  ####

  crosswalk <- st_drop_geometry(read_sf(gpkg, flowline_name))

  if("LevelpathID" %in% names(crosswalk)){
    crosswalk = rename(crosswalk, main_id = LevelPathID)
  }

  crosswalk = select(crosswalk, ID, member_COMID, main_id) %>%
    mutate(comid = strsplit(member_COMID, ",")) %>%
    tidyr::unnest_longer(col = c("comid")) %>%
    mutate(ID = ID, comid = as.numeric(comid))

  gwnc = RNetCDF::open.nc(file.path(nwm_dir, 'GWBUCKPARM_CONUS_FullRouting.nc'))

  vars = c("Area_sqkm", "ComID", "Coeff", "Expon", "Zinit", "Zmax")

  gwparams = suppressMessages({
    lapply(vars, function(x) RNetCDF::var.get.nc(gwnc, x)) %>%
      bind_cols() %>%
      setnames(vars) %>%
      rename(comid = ComID) %>%
      right_join(select(crosswalk, ID, comid), by = 'comid') %>%
      group_by(ID) %>%
      summarise(across(everything(), ~ round(
        weighted.mean(.x,
                      w = .data$Area_sqkm,
                      na.rm = TRUE), precision))) %>%
      select(-comid, -Area_sqkm)
  })

  data.table::fwrite(gwparams, file.path(out_dir, "GWBUCKPARM_FullRouting.csv"), row.names = FALSE)

}
