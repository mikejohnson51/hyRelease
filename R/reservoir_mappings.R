add_slope = function(flowpaths = agg$flowpaths){

  flowpaths = flowpaths %>%
    hyRefactor::add_lengthmap(length_table = nhdplusTools::get_vaa(c("lengthkm")))

  net_map  <- dplyr::select(st_drop_geometry(flowpaths), ID, lengthMap) %>%
    mutate(comid = strsplit(.data$lengthMap, ",")) %>%
    tidyr::unnest(cols = comid) %>%
    mutate(full_comids = floor(as.numeric(comid)),
           w = 10 * (as.numeric(comid) - full_comids),
           comid = NULL)

  df = nhdplusTools::get_vaa(c('lengthkm', "slope"))

  df = df %>%
    right_join(net_map, by = c('comid' = 'full_comids')) %>%
    select(-.data$lengthMap) %>%
    mutate(w = w * lengthkm) %>%
    group_by(ID) %>%
    summarise(across(everything(), ~ round(
      weighted.mean(x = .,
                    w = w,
                    na.rm = TRUE), 8))) %>%
    dplyr::select(-comid, -lengthkm, -w)

  left_join(flowpaths, df, by = "ID")
}

#' @title NetCDF Variable Names from a path
#' @description Given a path to a NetCDF file, list the variables
#' @param path a path to a NetCDF file
#' @return a vectore of variable names
#' @export
#' @importFrom RNetCDF open.nc file.inq.nc var.inq.nc

var_names = function(path) {
  nc   <- RNetCDF::open.nc(path)
  nvar <- RNetCDF::file.inq.nc(nc)$nvar
  varnames <- character(nvar)
  for (i in seq_len(nvar)) {
    varnames[i] <- RNetCDF::var.inq.nc(nc, i - 1)$name
  }
  varnames
}


# sf = read_sf('/Users/mjohnson/Downloads/nhdplus_waterbodies.gpkg')
# rl_vars = c("link", "NHDWaterbodyComID")
#
# lake_nc = RNetCDF::open.nc(lake)
# df2 = data.frame(do.call(cbind, lapply(vars, function(x)
#   RNetCDF::var.get.nc(lake_nc, x)))) %>%
#   setnames(vars) %>%
#   filter(lake_id > 0) |>
#   mutate(lake_id = as.numeric(lake_id))
#
# lakes = st_as_sf(df2, coords = c("lon", "lat"), crs = 4326)
#
# inner_join(df2, sf, by = c("lake_id" = "COMID"))
#
# dim(df2)
#
# nwm_lakes = sf |>
#   filter(COMID %in% df2$lake_id)

#' @title Build Lake Params for NGen
#' @description Map the NWM LAKEPARM_CONUS.nc file to an aggregated NGen Network
#' @param gpkg a geopackage with aggregation units
#' @param catchment_name the layer name of the aggregation units
#' @param flowline_name the layer name of the flowpath units
#' @param waterbody_gpkg  geopackage containing CONUS NHD waterbodies
#' @param nwm_dir a directory to NWM files
#' @param outfile optional. Path to JSON file to write outputs using `write_lakes_json`
#' @importFrom RNetCDF open.nc var.get.nc
#' @importFrom data.table setnames
#' @importFrom dplyr select filter mutate left_join right_join group_by summarise slice ungroup inner_join
#' @importFrom nhdplusTools get_vaa
#' @importFrom sf read_sf st_drop_geometry st_transform st_bbox st_as_sfc st_as_text st_intersection st_collection_extract
#' @importFrom hyAggregate add_lengthkm
#' @importFrom tidyr unnest
#' @export

build_lake_params = function(flowpaths,
                             waterbody_gpkg = '/Users/mjohnson/Downloads/nhdplus_waterbodies.gpkg',
                             nwm_dir,
                             outfile = NULL) {

  lake = file.path(nwm_dir, 'LAKEPARM_CONUS.nc')

  if(!file.exists(lake)){
    stop("This workflow cannot be run without ", lake)
  }

  rl   = file.path(nwm_dir, 'RouteLink_CONUS.nc')

  if(!file.exists(rl)){
    stop("This workflow cannot be run without ", rl)
  }

  vars = var_names(lake)

  rl_vars = c("link", "NHDWaterbodyComID")

  lake_nc = RNetCDF::open.nc(lake)
  rl_nc   = RNetCDF::open.nc(rl)

  df2 = data.frame(do.call(cbind, lapply(rl_vars, function(x)
    RNetCDF::var.get.nc(rl_nc, x)))) %>%
    setnames(c('comid', 'lake_id')) %>%
    filter(lake_id > 0)

  order = nhdplusTools::get_vaa("hydroseq")

  if(!is.null(flowpaths)){
    fl = flowpaths
  } else {
    fl = read_sf(gpkg, 'flowpaths')
  }

  flowpaths = fl %>%
    st_drop_geometry() %>%
    select(.data$ID, .data$member_COMID) %>%
    mutate(comid = strsplit(.data$member_COMID, ",")) %>%
    tidyr::unnest(cols = .data$comid) %>%
    mutate(comid = floor(as.numeric(.data$comid)),
           member_COMID = NULL) %>%
    left_join(order, by = "comid") %>%
    left_join(df2,   by = "comid")

  flowpaths2 = flowpaths %>%
    group_by(.data$lake_id) %>%
    summarise(outletId = .data$ID[which.min(.data$hydroseq)]) %>%
    slice(1) %>%
    ungroup() %>%
    inner_join(flowpaths, by = "lake_id") %>%
    filter(!is.na(.data$ID)) %>%
    group_by(.data$lake_id, .data$ID) %>%
    slice(1) %>%
    ungroup()

  df = data.frame(do.call(cbind, lapply(vars, function(x)
    RNetCDF::var.get.nc(lake_nc, x)))) %>%
    setnames(vars) %>%
    mutate(lake_id = as.numeric(.data$lake_id)) %>%
    inner_join(flowpaths2, by = "lake_id") %>%
    dplyr::select(-.data$comid, -.data$hydroseq, -.data$crs)

  wkt = st_bbox(fl) %>%
    st_as_sfc() %>%
    st_transform(4326) %>%
    st_as_text()

  wbs = read_sf(waterbody_gpkg, wkt_filter = wkt) %>%
    st_transform(5070) %>%
    filter(.data$COMID  %in% df$lake_id)

  ty = suppressWarnings({
    st_intersection(fl, wbs) %>%
      select(.data$ID, .data$COMID, .data$length_km) %>%
      st_collection_extract("LINESTRING") %>%
      mutate(partial_length_km = hyAggregate::add_lengthkm(.),
             partial_length_percent = .data$partial_length_km / .data$length_km) %>%
      st_drop_geometry() %>%
      mutate(
        lake_id = .data$COMID,
        COMID = NULL
      ) %>%
      right_join(df, by = c("ID", "lake_id"))
  })

  if(!is.null(write_json)){
    write_lakes_json(ty, outfile)
  } else {
    ty
  }
}


#' Quantify area of waterbody/land per catchment
#' Ngen needs a land area percentages to drive rainfall-runoff processes. Provided a Ngen geopackage,
#' and NHD waterbody geopackage, these are computed.
#' @param gpkg Geopackage containing catchment data.
#' @param waterbody_gpkg Geopackage of NHD waterbodies
#' @param nwm_dir if supplied, only waterbodies included in the NWM Routelink file will be
#' included in the area calculations. From the NWM directory, the Routelink file is searched for.
#' @return
#' @export
#' @importFrom RNetCDF open.nc var.get.nc
#' @importFrom data.table setnames
#' @importFrom dplyr filter mutate select group_by summarise ungroup ungroup right_join everything
#' @importFrom sf read_sf st_transform st_bbox st_as_sfc st_as_text st_intersection st_drop_geometry st_as_sf
#' @importFrom hyRefactor add_areasqkm
#' @importFrom nhdplusTools rename_geometry

catchment_waterbody_interaction = function(
    gpkg           = '/Users/mjohnson/Downloads/2021-10-28/spatial/hydrofabric.gpkg',
    catchment_name = "catchments",
    catchments     = NULL,
    waterbody_gpkg = NULL,
    nwm_dir        = NULL){

  if(is.null(catchments)){
    cats = read_sf(gpkg, catchment_name) %>%
      st_transform(5070)
  } else {
    cats = catchments
  }

  wkt = st_bbox(cats) %>%
    st_as_sfc() %>%
    st_transform(4326) %>%
    st_as_text()

  wbs = read_sf(waterbody_gpkg, wkt_filter = wkt) %>%
    st_transform(5070)

  rl   = file.path(nwm_dir, 'RouteLink_CONUS.nc')

 if(file.exists(rl)){

    rl_vars = c("link", "NHDWaterbodyComID")
    rl_nc   = RNetCDF::open.nc(rl)

    df2 = data.frame(do.call(cbind, lapply(rl_vars, function(x)
      RNetCDF::var.get.nc(rl_nc, x)))) %>%
      setnames(c('comid', 'lake_id')) %>%
      filter(lake_id > 0)

    wbs =  filter(wbs, .data$COMID  %in% df2$lake_id)
   }


  ty = suppressWarnings({
    st_intersection(cats, wbs) %>%
        dplyr::select(.data$ID) %>%
        dplyr::mutate(wb_area_sqkm = hyRefactor::add_areasqkm(.)) %>%
        st_drop_geometry() %>%
        dplyr::group_by(ID) %>%
        dplyr::summarise(wb_area_sqkm = sum(.data$wb_area_sqkm)) %>%
        dplyr::ungroup() %>%
        right_join(cats, by = "ID") %>%
    mutate(wb_area_sqkm = ifelse(is.na(.data$wb_area_sqkm), 0, .data$wb_area_sqkm),
           land_area_sqkm = .data$area_sqkm - .data$wb_area_sqkm) %>%
    st_as_sf() %>%
    nhdplusTools::rename_geometry("geometry") %>%
    select(.data$ID, .data$area_sqkm, .data$land_area_sqkm, .data$wb_area_sqkm, dplyr::everything())
  })
}

