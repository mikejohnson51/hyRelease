#' @title Add length averaged slope to flowpaths
#' @description Add length weighted slopes to refactored/aggregated flowpaths
#' using the base NHDPlus slope estimates
#' @param flowpaths sf LINESTRING
#' @return sf LINSETRING
#' @export
#' @importFrom hyRefactor add_lengthmap
#' @importFrom nhdplusTools get_vaa
#' @importFrom dplyr select mutate group_by summarise across everything left_join
#' @importFrom stats weighted.mean
#' @importFrom tidyr unnest

add_slope = function(flowpaths){

  .data <- NULL

  flowpaths = flowpaths %>%
    hyRefactor::add_lengthmap(length_table = nhdplusTools::get_vaa(c("lengthkm")))

  net_map  <- dplyr::select(st_drop_geometry(flowpaths), .data$ID, .data$lengthMap) %>%
    mutate(comid = strsplit(.data$lengthMap, ",")) %>%
    tidyr::unnest(cols = .data$comid) %>%
    mutate(full_comids = floor(as.numeric(.data$comid)),
           w = 10 * (as.numeric(.data$comid) - .data$full_comids),
           comid = NULL)

  df = nhdplusTools::get_vaa(c('lengthkm', "slope"))

  df = df %>%
    right_join(net_map, by = c('comid' = 'full_comids')) %>%
    select(-.data$lengthMap) %>%
    mutate(w = .data$w * .data$lengthkm) %>%
    group_by(.data$ID) %>%
    summarise(across(everything(), ~ round(
      weighted.mean(x = .,
                    w = .data$w,
                    na.rm = TRUE), 8))) %>%
    dplyr::select(-.data$comid, -.data$lengthkm, -.data$w)

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
#' @param out_file optional. Path to JSON file to write outputs using `write_lakes_json`
#' @importFrom RNetCDF open.nc var.get.nc
#' @importFrom data.table setnames
#' @importFrom dplyr select filter mutate left_join right_join group_by summarise slice ungroup inner_join
#' @importFrom nhdplusTools get_vaa
#' @importFrom sf read_sf st_drop_geometry st_transform st_bbox st_as_sfc st_as_text st_intersection st_collection_extract
#' @importFrom hyAggregate add_lengthkm
#' @importFrom tidyr unnest
#' @export

build_lake_params = function(gpkg,
                             catchment_name,
                             flowline_name,
                             waterbody_gpkg = '/Users/mjohnson/Downloads/nhdplus_waterbodies.gpkg',
                             nwm_dir,
                             out_file = NULL) {

    .data <- . <-  NULL

    lake = file.path(nwm_dir, 'LAKEPARM_CONUS.nc')

    if(!file.exists(lake)){
      stop("This workflow cannot be run without ", lake)
    }

    rl   = file.path(nwm_dir, 'RouteLink_CONUS.nc')

    if(!file.exists(rl)){
      stop("This workflow cannot be run without ", rl)
    }

    fl      = read_sf(gpkg, flowline_name)
    vars    = var_names(lake)
    rl_vars = c("link", "NHDWaterbodyComID")

    lake_nc = RNetCDF::open.nc(lake)
    rl_nc   = RNetCDF::open.nc(rl)

    rl_data = data.frame(do.call(cbind, lapply(rl_vars, function(x)
      RNetCDF::var.get.nc(rl_nc, x)))) %>%
      setnames(c('comid', 'lake_id')) %>%
      filter(.data$lake_id > 0) |>
      left_join(nhdplusTools::get_vaa("hydroseq"), by = "comid")

    map_file = fl  %>%
      st_drop_geometry()  %>%
      select(.data$ID, .data$member_COMID) %>%
      mutate(comid = strsplit(.data$member_COMID, ",")) %>%
      tidyr::unnest(cols = .data$comid) %>%
      mutate(comid = floor(as.numeric(.data$comid)), member_COMID = NULL) %>%
      inner_join(rl_data,   by = "comid")   %>%
      group_by(.data$lake_id) %>%
      slice_min(hydroseq) %>%
      ungroup()


    if(nrow(map_file) == 0 ){
      jsonlite::write_json(NULL, out_file)
    } else {

    lake_data = data.frame(do.call(cbind, lapply(vars, function(x)
      RNetCDF::var.get.nc(lake_nc, x)))) %>%
      setnames(vars) %>%
      mutate(lake_id = as.numeric(.data$lake_id)) %>%
      inner_join(map_file, by = "lake_id") %>%
      dplyr::select(-.data$comid, -.data$hydroseq, -.data$crs)

    wkt = st_bbox(fl) %>%
      st_as_sfc() %>%
      st_transform(4326) %>%
      st_as_text()

    wbs = read_sf(waterbody_gpkg, wkt_filter = wkt) %>%
      st_transform(5070) %>%
      filter(.data$COMID  %in% lake_data$lake_id)

    ty = suppressWarnings({
      st_intersection(fl, wbs) %>%
        select(.data$ID, .data$COMID, .data$lengthkm, .data$GNIS_NAME) %>%
        st_collection_extract("LINESTRING") %>%
        mutate(partial_length_km = hyAggregate::add_lengthkm(.),
               partial_length_percent = 100* (.data$partial_length_km / .data$lengthkm)) %>%
        st_drop_geometry() %>%
        filter(partial_length_km > .1) |>
        mutate(
          lake_id = .data$COMID,
          COMID = NULL
        ) %>%
        right_join(lake_data, by = c("ID", "lake_id"))
    })

  ind = data.frame(lake_id = unique(lake_data)) %>%
    left_join(select(lake_data, lake_id, outletID)) %>%
    dplyr::distinct()

  lake_crosswalk_list <- lapply(unique(lake_data),
                                function(x, df) {
                                  df_sub <- df[df$lake_id == x,]
                                  out <- list(member_wbs = df_sub$ID)
                                  out$partial_length_percent <- df_sub$partial_length_percent

                                  out$lake_id    <- unique(df_sub$lake_id)
                                  out$Dam_Length <- as.numeric(unique(df_sub$Dam_Length))
                                  out$LkMxE      <- as.numeric(unique(df_sub$LkMxE))
                                  out$OrificeE   <- as.numeric(unique(df_sub$OrificeE))
                                  out$WeirE      <- as.numeric(unique(df_sub$WeirE))
                                  out$LkArea     <- as.numeric(unique(df_sub$LkArea))
                                  out$WeirC      <- as.numeric(unique(df_sub$WeirC))
                                  out$WeirL      <- as.numeric(unique(df_sub$WeirL))
                                  out$OrificeC   <- as.numeric(unique(df_sub$OrificeC))
                                  out$OrificeA   <- as.numeric(unique(df_sub$OrificeA))
                                  out$lat        <- as.numeric(unique(df_sub$lat))
                                  out$lon        <- as.numeric(unique(df_sub$lon))
                                  out$time       <- as.numeric(unique(df_sub$time))

                                  out$ascendingIndex <- as.numeric(unique(df_sub$ascendingIndex))
                                  out$ifd <- as.numeric(unique(df_sub$ifd))

                                  out

                                }, df = lake_data)


  names(lake_crosswalk_list) <- ind$outletID

  out = jsonlite::toJSON(lake_crosswalk_list, pretty = TRUE)

  jsonlite::write_json(out, out_file)
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
    waterbody_gpkg = NULL,
    nwm_dir        = NULL){

    cats = read_sf(gpkg, catchment_name) %>% st_transform(5070)

    wkt = st_bbox(cats) %>%
      st_as_sfc() %>%
      st_transform(4326) %>%
      st_as_text()

    wbs = read_sf(waterbody_gpkg, wkt_filter = wkt) %>% st_transform(5070)

    rl   = file.path(nwm_dir, 'RouteLink_CONUS.nc')

 if(file.exists(rl)){

    rl_vars = c("link", "NHDWaterbodyComID")
    rl_nc   = RNetCDF::open.nc(rl)

    rl_data = data.frame(do.call(cbind, lapply(rl_vars, function(x)
      RNetCDF::var.get.nc(rl_nc, x)))) %>%
      setnames(c('comid', 'lake_id')) %>%
      filter(lake_id > 0)

    wbs =  filter(wbs, .data$COMID  %in% rl_data$lake_id)
   }


  if(nrow(wbs) > 0)
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

