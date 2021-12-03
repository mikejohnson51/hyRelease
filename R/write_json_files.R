#' Write lake params to JSON
#' @param df output of `build_lake_params`
#' @param outfile file path to write to
#' @return file path
#' @export
#' @importFrom dplyr left_join select distinct
#' @importFrom jsonlite write_json

write_lakes_json = function(df, outfile){

  ind = data.frame(lake_id = unique(df$lake_id)) %>%
    left_join(select(df, lake_id, outletId)) %>%
    dplyr::distinct()


  filter(df, outletId == 'wb-16278.1')

  lake_crosswalk_list <- lapply(unique(df$lake_id),
                                function(x, df) {
                                  df_sub <- df[df$lake_id == x,]
                                  out <-
                                    list(member_wbs = df_sub$ID)
                                  out$partial_length_percent <-df_sub$partial_length_percent

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

                                }, df = df)


  names(lake_crosswalk_list) <- ind$outletId


  grep('wb-16278.1', names(lake_crosswalk_list), value = TRUE)


  jsonlite::write_json(
    lake_crosswalk_list,
    outfile,
    pretty = TRUE,
    auto_unbox = TRUE
  )
}

#' Write routing params to JSON
#' Takes a geopackage, extracts the flowpath features, and builds a lengthweighted average of
#' all RouteLink variables using `length_average_routlink`. The output data as written as JSON
#' to the supplied outfile path.
#' @param gpkg geopackage with flowline features
#' @param fl_name name of flowline feature layer
#' @param nwm_dir directory to NWM base files. Searched for a file called `Routelink_CONUS.nc`
#' @param outfile file path to write to
#' @return file path
#' @export
#' @importFrom sf read_sf st_drop_geometry
#' @importFrom jsonlite write_json

write_waterbody_json = function(gpkg,
                                flowline_name,
                                flowpaths = NULL,
                                fl_name = "flowpaths",
                                nwm_dir = "/Users/mjohnson/Downloads/",
                                outfile){

  if(!is.null(flowpaths)){
    fl = flowpaths
  } else {
    fl = read_sf(gpkg, 'flowpaths')
  }


  fl = fl %>%
    length_average_routlink(
      rl_vars = c(
        "link", "Qi", "MusK", "MusX", "n", "So",  "ChSlp",
        "BtmWdth", "time", "Kchan", "nCC", "TopWdthCC", "TopWdth"),
      rl_path  = file.path(nwm_dir, "RouteLink_CONUS.nc"))  %>%
     st_drop_geometry()

 wb_feilds <- lapply(unique(fl$ID),
                                function(x, df) {
                                  df_sub <- df[df$ID == x,]
                                  out = list()
                                  out$member_COMID      <- strsplit(df_sub$member_COMID, ",")
                                  out$gages             <- strsplit(df_sub$gages, ",")
                                  out$NHDWaterbodyComID <- strsplit(df_sub$NHDWaterbodyComID, ",")

                                  out$toID       <- unique(df_sub$toID)
                                  out$main_id    <- unique(df_sub$main_id)
                                  out$realized_catchment <- unique(df_sub$realized_catchment)

                                  out$length_km  <- as.numeric(df_sub$length_km)
                                  out$slope_percentage   <- as.numeric(unique(df_sub$slope_percent))
                                  out$Qi      <- as.numeric(unique(df_sub$Qi))
                                  out$MusK    <- as.numeric(unique(df_sub$MusK))
                                  out$MusX    <- as.numeric(unique(df_sub$MusX))
                                  out$n       <- as.numeric(unique(df_sub$n))
                                  out$ChSlp   <- as.numeric(unique(df_sub$ChSlp))
                                  out$BtmWdth <- as.numeric(unique(df_sub$BtmWdth))
                                  out$time    <- as.numeric(unique(df_sub$time))
                                  out$Kchan   <- as.numeric(unique(df_sub$Kchan))
                                  out$nCC     <- as.numeric(unique(df_sub$nCC))
                                  out$TopWdthCC  <- as.numeric(unique(df_sub$TopWdthCC))
                                  out$TopWdth    <- as.numeric(unique(df_sub$TopWdth))
                                  out$Length_m   <- as.numeric(unique(df_sub$Length_m))

                                  ind = unlist(lapply(1:length(out), function(i) !is.na(out[[i]])))
                                  out[ind]

                                }, df = fl)


  names(wb_feilds) = fl$ID

  jsonlite::write_json(
    wb_feilds,
    outfile,
    pretty = TRUE,
    auto_unbox = TRUE
  )
}





#' Write NWIS mapping to JSON
#' Takes a geopackage or flowline layer, extracts the flowpath features,
#' and builds a lengthweighted average of all RouteLink variables using `length_average_routlink`.
#' The output data as written as JSON to the supplied outfile path.
#' @param gpkg geopackage with flowline features
#' @param fl_name name of flowline feature layer
#' @param nwm_dir directory to NWM base files. Searched for a file called `Routelink_CONUS.nc`
#' @param outfile file path to write to
#' @return file path
#' @export
#' @importFrom sf read_sf st_drop_geometry
#' @importFrom jsonlite write_json

write_nwis_crosswalk = function(gpkg,
                                flowline_name,
                                flowpaths = NULL,
                                fl_name = "flowpaths",
                                nwm_dir = "/Users/mjohnson/Downloads/",
                                outfile){

  if(!is.null(flowpaths)){
    fl_raw = flowpaths
  } else {
    fl_raw = read_sf(gpkg, 'flowpaths')
  }

  fl = fl_raw %>%
    st_drop_geometry() %>%
    select(.data$ID, .data$member_COMID) %>%
    mutate(comid = strsplit(.data$member_COMID, ",")) %>%
    tidyr::unnest(cols = .data$comid) %>%
    mutate(comid = floor(as.numeric(.data$comid)),
           member_COMID = NULL)

  nc = RNetCDF::open.nc(file.path(nwm_dir, "RouteLink_CONUS.nc"))

  df = suppressMessages({
    lapply(c("link", "gages", 'NHDWaterbodyComID'), function(x) x = RNetCDF::var.get.nc(nc, x)) %>%
      bind_cols()
  })

  names(df) = c("comid", "gages", 'NHDWaterbodyComID')

  df2 = df %>%
    right_join(fl, by = 'comid') %>%
    mutate(gages = trimws(gages),
           gages = ifelse(gages == "", NA, gages),
           NHDWaterbodyComID = ifelse(NHDWaterbodyComID == -9999, NA, NHDWaterbodyComID)
    ) %>%
    group_by(ID) %>%
    summarise(gages = list(gages[!is.na(gages)]),
              NHDWaterbodyComID = list(unique(NHDWaterbodyComID[!is.na(NHDWaterbodyComID)]))) %>%
    mutate(gages = ifelse(gages == "", NA, gages),
           NHDWaterbodyComID = ifelse(NHDWaterbodyComID == "", NA, NHDWaterbodyComID))


nhd_crosswalk <- st_drop_geometry(fl_raw) %>%
  select(ID, member_COMID, main_id) %>%
  mutate(comid = strsplit(member_COMID, ",")) %>%
  tidyr::unnest_longer(col = c("comid")) %>%
  mutate(ID = ID,
         comid = as.numeric(comid)) %>%
  select(ID, comid, main_id) %>%
  left_join(nhdplusTools::get_vaa("hydroseq"), by = "comid") %>%
  group_by(ID) %>%
  arrange(hydroseq) %>%
  mutate(outlet_comid = dplyr::first(comid)) %>%
  left_join(df2, by = "ID") %>%
  ungroup()



nhd_crosswalk_list <- lapply(unique(nhd_crosswalk$ID),
                             function(x, df) {
                               df_sub <- df[df$ID == x, ]
                               out <- list(member_comids = df_sub$comid)

                               out$gages <- unique(unlist(df_sub$gages))

                               out$waterbodies = unique(unlist(df_sub$NHDWaterbodyComID))

                               out$outlet_comid <- unique(df_sub$outlet_comid)
                               out$main = unique(df_sub$main_id)

                               ind = lapply(1:length(out), function(x){ sum(lengths(out[[x]]) != 0) } )

                               out[ind > 0]

                             }, df = nhd_crosswalk)


names(nhd_crosswalk_list) <- unique(nhd_crosswalk$ID)

jsonlite::write_json(
  nhd_crosswalk_list,
  outfile,
  pretty = TRUE,
  auto_unbox = TRUE
)

}
