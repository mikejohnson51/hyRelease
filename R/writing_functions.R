#' Write Core Hydrofabric Files
#' The core hydrofabric files for NextGen include (A) a flowpath edgelist, (B) a catchment geojson,
#' (C) a nexus geojson, and a geopackage which includes the catchments POLYGONs,
#' flowpaths LINESTRINGs, and nexus POINTs
#' @param gpkg a filename of a gpkg to process
#' @param catchment_name the name of the catchment layer
#' @param flowpath_name the name of the flowpath layer
#' @param ngen_pu the name of the current processing unit
#' @param spatial_path path to spat
#' @param parameter_path path to parameter set
#' @param waterbody_prefix waterbody prefix
#' @param catchment_prefix catchment prefix
#' @param nexus_prefix nexus prefix
#' @param terminal_nexus_prefix terminal nexus prefix
#' @param overwrite logical. overwrite existing files?
#' @return path to the resulting hydrofabric
#' @export
#' @importFrom sf read_sf
#' @importFrom dplyr mutate left_join rename
#' @importFrom hyAggregate get_nexus_locations

write_hydrofabric = function(gpkg,
                             catchment_name, flowpath_name,
                             ngen_pu,
                             spatial_path, parameter_path,
                             waterbody_prefix  = "wb-",
                             catchment_prefix = "cat-",
                             nexus_prefix = "nex-",
                             terminal_nexus_prefix = "tnx-",
                             overwrite = FALSE
){


  hyfab = file.path(spatial_path, "hydrofabric.gpkg")
  fp_el = file.path(parameter_path, "flowpath_edge_list.json")

  if(any(check_file(fp_el, overwrite), check_file(hyfab, overwrite))){

    fps   = sf::read_sf(gpkg, flowpath_name)
    ########## GRAPH DATA
    catchment_edge_list <- get_catchment_edges_terms(fps) |>
      mutate(ngen_pu = ngen_pu)

    fp_edge_list        <- get_catchment_edges_terms(fps, catchment_prefix = waterbody_prefix) |>
      mutate(ngen_pu = ngen_pu)

    write_json(fp_edge_list, fp_el, pretty = TRUE)

    ########## SPATIAL DATA

    nexus_data = hyAggregate::get_nexus_locations(fps) %>%
      mutate(prefix = ifelse(ID > 100000000, terminal_nexus_prefix, nexus_prefix)) %>%
      mutate(ID = paste0(prefix, ID), prefix = NULL) %>%
      left_join(catchment_edge_list, by = c("ID")) %>%
      mutate(toID  = ifelse(toID == "cat-0", NA, toID))

    catchment_data  = sf::read_sf(ref, catchment_name) %>%
      rename(area_sqkm = areasqkm) %>%
      get_catchment_data(catchment_edge_list, catchment_prefix = catchment_prefix)

    flowpath_data = get_flowpath_data(add_slope(fps), catchment_edge_list)

    hyfab = write_nextgen_spatial(flowpath_data, catchment_data, nexus_data, spatial_path)
  }

  return(hyfab)
}


#' Write routing params to JSON
#' Takes a geopackage, extracts the flowpath features, and builds a length weighted average of
#' all RouteLink variables using `length_average_routlink`. The output data ss written as JSON
#' to the supplied outfile path.
#' @param gpkg geopackage with flowline features
#' @param flowpath_name name of flowline feature layer
#' @param nwm_dir directory to NWM base files. Searches for a file called `Routelink_CONUS.nc`
#' @param outfile file path to write to
#' @return file path
#' @export
#' @importFrom sf read_sf st_drop_geometry
#' @importFrom jsonlite write_json

write_waterbody_json = function(gpkg, flowpath_name, nwm_dir, outfile){

  check_file(outfile)

   fps =  read_sf(gpkg, flowpath_name) #|>
  #  length_average_routlink(rl_path  = file.path(nwm_dir, "RouteLink_CONUS.nc"))

  req.names = c("ID", "member_COMID", "gages", "NHDWaterbodyComID",
                "toID", "main_id", "realized_catchment", "lengthkm",
                "slope_percent", "Qi", "MusK", "MusX", "n",
                "ChSlp", "BtmWdth", "time", "Kchan", "nCC", "TopWdthCC", "TopWdth", "Length_m")

  if(all(req.names %in% names(fps))){
    missing.names = req.names[!req.names %in% names(fps)]
  } else {
    stop("Need: ", paste(missing.names, collapse  = ", "))
  }

  wb_feilds <- lapply(unique(fps$ID),
                      function(x, df) {
                        df_sub <- df[df$ID == x,]
                        out = list()
                        # out$member_COMID      <- strsplit(df_sub$member_COMID, ",")[[1]]

                        # Only use the most downstream gage! Pick one...
                        out$gages               <- unique(strsplit(df_sub$gages, ",")[[1]])[1]
                        # out$NHDWaterbodyComID <- strsplit(df_sub$NHDWaterbodyComID, ",")[[1]]

                        out$toID                 <- unique(df_sub$toID)
                        # out$main_id            <- unique(df_sub$main_id)
                        # out$realized_catchment <- unique(df_sub$realized_catchment)

                        # out$length_km  <- as.numeric(df_sub$lengthkm)
                        out$So         <- as.numeric(unique(df_sub$slope_percent))
                        out$Qi         <- as.numeric(unique(df_sub$Qi))
                        out$MusK       <- as.numeric(unique(df_sub$MusK))
                        out$MusX       <- as.numeric(unique(df_sub$MusX))
                        out$n          <- as.numeric(unique(df_sub$n))
                        out$ChSlp      <- as.numeric(unique(df_sub$ChSlp))
                        out$BtmWdth    <- as.numeric(unique(df_sub$BtmWdth))

                        # out$time     <- as.numeric(unique(df_sub$time))

                        out$Kchan      <- as.numeric(unique(df_sub$Kchan))
                        out$nCC        <- as.numeric(unique(df_sub$nCC))
                        out$TopWdthCC  <- as.numeric(unique(df_sub$TopWdthCC))
                        out$TopWdth    <- as.numeric(unique(df_sub$TopWdth))
                        out$length_m   <- as.numeric(unique(df_sub$Length_m))

                        ind = unlist(lapply(1:length(out), function(i) !all(is.na(out[[i]]))))
                        out[ind]

                      }, df = fps)


  names(wb_feilds) = unique(fps$ID)

  jsonlite::write_json(
    wb_feilds,
    outfile,
    auto_unbox = TRUE,
    pretty = TRUE
  )
}

