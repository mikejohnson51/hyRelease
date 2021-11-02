#' Get catchment edge list
#' @description get a edge list for catchments
#' @param flowpaths  sf data.frame containing hyRefactor or hyAggregate output.
#' @param nexus_prefix  character prefix for nexus IDs
#' @param terminal_nexus_prefix character prefix for terminal nexus IDs
#' @param catchment_prefix character prefix for catchment IDs
#' @param cutoff terminal IDs begin above a defiend threshold
#' @return data.frame
#' @export
#' @importFrom sf st_drop_geometry
#' @importFrom dplyr select mutate left_join bind_rows

get_catchment_edges_terms = function(flowpaths,
                                     nexus_prefix = "nex-",
                                     terminal_nexus_prefix = "tnx-",
                                     catchment_prefix = "cat-",
                                     cutoff = 100000000){
  fline = select(st_drop_geometry(flowpaths), ID, toID)

  obj1 = fline %>%
    select(ID = .data$ID, toID = .data$toID) %>%
    mutate(ID = paste0(catchment_prefix, .data$ID),
           toID = paste0(ifelse(.data$toID > cutoff, terminal_nexus_prefix, nexus_prefix), .data$toID))

  obj2 =  data.frame(ID = unique(fline$toID)) %>%
    left_join(mutate(select(fline, ID), toID = ID), by = "ID") %>%
    mutate(toID = ifelse(is.na(.data$toID), 0, .data$toID)) %>%
    mutate(ID =  paste0(ifelse(.data$ID > cutoff, terminal_nexus_prefix, nexus_prefix), .data$ID),
           toID = paste0(catchment_prefix, .data$toID))

  bind_rows(obj1, obj2)
}

#' Get waterbody edge list
#' @description get a edge list for waterbodies
#' @param flowpaths  sf data.frame containing hyRefactor or hyAggregate output.
#' @param wb_prefix  character prefix for waterbody IDs
#' @param terminal_wb_prefix character prefix for terminal waterbody IDs
#' @param cutoff terminal IDs begin above a defined threshold
#' @return data.frame
#' @export
#' @importFrom sf st_drop_geometry
#' @importFrom dplyr select mutate

get_waterbody_edges_terms = function(flowpaths,
                                     wb_prefix = "nex-",
                                     terminal_wb_prefix = "tnx-",
                                     cutoff = 100000000){

  fline = select(st_drop_geometry(flowpaths), ID, toID)

  fline %>% select(.data$ID, .data$toID) %>%
    mutate(ID = paste0(ifelse(.data$ID > cutoff, terminal_wb_prefix, wb_prefix), .data$ID),
           toID = paste0(ifelse(.data$toID > cutoff, terminal_wb_prefix, wb_prefix), .data$toID))
}


get_catchment_data = function(catchment, catchment_edge_list, catchment_prefix = "catchment_") {

  select(catchment, ID = .data$ID, area_sqkm = .data$area_sqkm) %>%
    mutate(ID = paste0(catchment_prefix, .data$ID)) %>%
    left_join(catchment_edge_list,  by = "ID")
}


get_flowpath_data = function(fline, catchment_edge_list, catchment_prefix = "catchment_") {

  select(fline, ID = .data$ID, length_km = .data$length_km,
         slope_percent = .data$slope, main_id = .data$LevelPathID, member_COMID = .data$member_COMID) %>%
    mutate(ID = paste0(catchment_prefix, .data$ID)) %>%
    left_join(catchment_edge_list, by = "ID")
}
