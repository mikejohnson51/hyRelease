#' Flush existing ID prefixes
#' Given a data object and column, remove a prefix and adjoining "-"
#' @param input input data object
#' @param col column to remove prefix from
#' @return data object with updated column
#' @export

flush_prefix = function(input, col){
  for(i in col){
    input[[i]] = as.numeric(gsub(".*-","",input[[i]]))
  }
  input
}




#' Build catchment, flowpath and waterbody edgelists
#' Returns catchment, flowpath and waterbody edgelists. If `outdir` is supplied,
#' then JSON files are written there.
#' @param flowpaths  sf data.frame containing hyRefactor or hyAggregate output.
#' @param nexus_prefix  character prefix for nexus IDs
#' @param terminal_nexus_prefix character prefix for terminal nexus IDs
#' @param wb_prefix  character prefix for waterbody IDs
#' @param terminal_wb_prefix character prefix for terminal waterbody IDs
#' @param catchment_prefix character prefix for catchment IDs
#' @param cutoff terminal IDs begin above a defiend threshold
#' @param outdir optional. Directory to write JSON files
#' @return list of edgelist
#' @export
#' @importFrom jsonlite write_json

build_edge_lists = function(flowpaths,
                            nexus_prefix = "nex-",
                            terminal_nexus_prefix = "tnx-",
                            wb_prefix = "wb-",
                            terminal_wb_prefix = "twb-",
                            catchment_prefix = "cat-",
                            cutoff = 100000000,
                            outdir) {
  edges = list(
    catchment_edge_list  =
      get_catchment_edges_terms(
        flowpaths,
        nexus_prefix = nexus_prefix,
        terminal_nexus_prefix = terminal_nexus_prefix,
        catchment_prefix = catchment_prefix,
        cutoff = cutoff
      ),
    fp_edge_list         =
      get_catchment_edges_terms(
        flowpaths,
        nexus_prefix = nexus_prefix,
        terminal_nexus_prefix = terminal_nexus_prefix,
        catchment_prefix = wb_prefix,
        cutoff = cutoff
      ),
    waterbody_edge_list  = get_waterbody_edges_terms(
      flowpaths,
      wb_prefix = wb_prefix,
      terminal_wb_prefix = terminal_wb_prefix,
      cutoff = cutoff
    )
  )

  if (!is.null(outdir)) {
    message("Writing JSON files to: ", outdir)
    jsonlite::write_json(
      edges$catchment_edge_list,
      file.path(outdir,  "catchment_edge_list.json"),
      pretty = TRUE
    )
    jsonlite::write_json(
      edges$waterbody_edge_list,
      file.path(outdir, "waterbody_edge_list.json"),
      pretty = TRUE
    )
    jsonlite::write_json(edges$fp_edge_list,
                         file.path(outdir, "flowpath_edge_list.json"),
                         pretty = TRUE)
  }

  return(edges)

}

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

  fline = flush_prefix(fline, "ID")
  fline = flush_prefix(fline, "toID")

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
                                     wb_prefix = "wb-",
                                     terminal_wb_prefix = "twb-",
                                     cutoff = 100000000){

  fline = select(st_drop_geometry(flowpaths), ID, toID)

  fline = flush_prefix(fline, "ID")
  fline = flush_prefix(fline, "toID")

  fline %>% select(.data$ID, .data$toID) %>%
    mutate(ID = paste0(ifelse(.data$ID > cutoff, terminal_wb_prefix, wb_prefix), .data$ID),
           toID = paste0(ifelse(.data$toID > cutoff, terminal_wb_prefix, wb_prefix), .data$toID))
}



get_catchment_data = function(catchment, catchment_edge_list, catchment_prefix = "cat-") {

  catchment %>%
    mutate(ID = paste0(catchment_prefix, .data$ID)) %>%
    left_join(catchment_edge_list,  by = "ID")
}


get_flowpath_data = function(fline, catchment_edge_list, catchment_prefix = "cat-") {

  if("main_id" %in% names(fline)){
    fline = rename(fline, LevelPathID = main_id)
  }

  select(fline,
         ID = .data$ID,
         length_km = .data$length_km,
         slope_percent = .data$slope,
         main_id = .data$LevelPathID,
         member_COMID = .data$member_COMID) %>%
    mutate(ID = paste0(catchment_prefix, .data$ID)) %>%
    left_join(catchment_edge_list, by = "ID")
}
