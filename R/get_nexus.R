#' Get Nexus Locations
#'
#' @param fp a sf flowpath object
#' @return
#' @export
#' @importFrom dplyr left_join select mutate filter select bind_rows group_by ungroup arrange rename slice_head
#' @importFrom nhdplusTools get_node rename_geometry
#' @importFrom sf st_intersects st_drop_geometry st_as_sf


get_nex_locations = function(fp){

  nex = fp %>%
    left_join(st_drop_geometry(select(., toID = ID, ds_toID = toID)), by = c("toID")) %>%
    filter(ID %in% unique(toID)) %>%
    rename_geometry("geometry") %>%
    mutate(geometry = get_node(., "start")$geometry) %>%
    select(ID, toID)

  imap = st_intersects(nex, fp)

  df = data.frame(
    ID       = rep(nex$ID, times = lengths(imap)),
    touches  = fp$ID[unlist(imap)]) %>%
    mutate(cond = ifelse(ID == touches, "move","aaa")) %>%
    group_by(ID) %>%
    arrange(cond) %>%
    slice_head(n = 1) %>%
    ungroup()

  to_move = filter(fp, ID %in% filter(df, cond == "move")$ID) %>%
    select(ID) %>%
    rename_geometry("geometry")

  to_keep = filter(nex, ID %in% filter(df, cond != "move")$ID) %>%
    select(ID) %>%
    rename_geometry("geometry")

  fp_ends = bind_rows(
    select(fp, ID, Hydroseq, toID) %>%
      rename_geometry("geometry") %>%
      mutate(geometry = get_node(., "end")$geometry, pos = "end"),
    select(fp, ID, Hydroseq, toID) %>%
      rename_geometry("geometry") %>%
      mutate(geometry = get_node(., "start")$geometry, pos = "start")
  )

  imap = st_intersects(to_move, fp_ends)

  df = data.frame(
    ID       = rep(to_move$ID, times = lengths(imap)),
    touches  = fp_ends$ID[unlist(imap)],
    pos      = fp_ends$pos[unlist(imap)]) %>%
    filter(ID != touches) %>%
    left_join(rename(fp_ends, touches = ID), by = c('touches', 'pos')) %>%
    group_by(ID) %>%
    filter(ID == toID) %>%
    ungroup() %>%
    st_as_sf() %>%
    select(ID) %>%
    bind_rows(to_keep)
}

