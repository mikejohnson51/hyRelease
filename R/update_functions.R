omit.na = function(x){ x[!is.na(x)] }

write_nextgen_spatial = function(flowpaths, catchments, nexus, spatial_path){
  hyfab = file.path(spatial_path, "hydrofabric.gpkg")
  sf::write_sf(catchments, hyfab, 'catchments', overwrite = TRUE)
  sf::write_sf(flowpaths, hyfab, 'flowpaths', overwrite = TRUE)
  sf::write_sf(nexus, hyfab, 'nexus', overwrite = TRUE)

  write_geojson(catchment_data, file.path(spatial_path, "catchment_data.geojson"))
  write_geojson(nexus_data,     file.path(spatial_path, "nexus_data.geojson"))

  return(hyfab)
}

gpkg_layers = function(path, goal, pattern = NULL){

  if(!file.exists(path)){ return( FALSE )}

  ll = tryCatch({sf::st_layers(path)}, error = function(e){ NULL })
  name = if(!is.null(pattern)){
    grep(pattern, ll$name, value = TRUE)
  } else {
    ll$name
  }

  if(!is.null(goal)){
    return(length(name) == goal)
  } else {
    return(name)
  }
}

file.path.build = function(..., fsep = .Platform$file.sep){
  tmp.dir = .Internal(file.path(list(...), fsep))
  dir.create(tmp.dir, recursive = TRUE, showWarnings = FALSE)
  tmp.dir
}

get_UT_reference = function (network,
                             reference_fabric_dir,
                             comid,
                             outfile  = NULL) {

  reference_fabric = list.files(reference_fabric_dir, full.names = TRUE, pattern = "gpkg$")

  if(length(reference_fabric) == 0 ){
    stop("Reference Fabric not found", call. = FALSE)
  }

  master = dplyr::filter(nhdplusTools::get_vaa("rpuid"),
                         comid %in% nhdplusTools::get_UT(network, !!comid))

  RPUS = omit.na(unique(master$rpuid))

  cats <- fps <-  list()

  for (i in 1:length(RPUS)) {
    fin = grep(RPUS[i], reference_fabric, value = TRUE)


    db <- DBI::dbConnect(RSQLite::SQLite(), fin)
    mas = dplyr::filter(master, .data$rpuid == RPUS[i])$comid

    fps[[i]] = tbl(db, "flowpaths") %>%
      dplyr::filter(.data$COMID %in% mas) %>%
      dplyr::collect() %>%
      sf::st_as_sf(crs = 4326)

    cats[[i]] = dplyr::tbl(db, "nhd_catchment") %>%
      dplyr::filter(.data$featureid %in% mas) %>%
      dplyr::collect() %>%
      sf::st_as_sf(crs = 4326)

    DBI::dbDisconnect(db)
  }

  if(!is.null(outfile)){
    write_sf(dplyr::bind_rows(fps),  outfile, "reference_flowpaths")
    write_sf(dplyr::bind_rows(cats), outfile, "reference_catchments")
    return(outfile)
  } else {
    return(list(fps = dplyr::bind_rows(fps), cats = dplyr::bind_rows(cats)))
  }
}

refactor_wrapper = function (flowpaths, catchments,
                             events = NULL, avoid = NULL,
                             split_flines_meters = 10000,
                             collapse_flines_meters = 1000,
                             collapse_flines_main_meters = 1000,
                             cores = 1,
                             facfdr = NULL,
                             routing = NULL, keep = 0.9, outfile) {

  tf <- tempfile(pattern = "refactored", fileext = ".gpkg")
  tr <- tempfile(pattern = "reconciled", fileext = ".gpkg")

  if (!is.null(events)) {
    events = dplyr::filter(events, .data$COMID %in% flowpaths$COMID)
  }

  if (!is.null(avoid)) {
    avoid = avoid[avoid %in% flowpaths$COMID]
  }

  hyRefactor::refactor_nhdplus(nhdplus_flines = flowpaths, split_flines_meters = split_flines_meters,
                   split_flines_cores = 1, collapse_flines_meters = collapse_flines_meters,
                   collapse_flines_main_meters = collapse_flines_main_meters,
                   out_refactored = tf, out_reconciled = tr, three_pass = TRUE,
                   purge_non_dendritic = FALSE, events = events, exclude_cats = avoid,
                   warn = FALSE)

  rec = st_transform(read_sf(tr), 5070)

  if (!is.null(routing)) {
    rec$order = nhdplusTools::get_streamorder(st_drop_geometry(select(rec, .data$ID, .data$toID)), status = FALSE)
    rec = rec %>%
      rename(length_km = .data$lengthkm) %>%
      length_average_routlink(rl_vars = c("link", "Qi",
                                          "MusK", "MusX", "n", "So", "ChSlp", "BtmWdth",
                                          "time", "Kchan", "nCC", "TopWdthCC", "TopWdth"),
                              rl_path = routing)
  }

  if (!is.null(facfdr)) {

    rpus         = omit.na(unique(flowpaths$RPUID))
    fdrfac_files = list.files(facfdr, pattern = rpus, full.names = TRUE)

    if ("featureid" %in% names(catchments)) {
      catchments = dplyr::rename(catchments, FEATUREID = .data$featureid)
    }

    divides <- hyRefactor::reconcile_catchment_divides(catchment = catchments,
                                                       fline_ref = sf::read_sf(tf),
                                                       fline_rec = rec,
                                                       fdr = terra::rast(grep("_fdr.tif$", fdrfac_files, value = TRUE)),
                                                       fac = terra::rast(grep("_fac.tif$", fdrfac_files, value = TRUE)),
                                                       para = cores,
                                                       cache = NULL,
                                                       fix_catchments = TRUE)
  }

  unlink(list(tr, tf))

  if(!is.null(outfile)){
    write_sf(st_transform(rec, 5070),     outfile, "refactored_flowpaths",  overwrite = TRUE)
    write_sf(st_transform(divides, 5070), outfile, "refactored_catchments", overwrite = TRUE)
  } else {
    list(fps  = st_transform(rec, 5070), cats = st_transform(divides, 5070))
  }

}

length_average_routlink = function (flowpaths,
                                    rl_vars = c("link", "Qi", "MusK", "MusX", "n", "So", "ChSlp", "BtmWdth",
                                                "time", "Kchan", "nCC", "TopWdthCC", "TopWdth"),
                                    rl_path) {

  flowpaths =  hyRefactor::add_lengthmap(flowpaths, length_table = nhdplusTools::get_vaa("lengthkm"))

  if (!"Length" %in% rl_vars) {
    rl_vars = c("Length", rl_vars)
  }

  net_map <- dplyr::select(st_drop_geometry(flowpaths), .data$ID, .data$lengthMap) %>%
    mutate(comid = strsplit(.data$lengthMap, ",")) %>%
    tidyr::unnest(cols = .data$comid) %>%
    mutate(full_comids = floor(as.numeric(.data$comid)),
           w = 10 * (as.numeric(comid) - full_comids),
           comid = NULL)


  nc = RNetCDF::open.nc(rl_path)
  on.exit(RNetCDF::close.nc(nc))
  df = data.frame(do.call(cbind, lapply(rl_vars, function(x) RNetCDF::var.get.nc(nc, x))))
  names(df) = rl_vars

  df = df %>%
    rename(comid = .data$link) %>%
    right_join(net_map, by = c(comid = "full_comids")) %>%
    select(-.data$lengthMap) %>%
    mutate(w = .data$w * .data$Length) %>%
    group_by(.data$ID) %>%
    summarise(across(everything(), ~ round(
      weighted.mean(x = .,
                    w = .data$w, na.rm = TRUE), 3
    ))) %>%
    dplyr::select(-.data$comid, -.data$Length, -.data$w)

  df2 = suppressMessages({
    lapply(c("link", "gages", "NHDWaterbodyComID"), function(x)
      x = RNetCDF::var.get.nc(nc,
                              x)) %>%
      bind_cols()
  })

  names(df2) = c("comid", "gages", "NHDWaterbodyComID")

  df2 = df2 %>%
    right_join(net_map, by = c('comid' = "full_comids")) %>%
    mutate(gages = trimws(.data$gages),
           gages = ifelse(.data$gages == "", NA, .data$gages),
           NHDWaterbodyComID = ifelse(NHDWaterbodyComID == -9999, NA, NHDWaterbodyComID)) %>%
    mutate(NHDWaterbodyComID = as.numeric(.data$NHDWaterbodyComID)) %>%
    group_by(.data$ID) %>%
    summarise(gages = paste(.data$gages[!is.na(.data$gages)], collapse = ","),
              NHDWaterbodyComID = paste(unique(.data$NHDWaterbodyComID[!is.na(.data$NHDWaterbodyComID)]), collapse = ",")) %>%
  left_join(df) %>%
    mutate(gages = ifelse(.data$gages == "", NA, .data$gages),
           NHDWaterbodyComID = ifelse(.data$NHDWaterbodyComID == "", NA, .data$NHDWaterbodyComID)) |>
    mutate(gages = as.character(.data$gages),
           NHDWaterbodyComID = as.character(.data$NHDWaterbodyComID))

  left_join(flowpaths, df2, by = "ID") %>%
    mutate(Length_m = sf::st_length(.))
}


write_nwis_crosswalk2 = function(flowpaths, catchments, gages_iii, outfile){

  nwis_sites = read_sf(gages_iii) |>
    select(Gage_no, COMID) |>
    st_transform(st_crs(flowpaths))

  nwis_sites = sf::st_join(catchment_data, nwis_sites) |>
    select(ID, Gage_no) |>
    st_drop_geometry() |>
    filter(!is.na(Gage_no)) |>
    mutate(ID = as.numeric(gsub(".*-","", ID)))


  nhd_crosswalk <- st_drop_geometry(flowpaths) %>%
    select(ID, member_COMID, main_id) %>%
    mutate(comid = strsplit(member_COMID, ",")) %>%
    tidyr::unnest_longer(col = c("comid")) %>%
    mutate(comid = as.numeric(comid),
           ID = as.numeric(gsub(".*-","", ID))) %>%
    select(ID, comid, main_id) %>%
    left_join(nhdplusTools::get_vaa("hydroseq"), by = "comid") %>%
    group_by(ID) %>%
    arrange(hydroseq) %>%
    mutate(outlet_comid = dplyr::first(comid)) %>%
    left_join(nwis_sites, by = "ID") %>%
    ungroup()

  nhd_crosswalk_list <- lapply(unique(nhd_crosswalk$ID),
                               function(x, df) {
                                 df_sub <- df[df$ID == x,]
                                 out <-
                                   list(member_comids = df_sub$comid)
                                 if (any(!is.na(df_sub$Gage_no))) {
                                   out$Gage_no <- unique(df_sub$Gage_no[!is.na(df_sub$Gage_no)])
                                 }
                                 out$outlet_comid <-
                                   unique(df_sub$outlet_comid)
                                 out$main_id = unique(df_sub$main_id)
                                 out
                               }, df = nhd_crosswalk)


  names(nhd_crosswalk_list) <- paste0("wb-", unique(nhd_crosswalk$ID))


  jsonlite::write_json(
    nhd_crosswalk_list,
    outfile,
    pretty = TRUE
  )
}
