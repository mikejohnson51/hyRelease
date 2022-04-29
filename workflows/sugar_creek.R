COMID = 9731454
name  = "sugar_creek"
basename = "/Volumes/Transcend/ngen/"

######

dir = file.path.build(base_path, name)
spatial_path   = file.path.build(dir, "spatial")
parameter_path = file.path.build(dir, "parameters")
ref = file.path(dir,"reference.gpkg")

if(!gpkg_layers(path = ref, 2, "reference")){
  get_UT_reference(network,
                   reference_fabric_dir = reference_fabric_dir,
                   comid = COMID,
                   outfile = ref)
}

if(!gpkg_layers(path = ref, 2, "refactored")){
  reference_fabric = list.files(reference_fabric_dir, full.names = TRUE, pattern = paste0(camels$rpus[i], ".gpkg$"))

  refactor_wrapper(
    flowpaths  = sf::read_sf(ref, "reference_flowpaths"),
    catchments = sf::read_sf(ref, "reference_catchments"),
    events     = sf::read_sf(reference_fabric,  "events"),
    avoid      = sf::read_sf(reference_fabric,  "avoid")$COMID,
    facfdr     = facfdr,
    outfile    = ref
  )
}

if(!gpkg_layers(path = ref, 2, "aggregate")){
  o =  hyAggregate::aggregate_by_thresholds(gpkg = ref,
                                            fl_name  = 'refactored_flowpaths',
                                            cat_name = 'refactored_catchments',
                                            write    = TRUE)

  rm(o)
}

fps = sf::read_sf(ref, "aggregated_flowpaths")

########## GRAPH DATA
catchment_edge_list <- get_catchment_edges_terms(fps)
fp_edge_list        <- get_catchment_edges_terms(fps, catchment_prefix = waterbody_prefix)

########## SPATIAL DATA

hyfab = file.path(spatial_path, "hydrofabric.gpkg")

if(!file.exists(hyfab)){
  nexus_data = hyAggregate::get_nexus_locations(fps) %>%
    mutate(prefix = ifelse(ID > 100000000, terminal_nexus_prefix, nexus_prefix)) %>%
    mutate(ID = paste0(prefix, ID), prefix = NULL) %>%
    left_join(catchment_edge_list, by = c("ID")) %>%
    mutate(toID  = ifelse(toID == "cat-0", NA, toID))

  catchment_data  = sf::read_sf(ref, "aggregated_catchments") %>%
    rename(area_sqkm = areasqkm) %>%
    get_catchment_data(catchment_edge_list, catchment_prefix = catchment_prefix)

  flowpath_data =  fps %>%
    add_slope() %>%
    get_flowpath_data(catchment_edge_list) %>%
    length_average_routlink(rl_path  = file.path(nwm_dir, "RouteLink_CONUS.nc"))

  hyfab = write_nextgen_spatial(flowpath_data, catchment_data, nexus_data, spatial_path)

  write_json(fp_edge_list, file.path(parameter_path, "flowpath_edge_list.json"), pretty = TRUE)
  write_waterbody_json(flowpath_data, outfile = file.path(parameter_path, "flowpath_parameters.json"))
  write_nwis_crosswalk2(flowpath_data, gages_iii = gages_iii, outfile = file.path(parameter_path, "cross-walk.json"))
}

build_lake_params(gpkg = hyfab,
                  catchment_name = "catchments",
                  flowline_name  = 'flowpaths',
                  waterbody_gpkg = wb_gpkg,
                  nwm_dir = nwm_dir,
                  out_file = '/Users/mjohnson/Downloads/lakes_test.json')

  ### PARAMTERS
cfe    = file.path(parameter_path, 'cfe.csv')
b_atts = file.path(parameter_path, 'basin_attributes.csv')
noaowp = file.path(parameter_path, 'noahowp.csv')
atts   = file.path(parameter_path, 'attributes.parquet')


if(!file.exists(cfe)){
  aggregate_cfe(
    gpkg = hyfab,
    catchment_name = "catchments",
    flowline_name = "flowpaths",
    nwm_dir = nwm_dir,
    single_layer = TRUE,
    out_file  = cfe
  )
}

if(!file.exists(b_atts)){
  aggregate_basin_attributes(
    gpkg = hyfab,
    catchment_name = "catchments",
    geo_dir = file.path(geogrids::geo_path(), "climate"),
    out_file =  b_atts
  )
}

if(!file.exists(noaowp)){
  aggregate_noahowp(
    gpkg = hyfab,
    nwm_dir = nwm_dir,
    catchment_name = "catchments",
    out_file =  noaowp
  )
}

if(!file.exists(atts)){

  dfs = lapply(list(cfe, b_atts, noaowp), data.table::fread)

  xxx = Reduce(function(dtf1, dtf2) merge(dtf1, dtf2, by = "ID", all.x = TRUE), dfs)

  arrow::write_parquet(xxx, atts)
}
