## Read Refactor

## Aggregate

# gpkg = '/Users/mjohnson/Downloads/2021-10-28/spatial/hydrofabric.gpkg'
#
# agg = list(catchments = read_sf(gpkg, "catchments"),
#            flowpaths  = read_sf(gpkg, "flowpaths"))
#
# agg$catchments = flush_prefix(input = agg$catchments, col = c('ID', 'toID'))
# agg$flowpaths  = flush_prefix(input = agg$flowpaths, col = c('ID', 'toID'))
#
# write_sf(agg$catchments, "/Users/mjohnson/Desktop/test_release/agg_fabric.gpkg", "catchments")
# write_sf(agg$flowpaths, "/Users/mjohnson/Desktop/test_release/agg_fabric.gpkg", "flowpaths")

gpkg = "/Users/mjohnson/Desktop/test_release/agg_fabric.gpkg"

agg = hyAggregate::build_network_list(gpkg, "flowpaths", "catchments")

rel = build_release_directory('/Users/mjohnson/Desktop/test_release')

## Cross Area Tabulation

agg$catchments = rename(agg$catchments, area_sqkm = areasqkm)

agg$catchments = catchment_waterbody_interaction(catchments = agg$catchments,
                                                 waterbody_gpkg = '/Users/mjohnson/Downloads/nhdplus_waterbodies.gpkg',
                                                 nwm_dir = nwm_dir )

## Write Edgeslist

edgelists = build_edge_lists(agg$flowpaths, outdir = rel$release_graph)

### Write Spatial Data

catchment_data  =  get_catchment_data(agg$catchments,
                                      edgelists$catchment_edge_list,
                                      catchment_prefix = "cat-")

write_geojson(catchment_data, file.path(rel$release_geo, "catchment_data.geojson"))

flowpath_data <- agg$flowpaths %>%
  flush_prefix(c('toID', "ID")) %>%
  add_slope() %>%
  get_flowpath_data(edgelists$catchment_edge_list, catchment_prefix = wb_prefix) %>%
  mutate(realized_catchment = gsub(wb_prefix, catchment_prefix, ID))

write_geojson(flowpath_data,  file.path(rel$release_geo, "flowpath_data.geojson"))

build_lake_params(flowpaths = flowpath_data,
                  nwm_dir = nwm_dir,
                  outfile = file.path(rel$release_param, "lakes.json"))

write_waterbody_json(flowpaths = flowpath_data,
                     nwm_dir   = nwm_dir,
                     outfile   = file.path(rel$release_param, "waterbody.json"))

write_nwis_crosswalk(flowpaths = flowpath_data,
                     nwm_dir = nwm_dir,
                     outfile = file.path(rel$release_cw, "nwis.json"))

####

gpkg = file.path(rel$release_geo, "hydrofabric.gpkg")
write_sf(catchment_data, gpkg, "catchments")
write_sf(flowpath_data, gpkg, "flowpaths")

aggregate_nwm_params_reduce(
  gpkg = gpkg,
  catchment_name = "catchments",
  flowline_name = "flowpaths",
  nwm_dir =  "/Volumes/Transcend/nwmCONUS-v216",
  single_layer = TRUE,
  out_file  = file.path(rel$release_param, 'nwm.csv')
)
aggregate_lstm_params(
  gpkg = gpkg,
  catchment_name = "catchments",
  geo_dir = file.path(geogrids::geo_path(), "climate"),
  out_file =  file.path(rel$release_param, 'camels.csv')
)
aggregate_nlcd(
  gpkg = gpkg,
  catchment_name = "catchments",
  imperv_path    = '/Volumes/Transcend/ngen/nlcd_2019_impervious_l48_20210604/nlcd_2019_impervious_l48_20210604.img',
  lc_path        = '/Volumes/Transcend/ngen/nlcd_2019_land_cover_l48_20210604/nlcd_2019_land_cover_l48_20210604.img',
  out_file  = file.path(rel$release_param, 'camels.csv')
)



