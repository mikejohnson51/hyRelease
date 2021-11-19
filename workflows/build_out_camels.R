### CAMELS Processing

library(hyAggregate)
devtools::load_all(".")

nwm_dir =  "/Volumes/Transcend/nwmCONUS-v216"
nexus_prefix        = "nex-"
catchment_prefix    = "cat-"
waterbody_prefix    = "wb-"
terminal_nexus_prefix = "tnx-"
ternimal_wb_prefix  = "twb-"


# Read from Jonathans repo
gageIDs  = read.csv('https://raw.githubusercontent.com/jmframe/lstm/master/data/camels_basin_list_516.txt', header = FALSE) %>%
  mutate(gageID = sprintf("%08d", V1)) %>%
  pull(gageID)

out_dir = '/Volumes/Transcend/ngen/CAMELS'
#####################################################

for (i in 410:516) {

  here = file.path(out_dir, paste0("gage_", gageIDs[i]))

  if (length(list.files(here, recursive = TRUE)) != 11) {
    rel = build_release_directory(here)

    refactored_gpkg = glue::glue(
      '/Volumes/Transcend/ngen/refactor-tests/CAMELS/base-runs/gage_{ID}.gpkg',
      ID = gageIDs[i]
    )

    agg = aggregate_by_thresholds(gpkg     = refactored_gpkg,
                                  fl_name  = 'refactored_flowpaths',
                                  cat_name = 'refactored_catchments')


    agg$flowpaths = agg$flowpaths %>%
      mutate(length_km = add_lengthkm(.)) %>%
      length_average_routlink(
        rl_vars = c("link", "Qi", "MusK", "MusX", "n", "So", "ChSlp", "BtmWdth",
                    "time", "Kchan", "nCC", "TopWdthCC", "TopWdth"),
        rl_path  = file.path(nwm_dir, "RouteLink_CONUS.nc"))

    ########## GRAPH DATA
    catchment_edge_list <- get_catchment_edges_terms(agg$flowpaths)
    fp_edge_list        <- get_catchment_edges_terms(agg$flowpaths, catchment_prefix = waterbody_prefix)
    waterbody_edge_list <- get_waterbody_edges_terms(agg$flowpaths)

    jsonlite::write_json(
      catchment_edge_list,
      file.path(rel$release_graph,  "catchment_edge_list.json"),
      pretty = TRUE
    )

    jsonlite::write_json(
      waterbody_edge_list,
      file.path(rel$release_graph, "waterbody_edge_list.json"),
      pretty = TRUE
    )

    jsonlite::write_json(
      fp_edge_list,
      file.path(rel$release_graph, "flowpath_edge_list.json"),
      pretty = TRUE
    )

    ########## SPATIAL DATA

    nexus_data = hyAggregate::get_nexus_locations(agg$flowpaths) %>%
      mutate(prefix = ifelse(ID > 100000000, terminal_nexus_prefix, nexus_prefix)) %>%
      mutate(ID = paste0(prefix, ID), prefix = NULL) %>%
      left_join(catchment_edge_list, by = c("ID")) %>%
      mutate(toID  = ifelse(toID == "cat-0", NA, toID))

    catchment_data  = agg$catchments %>%
      rename(area_sqkm = areasqkm) %>%
      get_catchment_data(catchment_edge_list, catchment_prefix = catchment_prefix)

    flowpath_data <- mutate(agg$flowpaths, slope = So) %>%
      get_flowpath_data(catchment_edge_list, catchment_prefix = catchment_prefix) %>%
      mutate(realized_catchment = gsub(waterbody_prefix, catchment_prefix, ID))

    sf::write_sf(
      catchment_data,
      file.path(rel$release_geo, "hydrofabric.gpkg"),
      'catchments',
      overwrite = TRUE
    )
    sf::write_sf(
      flowpath_data,
      file.path(rel$release_geo, "hydrofabric.gpkg"),
      'flowpaths',
      overwrite = TRUE
    )
    sf::write_sf(
      nexus_data,
      file.path(rel$release_geo, "hydrofabric.gpkg"),
      'nexus',
      overwrite = TRUE
    )

    write_geojson(catchment_data,
                  file.path(rel$release_geo, "catchment_data.geojson"))
    write_geojson(flowpath_data,
                  file.path(rel$release_geo, "flowpath_data.geojson"))
    write_geojson(nexus_data,
                  file.path(rel$release_geo,  "nexus_data.geojson"))

    ########## CROSS WALK

    nwis_sites = agg$flowpaths %>%
      st_drop_geometry() %>%
      filter(!is.na(gages)) %>%
      select(ID, gages)

    nhd_crosswalk <- agg$flowpaths %>%
      st_drop_geometry() %>%
      select(ID, member_COMID, LevelPathID) %>%
      mutate(comid = strsplit(member_COMID, ",")) %>%
      tidyr::unnest_longer(col = c("comid")) %>%
      mutate(ID = ID,
             comid = as.numeric(comid)) %>%
      select(ID, comid, main = LevelPathID) %>%
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
                                   if (any(!is.na(df_sub$gages))) {
                                     out$gages <- unique(df_sub$gages[!is.na(df_sub$gages)])
                                   }
                                   out$outlet_comid <-
                                     unique(df_sub$outlet_comid)
                                   out$main = unique(df_sub$main)
                                   out
                                 }, df = nhd_crosswalk)


    names(nhd_crosswalk_list) <-
      paste0(catchment_prefix, unique(nhd_crosswalk$ID))


    jsonlite::write_json(
      nhd_crosswalk_list,
      file.path(rel$release_cw, "crosswalk-mapping.json"),
      pretty = TRUE,
      auto_unbox = FALSE
    )

    ### PARAMTERS

    aggregate_nwm_params_reduce(
      gpkg = file.path(rel$release_geo, "hydrofabric.gpkg"),
      catchment_name = "catchments",
      flowline_name = "flowpaths",
      nwm_dir = nwm_dir,
      single_layer = TRUE,
      out_file  = file.path(rel$release_param, 'nwm.csv')
    )

    aggregate_lstm_params(
      gpkg = file.path(rel$release_geo, "hydrofabric.gpkg"),
      catchment_name = "catchments",
      geo_dir = file.path(geogrids::geo_path(), "climate"),
      out_file =  file.path(rel$release_param, 'camels.csv')
    )


    wb_feilds = agg$flowpaths %>%
      select(-ID) %>%
      st_drop_geometry() %>%
      mutate(Length_m = as.numeric(Length_m),
             member_COMID = strsplit(member_COMID, ","))

    wb_feilds <- split(wb_feilds, seq(nrow(wb_feilds)))

    names(wb_feilds) = paste0(waterbody_prefix, agg$flowpaths$ID)

    jsonlite::write_json(wb_feilds,
                         file.path(rel$release_param, "waterbody-params.json"),
                         pretty = TRUE)
  }

  cat(crayon::blue("\nFinished:", i))
}



#
# files = file.path(list.dirs(out_dir, recursive = F), "spatial/hydrofabric.gpkg")
# out_csv = file.path(list.dirs(out_dir, recursive = F), "parameters/nwm.csv")
#
# for(i in 298:length(files)){
#
#   tryCatch({
#   aggregate_nwm_params_reduce(
#     gpkg = files[i],
#     catchment_name = "catchments",
#     flowline_name  = "flowpaths",
#     nwm_dir = nwm_dir,
#     single_layer = TRUE,
#     out_file  = out_csv[i]
#   )
#   }, error = function(e){ NULL})
#   message(i)
# }
#
# nlcd_csv = file.path(list.dirs(out_dir, recursive = F), "parameters/nlcd.csv")
#
# for(i in 1:length(files)){
#
#   if(!file.exists(nlcd_csv[i])){
#     tryCatch({
#       aggregate_nlcd(
#         gpkg = files[i],
#         catchment_name = "catchments",
#         imperv_path    = '/Volumes/Transcend/ngen/nlcd_2019_impervious_l48_20210604/nlcd_2019_impervious_l48_20210604.img',
#         lc_path        = '/Volumes/Transcend/ngen/nlcd_2019_land_cover_l48_20210604/nlcd_2019_land_cover_l48_20210604.img',
#         out_file  = nlcd_csv[i]
#       )
#     }, error = function(e){ NULL})
#   }
#   message(i)
# }
