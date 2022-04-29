### CAMELS Processing
# Read from Jonathans repo
# camels  = read.csv('https://raw.githubusercontent.com/jmframe/lstm/master/data/camels_basin_list_516.txt', header = FALSE) %>%
#   mutate(gageID = sprintf("%08d", V1)) %>%
#   mutate(V1 = NULL)
#
# nc = RNetCDF::open.nc(file.path(nwm_dir, "RouteLink_CONUS.nc"))
# df2 = suppressMessages({
#   lapply(c("link", "gages"), function(x) x = RNetCDF::var.get.nc(nc, x)) %>%
#     bind_cols()
# })
#
# names(df2) = c("comid", "gageID")
#
# camels_refine = df2 %>%
#   mutate(gageID = trimws(gageID),
#          gageID = ifelse(gageID == "", NA, gageID)
#   ) %>%
#   right_join(camels)
#
# network = arrow::read_parquet("/Volumes/Transcend/ngen/enhd_nhdplusatts.parquet")
#
# for(i in 1:nrow(camels_refine)){
#   if(is.na(camels_refine$comid[i])){
#     camels_refine$comid[i] =  dataRetrieval::findNLDI(nwis = camels_refine$gageID[i])$comid
#   }
#
#   master = left_join(data.frame(comid = nhdplusTools::get_UT(network, camels_refine$comid[i])),
#                      nhdplusTools::get_vaa("rpuid"),
#                      by = "comid")
#   RPUS = unique(master$rpuid)
#   camels_refine$rpus[i] = paste(RPUS[!is.na(RPUS)], collapse = ",")
# }
#

#data.table::fwrite(camels_refine, "camels.csv")
{
library(hyAggregate)
library(dplyr)
library(data.table)
library(arrow)
library(sf)
devtools::load_all(".")

base_path            = '/Volumes/Transcend/ngen/CAMELS20/'
##############################################################################
network_parquet      = "/Volumes/Transcend/ngen/enhd_nhdplusatts.parquet"
reference_fabric_dir = '/Volumes/Transcend/ngen/CONUS-hydrofabric/ngen-reference'
facfdr               = '/Volumes/Transcend/ngen/fdrfac_cog'
nwm_dir              = '/Volumes/Transcend/nwmCONUS-v216'
gages_iii            = '/Volumes/Transcend/ngen/gages_iii.gpkg'
wb_gpkg              = '/Users/mjohnson/Downloads/nhdplus_waterbodies.gpkg'

nexus_prefix          = "nex-"
catchment_prefix      = "cat-"
waterbody_prefix      = "wb-"
terminal_nexus_prefix = "tnx-"
ternimal_wb_prefix    = "twb-"
##############################################################################

camels    = mutate(fread("camels.csv"), gageID = sprintf("%08d", gageID))

network   = read_parquet(network_parquet,
                         col_select = c('comid', 'pathlength',
                                        'lengthkm', 'hydroseq',
                                         'levelpathi','dnhydroseq'))
}


#163 is a disconnected network problem!!

for(i in 164:nrow(camels)){

  COMID = camels$comid[i]
  name  = camels$name[i]

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

  write_hydrofabric = function(gpkg,
                               catchment_name, flowpath_name,
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


  hyfab = write_hydrofabric(ref,
                            catchment_name = "aggregated_catchments",
                            flowpath_name = "aggregated_flowpaths",
                            spatial_path, parameter_path )


  write_waterbody_json(hyfab,
                       flowpath_name = 'flowpaths',
                       outfile = file.path(parameter_path, "flowpath_parameters.json"))

  write_nwis_crosswalk2(flowpath_data,
                        gages_iii = gages_iii,
                        outfile = file.path(parameter_path, "cross-walk.json"))

  build_lake_params(gpkg = hyfab,
                    catchment_name = "catchments",
                    flowline_name  = 'flowpaths',
                    waterbody_gpkg = wb_gpkg,
                    nwm_dir = nwm_dir,
                    out_file = NULL)

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

  cat(crayon::blue("\nFinished:", i))
}



