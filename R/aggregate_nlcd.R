#' Aggregated NLCD variables
#' @param gpkg a geopackage with aggregation units
#' @param catchment_name the layer name of the aggregation units
#' @param imperv_path path to NLCD impervious IMG file
#' @param lc_path path to NLCD landcover IMG file
#' @param precision the precision of the computations
#' @param out_file the file to write to (CSV using data.table fwrite)
#' @return NULL
#' @export
#' @importFrom sf read_sf
#' @importFrom zonal weighting_grid execute_zonal execute_zonal_cat
#' @importFrom dplyr full_join mutate
#' @importFrom data.table fwrite setnames
#' @importFrom tidyr pivot_wider
#' @importFrom data.table fwrite

# imperv_path    = '/Volumes/Transcend/ngen/nlcd_2019_impervious_l48_20210604/nlcd_2019_impervious_l48_20210604.img'
# lc_path        = '/Volumes/Transcend/ngen/nlcd_2019_land_cover_l48_20210604/nlcd_2019_land_cover_l48_20210604.img'
# catchment_name = "refactored_catchments"
# gpkg           = '/Volumes/Transcend/ngen/refactor-tests/CAMELS/base-runs/gage_01022500.gpkg'


# system.time({
#   aggregate_nlcd(gpkg,
#                catchment_name = "refactored_catchments",
#                imperv_path    = '/Volumes/Transcend/ngen/nlcd_2019_impervious_l48_20210604/nlcd_2019_impervious_l48_20210604.img',
#                lc_path        = '/Volumes/Transcend/ngen/nlcd_2019_land_cover_l48_20210604/nlcd_2019_land_cover_l48_20210604.img',
#                out_file = '/Users/mjohnson/Downloads/nlcd_test.csv'
#                )
#
# })

aggregate_nlcd = function(gpkg,
                         catchment_name,
                         imperv_path = TRUE,
                         lc_path     = TRUE,
                         precision = 9,
                         out_file = NULL){

  cats = sf::read_sf(gpkg, catchment_name)
  lc_w = zonal::weighting_grid(imperv_path, cats, "ID" )

  xx = zonal::execute_zonal(file = imperv_path, w = lc_w) %>%
    setnames(c("ID", "impervious_percentage"))

  zz =  zonal::execute_zonal_cat(file = lc_path, w = lc_w) %>%
    mutate(value = paste0("percent_nlcd_lc_", value),
            percentage = 100 * percentage) %>%
    tidyr::pivot_wider(id_cols = "ID", names_from = "value", values_from = "percentage")

  dt = dplyr::full_join(xx, zz, by = "ID")

  cols = names(dt)[names(dt) != "ID"]
  dt[,(cols) := round(.SD, precision), .SDcols=cols]

  data.table::fwrite(dt, out_file, row.names = FALSE)

}
