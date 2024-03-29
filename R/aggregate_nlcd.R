#' @title Aggregated NLCD variables
#' @param gpkg a geopackage with aggregation units
#' @param catchment_name the layer name of the aggregation units
#' @param imperv_path path to NLCD impervious IMG file
#' @param lc_path path to NLCD landcover IMG file
#' @param precision the precision of the computations
#' @param out_file the file to write to (CSV using data.table fwrite)
#' @return NULL
#' @export
#' @importFrom sf read_sf
#' @importFrom zonal weighting_grid execute_zonal
#' @importFrom dplyr full_join mutate
#' @importFrom data.table fwrite setnames
#' @importFrom tidyr pivot_wider
#' @importFrom data.table fwrite

aggregate_nlcd = function(gpkg,
                         catchment_name,
                         imperv_path = TRUE,
                         lc_path     = TRUE,
                         precision = 9,
                         out_file = NULL){

  .SD <- .data <- NULL

  cats = sf::read_sf(gpkg, catchment_name)
  lc_w = zonal::weighting_grid(imperv_path, cats, "ID" )

  xx = zonal::execute_zonal(file = imperv_path, w = lc_w, join = FALSE) %>%
    setnames(c("ID", "impervious_percentage"))

  zz =  zonal::execute_zonal(file = lc_path, w = lc_w, FUN = "freq", join = FALSE) %>%
    mutate(value = paste0("percent_nlcd_lc_", .data$value),
            percentage = 100 * .data$percentage) %>%
    tidyr::pivot_wider(id_cols = "ID",
                       names_from = "value",
                       values_from = "percentage")

  dt = dplyr::full_join(xx, zz, by = "ID")

  cols = names(dt)[names(dt) != "ID"]
  dt[,(cols) := round(.SD, precision), .SDcols=cols]

  data.table::fwrite(dt, out_file, row.names = FALSE)

}
