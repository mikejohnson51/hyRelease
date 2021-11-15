#' Build Release Directory Structure
#' @description in a root directory, create a graph, spatial, parameters,
#' and crosswalks sub directories
#' @param root the root path
#' @return list of directory paths
#' @export
#'

build_release_directory = function(root = "" ){
  dir.create(root, showWarnings = FALSE, recursive = TRUE)
  release_graph = file.path(root, "graph")
  dir.create(release_graph, showWarnings = FALSE, recursive = TRUE)
  release_geo   = file.path(root, "spatial")
  dir.create(release_geo, showWarnings = FALSE, recursive = TRUE)
  release_param = file.path(root, "parameters")
  dir.create(release_param, showWarnings = FALSE, recursive = TRUE)
  release_cw    = file.path(root, "crosswalks")
  dir.create(release_cw, showWarnings = FALSE, recursive = TRUE)
  return(list(release_graph = release_graph,
              release_geo = release_geo,
              release_param = release_param,
              release_cw = release_cw
  ))
}


#' Write GeoJSON
#' @param x sf object
#' @param file file path
#' @export
#' @importFrom sf write_sf st_make_valid st_transform
#'
write_geojson <- function(x, file) {
  names(x) <- tolower(names(x))
  unlink(file)
  write_sf(st_make_valid(st_transform(x, 4326)), file,
           layer_options = c("ID_FIELD=id", "ID_TYPE=String"))
}
