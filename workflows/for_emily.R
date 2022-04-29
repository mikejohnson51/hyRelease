library(sf)
library(dplyr)


wg_files = list.files('/Users/mjohnson/Downloads/WalnutGulch_subbasins', pattern = ".shp$", full.names = TRUE)

bind_rows(lapply(wg_files, read_sf)) %>%
  st_transform(5070) %>%
  mutate(ID = gsub(".shp", "", basename(wg_files))) %>%
  write_sf("/Users/mjohnson/Downloads/WalnutGulch.gpkg", "catchments")


aggregate_basin_attributes(gpkg = '/Users/mjohnson/Downloads/WalnutGulch.gpkg' ,
                      catchment_name = "catchments",
                      geo_dir = file.path(geogrids::geo_path(), "climate"),
                      years = 5,
                      precision = 9,
                      out_file = '/Users/mjohnson/Downloads/wg_lstm2.csv')

rc_files = data.frame(files = list.files('/Users/mjohnson/Downloads/ReynoldsCreek_subbasins',
                                         pattern = ".shp$",
                                         full.names = TRUE)) %>%
  dplyr::filter(!grepl("watershed_036", files))

rc = bind_rows(lapply(rc_files$files, read_sf)) %>%
  st_transform(5070) %>%
  mutate(ID = gsub(".shp", "", basename(rc_files$files)), Id = NULL) %>%
  sf::write_sf("/Users/mjohnson/Downloads/ReynoldsCreek.gpkg", "catchments")

aggregate_lstm_params(gpkg = '/Users/mjohnson/Downloads/ReynoldsCreek.gpkg' ,
                      catchment_name = "catchments",
                      geo_dir = file.path(geogrids::geo_path(), "climate"),
                      years = 5,
                      precision = 9,
                      out_file = '/Users/mjohnson/Downloads/rc_lstm3.csv')
