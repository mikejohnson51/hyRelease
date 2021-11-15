#' Assign Length Weighted Flowpath Attributes from Routelink file
#' @param flowpaths flowpaths with a`member_COMID` attribute
#' @param rl_vars RouteLink variable to append
#' @param rl_path Path to RouteLink
#' @return flowpaths with length weighted RouteLink Attributes
#' @export
#' @importFrom sf st_drop_geometry st_length
#' @importFrom dplyr select mutate rename right_join group_by summarise across everything left_join
#' @importFrom tidyr unnest
#' @importFrom RNetCDF open.nc var.get.nc
#' @importFrom stats weighted.mean

length_average_routlink = function(flowpaths,
         rl_vars,
         rl_path){

  flowpaths = flowpaths %>%
    hyRefactor::add_lengthmap(length_table = nhdplusTools::get_vaa("lengthkm"))

  if(!"Length" %in% rl_vars){ rl_vars = c("Length", rl_vars) }

  net_map  <- dplyr::select(st_drop_geometry(flowpaths), ID, lengthMap) %>%
    mutate(comid = strsplit(.data$lengthMap, ",")) %>%
    tidyr::unnest(cols = comid) %>%
    mutate(full_comids = floor(as.numeric(comid)),
           w = 10 * (as.numeric(comid) - full_comids),
           #w = ifelse(rep(length_weight, dplyr::n()), w, 1),
           comid = NULL)

  nc = RNetCDF::open.nc(rl_path)

  df = data.frame(do.call(cbind, lapply(rl_vars, function(x) RNetCDF::var.get.nc(nc, x))))

  names(df) = rl_vars

  df = df %>%
    rename(comid = link) %>%
    right_join(net_map, by = c('comid' = 'full_comids')) %>%
    select(-.data$lengthMap) %>%
    mutate(w = w * Length) %>%
    group_by(ID) %>%
    summarise(across(everything(), ~ round(
      weighted.mean(x = .,
                    w = w,
                    na.rm = TRUE), 3))) %>%
    dplyr::select(-comid, -Length, -w)

  df2 = suppressMessages({
    lapply(c("link", "gages", 'NHDWaterbodyComID'), function(x) x = RNetCDF::var.get.nc(nc, x)) %>%
    bind_cols()
  })

  names(df2) = c("link", "gages", 'NHDWaterbodyComID')

  df2 = df2 %>%
    rename(comid = link) %>%
    right_join(net_map, by = c('comid' = 'full_comids')) %>%
    mutate(gages = trimws(gages),
           gages = ifelse(gages == "", NA, gages),
           NHDWaterbodyComID = ifelse(NHDWaterbodyComID == -9999, NA, NHDWaterbodyComID)
    ) %>%
    group_by(ID) %>%
    summarise(gages = paste(gages[!is.na(gages)], collapse = ","),
              NHDWaterbodyComID = paste(unique(NHDWaterbodyComID[!is.na(NHDWaterbodyComID)]), collapse = ",")) %>%
    left_join(df, by = "ID") %>%
    mutate(gages = ifelse(gages == "", NA, gages),
           NHDWaterbodyComID = ifelse(NHDWaterbodyComID == "", NA, NHDWaterbodyComID))

  left_join(flowpaths, df2, by = "ID") %>%
    mutate(Length_m = sf::st_length(.))
}
