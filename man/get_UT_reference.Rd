% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/update_functions.R
\name{get_UT_reference}
\alias{get_UT_reference}
\title{Extract upstream network from reference fabric}
\usage{
get_UT_reference(
  network,
  reference_flowlines,
  reference_catchments = NULL,
  comid,
  outfile = NULL
)
}
\arguments{
\item{network}{enhd_nhdplusatts data.frame}

\item{comid}{Origin of network trace (COMID)}

\item{reference_fabric_dir}{path to refernece network directory}
}
\value{
list containing catchment and flowpath networks
}
\description{
Extract upstream network from reference fabric
}
