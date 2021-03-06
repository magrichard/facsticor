% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/visualize_data.R
\name{plot_tc}
\alias{plot_tc}
\title{Plot the time course distributions accoording to color codes.}
\usage{
plot_tc(data, idx_data, colors = c("white", "darkblue"), ...)
}
\arguments{
\item{data}{The gated data generated by gate_facs_yeast_HTS_data}

\item{idx_data}{The index of the data that will be ploted}

\item{colors}{The colors to use to construct the ramp, default is c("white", "darkblue")}
}
\value{
A plot of distribution according to time
}
\description{
This function works on FL1H, FSC and SSC data acquired on 96-well plates
}
