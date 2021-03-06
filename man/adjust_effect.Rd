% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/normalize_data.R
\name{adjust_effect}
\alias{adjust_effect}
\title{Normalize the data according to the annova model}
\usage{
adjust_effect(data, offsets, offset_key)
}
\arguments{
\item{data}{The data to be corrected}

\item{offsets}{The offsets to be used}
}
\value{
A normalized data list
}
\description{
This function works on FL1H, FSC and SSC data acquired on 96-well plates
}
\examples{
à finaliser
load("data/data.rda")
final_data = filter_data(data, data$anno$rep != 1)
offsets_plate = model.tables(m_annov)$tables$plate
final_data = adjust_effect(final_data, offsets_plate, "plate")
}
