% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/coordinate-mapping.R
\name{locate_single_point}
\alias{locate_single_point}
\title{Create anatomical location for a single gaze point}
\usage{
locate_single_point(
  x,
  y,
  z,
  nifti_data,
  plane = "AXIAL",
  timestamp_us = NULL,
  adjust_coordinates = TRUE
)
}
\arguments{
\item{x}{Gaze x coordinate (0-1 normalized to full display)}

\item{y}{Gaze y coordinate (0-1 normalized to full display)}

\item{z}{Normalized slice position (0-1)}

\item{nifti_data}{List from preload_nifti_data()}

\item{plane}{Imaging plane}

\item{timestamp_us}{Timestamp in microseconds}

\item{adjust_coordinates}{Whether to adjust for display frame (default TRUE)}
}
\value{
List with location information
}
\description{
This is equivalent to the JavaScript locate() function for a single point.
Takes gaze coordinates and returns anatomical location information.
}
