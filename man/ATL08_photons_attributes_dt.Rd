% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ATL08_photons_attributes_dt.R
\name{ATL08_photons_attributes_dt}
\alias{ATL08_photons_attributes_dt}
\title{ATL08 computed photons attributes}
\usage{
ATL08_photons_attributes_dt(
  atl08_h5,
  beam = c("gt1l", "gt1r", "gt2l", "gt2r", "gt3l", "gt3r"),
  photon_attributes = c("ph_segment_id", "classed_pc_indx", "classed_pc_flag", "ph_h",
    "d_flag", "delta_time")
)
}
\arguments{
\item{atl08_h5}{A ICESat-2 ATL08 object (output of \code{\link[=ATL08_read]{ATL08_read()}} function).
An S4 object of class \code{\linkS4class{icesat2.atl08_dt}}.}

\item{beam}{Character vector indicating beams to process (e.g. "gt1l", "gt1r", "gt2l", "gt2r", "gt3l", "gt3r")}

\item{photon_attributes}{\linkS4class{character} vector indicating the attributes to extract from the ATL08 photons.
Default all \code{c("ph_segment_id", "classed_pc_indx", "classed_pc_flag", "ph_h", "d_flag", "delta_time")}.}
}
\value{
Returns an S4 object of class \link[data.table:data.table]{data.table::data.table}
containing the ATL08 computed photons attributes.
}
\description{
This function extracts computed photons attributes from ICESat-2 ATL08 data
}
\details{
These are the photons attributes extracted by default:
\itemize{
\item \strong{ph_segment_id}: Georeferenced bin number (20-m) associated with each photon
\item \strong{classed_pc_indx}: Indices of photons tracking back to ATL03 that surface finding
software identified and used within the creation of the data products.
\item \strong{classed_pc_flag}: The L2B algorithm is run if this flag is set to 1 indicating data have sufficient waveform fidelity for L2B to run
\item \strong{ph_h}: Height of photon above interpolated ground surface
\item \strong{d_flag}: Flag indicating whether DRAGANN labeled the photon as noise or signal
\item \strong{delta_time}: Mid-segment GPS time in seconds past an epoch. The epoch is provided in the metadata at the file level
}
}
\examples{
# Specifying the path to ATL08 file
atl08_path <- system.file("extdata",
  "atl08_clip.h5",
  package = "ICESat2VegR"
)

# Reading ATL08 data (h5 file)
atl08_h5 <- ATL08_read(atl08_path)

# Extracting ATL08 classified photons and heights
atl08_photons <- ATL08_photons_attributes_dt(atl08_h5 = atl08_h5)

head(atl08_photons)
close(atl08_h5)
}
\seealso{
\url{https://icesat-2.gsfc.nasa.gov/sites/default/files/page_files/ICESat2_ATL08_ATBD_r006.pdf}
}
