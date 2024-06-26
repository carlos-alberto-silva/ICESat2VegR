% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ATL03_ATL08_photons_attributes_dt_join.R
\name{ATL03_ATL08_photons_attributes_dt_join}
\alias{ATL03_ATL08_photons_attributes_dt_join}
\title{Join ATL03 and ATL08 photons attributes}
\usage{
ATL03_ATL08_photons_attributes_dt_join(
  atl03_h5,
  atl08_h5,
  beam = c("gt1l", "gt1r", "gt2l", "gt2r", "gt3l", "gt3r")
)
}
\arguments{
\item{atl03_h5}{A ICESat-2 ATL03 object (output of \code{\link[=ATL03_read]{ATL03_read()}} function).
An S4 object of class \code{\linkS4class{icesat2.atl03_dt}}.}

\item{atl08_h5}{A ICESat-2 ATL08 object (output of \code{\link[=ATL08_read]{ATL08_read()}} function).
An S4 object of class \code{\linkS4class{icesat2.atl08_dt}}.}

\item{beam}{Character vector indicating beams to process
(e.g. "gt1l", "gt1r", "gt2l", "gt2r", "gt3l", "gt3r")}
}
\value{
Returns an S4 object of class \code{\linkS4class{icesat2.atl03atl08_dt}}
containing the ATL08 computed photons attributes.
}
\description{
This function joins ATL03 and ATL08 computed photons attributes
}
\details{
These are the photons attributes extracted by default:
\itemize{
\item \code{ph_segment_id}: Georeferenced segment id (20-m) associated with each photon.
\item \code{lon_ph}: Longitude of each received photon. Computed from the ECEF Cartesian coordinates
of the bounce point.
\item \code{lat_ph}: Latitude of each received photon. Computed from the ECEF Cartesian coordinates
of the bounce point.
\item \code{h_ph}: Height of each received photon, relative to the WGS-84 ellipsoid including the
geophysical corrections noted in section 6.0. Please note that neither the geoid,
ocean tide nor the dynamic atmospheric corrections (DAC) are applied to the
ellipsoidal heights.
\item \code{quality_ph}: Indicates the quality of the associated photon. 0=nominal,
1=possible_afterpulse, 2=possible_impulse_response_effect, 3=possible_tep.
Use this flag in conjunction with \code{signal_conf_ph} to identify those photons
that are likely noise or likely signal.
\item \code{solar_elevation}: Elevation of the sun above the horizon at the photon bounce point.
\item \code{dist_ph_along}: Along-track distance of the photon from the beginning of the segment.
\item \code{dist_ph_across}: Across-track distance of the photon from the center of the segment.
\item \code{night_flag}: Flag indicating the data were acquired in night conditions: 0=day, 1=night.
Night flag is set when solar elevation is below 0.0 degrees.
\item \code{classed_pc_indx}: Indices of photons tracking back to ATL03 that surface finding software
identified and used within the creation of the data products.
\item \code{classed_pc_flag}: The L2B algorithm is run if this flag is set to 1 indicating data have
sufficient waveform fidelity for L2B to run.
\item \code{ph_h}: Height of photon above interpolated ground surface.
\item \code{d_flag}: Flag indicating whether DRAGANN labeled the photon as noise or signal.
\item \code{delta_time}: Mid-segment GPS time in seconds past an epoch. The epoch is provided in
the metadata at the file level.
\item \code{orbit_number}: Orbit number identifier to identify data from different orbits.
\item \code{beam}: Beam identifier.
\item \code{strong_beam}: Logical indicating if the beam is a strong beam.
}
}
\examples{
# Specifying the path to ATL03 file
atl03_path <- system.file("extdata",
  "atl03_clip.h5",
  package = "ICESat2VegR"
)

# Specifying the path to ATL08 file
atl08_path <- system.file("extdata",
  "atl08_clip.h5",
  package = "ICESat2VegR"
)

# Reading ATL03 data (h5 file)
atl03_h5 <- ATL03_read(atl03_path = atl03_path)

# Reading ATL08 data (h5 file)
atl08_h5 <- ATL08_read(atl08_path = atl08_path)

# # Extracting ATL03 and ATL08 photons and heights
atl03_atl08_dt <- ATL03_ATL08_photons_attributes_dt_join(atl03_h5, atl08_h5)
head(atl03_atl08_dt)

close(atl03_h5)
close(atl08_h5)
}
\seealso{
\url{https://icesat-2.gsfc.nasa.gov/sites/default/files/page_files/ICESat2_ATL08_ATBD_r006.pdf}
}
