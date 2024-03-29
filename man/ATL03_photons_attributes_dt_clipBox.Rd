% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ATL03_photons_attributes_dt_clipBox.R
\name{ATL03_photons_attributes_dt_clipBox}
\alias{ATL03_photons_attributes_dt_clipBox}
\title{Clip ATL03 photons by Coordinates}
\usage{
ATL03_photons_attributes_dt_clipBox(atl03_photons_dt, xmin, xmax, ymin, ymax)
}
\arguments{
\item{atl03_photons_dt}{A atl03_photons_dt object (output of \code{\link[=atl03_photons_attributes_dt]{atl03_photons_attributes_dt()}} function).
An S4 object of class \link{icesat2.atl03_dt}}

\item{xmin}{Numeric. West longitude (x) coordinate of bounding rectangle, in decimal degrees.}

\item{xmax}{Numeric. East longitude (x) coordinate of bounding rectangle, in decimal degrees.}

\item{ymin}{Numeric. South latitude (y) coordinate of bounding rectangle, in decimal degrees.}

\item{ymax}{Numeric. North latitude (y) coordinate of bounding rectangle, in decimal degrees.}
}
\value{
Returns an S4 object of class \link{icesat2.atl03_dt}
containing the ATL03 photons attributes.
}
\description{
This function clips ATL03 photons attributes within a given bounding coordinates
}
\examples{
# Specifying the path to ATL03 file (zip file)
outdir <- tempdir()
atl03_zip <- system.file("extdata",
  "atl03_20220401221822_01501506_005_01.zip",
  package = "ICESat2VegR"
)

# Unzipping ATL03 file
atl03_path <- unzip(atl03_zip, exdir = outdir)

# Reading ATL03 data (h5 file)

# Extracting ATL03 photons attributes
atl03_photons_dt <- ATL03_photons_attributes_dt(atl03_h5 = atl03_h5)

# Bounding rectangle coordinates
xmin <- -107.7
xmax <- -106.5
ymin <- 32.75
ymax <- 42.75

# Clipping ATL03 photons  by boundary box extent
atl03_photons_dt_clip <- ATL03_photons_attributes_dt_clipBox(atl03_photons_dt, xmin, xmax, ymin, ymax)

close(atl03_h5)
}
\seealso{
\url{https://icesat-2.gsfc.nasa.gov/sites/default/files/page_files/ICESat2_atl03_ATBD_r006.pdf}
}
