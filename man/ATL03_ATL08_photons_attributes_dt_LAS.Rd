% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ATL03_ATL08_photons_attributes_dt_LAS.R
\name{ATL03_ATL08_photons_attributes_dt_LAS}
\alias{ATL03_ATL08_photons_attributes_dt_LAS}
\title{Converts ATL03/ATL08 classified photon cloud to LAS}
\usage{
ATL03_ATL08_photons_attributes_dt_LAS(
  atl03_atl08_dt,
  output,
  normalized = TRUE
)
}
\arguments{
\item{atl03_atl08_dt}{An S4 object of class \link{icesat2.atl08_dt} containing ATL03 and ATL08 data
(output of \code{\link[=ATL03_ATL08_photons_attributes_dt_join]{ATL03_ATL08_photons_attributes_dt_join()}} function).}

\item{output}{character. The output path of for the LAS(Z) file(s)
The function will create one LAS file per UTM Zone in WGS84 datum.}

\item{normalized}{logical, default TRUE. Whether the output should be normalized LAS or raw altitude.}
}
\value{
Nothing, it just saves outputs as LAS file in disk
}
\description{
Converts ATL03/ATL08 classified photon cloud to LAS
}
\examples{

# Specifying the path to ATL03 and ATL08 file (zip file)
outdir <- tempdir()
atl03_zip <- system.file("extdata",
  "ATL03_20220401221822_01501506_005_01.zip",
  package = "ICESat2VegR"
)

atl08_zip <- system.file("extdata",
  "ATL08_20220401221822_01501506_005_01.zip",
  package = "ICESat2VegR"
)

# Unzipping ATL03 file
atl03_path <- unzip(atl03_zip, exdir = outdir)

# Unzipping ATL08 file
atl08_path <- unzip(atl08_zip, exdir = outdir)

# Reading ATL03 data (h5 file)

# Reading ATL08 data (h5 file)

# # Extracting ATL03 and ATL08 photons and heights
atl03_atl08_dt <- ATL03_ATL08_photons_attributes_dt_join(atl03_h5, atl08_h5)

ATL03_ATL08_photons_attributes_dt_LAS(
  atl03_atl08_dt,
  file.path(outdir, "output.laz"),
  normalized = TRUE
)

close(atl03_h5)
close(atl08_h5)
}
