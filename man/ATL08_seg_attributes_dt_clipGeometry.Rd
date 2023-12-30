% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ATL08_seg_attributes_dt_clipGeometry.R
\name{ATL08_seg_attributes_dt_clipGeometry}
\alias{ATL08_seg_attributes_dt_clipGeometry}
\title{Clip ATL08 Terrain and Canopy Attributes by Geometry}
\usage{
ATL08_seg_attributes_dt_clipGeometry(atl08_seg_att_dt, xmin, xmax, ymin, ymax)
}
\arguments{
\item{atl08_seg_att_dt}{A atl08_seg_att_dt object (output of \code{\link[=ATL08_seg_attribute_dt]{ATL08_seg_attribute_dt()}} function).
An S4 object of class \link{icesat2.atl08_dt}}

\item{polygon}{Polygon. An object of class \code{\link[terra:SpatVector-class]{terra::SpatVector}},
which can be loaded as an ESRI shapefile using \link[terra:vect]{terra::vect} function in the
\emph{sf} package.}

\item{split_by}{Polygon id. If defined, ATL08 data will be clipped by each polygon using
the polygon id from table of attribute defined by the user}
}
\value{
Returns an S4 object of class \link{icesat2.atl08_dt}
containing the clipped ATL08 Terrain and Canopy Attributes.
}
\description{
This function clips ATL08 Terrain and Canopy Attributes within a given geometry
}
\examples{
# Specifying the path to ATL08 file (zip file)
outdir <- tempdir()
atl08_zip <- system.file("extdata",
  "ATL08_20220401221822_01501506_005_01.zip",
  package = "rICESat2Veg"
)

# Unzipping ATL08 file
atl08_path <- unzip(atl08_zip, exdir = outdir)

# Reading ATL08 data (h5 file)

# Extracting ATL08-derived Canopy Metrics
atl08_seg_att_dt <- ATL08_seg_attributes_dt(atl08_h5 = atl08_h5)

# Specifying the path to shapefile
polygon_filepath <- system.file("extdata", "polygon.shp", package = "rICESat2Veg")

# Reading shapefile as sf object
polygon <- terra::vect(polygon_filepath)

# Clipping ATL08 Terrain and Canopy Attributes by Geometry
atl08_seg_att_dt_clip <- ATL08_seg_attributes_dt_clipGeometry(atl08_seg_att_dt, polygon, split_by = "FID")

hasLeaflet <- require(leaflet)

if (hasLeaflet) {
  leaflet() \%>\%
    addCircleMarkers(atl08_seg_att_dt_clip$longitude,
      atl08_seg_att_dt_clip$latitude,
      radius = 1,
      opacity = 1,
      color = "red"
    ) \%>\%
    addScaleBar(options = list(imperial = FALSE)) \%>\%
    addPolygons(
      data = polygon, weight = 1, col = "white",
      opacity = 1, fillOpacity = 0
    ) \%>\%
    addProviderTiles(providers$Esri.WorldImagery)
}
close(atl08_h5)
}