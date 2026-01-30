require(ICESat2VegR)
require(data.table)
# Specifying the path to ATL03 file
atl03_path <- system.file("extdata",
  "atl03_clip.h5",
  package = "ICESat2VegR"
)

# Reading ICESat-2 ATL03 data (h5 file)
ATL03 <- ATL03_read(atl03_path = atl03_path)
close(ATL03$h5)

# Specifying the path to ATL08 file
atl08_path <- system.file("extdata",
  "atl08_clip.h5",
  package = "ICESat2VegR"
)

# Reading ICESat-2 ATL08 data (h5 file)
atl08 <- ATL08_read(atl08_path = atl08_path)
close(atl08)


atl03_path <- system.file(
  "extdata", "atl03_clip.h5",
  package = "ICESat2VegR"
)

atl03_h5 <- ATL03_read(atl03_path)

# Extract ATL03 geolocation segment metadata
atl03_segment_dt <- ATL03_seg_metadata_dt(atl03_h5)

head(atl03_segment_dt)
close(atl03_h5)

atl03_photons_dt <- ATL03_photons_attributes_dt(atl03_h5)
head(atl03_photons_dt)


# ATL03 file path
atl03_path <- system.file("extdata",
  "atl03_clip.h5",
  package = "ICESat2VegR"
)

# Reading ATL03 data (h5 file)
atl03_h5 <- ATL03_read(atl03_path = atl03_path)

# Extracting ATL03 and ATL08 photons and heights
atl03_dt <- ATL03_photons_attributes_dt(atl03_h5, beam = "gt1r")

outdir <- tempdir()
ATL03_photons_attributes_dt_LAS(
  atl03_dt,
  file.path(outdir, "output.laz")
)


# ATL03 file path
atl03_path <- system.file("extdata",
  "atl03_clip.h5",
  package = "ICESat2VegR"
)

# Reading ATL03 data (h5 file)
atl03_h5 <- ATL03_read(atl03_path = atl03_path)

# Extracting ATL03 photon attributes
atl03_photons_dt <- ATL03_photons_attributes_dt(atl03_h5 = atl03_h5)

# Specifying the path to shapefile
polygon_filepath <-
  system.file(
    "extdata",
    "clip_geom.shp",
    package = "ICESat2VegR"
  )

# Reading shapefile as sf object
sppoly <- terra::vect(polygon_filepath)
plot(sppoly, add=T)
windows()
plot(atl03_photons_dt$lon_ph,atl03_photons_dt$lat_ph)

# Clipping ATL03 photon attributes by Geometry
atl03_photons_dt_clip <-
  ATL03_photons_attributes_dt_clipGeometry(atl03_photons_dt, sppoly, split_by = "id")

atl03_photons_dt_clip <-clip(atl03_photons_dt, sppoly)

head(atl03_photons_dt_clip)


points(atl03_photons_dt$lon_ph,atl03_photons_dt$lat_ph, col=atl03_photons_dt_clip$poly_id)

library(terra)
library(sf)
library(data.table)

# --- 1. Read polygon shapefile ---
clip_obj <- vect(polygon_filepath)
plot(clip_obj)   # initial plot

# --- 2. Convert points to spatial (SpatVector) ---
pts <- vect(data.table(atl03_photons_dt_clip[,c("lon_ph", "lat_ph")]),
            geom = c("lon_ph", "lat_ph"))   # IMPORTANT: match CRS

# Optional: view
plot(pts, add = TRUE, pch = 20, col = "blue")


# --- 4. Plot points colored by polygon ---
windows()
plot(clip_obj)
points(pts, col = atl03_photons_dt_clip$poly_id, pch = 20)

#'
#' close(atl03_h5)
