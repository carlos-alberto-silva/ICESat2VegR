test_that("bundled ATL03/ATL08 fixtures support the local README workflow primitives", {
  skip_if_not_installed("terra")

  atl03_path <- system.file("extdata", "atl03_clip.h5", package = "ICESat2VegR")
  atl08_path <- system.file("extdata", "atl08_clip.h5", package = "ICESat2VegR")
  skip_if(
    atl03_path == "" || atl08_path == "",
    "Bundled HDF5 fixtures are unavailable."
  )

  atl03 <- ATL03_read(atl03_path)
  atl08 <- ATL08_read(atl08_path)
  on.exit({ close(atl03); close(atl08) }, add = TRUE)

  expect_true(inherits(atl03, "icesat2.atl03_h5"))
  expect_true(inherits(atl08, "icesat2.atl08_h5"))

  atl03_seg_dt <- ATL03_seg_metadata_dt(
    atl03,
    attributes = c("delta_time", "solar_elevation", "pitch", "h_ph", "ref_elev")
  )
  atl08_seg_dt <- ATL08_seg_attributes_dt(
    atl08,
    attributes = c("h_canopy", "h_te_mean", "terrain_slope", "canopy_openness", "night_flag")
  )

  expect_true(inherits(atl03_seg_dt, "icesat2.atl03_seg_dt"))
  expect_true(inherits(atl08_seg_dt, "icesat2.atl08_dt"))
  expect_true(
    all(c("reference_photon_lon", "reference_photon_lat", "beam", "h_ph") %in% names(atl03_seg_dt))
  )
  expect_true(
    all(c("longitude", "latitude", "beam", "h_canopy", "terrain_slope") %in% names(atl08_seg_dt))
  )
  expect_gt(nrow(atl03_seg_dt), 0)
  expect_gt(nrow(atl08_seg_dt), 0)

  atl03_vect <- to_vect(atl03_seg_dt)
  atl08_vect <- to_vect(atl08_seg_dt)
  max_h_canopy <- ATL08_seg_attributes_dt_gridStat(
    atl08_seg_dt,
    func = max(h_canopy),
    res = 0.02
  )

  expect_s4_class(atl03_vect, "SpatVector")
  expect_s4_class(atl08_vect, "SpatVector")
  expect_s4_class(max_h_canopy, "SpatRaster")
  expect_equal(terra::nlyr(max_h_canopy), 1)
})

test_that("joined ATL03/ATL08 photons can be segmented and summarized", {
  skip_if_not_installed("sf")
  skip_if_not_installed("terra")

  atl03_path <- system.file("extdata", "atl03_clip.h5", package = "ICESat2VegR")
  atl08_path <- system.file("extdata", "atl08_clip.h5", package = "ICESat2VegR")
  aoi_path <- system.file("extdata", "aoi_4326.geojson", package = "ICESat2VegR")
  skip_if(
    atl03_path == "" || atl08_path == "" || aoi_path == "",
    "Bundled workflow fixtures are unavailable."
  )

  atl03 <- ATL03_read(atl03_path)
  atl08 <- ATL08_read(atl08_path)
  on.exit({ close(atl03); close(atl08) }, add = TRUE)

  joined <- ATL03_ATL08_photons_attributes_dt_join(atl03, atl08)
  seg20 <- ATL03_ATL08_segment_create(joined, 20, centroid = "mean", output = NA)
  stats20 <- ATL03_ATL08_compute_seg_attributes_dt_segStat(
    seg20,
    list(
      rh98 = stats::quantile(
        ph_h[ph_h > 0 & classed_pc_flag %in% c(2, 3)],
        0.98,
        na.rm = TRUE
      ),
      n_canopy_total = sum(classed_pc_flag >= 2),
      mean_solar     = mean(solar_elevation, na.rm = TRUE)
    ),
    ph_class = c(2, 3)
  )

  expect_true(inherits(joined, "icesat2.atl03atl08_dt"))
  expect_true(inherits(seg20, "icesat2.atl03_atl08_seg_dt"))
  expect_true(inherits(stats20, "icesat2.atl03_atl08_seg_dt"))
  expect_true(
    all(c("segment_id", "beam", "longitude", "latitude", "rh98") %in% names(stats20))
  )
  expect_gt(nrow(stats20), 0)

  boundary <- sf::st_read(aoi_path, quiet = TRUE)
  clipped <- suppressWarnings(ATL03_ATL08_seg_attributes_dt_clipGeometry(stats20, boundary))

  expect_true(inherits(clipped, "icesat2.atl03_atl08_seg_dt"))
  expect_true(nrow(clipped) <= nrow(stats20))
})
