test_that("to_vect converts common coordinate table shapes", {
  skip_if_not_installed("terra")

  dt <- data.table::data.table(
    longitude = c(-83.1, -83.2, NA),
    latitude = c(32.1, 32.2, 32.3),
    h_canopy = c(10, 15, 20)
  )

  vect_obj <- to_vect(dt)

  expect_s4_class(vect_obj, "SpatVector")
  expect_equal(nrow(vect_obj), 2)
  expect_true("h_canopy" %in% names(vect_obj))
  expect_match(terra::crs(vect_obj), "4326")
})

test_that("to_vect infers ICESat-2 photon coordinate columns", {
  skip_if_not_installed("terra")

  photons <- data.frame(
    lon_ph = c(-106.1, -106.2),
    lat_ph = c(41.1, 41.2),
    h_ph = c(30, 35)
  )

  vect_obj <- to_vect(photons)

  expect_s4_class(vect_obj, "SpatVector")
  expect_equal(nrow(vect_obj), 2)
  expect_true("h_ph" %in% names(vect_obj))
})

test_that("to_vect supports explicit coordinate names and helpful failures", {
  skip_if_not_installed("terra")

  points <- data.frame(
    x_coord = c(-84, -84.1),
    y_coord = c(29.6, 29.7),
    id = 1:2
  )

  expect_s4_class(to_vect(points, lon = "x_coord", lat = "y_coord"), "SpatVector")
  expect_error(to_vect(data.frame(a = 1, b = 2)), "Could not infer coordinate columns")
  expect_error(to_vect(data.frame(longitude = NA_real_, latitude = NA_real_)), "No finite coordinates")
})

test_that("sampling methods return stable data.table-compatible subsets", {
  dt <- data.table::data.table(
    longitude = seq(-84, -83.91, length.out = 10),
    latitude = seq(29.6, 29.69, length.out = 10),
    h_canopy = rep(c(5, 15), each = 5)
  )
  class(dt) <- c("icesat2.atl08_dt", "data.table", "data.frame")

  set.seed(1)
  sampled_abs <- ICESat2VegR::sample(dt, method = randomSampling(4))
  sampled_frac <- ICESat2VegR::sample(dt, method = randomSampling(0.5))
  sampled_grid <- ICESat2VegR::sample(dt, method = gridSampling(size = 1, grid_size = 0.05))
  sampled_strata <- ICESat2VegR::sample(
    dt,
    method = stratifiedSampling(size = 1, variable = "h_canopy")
  )

  expect_equal(nrow(sampled_abs), 4)
  expect_equal(nrow(sampled_frac), 5)
  expect_true(nrow(sampled_grid) >= 1)
  expect_true(all(c("x_grid", "y_grid") %in% names(sampled_grid)))
  expect_true(nrow(sampled_strata) >= 1)
  expect_true("breaks" %in% names(sampled_strata))
})
