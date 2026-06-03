test_that("fit_metrics computes known error statistics", {
  observed <- c(178, 33, 60, 80, 104, 204, 146)
  predicted <- c(184, 28.5, 55, 85, 105, 210, 155)

  metrics <- fit_metrics(observed, predicted)

  expect_s3_class(metrics, "data.frame")
  expect_named(metrics, c("stat", "value", "unit"))
  expect_equal(metrics$stat, c("rmse", "rmseR", "mae", "maeR", "bias", "biasR", "r", "adj_r2"))
  expect_equal(metrics$value[metrics$stat == "rmse"], 5.6600101, tolerance = 1e-6)
  expect_equal(metrics$value[metrics$stat == "mae"], 5.2142857, tolerance = 1e-6)
  expect_equal(metrics$value[metrics$stat == "bias"], 2.5, tolerance = 1e-6)
  expect_equal(metrics$unit[metrics$stat == "rmseR"], "%")
})

test_that("fit_metrics validates input shape and complete pairs", {
  expect_error(fit_metrics(c(1, 2), c(1, 2, 3)), "same length")
  expect_error(fit_metrics(c(1, NA), c(1, 2)), "at least 2 complete pairs")

  metrics <- fit_metrics(c(-1, 1, NA), c(-2, 2, 10))
  expect_true(is.na(metrics$value[metrics$stat == "rmseR"]))
  expect_true(is.na(metrics$value[metrics$stat == "maeR"]))
  expect_true(is.na(metrics$value[metrics$stat == "biasR"]))
})

test_that("ATL08 grid statistics summarize canopy values from local tables", {
  skip_if_not_installed("terra")

  dt <- data.table::data.table(
    longitude = c(-83.101, -83.102, -83.201),
    latitude = c(32.101, 32.102, 32.201),
    beam = c("gt1r", "gt1r", "gt2r"),
    strong_beam = c(TRUE, TRUE, FALSE),
    h_canopy = c(10, 20, 30)
  )
  class(dt) <- c("icesat2.atl08_dt", "data.table", "data.frame")

  grid <- ATL08_seg_attributes_dt_gridStat(
    dt,
    func = max(h_canopy),
    res = 0.05
  )

  expect_s4_class(grid, "SpatRaster")
  expect_equal(terra::nlyr(grid), 1)
  expect_true(any(terra::values(grid, mat = FALSE) %in% c(20, 30), na.rm = TRUE))
})
