test_that("ANN spatial indexes survive repeated creation and collection", {
  x <- c(0, 1, 2, 3)
  y <- c(0, 1, 2, 3)

  for (i in seq_len(5L)) {
    index <- ANNIndex$new(x, y)
    nearby <- index$searchFixedRadius(0, 0, 0.01)
    expect_type(nearby, "integer")
    expect_true(0L %in% nearby)
    rm(index)
    gc()
  }
})

test_that("ANN radius and spaced sampling use the requested distance", {
  dt <- data.table::data.table(
    longitude = c(0, 0.05),
    latitude = c(0, 0),
    h_canopy = c(1, 2)
  )
  class(dt) <- c("icesat2.atl08_dt", "data.table", "data.frame")

  index <- ANNIndex$new(dt$longitude, dt$latitude)
  expect_setequal(index$searchFixedRadius(0, 0, 0.01), 0L)
  expect_setequal(index$searchFixedRadius(0, 0, 0.06), c(0L, 1L))

  set.seed(1)
  spaced <- ICESat2VegR::sample(dt, method = spacedSampling(size = 2, radius = 0.01))
  expect_equal(nrow(spaced), 2L)
})
