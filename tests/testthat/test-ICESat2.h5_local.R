h5_path <- system.file("extdata", "atl08_clip.h5", package = "ICESat2VegR")

test_that("Can create a ICESat2.h5_local", {
  atl08_h5 <- ICESat2.h5_local$new(h5 = h5_path)
  expect_true(inherits(atl08_h5, "ICESat2.h5_local"))
})

test_that("Can list groups", {
  atl08_h5 <- ICESat2.h5_local$new(h5 = h5_path)
  expect_true("gt1r" %in% atl08_h5$ls())
})

test_that("Can read attribute", {
  atl08_h5 <- ICESat2.h5_local$new(h5 = h5_path)
  message(atl08_h5$attr("short_name"))
  expect_equal(atl08_h5$attr("short_name"), "ATL08")
})

test_that("Beams are fetched", {
  atl08_h5 <- ICESat2.h5_local$new(h5 = h5_path)
  expect_true("gt1r" %in% atl08_h5$beams)
})

test_that("Can open group", {
  atl08_h5 <- ICESat2.h5_local$new(h5 = h5_path)
  group <- atl08_h5[["gt1r"]]
  expect_true(inherits(group, "ICESat2.h5_local"))
})

test_that("Exists work", {
  atl08_h5 <- ICESat2.h5_local$new(h5 = h5_path)
  expect_true(atl08_h5$exists("gt1r"))
  expect_false(atl08_h5$exists("gt1s"))
})

test_that("Group can also list groups", {
  atl08_h5 <- ICESat2.h5_local$new(h5 = h5_path)
  group <- atl08_h5[["gt1r"]]
  expect_true("land_segments" %in% group$ls())
})

test_that("Can open dataset from file", {
  atl08_h5 <- ICESat2.h5_local$new(h5 = h5_path)
  ds <- atl08_h5[["gt1r/land_segments/latitude"]]
  expect_true(inherits(ds, "ICESat2.h5ds_local"))
})

test_that("Can open dataset from group", {
  atl08_h5 <- ICESat2.h5_local$new(h5 = h5_path)
  group <- atl08_h5[["gt1r"]]
  ds <- group[["land_segments/latitude"]]
  expect_true(inherits(ds, "ICESat2.h5ds_local"))
})

test_that("Can find dataset dimension", {
  atl08_h5 <- ICESat2.h5_local$new(h5 = h5_path)
  group <- atl08_h5[["gt1r"]]
  ds <- group[["land_segments/latitude"]]
  expect_true(all(ds$dims > 0))
})

test_that("Can find open one value from dataset", {
  atl08_h5 <- ICESat2.h5_local$new(h5 = h5_path)
  group <- atl08_h5[["gt1r"]]
  ds <- group[["land_segments/latitude"]]
  val <- ds[1]
  message(val)
  expect_true(inherits(val, "numeric"))
})

test_that("Can find open five values from dataset", {
  atl08_h5 <- ICESat2.h5_local$new(h5 = h5_path)
  group <- atl08_h5[["gt1r"]]
  ds <- group[["land_segments/latitude"]]
  val <- ds[1:5]
  message(paste(val, collapse = " "))
  expect_true(length(val) == 5)
})

test_that("Can open entire datase", {
  atl08_h5 <- ICESat2.h5_local$new(h5 = h5_path)
  group <- atl08_h5[["gt1r"]]
  ds <- group[["land_segments/latitude"]]
  vals <- ds[]
  expect_true(length(vals) == ds$dims)
})

test_that("Can open by mask", {
  atl08_h5 <- ICESat2.h5_local$new(h5 = h5_path)
  group <- atl08_h5[["gt1r"]]
  ds <- group[["land_segments/latitude"]]
  vals <- ds[]
  vals2 <- ds[vals < mean(vals)]
  expect_true(length(vals2) < length(vals))
})
