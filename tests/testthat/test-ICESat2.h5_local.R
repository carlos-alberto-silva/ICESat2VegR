h5_path <- Sys.getenv("ATL08_PATH")
atl08_h5 <- NA
group <- NA
ds <- NA
vals <- NA

test_that("Can create a ICESat2.h5_local", {
  atl08_h5 <<- ICESat2.h5_local$new(h5 = h5_path)
  expect_true(inherits(atl08_h5, "ICESat2.h5_local"))
})

test_that("Can list groups", {
  expect_true("gt1l" %in% atl08_h5$ls())
})

test_that("Can open group", {
  group <<- atl08_h5[["gt1l"]]
  expect_true(inherits(group, "ICESat2.h5_local"))
})

test_that("Exists work", {
  expect_true(atl08_h5$exists("gt1l"))
  expect_false(atl08_h5$exists("gt1s"))
})

test_that("Group can also list groups", {
  expect_true("land_segments" %in% group$ls())
})

test_that("Can open dataset from file", {
  ds <<- atl08_h5[["gt1l/land_segments/latitude"]]
  expect_true(inherits(ds, "ICESat2.h5_ds_local"))
})

test_that("Can open dataset from group", {
  ds <<- group[["land_segments/latitude"]]
  expect_true(inherits(ds, "ICESat2.h5_ds_local"))
})

test_that("Can find dataset dimension", {
  expect_true(all(ds$dims > 0))
})

test_that("Can find open one value from dataset", {
  val <- ds[1]
  message(val)
  expect_true(inherits(val, "numeric"))
})

test_that("Can find open five values from dataset", {
  val <- ds[1:5]
  message(paste(val, collapse = " "))
  expect_true(length(val) == 5)
})

test_that("Can open entire datase", {
  vals <<- ds[]
  expect_true(length(vals) == ds$dims)
})

test_that("Can open by mask", {
  vals2 <- ds[vals < mean(vals)]
  expect_true(length(vals2) < length(vals))
})
