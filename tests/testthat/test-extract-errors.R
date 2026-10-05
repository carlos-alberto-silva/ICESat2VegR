test_that("Earth Engine sampling errors are not reported as empty results", {
  failed_sample <- list(size = function() {
    list(getInfo = function() stop("server failure"))
  })
  testthat::local_mocked_bindings(
    extract = function(...) failed_sample,
    .package = "ICESat2VegR"
  )

  expect_error(
    seg_ancillary_extract(stack = NULL, geom = matrix(1, nrow = 1),
                          chunk_size = 1),
    "Earth Engine sampling failed for rows 1-1: server failure"
  )
  expect_error(ee_to_dt(failed_sample), "server failure")
})
