library(testthat)

test_that("predict_h5 exists", {
  expect_true(isGeneric("predict_h5"))
})


test_that("predict_h5 accepts model, seg_dt and output", {
  expect_true(existsMethod("predict_h5", signature = c("ANY", "icesat2.atl03_seg_dt", "character")))
})


test_that("predict_h5 accepts model, atl08_dt and output", {
  expect_true(existsMethod("predict_h5", signature = c("ANY", "icesat2.atl08_dt", "character")))
})


test_that("predict_h5 can be called with icesat2.atl03_seg_dt", {
  atl03_path <- system.file(
    "extdata",
    "atl03_clip.h5",
    package = "ICESat2VegR"
  )

  atl03_h5 <- ATL03_read(atl03_path = atl03_path)
  atl03_seg_dt <- ATL03_seg_attributes_dt(atl03_h5)

  linear_model <- stats::lm(h_ph ~ segment_ph_cnt, data = atl03_seg_dt)

  output_h5 <- tempfile(fileext = ".h5")

  predict_h5(linear_model, atl03_seg_dt, output_h5)
  expect_true(TRUE)
})

test_that("predict_h5 on icesat2.atl03_seg_dt will output longitude latitude and predicted groups", {
  atl03_path <- system.file(
    "extdata",
    "atl03_clip.h5",
    package = "ICESat2VegR"
  )

  atl03_h5 <- ATL03_read(atl03_path = atl03_path)
  atl03_seg_dt <- ATL03_seg_attributes_dt(atl03_h5)

  linear_model <- stats::lm(h_ph ~ segment_ph_cnt, data = atl03_seg_dt)
  output_h5 <- tempfile(fileext = ".h5")
  predicted_h5 <- predict_h5(linear_model, atl03_seg_dt, output_h5)

  expect_true(all(predicted_h5$ls()$name %in% c("latitude", "longitude", "prediction")))
})

test_that("predict_h5 on icesat2.atl03_seg_dt latitude column is the same", {
  atl03_path <- system.file(
    "extdata",
    "atl03_clip.h5",
    package = "ICESat2VegR"
  )

  atl03_h5 <- ATL03_read(atl03_path = atl03_path)
  atl03_seg_dt <- ATL03_seg_attributes_dt(atl03_h5)

  linear_model <- stats::lm(h_ph ~ segment_ph_cnt, data = atl03_seg_dt)
  output_h5 <- tempfile(fileext = ".h5")
  predicted_h5 <- predict_h5(linear_model, atl03_seg_dt, output_h5)

  expect_true(all(predicted_h5[["latitude"]][] == atl03_seg_dt$reference_photon_lat))
  close(predicted_h5)
})

test_that("predict_h5 can be called with icesat2.atl08_dt", {
  atl08_path <- system.file(
    "extdata",
    "atl08_clip.h5",
    package = "ICESat2VegR"
  )

  atl08_h5 <- ATL08_read(atl08_path = atl08_path)
  atl08_dt <- ATL08_seg_attributes_dt(atl08_h5)

  linear_model <- stats::lm(h_canopy ~ canopy_openness, data = atl08_dt)

  output_h5 <- tempfile(fileext = ".h5")

  predicted_h5 <- predict_h5(linear_model, atl08_dt, output_h5)
  expect_true(TRUE)
  
  close(atl08_h5)
  close(predicted_h5)
})

test_that("predict_h5 on icesat2.atl08_dt will output longitude latitude and predicted groups", {
  atl08_path <- system.file(
    "extdata",
    "atl08_clip.h5",
    package = "ICESat2VegR"
  )

  atl08_h5 <- ATL08_read(atl08_path = atl08_path)
  atl08_dt <- ATL08_seg_attributes_dt(atl08_h5)

  linear_model <- stats::lm(h_canopy ~ canopy_openness, data = atl08_dt)
  output_h5 <- tempfile(fileext = ".h5")
  predicted_h5 <- predict_h5(linear_model, atl08_dt, output_h5)

  expect_true(all(predicted_h5$ls()$name %in% c("latitude", "longitude", "prediction")))
  close(predicted_h5)
})

test_that("predict_h5 on icesat2.atl08_dt latitude column is the same", {
  atl08_path <- system.file(
    "extdata",
    "atl08_clip.h5",
    package = "ICESat2VegR"
  )

  atl08_h5 <- ATL08_read(atl08_path = atl08_path)
  atl08_dt <- ATL08_seg_attributes_dt(atl08_h5)

  linear_model <- stats::lm(h_canopy ~ canopy_openness, data = atl08_dt)
  output_h5 <- tempfile(fileext = ".h5")
  predicted_h5 <- predict_h5(linear_model, atl08_dt, output_h5)

  expect_true(all(predicted_h5[["latitude"]][] == atl08_dt$reference_photon_lat))
  close(predicted_h5)
})


# testthat::test_file("tests//testthat/test-predict_h5.R")

# test_that("predict_h5 accepts model, seg_dt and outdir", {
#   atl08_h5_path <- system.file(
#     "extdata",
#     "atl08_clip.h5",
#     package = "ICESat2VegR"
#   )

#   atl08_h5 <- ATL08_read(atl08_path = atl08_h5_path)

#   exists("predict_h5")
# })
