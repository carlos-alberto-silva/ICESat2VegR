library(testthat)


atl03_path <- system.file(
  "extdata",
  "atl03_clip.h5",
  package = "ICESat2VegR"
)
atl03_h5 <- ATL03_read(atl03_path)
atl03_seg_dt <- ATL03_seg_attributes_dt(atl03_h5)

atl03_seg_vec <- to_vect(atl03_seg_dt)

linear_model <- stats::lm(h_ph ~ segment_ph_cnt, data = atl03_seg_dt)
atl03_seg_vec$predict <- predict(linear_model, atl03_seg_dt[, .(segment_ph_cnt)])
writeVector(atl03_seg_vec, "../atl03_seg_vect.gpkg", overwrite = TRUE)
output_h5 <- tempfile(fileext = ".h5")
predicted_h5 <- predict_h5(linear_model, atl03_seg_dt, output_h5)

h5_input <- predicted_h5
res <- 0.0005
bbox <- terra::ext(min(h5_input[["longitude"]][]), max(h5_input[["longitude"]][]), min(h5_input[["latitude"]][]), max(h5_input[["latitude"]][]))


test_that("rasterize_h5 exists", {
  expect_true(isGeneric("rasterize_h5"))
})

test_that("rasterize_h5 accepts predict_h5, output, bbox, res", {
  expect_true(existsMethod("rasterize_h5", signature = c("icesat2.predict_h5", output = "character", bbox = "SpatExtent", res = "numeric")))
})


test_that("rasterize_h5 can be called with icesat2.predict_h5", {
  output <- tempfile(fileext = ".tif")

  expect_error(rasterize_h5(predicted_h5))
  expect_error(rasterize_h5(predicted_h5, output))
  expect_error(rasterize_h5(predicted_h5, output, bbox = bbox))
  rasterize_h5(predicted_h5, output, bbox = bbox, res = res)

  close(predicted_h5)
})
