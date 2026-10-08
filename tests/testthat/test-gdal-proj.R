test_that("GDAL uses the installed PROJ data", {
  expect_identical(system.file("proj", package = "ICESat2VegR"), "")

  raster_path <- tempfile(fileext = ".tif")
  on.exit(unlink(raster_path), add = TRUE)

  ds <- createDataset(
    raster_path = raster_path,
    nbands = 1,
    datatype = GDALDataType$GDT_Int32,
    projstring = "EPSG:4326",
    ul_lat = 1, ul_lon = 0, lr_lat = 0, lr_lon = 1,
    res = c(0.1, -0.1),
    nodata = -1,
    co = character()
  )
  expect_true(file.exists(raster_path))
  ds$Close()

  reopened <- GDALOpen(raster_path)
  expect_equal(reopened$GetRasterXSize(), 10)
  expect_equal(reopened$GetRasterYSize(), 10)
  reopened$Close()
})
