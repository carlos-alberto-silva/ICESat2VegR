h5_path <- system.file("extdata", "atl08_clip.h5", package = "ICESat2VegR")
atl08_h5 <- NA
group <- NA
ds <- NA
vals <- NA
atl08_h5 <- ICESat2.h5_local$new(h5 = h5_path)


test_that("Can retrieve dataset dims", {
  expect_equal(atl08_h5[["orbit_info/rgt"]]$dims, 1)

  land_segments <- atl08_h5[["/gt1r/land_segments"]]
  delta_time_dims <- land_segments[["delta_time"]]$dims
  latitude_dims <- land_segments[["latitude"]]$dims
  expect_gt(delta_time_dims, 1)
  expect_equal(delta_time_dims, latitude_dims)
})


test_that("Can get type", {
  h5dtypes <- list()
  h5dtypes[["gt1r/land_segments/rgt"]] <- hdf5r::h5types$H5T_STD_I16LE
  h5dtypes[["gt1r/land_segments/ph_ndx_beg"]] <- hdf5r::h5types$H5T_STD_I64LE
  h5dtypes[["orbit_info/orbit_number"]] <- hdf5r::h5types$H5T_STD_U16LE
  h5dtypes[["gt1r/land_segments/canopy/canopy_h_metrics"]] <- hdf5r::h5types$H5T_IEEE_F32LE
  h5dtypes[["orbit_info/crossing_time"]] <- hdf5r::h5types$H5T_IEEE_F64LE

  for (path in names(h5dtypes)) {
    expect_equal(atl08_h5[[path]]$get_type(), h5dtypes[[path]])
  }
})

test_that("Can get chunk dimensions", {
  chunk_dims_latitude <- atl08_h5[["gt1r/land_segments/latitude"]]$chunk_dims
  chunk_dims_canopy <- atl08_h5[["gt1r/land_segments/canopy/canopy_h_metrics"]]$chunk_dims

  expect_equal(chunk_dims_canopy[1], 18)
})

test_that("Can get fill value", {
  terrain_flg_latitude <- atl08_h5[["gt1r/land_segments/terrain_flg"]]$get_fill_value()
  fill_value_canopy <- atl08_h5[["gt1r/land_segments/canopy/canopy_h_metrics"]]$get_fill_value()

  expect_true(TRUE)
})