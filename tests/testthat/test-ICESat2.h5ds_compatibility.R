lower_left_lon <- -85
lower_left_lat <- 30.0
upper_right_lon <- -84.0
upper_right_lat <- 31.0
version <- "006"
short_name <- "ATL08"
persist <- TRUE
daterange <- c("2019-08-19", "2019-08-20")
res <- NA
atl08_h5 <- NA
group <- NA
ds <- NA
vals <- NA


res <- ATLAS_dataFinder_cloud(
  short_name,
  lower_left_lon,
  lower_left_lat,
  upper_right_lon,
  upper_right_lat,
  version,
  daterange,
  persist
)

atl08_h5_cloud <- ICESat2.h5_cloud$new(h5 = res[1])
h5_path <- Sys.getenv("ATL08_PATH")
atl08_h5 <- ICESat2.h5_local$new(h5 = h5_path)

test_that("Retrieve dataset dims", {
  expect_true(
    all.equal(
      atl08_h5[["ancillary_data/control"]]$dims,
      atl08_h5_cloud[["ancillary_data/control"]]$dims
    )
  )

  land_segments <- atl08_h5[["/gt1l/land_segments"]]
  land_segments_cloud <- atl08_h5_cloud[["/gt1l/land_segments"]]
  expect_true(
    all.equal(
      land_segments[["delta_time"]]$dims,
      land_segments_cloud[["delta_time"]]$dims
    )
  )
})



test_that("Get type are same", {
  h5dtypes <- list()
  h5dtypes[["ds_metrics"]] <- hdf5r::h5types$H5T_STD_I8LE
  h5dtypes[["ds_surf_type"]] <- hdf5r::h5types$H5T_STD_I32LE
  h5dtypes[["gt1l/land_segments/rgt"]] <- hdf5r::h5types$H5T_STD_I16LE
  h5dtypes[["gt1l/land_segments/ph_ndx_beg"]] <- hdf5r::h5types$H5T_STD_I64LE
  h5dtypes[["orbit_info/orbit_number"]] <- hdf5r::h5types$H5T_STD_U16LE
  h5dtypes[["gt1l/land_segments/canopy/canopy_h_metrics"]] <- hdf5r::h5types$H5T_IEEE_F32LE
  h5dtypes[["orbit_info/crossing_time"]] <- hdf5r::h5types$H5T_IEEE_F64LE

  for (path in names(h5dtypes)) {
    expect_equal(atl08_h5[[path]]$get_type(), atl08_h5_cloud[[path]]$get_type())
  }
})

test_that("Chunk dimensions are same", {
  chunk_dims_latitude <- atl08_h5[["gt2r/land_segments/latitude"]]$chunk_dims
  chunk_dims_canopy <- atl08_h5[["gt2r/land_segments/canopy/canopy_h_metrics"]]$chunk_dims


  expect_true(
    all.equal(
      atl08_h5[["gt2r/land_segments/latitude"]]$chunk_dims,
      atl08_h5_cloud[["gt2r/land_segments/latitude"]]$chunk_dims
    )
  )

  expect_true(
    all.equal(
      atl08_h5[["gt2r/land_segments/canopy/canopy_h_metrics"]]$chunk_dims,
      atl08_h5_cloud[["gt2r/land_segments/canopy/canopy_h_metrics"]]$chunk_dims
    )
  )
})

test_that("Can get fill value", {
  
  expect_equal(
    atl08_h5[["gt2r/land_segments/terrain_flg"]]$get_fill_value(),
    atl08_h5_cloud[["gt2r/land_segments/terrain_flg"]]$get_fill_value()
  )

  expect_true(
    all.equal(
      atl08_h5[["gt2r/land_segments/canopy/canopy_h_metrics"]]$get_fill_value(),
      atl08_h5_cloud[["gt2r/land_segments/canopy/canopy_h_metrics"]]$get_fill_value()
    )
  )
})
