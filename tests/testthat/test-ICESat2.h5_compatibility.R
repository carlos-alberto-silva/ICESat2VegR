h5_path <- Sys.getenv("ATL08_PATH")
atl08_h5 <- NA
group <- NA
ds <- NA
vals <- NA
lower_left_lon <- -85
lower_left_lat <- 30.0
upper_right_lon <- -84.0
upper_right_lat <- 31.0
version <- "006"
short_name <- "ATL08"
persist <- TRUE
daterange <- c("2019-08-19", "2019-08-20")
res2 <- NA
group2 <- NA
ds2 <- NA
vals2 <- NA

atl08_h5 <- ICESat2.h5_local$new(h5 = h5_path)
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

test_that("Groups are equal", {
    list1 <- atl08_h5$ls()
    list2 <- atl08_h5_cloud$ls()
    list3 <- intersect(list1, list2)
    expect_equal(length(list1), length(list2))
    expect_equal(length(list1), length(list3))
})

test_that("Datasets dims are equal", {
    ds1 <- atl08_h5[["gt1l/land_segments/latitude"]]
    ds2 <- atl08_h5_cloud[["gt1l/land_segments/latitude"]]
    expect_true(all.equal(ds1$dims, ds2$dims))
})

test_that("Datasets complex dims are equal", {
    ds1 <- atl08_h5[["gt1l/land_segments/canopy/canopy_h_metrics_abs"]]
    ds2 <- atl08_h5_cloud[["gt1l/land_segments/canopy/canopy_h_metrics_abs"]]
    expect_true(all.equal(ds1$dims, ds2$dims))
})
