test_that("Earthdata E2E search is available as an opt-in credentialed workflow", {
  skip_if_offline()
  skip_if(
    Sys.getenv("RUN_EARTHDATA_E2E") != "true",
    "Set RUN_EARTHDATA_E2E=true to run this credentialed E2E test."
  )
  skip_if(
    Sys.getenv("EARTHDATA_USERNAME") == "" || Sys.getenv("EARTHDATA_PASSWORD") == "",
    "Earthdata credentials missing."
  )

  lower_left_lon <- -96.0
  lower_left_lat <- 40.0
  upper_right_lon <- -100.0
  upper_right_lat <- 42.0
  daterange <- c("2021-10-02", "2021-10-03")

  atl03_granules <- ATLAS_dataFinder(
    short_name = "ATL03",
    lower_left_lon,
    lower_left_lat,
    upper_right_lon,
    upper_right_lat,
    version = "007",
    daterange = daterange,
    cloud_computing = FALSE
  )
  atl08_granules <- ATLAS_dataFinder(
    short_name = "ATL08",
    lower_left_lon,
    lower_left_lat,
    upper_right_lon,
    upper_right_lat,
    version = "007",
    daterange = daterange,
    cloud_computing = FALSE
  )

  expect_type(atl03_granules, "character")
  expect_type(atl08_granules, "character")
})
