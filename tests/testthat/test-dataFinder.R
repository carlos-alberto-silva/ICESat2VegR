test_that("ATLAS_dataDownload validates arguments before prompting for credentials", {
  expect_error(ATLAS_dataDownload(character(), overwrite = "yes"), "`overwrite`")
  expect_error(ATLAS_dataDownload(character(), buffer_size = 0), "`buffer_size`")
  expect_error(ATLAS_dataDownload(character(), timeout = -1), "`timeout`")
  expect_error(ATLAS_dataDownload(character(), retries = 0), "`retries`")
  expect_error(ATLAS_dataDownload(character(), backoff = 0), "`backoff`")
  expect_error(ATLAS_dataDownload(character(), quiet = "yes"), "`quiet`")
})

test_that("ATLAS_dataDownload handles empty and skipped URLs deterministically", {
  skip_if_not_installed("reticulate")

  withr::local_envvar(c(
    RETICULATE_PYTHON = NA,
    EARTHDATA_USERNAME = "unit_user",
    EARTHDATA_PASSWORD = "unit_password"
  ))

  outdir <- withr::local_tempdir()
  existing <- file.path(outdir, "ATL08_existing.h5")
  file.create(existing)

  status <- suppressWarnings(
    suppressMessages(
      ATLAS_dataDownload(
        c(NA_character_, "", "https://example.com/ATL08_existing.h5"),
        outdir = outdir,
        quiet = TRUE
      )
    )
  )

  expect_equal(status$status, c("skipped", "skipped", "skipped"))
  expect_equal(status$attempts, c(0L, 0L, 0L))
  expect_equal(status$file[3], "ATL08_existing.h5")
})

test_that("Earthdata search can run when credentials and network are available", {
  skip_if_offline()
  skip_if(
    Sys.getenv("RUN_EARTHDATA_E2E") != "true",
    "Set RUN_EARTHDATA_E2E=true to run Earthdata search tests."
  )
  skip_if(
    Sys.getenv("EARTHDATA_USERNAME") == "" || Sys.getenv("EARTHDATA_PASSWORD") == "",
    "Earthdata credentials missing."
  )

  granules <- ATLAS_dataFinder(
    short_name = "ATL08",
    lower_left_lon = -96.0,
    lower_left_lat = 40.0,
    upper_right_lon = -100.0,
    upper_right_lat = 42.0,
    version = "007",
    daterange = c("2021-10-02", "2021-10-03"),
    cloud_computing = FALSE
  )

  expect_type(granules, "character")
})
