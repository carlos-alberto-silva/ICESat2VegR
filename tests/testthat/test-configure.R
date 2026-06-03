test_that("configuration status probe can be assembled without side effects", {
  safely <- function(expr, default = NA) {
    tryCatch(expr, error = function(e) default)
  }

  have_reticulate <- requireNamespace("reticulate", quietly = TRUE)
  py_path <- if (have_reticulate) {
    safely(reticulate::py_discover_config(required_module = NULL)$python, NA_character_)
  } else {
    NA_character_
  }

  status <- list(
    ICESat2VegR_loaded = "package:ICESat2VegR" %in% search(),
    python_used        = py_path,
    h5py               = FALSE,
    earthaccess        = FALSE,
    ee                 = FALSE
  )

  expect_named(
    status,
    c("ICESat2VegR_loaded", "python_used", "h5py", "earthaccess", "ee")
  )
  expect_true(status$ICESat2VegR_loaded)
  expect_type(status$h5py, "logical")
})

test_that("earthdata_login writes a netrc from environment credentials", {
  skip_if_not_installed("reticulate")

  withr::local_envvar(c(
    RETICULATE_PYTHON = NA,
    EARTHDATA_USERNAME = "unit_user",
    EARTHDATA_PASSWORD = "unit_password"
  ))

  outdir <- withr::local_tempdir()
  netrc <- suppressWarnings(earthdata_login(outdir))
  netrc_lines <- readLines(netrc)

  expect_true(file.exists(netrc))
  expect_true(any(grepl("urs\\.earthdata\\.nasa\\.gov", netrc_lines)))
  expect_true(any(grepl("login unit_user", netrc_lines, fixed = TRUE)))
  expect_equal(
    normalizePath(Sys.getenv("NETRC"), winslash = "/"),
    normalizePath(netrc, winslash = "/")
  )
})
