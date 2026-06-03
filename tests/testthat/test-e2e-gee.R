test_that("GEE E2E configures an isolated Python environment and initializes Earth Engine", {
  skip_if(
    Sys.getenv("RUN_GEE_E2E") != "true",
    "Set RUN_GEE_E2E=true to run Google Earth Engine E2E tests."
  )
  skip_if_not_installed("reticulate")
  skip_if(
    Sys.getenv("EE_PROJECT") == "",
    "Set EE_PROJECT to a Google Earth Engine-enabled project."
  )

  repo_root <- normalizePath(testthat::test_path("..", ".."), winslash = "/")
  project <- Sys.getenv("EE_PROJECT")
  env_root <- file.path(tempdir(), "icesat2vegr-reticulate-envs")
  env_name <- paste0("icesat2vegr-gee-", as.integer(Sys.time()))
  rscript <- file.path(
    R.home("bin"),
    if (.Platform$OS.type == "windows") "Rscript.exe" else "Rscript"
  )

  # Run in a fresh R process because reticulate binds to one Python per session.
  script <- tempfile(fileext = ".R")
  writeLines(c(
    "setwd(normalizePath(commandArgs(TRUE)[[1]], winslash = '/'))",
    "project <- commandArgs(TRUE)[[2]]",
    "env_root <- commandArgs(TRUE)[[3]]",
    "env_name <- commandArgs(TRUE)[[4]]",
    "Sys.setenv(EE_PROJECT = project)",
    "Sys.setenv(RETICULATE_VIRTUALENV_ROOT = env_root)",
    "Sys.unsetenv('RETICULATE_PYTHON')",
    "Sys.unsetenv('RETICULATE_PYTHON_ENV')",
    "suppressPackageStartupMessages(devtools::load_all('.', quiet = TRUE))",
    "ok <- ICESat2VegR_configure(",
    "  envname = env_name,",
    "  prefer_conda = FALSE,",
    "  allow_session_installs = FALSE,",
    "  manage_aiobotocore = FALSE,",
    "  auto_restart = FALSE,",
    "  verbose = TRUE",
    ")",
    "if (!isTRUE(ok)) stop('ICESat2VegR_configure() did not return TRUE')",
    "cfg <- reticulate::py_config()",
    "message('Configured Python: ', cfg$python)",
    "mods <- c('h5py', 'earthaccess', 'ee')",
    "missing <- mods[!vapply(mods, reticulate::py_module_available, logical(1))]",
    "if (length(missing)) stop('Missing configured Python modules: ', paste(missing, collapse = ', '))",
    "if (!isTRUE(ee_initialize(project = project, quiet = TRUE))) stop('ee_initialize() did not return TRUE')",
    "ee_mod <- reticulate::import('ee', delay_load = FALSE)",
    "if (!isTRUE(.ee_ping(ee_mod))) stop('.ee_ping() did not return TRUE')"
  ), script)

  out <- system2(
    rscript,
    args = c(
      shQuote(script),
      shQuote(repo_root),
      shQuote(project),
      shQuote(env_root),
      shQuote(env_name)
    ),
    stdout = TRUE,
    stderr = TRUE
  )
  status <- attr(out, "status")
  if (is.null(status)) status <- 0L

  expect_equal(
    status,
    0L,
    info = paste(out, collapse = "\n")
  )
  expect_true(any(grepl("Configured Python:", out, fixed = TRUE)))
})

test_that("README documents the manual GEE workflow", {
  readme_path <- testthat::test_path("..", "..", "README.md")

  expect_true(file.exists(readme_path))
  readme <- readLines(readme_path, warn = FALSE)
  expect_true(any(grepl("ee_build_AlphaEarth_embedding_terrain_stack", readme, fixed = TRUE)))
  expect_true(any(grepl("map_download", readme, fixed = TRUE)))
})
