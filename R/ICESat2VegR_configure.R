# ==============================================================================
# Earth Engine project normalization and initialization
# ==============================================================================

#' Normalize and validate a Google Cloud project ID for Earth Engine
#'
#' This is an internal helper used by `tryInitializeEarthEngine()`. It
#' accepts explicit project IDs, numeric project numbers, or pulls from the
#' environment variable `EE_PROJECT`. It:
#'
#' \itemize{
#'   \item Accepts accidental `projects/<id>` strings and strips the prefix.
#'   \item Allows numeric project numbers as valid.
#'   \item Rejects IDs with underscores (suggesting replacement with hyphens).
#'   \item Enforces basic format rules for project IDs.
#' }
#'
#' @param project Character project ID or number (may be `NULL` or empty).
#' @param quiet Logical. If `FALSE`, print informational messages.
#'
#' @return A normalized project ID string.
#' @keywords internal
.ee_normalize_project <- function(project, quiet = FALSE) {
  msg <- function(...) if (!quiet) message("[EE] ", ...)

  if (is.null(project) || is.na(project) || !nzchar(project)) {
    project <- Sys.getenv("EE_PROJECT", unset = "")
    if (!nzchar(project)) {
      stop("No GCP project provided. Pass project= or set EE_PROJECT env var.")
    }
    msg("Using EE_PROJECT from environment: ", project)
  }

  # Accept accidental "projects/<id>" and strip prefix
  if (grepl("^projects\\/", project)) {
    msg("Stripping 'projects/' prefix from project: ", project)
    project <- sub("^projects\\/", "", project)
  }

  # If it's all digits, it's a project number; allowed, return as-is
  if (grepl("^[0-9]+$", project)) {
    return(project)
  }

  # Block underscores (not allowed in GCP project IDs)
  if (grepl("_", project)) {
    suggestion <- gsub("_", "-", project)
    stop(
      "Invalid project ID '", project, "': underscores are not allowed.\n",
      "Use hyphens instead, e.g. project = '", suggestion, "'."
    )
  }

  # Validate format: start letter; 6-30 chars; letters/digits/hyphens; not end with hyphen
  if (!grepl("^[a-z][a-z0-9-]{4,29}$", project) || grepl("-$", project)) {
    stop(
      "Invalid project ID '", project, "'. A valid ID must:\n",
      " - start with a lowercase letter\n",
      " - be 6-30 characters long\n",
      " - contain only lowercase letters, digits, and hyphens\n",
      " - not end with a hyphen"
    )
  }

  project
}

# Lightweight Earth Engine initialization helper
#
# This function is a low-level helper to initialize the Python
# `earthengine-api` client through `reticulate`. It supports both:
#
# \itemize{
#   \item Service Account (SA) authentication via JSON key file.
#   \item User OAuth authentication (browser flow).
# }
#
# It validates and normalizes the Google Cloud project ID, initializes the
# Earth Engine client, and sets `EE_PROJECT` in the environment for
# downstream use (e.g., in Python code).
#
# @param project Character. GCP Project ID or numeric project number. If
#   `NULL` or empty, the environment variable `EE_PROJECT` is used
#   (and must be set).
# @param service_account Optional service account email. If provided, a
#   `keyfile` must also be supplied.
# @param keyfile Optional path to the service account JSON key file (required
#   if `service_account` is provided).
# @param quiet Logical. Suppress messages if `TRUE`.
# @param force_auth Logical. If `TRUE`, force a user OAuth flow prior to
#   calling `ee$Initialize()`.
#
# @return Logical `TRUE` on success; otherwise an error is raised.
# @keywords internal
#
# @examples
# \dontrun{
#   # Using user OAuth and EE_PROJECT from environment:
#   tryInitializeEarthEngine()
#
#   # Using explicit project id and OAuth
#   tryInitializeEarthEngine(project = "my-ee-project")
#
#   # Using a service account
#   tryInitializeEarthEngine(
#     project         = "my-ee-project",
#     service_account = "my-sa@my-ee-project.iam.gserviceaccount.com",
#     keyfile         = "path/to/key.json"
#   )
# }
tryInitializeEarthEngine <- function(
    project = Sys.getenv("EE_PROJECT", unset = NA),
    service_account = NULL,
    keyfile = NULL,
    quiet = FALSE,
    force_auth = FALSE
) {
  if (!requireNamespace("reticulate", quietly = TRUE))
    stop("Package 'reticulate' is required to initialize Earth Engine.")

  ee <- reticulate::import("ee", delay_load = FALSE)
  msg <- function(...) if (!quiet) message("[EE] ", ...)

  # Validate/normalize project
  project <- .ee_normalize_project(project, quiet = quiet)
  Sys.setenv(EE_PROJECT = project)

  # --- Service Account branch -------------------------------------------------
  if (!is.null(service_account) && nzchar(service_account)) {
    if (is.null(keyfile) || !file.exists(keyfile)) {
      stop("Service account specified but 'keyfile' is missing or not found.")
    }
    msg("Initializing Earth Engine with Service Account: ", service_account)
    creds <- ee$ServiceAccountCredentials(service_account, keyfile)
    ok <- tryCatch(
      {
        ee$Initialize(credentials = creds, project = project)
        TRUE
      },
      error = function(e) {
        msg("SA Initialize failed: ", conditionMessage(e))
        FALSE
      }
    )
    if (!ok) {
      stop(
        "Could not initialize Earth Engine with service account.\n",
        "Check that the SA has:\n",
        " - roles/serviceusage.serviceUsageConsumer on project '", project, "'\n",
        " - Earth Engine API is enabled in that project"
      )
    }
    msg("Initialized Earth Engine with project: ", project)
    return(TRUE)
  }

  # --- User OAuth branch ------------------------------------------------------
  if (isTRUE(force_auth)) {
    msg("Forcing OAuth flow...")
    try(ee$Authenticate(), silent = TRUE)
  }

  ok <- tryCatch(
    {
      ee$Initialize(project = project)
      TRUE
    },
    error = function(e) {
      msg("Initialize() failed (", conditionMessage(e), "). Attempting OAuth...")
      FALSE
    }
  )

  if (!ok) {
    auth_ok <- tryCatch(
      {
        ee$Authenticate(); TRUE
      },
      error = function(e) {
        msg("Authenticate failed: ", conditionMessage(e))
        FALSE
      }
    )
    if (!auth_ok) stop("Could not authenticate with Earth Engine (OAuth).")

    ok <- tryCatch(
      {
        ee$Initialize(project = project)
        TRUE
      },
      error = function(e) {
        msg("Initialize after auth failed: ", conditionMessage(e))
        FALSE
      }
    )
    if (!ok) {
      stop(
        "Could not initialize Earth Engine with project '", project, "'.\n",
        "Double-check that:\n",
        " - The project exists and you have access\n",
        " - The Earth Engine API is enabled\n",
        " - Your identity has roles/serviceusage.serviceUsageConsumer on the project"
      )
    }
  }

  msg("Initialized Earth Engine with project: ", project)
  TRUE
}

#' Reload and verify ICESat2VegR cloud capabilities
#'
#' This helper attempts to unload/reload the \pkg{ICESat2VegR} namespace and
#' then checks whether critical Python modules (e.g., `h5py`,
#' `earthaccess`) have been correctly imported into the package-level
#' environment. It is intended primarily for diagnostics after calling
#' `ICESat2VegR_configure()`.
#'
#' The function is conservative and does not throw errors if the reload fails;
#' it simply emits warnings if critical Python bindings appear to be missing.
#'
#' @return Invisibly returns `TRUE`; warnings are issued if something
#'   looks misconfigured.
#' @keywords internal
tryInitializeCloudCapabilities <- function() {
  if ("ICESat2VegR" %in% loadedNamespaces()) {
    try(unloadNamespace("ICESat2VegR"), silent = TRUE)
    try(loadNamespace("ICESat2VegR"),  silent = TRUE)
  }

  # These objects are assumed to be defined in the package's Python bindings;
  # if they are NULL, it suggests imports did not succeed.
  if (is.null(get0("h5py",       envir = asNamespace("ICESat2VegR"), inherits = FALSE)) ||
      is.null(get0("earthaccess", envir = asNamespace("ICESat2VegR"), inherits = FALSE))) {
    warning("The package's earthaccess / h5py functions were not loaded properly.")
    warning("Use ICESat2VegR_configure() and reload the package.")
  }

  invisible(TRUE)
}
