# ==============================================================================
# ICESat2VegR — Python Environment Configuration + Earth Engine Initialization
# ==============================================================================

#' Configure Python environment for ICESat2VegR cloud features
#'
#' This function sets up a robust, persistent Python environment for the
#' cloud-related features of **ICESat2VegR**, including:
#'
#' * Downloading ICESat-2 data via \code{earthaccess}
#' * Interfacing with Google Earth Engine via \code{earthengine-api}
#' * Handling S3 access via \code{boto3}/\code{botocore}
#'
#' The function:
#'
#' \itemize{
#'   \item Detects whether Python is already initialized in the R session.
#'   \item Performs a fast "quick probe" to see if required modules are present
#'         (\code{h5py}, \code{earthaccess}, \code{ee} by default).
#'   \item If everything is already available, it returns immediately.
#'   \item Otherwise, it:
#'     \itemize{
#'       \item Ensures Miniconda is installed (if \code{prefer_conda = TRUE}).
#'       \item Creates or reuses a conda environment (default: \code{"icesat2-env"}).
#'       \item Installs required Python modules using conda-forge (and pip fallback).
#'       \item Optionally tries in-session installs if Python is already initialized.
#'       \item Optionally resolves common \code{aiobotocore}/\code{botocore} conflicts.
#'       \item Optionally restarts the R session in RStudio to bind the new env.
#'     }
#' }
#'
#' Typical usage is simply:
#'
#' \preformatted{
#'   ICESat2VegR_configure()
#' }
#'
#' You usually only need to run this:
#' \itemize{
#'   \item Once per machine, or
#'   \item After major Python / Miniconda changes, or
#'   \item If the package reports missing Python modules.
#' }
#'
#' @section Checked/installed Python modules:
#'
#' The function checks for the following import names and installs the
#' corresponding packages when missing:
#'
#' \tabular{ll}{
#'   \strong{Import name}      \tab \strong{Install name}           \cr
#'   \code{numpy}              \tab \code{numpy}                    \cr
#'   \code{h5py}               \tab \code{h5py}                     \cr
#'   \code{packaging}          \tab \code{packaging}                \cr
#'   \code{ee}                 \tab \code{earthengine-api}          \cr
#'   \code{earthaccess}        \tab \code{earthaccess}              \cr
#'   \code{googleapiclient}    \tab \code{google-api-python-client} \cr
#'   \code{boto3}              \tab \code{boto3}                    \cr
#'   \code{botocore}           \tab \code{botocore}                 \cr
#'   \code{s3transfer}         \tab \code{s3transfer}               \cr
#' }
#'
#' @param envname Character. Name of the conda/virtualenv environment to use or
#'   create. Default is \code{"icesat2-env"}. If this is the default value and
#'   the environment variable \code{CONDA_DEFAULT_ENV} is set, that value is
#'   used instead (useful when running inside an existing conda environment).
#' @param py_ver Character. Python version to request when creating a new env
#'   (default \code{"3.10"}).
#' @param prefer_conda Logical. If \code{TRUE} (default), use conda/Miniconda.
#'   If \code{FALSE}, a \code{virtualenv} is used instead.
#' @param force_conda_forge Logical. If \code{TRUE} (default), configure conda
#'   to use \code{conda-forge} with strict channel priority (recommended).
#' @param allow_session_installs Logical. If \code{TRUE} (default) and Python
#'   is already initialized in the current R session, attempt session-only
#'   installs first via \code{reticulate::py_require()} and, if needed, a
#'   pip call from the active interpreter. If that fails, a persistent env
#'   is prepared instead.
#' @param manage_aiobotocore Logical. If \code{TRUE} (default), attempt to
#'   detect and resolve common \code{aiobotocore}/\code{botocore}/\code{boto3}
#'   version conflicts by pinning compatible versions and optionally
#'   uninstalling \code{aiobotocore} (which is not required for ICESat2VegR).
#' @param auto_restart Logical. If \code{TRUE} (default), attempts to restart
#'   the R session in RStudio via \code{rstudioapi::restartSession()} after
#'   preparing the conda environment, so that the new Python binding is cleanly
#'   established.
#' @param verbose Logical. If \code{TRUE} (default), print progress messages.
#' @param quick_probe Character vector of Python import names to test for a
#'   fast early exit (default: \code{c("h5py", "earthaccess", "ee")}).
#'
#' @return Logical \code{TRUE} on success, \code{FALSE} otherwise.
#' @export
#'
#' @examples
#' \dontrun{
#'   # Basic configuration using defaults
#'   ICESat2VegR_configure()
#'
#'   # Use existing conda env named "icesat2"
#'   ICESat2VegR_configure(envname = "icesat2")
#' }
ICESat2VegR_configure <- function(
    envname = "icesat2-env",
    py_ver = "3.10",
    prefer_conda = TRUE,
    force_conda_forge = TRUE,
    allow_session_installs = TRUE,
    manage_aiobotocore = TRUE,
    auto_restart = TRUE,
    verbose = TRUE,
    quick_probe = c("h5py", "earthaccess", "ee")
) {

  # ---------------------------------------------------------------------------
  # Internal helper: conditional message
  # ---------------------------------------------------------------------------
  say <- function(...) if (isTRUE(verbose)) message("[ICESat2VegR] ", ...)

  # ---------------------------------------------------------------------------
  # Internal helper: quick probe to see if environment is already OK
  # ---------------------------------------------------------------------------
  .quick_probe <- function(mods = quick_probe) {
    if (!requireNamespace("reticulate", quietly = TRUE))
      return(FALSE)

    cfg_ok <- !inherits(try(reticulate::py_config(), silent = TRUE), "try-error")
    if (!cfg_ok) return(FALSE)

    all(vapply(mods, reticulate::py_module_available, logical(1)))
  }

  # ---------------------------------------------------------------------------
  # FAST PATH: If everything is already configured, exit early and quietly
  # ---------------------------------------------------------------------------
  if (.quick_probe()) {
    say("Python environment already has required modules (quick probe).")
    return(TRUE)
  }

  # ---------------------------------------------------------------------------
  # Local helpers for the full configuration flow
  # ---------------------------------------------------------------------------

  # Detect whether Python is initialized in the current session
  py_initialized <- function() {
    out <- FALSE
    try({ reticulate::py_config(); out <- TRUE }, silent = TRUE)
    out
  }

  # Check if conda is available
  has_conda <- function() {
    tryCatch(!is.null(reticulate::conda_version()), error = function(e) FALSE)
  }

  # Locate conda binary
  conda_bin <- function() {
    tryCatch(reticulate::conda_binary(), error = function(e) NULL)
  }

  # Run a conda command, returning status code (0 = success)
  conda_run <- function(args) {
    cb <- conda_bin()
    if (is.null(cb)) return(invisible(1L))
    res <- suppressWarnings(system2(cb, args, stdout = TRUE, stderr = TRUE))
    status <- attr(res, "status")
    if (is.null(status)) status <- 0L
    invisible(status)
  }

  # Missing Python imports checker
  missing_imports <- function(import_names) {
    import_names[!vapply(import_names, reticulate::py_module_available, logical(1))]
  }

  # pip check: look for dependency inconsistencies
  pip_check <- function() {
    tryCatch(
      reticulate::py_run_string(
        "import sys, subprocess; subprocess.check_call([sys.executable,'-m','pip','check'])"
      ),
      error = function(e) e
    )
  }

  # Summarize Python binding
  print_summary <- function(bound_env = NULL) {
    cfg <- tryCatch(reticulate::py_config(), error = function(e) NULL)
    if (is.null(cfg)) {
      say("Python not initialized.")
      return(invisible())
    }
    say(sprintf("Python:  %s", cfg$python))
    if (!is.null(bound_env)) say(sprintf("Env:     %s", bound_env))
    say(sprintf("Version: %s", cfg$version))
  }

  # Attempt to restart R session in RStudio
  try_restart <- function() {
    if (!isTRUE(auto_restart)) return(invisible(FALSE))
    if (requireNamespace("rstudioapi", quietly = TRUE)) {
      say("Restarting R session to bind the new environment…")
      try(rstudioapi::restartSession(), silent = TRUE)
      return(invisible(TRUE))
    }
    FALSE
  }

  # Map Python import names -> install names
  install_name <- c(
    numpy           = "numpy",
    h5py            = "h5py",
    packaging       = "packaging",
    ee              = "earthengine-api",
    earthaccess     = "earthaccess",
    googleapiclient = "google-api-python-client",
    boto3           = "boto3",
    botocore        = "botocore",
    s3transfer      = "s3transfer"
  )

  py_imports <- names(install_name)

  # Allow binding to current conda env if envname is default
  if (envname == "icesat2-env" && nzchar(Sys.getenv("CONDA_DEFAULT_ENV"))) {
    envname <- Sys.getenv("CONDA_DEFAULT_ENV")
  }

  # ---------------------------------------------------------------------------
  # MAIN FLOW
  # ---------------------------------------------------------------------------
  is_init <- py_initialized()

  # ----------------------------- FRESH SESSION --------------------------------
  if (!is_init) {
    # ---- Path 1: use conda/Miniconda ----------------------------------------
    if (prefer_conda) {

      # Ensure Miniconda is available
      if (!has_conda()) {
        say("Miniconda not found; installing Miniconda…")
        reticulate::install_miniconda()
        if (!has_conda()) {
          warning("Miniconda installation failed.")
          return(FALSE)
        }
      }

      # Configure conda-forge as the only/primary channel (strict priority)
      if (isTRUE(force_conda_forge)) {
        suppressWarnings(conda_run(c("config", "--remove-key", "default_channels")))
        suppressWarnings(conda_run(c("config", "--remove-key", "channels")))
        conda_run(c("config", "--add", "channels", "conda-forge"))
        conda_run(c("config", "--set", "channel_priority", "strict"))
      }

      # Create env if needed
      envs <- tryCatch(reticulate::conda_list()$name, error = function(e) character())
      if (!(envname %in% envs)) {
        say(sprintf("Creating conda env '%s' (Python %s)…", envname, py_ver))
        st <- conda_run(c("create", "--yes", "--name", envname, paste0("python=", py_ver)))
        if (!identical(st, 0L)) {
          warning(sprintf("Failed to create env '%s'.", envname))
          return(FALSE)
        }
      }

      # Bind to this environment
      reticulate::use_condaenv(envname, required = TRUE)
      say(sprintf("Bound to '%s'. Checking and installing missing dependencies…", envname))

      # Install missing imports via conda (pip fallback)
      miss <- missing_imports(py_imports)
      if (length(miss)) {
        pkgs <- unname(install_name[miss])
        ok <- TRUE
        tryCatch(
          reticulate::conda_install(envname = envname, packages = pkgs, channel = "conda-forge"),
          error = function(e) { ok <<- FALSE }
        )
        if (!ok) {
          say("conda solver failed; falling back to pip for missing packages…")
          reticulate::py_install(pkgs, method = "pip", pip = TRUE)
        }
      }

    } else {
      # ---- Path 2: virtualenv route -----------------------------------------
      if (!reticulate::virtualenv_exists(envname)) {
        say(sprintf("Creating virtualenv '%s'…", envname))
        reticulate::virtualenv_create(envname)
      }
      reticulate::use_virtualenv(envname, required = TRUE)
      say("Checking and installing missing dependencies into virtualenv via pip…")
      miss <- missing_imports(py_imports)
      if (length(miss)) {
        pkgs <- unname(install_name[miss])
        reticulate::virtualenv_install(envname, packages = pkgs, ignore_installed = FALSE)
      }
    }

    # ----------------------- PYTHON ALREADY INITIALIZED -----------------------
  } else {
    miss <- missing_imports(py_imports)
    if (length(miss) == 0) {
      say("Active Python already has all required modules.")
      print_summary()
      return(TRUE)
    }

    pkgs <- unname(install_name[miss])

    # 1) Try session-only install
    if (isTRUE(allow_session_installs)) {
      say(paste0(
        "Active Python is initialized (possibly ephemeral). Missing modules: ",
        paste(miss, collapse = ", "),
        "\nTrying py_require() for this session…"
      ))

      req_ok <- TRUE
      tryCatch(
        { reticulate::py_require(pkgs) },
        error = function(e) { req_ok <<- FALSE }
      )

      if (!req_ok) {
        say("py_require() failed; attempting in-session pip install via py_run_string()…")
        pip_ok <- TRUE
        pip_cmd <- sprintf(
          "import sys, subprocess; subprocess.check_call([sys.executable,'-m','pip','install','-U',%s])",
          paste(sprintf("'%s'", pkgs), collapse = ",")
        )
        tryCatch(
          { reticulate::py_run_string(pip_cmd) },
          error = function(e) {
            pip_ok <<- FALSE
            warning("Session pip install failed: ", conditionMessage(e))
          }
        )
        if (pip_ok) miss <- missing_imports(py_imports)
      }
    }

    # 2) If still missing, prepare persistent conda env
    if (length(miss)) {
      say("Preparing a persistent conda environment with required packages…")

      if (!has_conda()) {
        say("Miniconda not found; installing Miniconda…")
        reticulate::install_miniconda()
        if (!has_conda()) {
          warning("Miniconda installation failed; cannot prepare persistent environment.")
          return(FALSE)
        }
      }

      if (isTRUE(force_conda_forge)) {
        suppressWarnings(conda_run(c("config", "--remove-key", "default_channels")))
        suppressWarnings(conda_run(c("config", "--remove-key", "channels")))
        conda_run(c("config", "--add", "channels", "conda-forge"))
        conda_run(c("config", "--set", "channel_priority", "strict"))
      }

      envs <- tryCatch(reticulate::conda_list()$name, error = function(e) character())
      if (!(envname %in% envs)) {
        say(sprintf("Creating conda env '%s' (Python %s)…", envname, py_ver))
        st <- conda_run(c("create", "--yes", "--name", envname, paste0("python=", py_ver)))
        if (!identical(st, 0L)) {
          warning(sprintf("Failed to create env '%s'.", envname))
          return(FALSE)
        }
      }

      # Pre-install all required packages into the new env
      need_pkgs <- unname(install_name[names(install_name)])
      ok <- TRUE
      tryCatch(
        reticulate::conda_install(envname = envname, packages = need_pkgs, channel = "conda-forge"),
        error = function(e) { ok <<- FALSE }
      )
      if (!ok) {
        say("conda solver failed; falling back to pip for required packages…")
        # Let reticulate handle pip into target env
        try(
          reticulate::py_install(need_pkgs, envname = envname, method = "pip", pip = TRUE),
          silent = TRUE
        )
      }

      # Try to auto-restart to bind the new env
      if (try_restart()) return(invisible(TRUE))

      # Otherwise, guide the user
      say(paste0(
        "A persistent environment '", envname, "' is now prepared with required packages.\n",
        "Please restart R and bind it with:\n",
        "  reticulate::use_condaenv('", envname, "', required = TRUE)\n",
        "Optionally, add to ~/.Renviron for persistence, e.g.:\n",
        "  RETICULATE_PYTHON=", normalizePath(
          file.path(Sys.getenv("LOCALAPPDATA"), "r-miniconda", "envs", envname, "python.exe"),
          winslash = "/"
        ),
        "\nThen re-run ICESat2VegR_configure()."
      ))
      return(FALSE)
    }
  }

  # ---------------------------------------------------------------------------
  # Verify critical imports
  # ---------------------------------------------------------------------------
  core <- c("packaging", "h5py", "earthaccess", "ee")
  miss_core <- missing_imports(core)
  if (length(miss_core)) {
    warning("Missing required Python modules after installation: ",
            paste(miss_core, collapse = ", "))
    return(FALSE)
  }

  # ---------------------------------------------------------------------------
  # Optionally resolve aiobotocore/botocore incompatibilities
  # ---------------------------------------------------------------------------
  if (isTRUE(manage_aiobotocore)) {
    chk <- pip_check()
    if (inherits(chk, "error") &&
        grepl("aiobotocore .* requires botocore", conditionMessage(chk))) {

      say("Resolving aiobotocore/botocore conflict by pinning boto3/botocore/s3transfer…")
      reticulate::py_install(
        c("botocore==1.40.49", "boto3==1.40.49", "s3transfer==0.14.0"),
        method = "pip",
        pip = TRUE
      )

      chk2 <- pip_check()
      if (inherits(chk2, "error") &&
          grepl("aiobotocore .* requires botocore", conditionMessage(chk2))) {
        say("Uninstalling aiobotocore (not required by ICESat2VegR)…")
        reticulate::py_run_string(
          "import sys, subprocess; subprocess.call([sys.executable,'-m','pip','uninstall','-y','aiobotocore'])"
        )
        pip_check()
      }
    }
  }

  print_summary(bound_env = if (!is_init) envname else NULL)
  say("Configuration complete.")
  TRUE
}

# ==============================================================================
# Earth Engine project normalization and initialization
# ==============================================================================

#' Normalize and validate a Google Cloud project ID for Earth Engine
#'
#' This is an internal helper used by \code{tryInitializeEarthEngine()}. It
#' accepts explicit project IDs, numeric project numbers, or pulls from the
#' environment variable \code{EE_PROJECT}. It:
#'
#' \itemize{
#'   \item Accepts accidental \code{"projects/<id>"} strings and strips the prefix.
#'   \item Allows numeric project numbers as valid.
#'   \item Rejects IDs with underscores (suggesting replacement with hyphens).
#'   \item Enforces basic format rules for project IDs.
#' }
#'
#' @param project Character project ID or number (may be \code{NULL} or empty).
#' @param quiet Logical. If \code{FALSE}, print informational messages.
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

  # Validate format: start letter; 6–30 chars; letters/digits/hyphens; not end with hyphen
  if (!grepl("^[a-z][a-z0-9-]{4,29}$", project) || grepl("-$", project)) {
    stop(
      "Invalid project ID '", project, "'. A valid ID must:\n",
      " - start with a lowercase letter\n",
      " - be 6–30 characters long\n",
      " - contain only lowercase letters, digits, and hyphens\n",
      " - not end with a hyphen"
    )
  }

  project
}

#' Lightweight Earth Engine initialization helper
#'
#' This function is a low-level helper to initialize the Python
#' \code{earthengine-api} client through \code{reticulate}. It supports both:
#'
#' \itemize{
#'   \item Service Account (SA) authentication via JSON key file.
#'   \item User OAuth authentication (browser flow).
#' }
#'
#' It validates and normalizes the Google Cloud project ID, initializes the
#' Earth Engine client, and sets \code{EE_PROJECT} in the environment for
#' downstream use (e.g., in Python code).
#'
#' @param project Character. GCP Project ID or numeric project number. If
#'   \code{NULL} or empty, the environment variable \code{EE_PROJECT} is used
#'   (and must be set).
#' @param service_account Optional service account email. If provided, a
#'   \code{keyfile} must also be supplied.
#' @param keyfile Optional path to the service account JSON key file (required
#'   if \code{service_account} is provided).
#' @param quiet Logical. Suppress messages if \code{TRUE}.
#' @param force_auth Logical. If \code{TRUE}, force a user OAuth flow prior to
#'   calling \code{ee$Initialize()}.
#'
#' @return Logical \code{TRUE} on success; otherwise an error is raised.
#' @keywords internal
#'
#' @examples
#' \dontrun{
#'   # Using user OAuth and EE_PROJECT from environment:
#'   tryInitializeEarthEngine()
#'
#'   # Using explicit project id and OAuth
#'   tryInitializeEarthEngine(project = "my-ee-project")
#'
#'   # Using a service account
#'   tryInitializeEarthEngine(
#'     project         = "my-ee-project",
#'     service_account = "my-sa@my-ee-project.iam.gserviceaccount.com",
#'     keyfile         = "path/to/key.json"
#'   )
#' }
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
    msg("Forcing OAuth flow…")
    try(ee$Authenticate(), silent = TRUE)
  }

  ok <- tryCatch(
    {
      ee$Initialize(project = project)
      TRUE
    },
    error = function(e) {
      msg("Initialize() failed (", conditionMessage(e), "). Attempting OAuth…")
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
#' then checks whether critical Python modules (e.g., \code{h5py},
#' \code{earthaccess}) have been correctly imported into the package-level
#' environment. It is intended primarily for diagnostics after calling
#' \code{ICESat2VegR_configure()}.
#'
#' The function is conservative and does not throw errors if the reload fails;
#' it simply emits warnings if critical Python bindings appear to be missing.
#'
#' @return Invisibly returns \code{TRUE}; warnings are issued if something
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
