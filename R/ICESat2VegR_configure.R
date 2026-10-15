# ==============================================================================
# ICESat2VegR - Python Environment Configuration + Earth Engine Initialization
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
#'       \item Reuses the interpreter reticulate already selected when it has
#'             every required module.
#'       \item Declares the required modules with \code{reticulate::py_require()}.
#'             With reticulate >= 1.41 this provisions reticulate's own managed
#'             environment, so no conda and no Miniconda are involved.
#'       \item Falls back to a \code{pip} install into the active interpreter.
#'       \item Falls back to a persistent conda environment (or a virtualenv
#'             when \code{prefer_conda = FALSE}), installing Miniconda
#'             automatically when conda is missing or cannot create environments.
#'       \item Prints the command and the output of anything that fails, so a
#'             failure can be understood from this call alone.
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
#' @param py_ver Character. Python version to request when creating a new
#'   persistent (conda/virtualenv) environment (default \code{"3.10"}).
#' @param prefer_conda Logical. If \code{TRUE} (default), conda/Miniconda is
#'   used for a persistent environment. If \code{FALSE}, a \code{virtualenv} is
#'   used instead.
#' @param force_conda_forge Logical. If \code{TRUE} (default), conda is run with
#'   \code{--override-channels -c conda-forge}. The machine's global conda
#'   channel configuration is never modified.
#' @param allow_session_installs Logical. If \code{TRUE} (default), missing
#'   modules are installed into the Python that is already active first, through
#'   \code{reticulate::py_require()} and, if needed, a \code{pip} call from the
#'   active interpreter, before any persistent environment is created.
#' @param allow_managed_env Logical. If \code{TRUE} (default), the required
#'   modules are declared with \code{reticulate::py_require()}. With
#'   reticulate >= 1.41 this lets reticulate build and reuse its own managed
#'   environment, which needs neither conda nor Miniconda.
#' @param auto_install_miniconda Logical. If \code{TRUE} (default) and no usable
#'   conda is found, Miniconda is installed for the current user (about a 100 MB
#'   download). Set to \code{FALSE} to forbid that download.
#' @param manage_aiobotocore Logical. If \code{TRUE} (default), attempt to
#'   detect and resolve common \code{aiobotocore}/\code{botocore}/\code{boto3}
#'   version conflicts by pinning compatible versions and optionally
#'   uninstalling \code{aiobotocore} (which is not required for ICESat2VegR).
#' @param auto_restart Logical. If \code{TRUE} (default), attempts to restart
#'   the R session in RStudio via \code{rstudioapi::restartSession()} after
#'   preparing a persistent environment, so that the new Python binding is
#'   cleanly established.
#' @param verbose Logical. If \code{TRUE} (default), print progress messages.
#' @param quick_probe Character vector of Python import names to test for a
#'   fast early exit (default: \code{c("h5py", "earthaccess", "ee")}).
#'
#' @return Logical \code{TRUE} when every required module is available in the
#'   current R session, \code{FALSE} otherwise. Failures print the command that
#'   failed together with its output.
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
    allow_managed_env = TRUE,
    auto_install_miniconda = TRUE,
    manage_aiobotocore = TRUE,
    auto_restart = TRUE,
    verbose = TRUE,
    quick_probe = c("h5py", "earthaccess", "ee")
) {

  # ---------------------------------------------------------------------------
  # Messaging. Anything that fails is printed together with the command that
  # produced it, so a failure can be understood from this call alone.
  # ---------------------------------------------------------------------------
  say  <- function(...) if (isTRUE(verbose)) message("[ICESat2VegR] ", ...)
  say0 <- function(...) message("[ICESat2VegR] ", ...)

  print_command_output <- function(out) {
    out <- as.character(out)
    if (!length(out)) return(invisible(FALSE))
    if (length(out) > 40L) {
      out <- c(utils::head(out, 20L), "... (output truncated) ...", utils::tail(out, 20L))
    }
    say0("--- output of the command above ---")
    for (line in out) message("    ", line)
    say0("-----------------------------------")
    invisible(TRUE)
  }

  # Run an external command and keep what it printed instead of discarding it.
  run_cmd <- function(command, args) {
    say("Running: ", command, " ", paste(args, collapse = " "))
    out <- tryCatch(
      suppressWarnings(system2(command, args, stdout = TRUE, stderr = TRUE)),
      error = function(e) structure(conditionMessage(e), status = 126L)
    )
    status <- attr(out, "status")
    if (is.null(status)) status <- 0L
    status <- as.integer(status)
    if (status != 0L) {
      say0("Command failed with exit code ", status, ".")
      print_command_output(out)
    }
    invisible(status)
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

  # An environment was requested by name: honour that request.
  named_env <- !identical(envname, "icesat2-env")

  # Is a Python interpreter already running in this session?
  py_started <- function() {
    isTRUE(tryCatch(reticulate::py_available(initialize = FALSE), error = function(e) FALSE))
  }

  # Does the installed reticulate provide `fn`? (py_require() needs >= 1.41)
  ret_has <- function(fn) {
    is.function(get0(fn, envir = asNamespace("reticulate"), inherits = FALSE))
  }

  python_info <- function() {
    cfg <- tryCatch(reticulate::py_config(), error = function(e) NULL)
    if (is.null(cfg) || is.null(cfg$python)) return("no Python interpreter available")
    sprintf("%s (%s)", cfg$python, cfg$version)
  }

  # Which required modules can be imported right now? Never throws.
  inspect_python <- function() {
    up <- isTRUE(tryCatch(reticulate::py_available(initialize = TRUE), error = function(e) FALSE))
    have <- rep(FALSE, length(py_imports))
    names(have) <- py_imports
    if (up) {
      have <- vapply(
        py_imports,
        function(m) isTRUE(tryCatch(reticulate::py_module_available(m), error = function(e) FALSE)),
        logical(1)
      )
      have <- as.logical(have)
      names(have) <- py_imports
    }
    list(have = have, missing = py_imports[!have])
  }

  # Declare the requirements, so that reticulate installs them while it starts
  # its own (uv managed) environment, before any Python is running.
  declare_requirements <- function(missing = py_imports) {
    if (!isTRUE(allow_managed_env) || !ret_has("py_require") || !length(missing)) return(FALSE)
    say("Declaring Python requirements: ", paste(unname(install_name[missing]), collapse = ", "))
    tryCatch(
      reticulate::py_require(unname(install_name[missing])),
      error   = function(e) say("py_require() reported: ", conditionMessage(e)),
      warning = function(w) say("py_require() reported: ", conditionMessage(w))
    )
    invisible(TRUE)
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

  # Attempt to restart R session in RStudio
  try_restart <- function() {
    if (!isTRUE(auto_restart)) return(invisible(FALSE))
    if (requireNamespace("rstudioapi", quietly = TRUE)) {
      say("Restarting the R session to bind the newly prepared environment...")
      try(rstudioapi::restartSession(), silent = TRUE)
      return(invisible(TRUE))
    }
    FALSE
  }


  # ---------------------------------------------------------------------------
  # Session level repair: needs neither conda nor Miniconda.
  # ---------------------------------------------------------------------------
  session_install <- function(missing) {
    if (!length(missing) || !isTRUE(allow_session_installs)) return(FALSE)

    # 1) reticulate's own environment (reticulate >= 1.41, backed by uv).
    declare_requirements(missing)
    missing <- inspect_python()$missing
    if (!length(missing)) {
      say("The running Python now provides every required module.")
      return(TRUE)
    }

    # 2) pip into the interpreter that is active right now.
    say("Trying pip into the active interpreter: ",
        paste(unname(install_name[missing]), collapse = ", "))
    pip_cmd <- sprintf(
      "import sys, subprocess; subprocess.check_call([sys.executable,'-m','pip','install','-U',%s])",
      paste(sprintf("'%s'", unname(install_name[missing])), collapse = ",")
    )
    tryCatch(
      reticulate::py_run_string(pip_cmd),
      error = function(e) say0("pip reported: ", conditionMessage(e))
    )
    !length(inspect_python()$missing)
  }

  # ---------------------------------------------------------------------------
  # Persistent environment (conda, or virtualenv when prefer_conda = FALSE)
  # ---------------------------------------------------------------------------
  channel_args <- function() {
    if (isTRUE(force_conda_forge)) c("--override-channels", "-c", "conda-forge") else character()
  }

  conda_bin <- function() {
    tryCatch(reticulate::conda_binary(), error = function(e) NULL)
  }

  # `--name` for a named environment, `--prefix` for a path.
  conda_target <- function(env) {
    if (grepl("[/\\\\]", env)) c("--prefix", env) else c("--name", env)
  }

  reticulate_miniconda_root <- function() {
    tryCatch(reticulate::miniconda_path(), error = function(e) NULL)
  }

  is_reticulate_conda <- function(cb) {
    root <- reticulate_miniconda_root()
    if (is.null(root)) return(FALSE)
    norm <- function(p) tryCatch(
      normalizePath(p, winslash = "/", mustWork = FALSE),
      error = function(e) p
    )
    startsWith(norm(cb), norm(root))
  }

  # A prefix inside the user's data directory, which is always writable.
  private_env_prefix <- function() {
    root <- tryCatch(
      tools::R_user_dir("ICESat2VegR", which = "data"),
      error = function(e) file.path(tempdir(), "ICESat2VegR")
    )
    file.path(root, "envs", gsub("[^A-Za-z0-9._-]", "_", envname))
  }

  # Make sure a conda is available, installing Miniconda when allowed.
  ensure_conda <- function() {
    if (!is.null(conda_bin())) return(TRUE)
    if (!isTRUE(auto_install_miniconda)) {
      say0("No conda found and auto_install_miniconda = FALSE.")
      return(FALSE)
    }
    say0("No conda found; installing Miniconda for the current user ",
         "(about 100 MB, one time)...")
    ok <- tryCatch(
      { reticulate::install_miniconda(update = FALSE); TRUE },
      error = function(e) { say0("install_miniconda() reported: ", conditionMessage(e)); FALSE }
    )
    if (!ok || is.null(conda_bin())) {
      root <- reticulate_miniconda_root()
      say0("Miniconda was not installed (expected under ",
           if (is.null(root)) "the reticulate cache" else root, ").")
      return(FALSE)
    }
    TRUE
  }

  # Create the conda environment, trying progressively harder.
  conda_prepare_env <- function() {
    cb <- conda_bin()
    if (is.null(cb)) return(NULL)

    envs <- tryCatch(as.character(reticulate::conda_list()$name), error = function(e) character())
    if (envname %in% envs) {
      say("Reusing the existing conda environment '", envname, "'.")
      return(envname)
    }

    python_arg <- paste0("python=", py_ver)

    say("Creating conda environment '", envname, "' (Python ", py_ver, ")...")
    if (run_cmd(cb, c("create", "--yes", "--name", envname, python_arg, channel_args())) == 0L) {
      return(envname)
    }

    # A named environment fails when the environments directory is not
    # writable, or when an earlier attempt left a directory behind. A private
    # prefix avoids both.
    prefix <- private_env_prefix()
    say("Retrying with a private prefix: ", prefix)
    dir.create(dirname(prefix), recursive = TRUE, showWarnings = FALSE)
    if (run_cmd(cb, c("create", "--yes", "--prefix", prefix, python_arg, channel_args())) == 0L) {
      envname <<- prefix
      return(prefix)
    }

    # The conda that was found may itself be unusable: reticulate can install
    # its own Miniconda next to the R library.
    if (isTRUE(auto_install_miniconda) && !is_reticulate_conda(cb)) {
      say0("The conda at '", cb, "' cannot create environments; trying reticulate's Miniconda...")
      installed <- tryCatch(
        { reticulate::install_miniconda(update = FALSE); TRUE },
        error = function(e) { say0("install_miniconda() reported: ", conditionMessage(e)); FALSE }
      )
      cb2 <- if (installed) conda_bin() else NULL
      if (!is.null(cb2) && !identical(cb2, cb) &&
          run_cmd(cb2, c("create", "--yes", "--name", envname, python_arg, channel_args())) == 0L) {
        return(envname)
      }
    }

    NULL
  }

  venv_prepare <- function(env) {
    if (reticulate::virtualenv_exists(env)) {
      say("Reusing the existing virtualenv '", env, "'.")
      return(TRUE)
    }
    say("Creating virtualenv '", env, "'...")
    tryCatch(
      { reticulate::virtualenv_create(env); TRUE },
      error = function(e) { say0("virtualenv_create() reported: ", conditionMessage(e)); FALSE }
    )
  }

  # Returns list(env = <identifier or character()>, bound = <TRUE/FALSE>).
  persistent_install <- function(missing) {
    if (!length(missing)) return(list(env = character(), bound = FALSE))

    target <- envname

    if (isTRUE(prefer_conda)) {
      if (!ensure_conda()) return(list(env = character(), bound = FALSE))
      cb <- conda_bin()
      target <- conda_prepare_env()
      if (is.null(target) || is.null(cb)) return(list(env = character(), bound = FALSE))
      say("Installing required modules into '", target, "' with conda...")
      status <- run_cmd(conda_bin(), c("install", "--yes", conda_target(target),
                                       unname(install_name[missing]), channel_args()))
      if (status != 0L) {
        say("conda install failed; falling back to pip inside '", target, "'...")
        try(
          reticulate::py_install(unname(install_name[missing]), envname = target,
                                 method = "pip", pip = TRUE),
          silent = TRUE
        )
      }
    } else {
      if (!venv_prepare(target)) return(list(env = character(), bound = FALSE))
      tryCatch(
        reticulate::virtualenv_install(target, packages = unname(install_name[missing]),
                                       ignore_installed = FALSE),
        error = function(e) say0("virtualenv_install() reported: ", conditionMessage(e))
      )
    }

    # A running interpreter cannot be replaced mid session: in that case the
    # environment becomes effective in the next session.
    if (py_started()) return(list(env = target, bound = FALSE))

    bound <- tryCatch({
      if (isTRUE(prefer_conda)) reticulate::use_condaenv(target, required = TRUE)
      else reticulate::use_virtualenv(target, required = TRUE)
      TRUE
    }, error = function(e) { say0("Could not bind '", target, "': ", conditionMessage(e)); FALSE })

    list(env = target, bound = isTRUE(bound) && !length(inspect_python()$missing))
  }

  # ---------------------------------------------------------------------------
  # MAIN FLOW
  # ---------------------------------------------------------------------------
  if (!requireNamespace("reticulate", quietly = TRUE)) {
    say0("The 'reticulate' package is required: install.packages('reticulate').")
    return(FALSE)
  }

  say("Checking the Python environment (requested envname: '", envname, "')...")

  # 0) Reuse an environment that already exists, before starting anything.
  bound_early <- FALSE
  if (!py_started()) {
    if (isTRUE(prefer_conda)) {
      envs <- tryCatch(as.character(reticulate::conda_list()$name), error = function(e) character())
      if (envname %in% envs) {
        say("Binding the existing conda environment '", envname, "'...")
        bound_early <- isTRUE(tryCatch(
          { reticulate::use_condaenv(envname, required = TRUE); TRUE },
          error = function(e) { say("Could not bind it: ", conditionMessage(e)); FALSE }
        ))
      }
    } else if (reticulate::virtualenv_exists(envname)) {
      say("Binding the existing virtualenv '", envname, "'...")
      bound_early <- isTRUE(tryCatch(
        { reticulate::use_virtualenv(envname, required = TRUE); TRUE },
        error = function(e) { say("Could not bind it: ", conditionMessage(e)); FALSE }
      ))
    }
  }

  # 1) Start from what this session already has. When Python is still down and
  #    no specific environment was requested, the requirements are declared
  #    first, so reticulate provisions them while starting its own environment
  #    (reticulate >= 1.41; no conda involved).
  if (!py_started() && !named_env && !bound_early) declare_requirements()

  state <- inspect_python()
  say("Active Python: ", python_info())
  say("Modules: ", paste(sprintf("%s=%s", py_imports, ifelse(state$have, "ok", "missing")),
                         collapse = "  "))

  if (!length(state$missing)) {
    say("Nothing to do: every required module is already available.")
    return(TRUE)
  }

  if (!named_env && length(quick_probe) && all(state$have[quick_probe])) {
    say("Core modules (", paste(quick_probe, collapse = ", "),
        ") are available; skipping installation.")
    return(TRUE)
  }

  # 2) Repair the Python that is already active. This needs no installation
  #    infrastructure at all.
  if (!named_env) {
    session_install(state$missing)
    state <- inspect_python()
  }

  # 3) Build a persistent environment when the session route was not enough.
  prepared <- list(env = character(), bound = FALSE)
  if (length(state$missing)) {
    say("Preparing a persistent environment for the missing modules...")
    prepared <- persistent_install(state$missing)
    if (length(prepared$env)) {
      say("Environment prepared: ", prepared$env,
          if (isTRUE(prepared$bound)) " (active in this session)" else " (for the next session)")
    }
    state <- inspect_python()
  }

  # 4) A named environment cannot replace a running interpreter: make sure the
  #    current session works anyway.
  if (length(state$missing) && named_env && isTRUE(allow_session_installs)) {
    session_install(state$missing)
    state <- inspect_python()
  }

  # 5) A prepared environment cannot replace an interpreter that is already
  #    running: it becomes active in the next session.
  if (length(prepared$env) && !isTRUE(prepared$bound)) {
    say0("The environment '", prepared$env, "' is ready. This R session keeps using ",
         "the previous Python: restart R (or set RETICULATE_PYTHON) to bind it.")
    if (try_restart()) return(invisible(TRUE))
  }

  # 6) Nothing left to try: report exactly what failed and how to finish by hand.
  if (length(state$missing)) {
    say0(
      "\nICESat2VegR could not provide every required Python module automatically.",
      "\n  Python:  ", python_info(),
      "\n  Missing: ", paste(state$missing, collapse = ", "),
      "\n\nSteps that reproduce what this function tried:",
      "\n  * reticulate's own environment, no conda needed:",
      "\n      reticulate::py_require(c(",
      paste(sprintf("'%s'", unname(install_name)), collapse = ", "), "))",
      "\n  * conda environment, with the conda errors visible:",
      "\n      reticulate::conda_create('", envname, "', python_version = '", py_ver,
      "', channel = 'conda-forge')",
      "\n      reticulate::conda_install('", envname, "', channel = 'conda-forge', packages = c(",
      paste(sprintf("'%s'", unname(install_name)), collapse = ", "), "))",
      "\n  * when there is no conda at all: reticulate::install_miniconda()",
      "\n  * or point R at a Python you prepared yourself, in .Renviron:",
      "\n      RETICULATE_PYTHON=<path to python executable>"
    )
    return(FALSE)
  }

  # ---------------------------------------------------------------------------
  # Optionally resolve aiobotocore/botocore incompatibilities
  # ---------------------------------------------------------------------------
  if (isTRUE(manage_aiobotocore)) {
    chk <- pip_check()
    if (inherits(chk, "error") &&
        grepl("aiobotocore .* requires botocore", conditionMessage(chk))) {

      say("Resolving aiobotocore/botocore conflict by pinning boto3/botocore/s3transfer...")
      reticulate::py_install(
        c("botocore==1.40.49", "boto3==1.40.49", "s3transfer==0.14.0"),
        method = "pip",
        pip = TRUE
      )

      chk2 <- pip_check()
      if (inherits(chk2, "error") &&
          grepl("aiobotocore .* requires botocore", conditionMessage(chk2))) {
        say("Uninstalling aiobotocore (not required by ICESat2VegR)...")
        reticulate::py_run_string(
          "import sys, subprocess; subprocess.call([sys.executable,'-m','pip','uninstall','-y','aiobotocore'])"
        )
        pip_check()
      }
    }
  }

  # Every module is present: state it plainly, with the interpreter in use.
  say0("Python:  ", python_info())
  say0("Modules: ", paste(sprintf("%s=%s", py_imports, ifelse(state$have, "ok", "MISSING")),
                          collapse = "  "))
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
