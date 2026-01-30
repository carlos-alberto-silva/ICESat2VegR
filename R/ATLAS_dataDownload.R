#' Download ICESat-2 ATL03/ATL08 data
#'
#' @description
#' Download ICESat-2 ATL03 and ATL08 data from the LP DAAC / NSIDC endpoints.
#' Handles connection drops by resuming partial files and retrying with
#' exponential backoff. Authentication is handled via a `.netrc` file:
#' if it does not exist, the user is prompted for Earthdata Login credentials.
#' If authentication fails (e.g., wrong username/password), the user is
#' informed and can choose to re-enter credentials or cancel.
#'
#' @param url Character vector; URLs to ATL03/ATL08 files
#'   (e.g., returned by \code{ATLAS_dataFinder()}).
#' @param outdir Character; output directory (default: [tempdir()]).
#' @param overwrite Logical; overwrite existing files? Default \code{FALSE}.
#' @param buffer_size Integer; chunk size in KB per \code{readBin} call
#'   (default 512).
#' @param timeout Numeric; connection timeout (seconds) for establishing
#'   the connection (default 10).
#' @param retries Integer; maximum attempts per file (default 3).
#' @param backoff Numeric; exponential backoff base used as
#'   \code{sleep = backoff^(attempt - 1)} (default 2).
#' @param quiet Logical; suppress per-file messages (default \code{FALSE}).
#' @param use_home_netrc Logical; if \code{TRUE} (default), use a single
#'   \file{~/.netrc} file in the user's home directory for Earthdata
#'   credentials so they can be reused across projects. If \code{FALSE},
#'   a \file{.netrc} file is created in \code{outdir}.
#'
#' @return (invisibly) a \code{data.frame} with columns:
#'   \itemize{
#'     \item \code{file}: file name
#'     \item \code{dest}: full destination path
#'     \item \code{status}: one of \code{"ok"}, \code{"failed"}, \code{"skipped"}
#'     \item \code{attempts}: number of attempts used
#'     \item \code{error}: error message (if any)
#'   }
#'
#' @references
#' Credits to Cole Krehbiel. Code adapted from:
#' \url{https://git.earthdata.nasa.gov/projects/LPDUR/repos/daac_data_download_r/browse/DAACDataDownload.R}
#'
#' @examples
#' \dontrun{
#' urls <- c(
#'   "https://example.com/path/to/ATL03_001.h5",
#'   "https://example.com/path/to/ATL08_001.h5"
#' )
#' status <- ATLAS_dataDownload(urls,
#'                              outdir        = "data",
#'                              retries       = 5,
#'                              backoff       = 2,
#'                              use_home_netrc = TRUE)
#' subset(status, status == "failed")
#' }
#'
#' @import curl
#' @export
ATLAS_dataDownload <- function(url,
                               outdir        = NULL,
                               overwrite     = FALSE,
                               buffer_size   = 512,
                               timeout       = 10,
                               retries       = 3,
                               backoff       = 2,
                               quiet         = FALSE,
                               use_home_netrc = TRUE) {
  
  # -------- argument checks --------
  if (is.null(outdir)) outdir <- tempdir()
  if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
  
  if (!is.logical(overwrite) || length(overwrite) != 1L)
    stop("`overwrite` must be TRUE/FALSE")
  if (!is.numeric(buffer_size) || length(buffer_size) != 1L || buffer_size <= 0)
    stop("`buffer_size` must be a positive number (KB)")
  if (!is.numeric(timeout) || length(timeout) != 1L || timeout < 0)
    stop("`timeout` must be a non-negative number (seconds)")
  if (!is.numeric(retries) || length(retries) != 1L || retries < 1)
    stop("`retries` must be a positive integer")
  if (!is.numeric(backoff) || length(backoff) != 1L || backoff <= 0)
    stop("`backoff` must be a positive number")
  if (!is.logical(quiet) || length(quiet) != 1L)
    stop("`quiet` must be TRUE/FALSE")
  if (!is.logical(use_home_netrc) || length(use_home_netrc) != 1L)
    stop("`use_home_netrc` must be TRUE/FALSE")
  
  buffer_size <- as.integer(buffer_size)
  retries     <- as.integer(retries)
  
  # Ensure we have a usable .netrc (or interactively create one if needed)
  netrc <- getNetRC(outdir, use_home_netrc = use_home_netrc)
  
  files   <- as.character(url)
  n_files <- length(files)
  results <- vector("list", n_files)
  
  for (i in seq_len(n_files)) {
    u <- files[i]
    if (is.na(u) || !nzchar(u)) {
      results[[i]] <- data.frame(
        file     = NA_character_,
        dest     = NA_character_,
        status   = "skipped",
        attempts = 0L,
        error    = "empty URL",
        stringsAsFactors = FALSE
      )
      next
    }
    
    dest <- file.path(outdir, basename(u))
    
    if (!quiet) {
      message("------------------------------")
      message(sprintf("Downloading file %d/%d: %s", i, n_files, basename(u)))
      message("------------------------------")
    }
    
    if (file.exists(dest) && !overwrite) {
      if (!quiet) message("Exists & overwrite=FALSE -> skipped.")
      results[[i]] <- data.frame(
        file     = basename(u),
        dest     = dest,
        status   = "skipped",
        attempts = 0L,
        error    = NA_character_,
        stringsAsFactors = FALSE
      )
      next
    }
    
    attempt   <- 1L
    success   <- FALSE
    last_err  <- NULL
    
    while (attempt <= retries && !success) {
      rc <- try(
        icesat2DownloadFile(
          url          = u,
          outdir       = outdir,
          overwrite    = overwrite,
          buffer_size  = buffer_size,
          netrc        = netrc,
          timeout      = timeout,
          quiet        = quiet
        ),
        silent = TRUE
      )
      
      if (!inherits(rc, "try-error") && is.numeric(rc) && rc == 0) {
        success <- TRUE
        if (!quiet) message("Finished successfully!")
      } else {
        if (inherits(rc, "try-error")) {
          last_err <- as.character(rc)
        } else {
          last_err <- sprintf("download failed (code %s)", as.character(rc))
        }
        if (!quiet) message(sprintf("Attempt %d failed: %s", attempt, last_err))
        if (attempt < retries) {
          sleep_s <- backoff^(attempt - 1)
          if (!quiet) message(sprintf("Retrying in %.1f s ...", sleep_s))
          Sys.sleep(sleep_s)
        }
        attempt <- attempt + 1L
      }
    }
    
    if (success) {
      results[[i]] <- data.frame(
        file     = basename(u),
        dest     = dest,
        status   = "ok",
        attempts = attempt - 1L,
        error    = NA_character_,
        stringsAsFactors = FALSE
      )
    } else {
      results[[i]] <- data.frame(
        file     = basename(u),
        dest     = dest,
        status   = "failed",
        attempts = retries,
        error    = last_err %||% NA_character_,
        stringsAsFactors = FALSE
      )
      message("⚠ Failed: ", basename(u), if (!is.null(last_err)) paste0(" | ", last_err))
    }
  }
  
  out <- do.call(rbind, results)
  if (!quiet) {
    message(sprintf(
      "Summary: ok=%d, failed=%d, skipped=%d (out of %d)",
      sum(out$status == "ok", na.rm = TRUE),
      sum(out$status == "failed", na.rm = TRUE),
      sum(out$status == "skipped", na.rm = TRUE),
      n_files
    ))
  }
  invisible(out)
}

# -------------------------------------------------------------------------
# Internal helpers
# -------------------------------------------------------------------------

#' @keywords internal
`%||%` <- function(a, b) if (!is.null(a)) a else b

# Re-write or create the .netrc at a given path after failed auth.
# Ask user if they want to retry or cancel.
# Returns netrc_path (invisible) if credentials updated,
# or NA if user chose to cancel.
#' @keywords internal
.update_netrc_file <- function(netrc_path, quiet = FALSE) {
  message("\nAuthentication failed (wrong username or password).")
  ans <- readline("Do you want to re-enter your Earthdata credentials? [y/N]: ")
  if (!nzchar(ans) || !tolower(substr(ans, 1, 1)) %in% c("y")) {
    message("User chose not to re-enter credentials. Cancelling this download.")
    return(NA_character_)
  }
  
  if (!quiet) {
    message("Please re-enter your Earthdata credentials.")
  }
  usr <- readline("Earthdata username: ")
  pwd <- readline("Earthdata password: ")
  
  lines <- c(
    "machine urs.earthdata.nasa.gov",
    paste0("  login ", usr),
    paste0("  password ", pwd)
  )
  
  dir.create(dirname(netrc_path), recursive = TRUE, showWarnings = FALSE)
  writeLines(lines, netrc_path)
  
  suppressWarnings(try(Sys.chmod(netrc_path, mode = "600"), silent = TRUE))
  
  if (!quiet) message("Updated credentials saved to: ", netrc_path)
  invisible(netrc_path)
}

# Single-file downloader with resume + progress + re-auth on 401/403.
# Returns:
#   0  on success
#  -1  on error (connection/auth/user-cancel) so caller can record status.
#' @keywords internal
icesat2DownloadFile <- function(url, outdir, overwrite, buffer_size, netrc, timeout, quiet = FALSE) {
  
  filename <- file.path(outdir, basename(url))
  if ((!overwrite) && file.exists(filename)) {
    if (!quiet) message("Skipping this file, already downloaded!")
    return(0)
  }
  
  # Allow one retry after re-entering credentials if 401/403 is detected.
  max_auth_attempts <- 2L
  auth_attempt      <- 1L
  
  repeat {
    
    resume      <- paste0(filename, ".curltmp")
    resume_from <- if (file.exists(resume)) file.info(resume)$size %||% 0 else 0
    
    h <- curl::new_handle()
    curl::handle_setopt(
      h,
      netrc          = 1L,
      netrc_file     = netrc,
      resume_from    = as.numeric(resume_from),
      connecttimeout = as.numeric(timeout),
      followlocation = 1L
    )
    
    fileHandle <- conn <- NULL
    fileHandle <- file(resume, open = "ab", raw = TRUE)
    
    if (!quiet) message("Connecting...")
    
    # ---- OPEN CONNECTION (catch 401/403 here too) ----
    conn_try <- suppressWarnings(try(curl::curl(url, handle = h, open = "rb"), silent = TRUE))
    if (inherits(conn_try, "try-error")) {
      err_msg <- paste0(conn_try)
      status  <- curl::handle_data(h)$status_code %||% NA_integer_
      
      auth_like <- (!is.na(status) && status %in% c(401L, 403L)) ||
        grepl("401", err_msg, fixed = TRUE) ||
        grepl("403", err_msg, fixed = TRUE)
      
      # Always close file handle before leaving / retrying
      try(close(fileHandle), silent = TRUE); fileHandle <- NULL
      
      if (auth_like) {
        message("\nAuthentication appears to have failed when opening connection.")
        
        if (file.exists(resume)) {
          unlink(resume)
        }
        
        if (auth_attempt < max_auth_attempts) {
          new_netrc <- .update_netrc_file(netrc, quiet = quiet)
          if (is.na(new_netrc)) {
            # user chose to cancel; no error thrown, just return code
            return(-1)
          }
          netrc        <- new_netrc
          auth_attempt <- auth_attempt + 1L
          next  # retry entire download with new credentials
        } else {
          message("Authentication still failing after re-entering credentials. Skipping this file.")
          return(-1)
        }
      }
      
      # Not clearly auth-related: generic connection error
      if (!quiet) {
        message("\nConnection error while opening the connection:")
        message(err_msg)
      }
      return(-1)
    }
    
    # If we got here, opening worked
    conn <- conn_try
    if (!quiet) message("Connected successfully, downloading...")
    
    # total size (may be NA if server doesn't send it)
    headers_raw <- curl::handle_data(h)$headers
    total_size <- NA_real_
    if (length(headers_raw)) {
      headers <- try(curl::parse_headers_list(rawToChar(headers_raw)), silent = TRUE)
      if (!inherits(headers, "try-error") && !is.null(headers[["content-length"]])) {
        total_size <- suppressWarnings(as.numeric(headers[["content-length"]]))
      }
    }
    
    downloaded <- as.numeric(resume_from)
    
    repeat {
      if (!quiet) {
        pct <- if (is.finite(total_size)) (100 * downloaded / total_size) else NA_real_
        msg <- if (is.finite(total_size)) {
          sprintf("\rDownloading... %.2f/%.2fMB (%.2f%%)    ",
                  downloaded/1048576, total_size/1048576, pct)
        } else {
          sprintf("\rDownloading... %.2fMB (size unknown)    ", downloaded/1048576)
        }
        message(msg, appendLF = FALSE)
      }
      
      chunk <- suppressWarnings(try(readBin(conn, what = raw(), n = 1024L * buffer_size), silent = TRUE))
      if (inherits(chunk, "try-error")) {
        err_msg <- paste0(chunk)
        status  <- curl::handle_data(h)$status_code %||% NA_integer_
        
        auth_like <- (!is.na(status) && status %in% c(401L, 403L)) ||
          grepl("401", err_msg, fixed = TRUE) ||
          grepl("403", err_msg, fixed = TRUE)
        
        # Close both before leaving / retrying
        try(close(conn),       silent = TRUE); conn <- NULL
        try(close(fileHandle), silent = TRUE); fileHandle <- NULL
        
        if (auth_like) {
          message("\nAuthentication appears to have failed while reading data.")
          
          if (file.exists(resume)) {
            unlink(resume)
          }
          
          if (auth_attempt < max_auth_attempts) {
            new_netrc <- .update_netrc_file(netrc, quiet = quiet)
            if (is.na(new_netrc)) {
              return(-1)
            }
            netrc        <- new_netrc
            auth_attempt <- auth_attempt + 1L
            next  # retry whole download
          } else {
            message("Authentication still failing after re-entering credentials. Skipping this file.")
            return(-1)
          }
        }
        
        if (!quiet) {
          message("\nConnection error while reading:")
          message(err_msg)
        }
        return(-1)
      }
      
      nbytes <- length(chunk)
      if (nbytes == 0L) break  # EOF
      
      writeBin(chunk, fileHandle, useBytes = TRUE)
      downloaded <- downloaded + nbytes
    }
    
    # normal close
    try(close(conn),       silent = TRUE); conn <- NULL
    try(close(fileHandle), silent = TRUE); fileHandle <- NULL
    
    # Final status check
    status <- curl::handle_data(h)$status_code %||% NA_integer_
    
    if (!is.na(status) && status %in% c(401L, 403L)) {
      message("\nHTTP status ", status, " (authentication failed).")
      
      if (file.exists(resume)) {
        unlink(resume)
      }
      
      if (auth_attempt < max_auth_attempts) {
        new_netrc <- .update_netrc_file(netrc, quiet = quiet)
        if (is.na(new_netrc)) {
          return(-1)
        }
        netrc        <- new_netrc
        auth_attempt <- auth_attempt + 1L
        next
      } else {
        message("Authentication still failing after re-entering credentials. Skipping this file.")
        return(-1)
      }
    }
    
    if (!quiet) {
      if (is.finite(total_size)) {
        message(sprintf("\rDownloading... %.2f/%.2fMB (100%%)    ",
                        total_size/1048576, total_size/1048576))
      } else {
        message(sprintf("\rDownloading... %.2fMB (completed)    ", downloaded/1048576))
      }
    }
    
    if (file.exists(resume)) {
      file.rename(resume, filename)
    }
    
    return(0)
  }
}

#' Get or create a .netrc file for Earthdata
#'
#' @keywords internal
getNetRC <- function(dl_dir, use_home_netrc = TRUE) {
  
  if (use_home_netrc) {
    home_dir <- normalizePath("~", mustWork = TRUE)
    netrc    <- file.path(home_dir, ".netrc")
  } else {
    netrc    <- file.path(dl_dir, ".netrc")
  }
  
  if (!file.exists(netrc) ||
      all(!grepl("urs.earthdata.nasa.gov", readLines(netrc), fixed = TRUE))) {
    
    netrc_conn <- file(netrc)
    on.exit(close(netrc_conn), add = TRUE)
    
    user <- Sys.getenv("NASA_USER")
    if (user == "") {
      if (!requireNamespace("getPass", quietly = TRUE)) {
        stop("Please install 'getPass' or set NASA_USER/NASA_PASSWORD env vars.")
      }
      user <- getPass::getPass(
        "Enter NASA Earthdata Login Username \n (or create an account at urs.earthdata.nasa.gov):"
      )
    }
    
    password <- Sys.getenv("NASA_PASSWORD")
    if (password == "") {
      if (!requireNamespace("getPass", quietly = TRUE)) {
        stop("Please install 'getPass' or set NASA_USER/NASA_PASSWORD env vars.")
      }
      password <- getPass::getPass("Enter NASA Earthdata Login Password:")
    }
    
    writeLines(c(
      "machine urs.earthdata.nasa.gov",
      sprintf("  login %s", user),
      sprintf("  password %s", password)
    ), netrc_conn)
    
    message("A .netrc file with your Earthdata Login credentials was stored at: ", netrc)
  }
  
  netrc
}
