icesat2DownloadFile <- function(url, outdir, overwrite, buffer_size, netrc, timeout) {
  filename <- file.path(outdir, basename(url)) # Keep original filename
  if ((!overwrite) && file.exists(filename)) {
    message("Skipping this file, already downloaded!")
    return(0)
  } # SKip if already downloaded

  # Temporary to file to resume to
  resume <- paste0(filename, ".curltmp")
  if (file.exists(resume)) {
    resume_from <- file.info(resume)$size # Get current size to resume from
  } else {
    resume_from <- 0
  }

  # Connection config
  h <- curl::new_handle()
  curl::handle_setopt(h, netrc = as.integer(1), netrc_file = netrc, resume_from = resume_from, connecttimeout = timeout, followlocation = TRUE)

  tryCatch({
    fileHandle <- file(resume, open = "ab", raw = TRUE)
    message("Connecting...")
    conn <- tryCatch(curl::curl(url, handle = h, open = "rb"), error = function(e) {
      tryCatch(curl::curl(url, handle = h, open = "rb"), error = function(e) {
        stop(e)
      })
    })

    message("Connected successfully, downloading...")
    headers <- rawToChar(curl::handle_data(h)$headers)
    total_size <- as.numeric(gsub("[^\u00e7]*Content-Length: ([0-9]+)[^\u00e7]*", "\\1", x = headers, perl = TRUE))
    while (TRUE) {
      message(
        sprintf(
          "\rDownloading... %.2f/%.2fMB (%.2f%%)    ",
          resume_from / 1024.0 / 1024.0,
          total_size / 1024.0 / 1024.0,
          100.0 * resume_from / total_size
        ),
        appendLF = FALSE
      )
      data <- readBin(conn, what = raw(), n = 1024 * buffer_size)
      size <- length(data)
      if (size == 0) {
        break
      }
      writeBin(data, fileHandle, useBytes = TRUE)
      resume_from <- resume_from + size
    }
    message(sprintf(
      "\rDownloading... %.2f/%.2fMB (100%%)    ",
      total_size / 1024.0 / 1024.0,
      total_size / 1024.0 / 1024.0
    ))
    close(fileHandle)
    close(conn)
    file.rename(resume, filename)
    return(0)
  }, interrupt = function(e) {
    warning("\nDownload interrupted!!!")
    try(close(conn), silent = TRUE)
    try(close(fileHandle), silent = TRUE)
  }, finally = {
    try(close(conn), silent = TRUE)
    try(close(fileHandle), silent = TRUE)
  })
  return(-1)
}

getNetRC <- function(dl_dir) {
  writeLines(
    sprintf(
      "machine %s login %s password %s",
      "urs.earthdata.nasa.gov",
      user,
      password
    ),
    ".netrc"
  )
}
