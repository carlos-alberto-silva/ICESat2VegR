#' LAS Public Header Description
#'
#' @description
#' Internal function that returns the ASPRS LAS public header structure.
#'
#' @return A data.frame describing the LAS public header fields.
#'
#' @keywords internal
.publicHeaderDescription <- function() {

  hd <- structure(list(
    Item = c(
      "File Signature (\"LASF\")",
      "(1.1) File Source ID",
      "(1.1) Global Encoding",
      "(1.1) Project ID - GUID data 1",
      "(1.1) Project ID - GUID data 2",
      "(1.1) Project ID - GUID data 3",
      "(1.1) Project ID - GUID data 4",
      "Version Major",
      "Version Minor",
      "(1.1) System Identifier",
      "Generating Software",
      "(1.1) File Creation Day of Year",
      "(1.1) File Creation Year",
      "Header Size",
      "Offset to point data",
      "Number of variable length records",
      "Point Data Format ID (0-99 for spec)",
      "Point Data Record Length",
      "Number of point records",
      "Number of points by return",
      "X scale factor",
      "Y scale factor",
      "Z scale factor",
      "X offset",
      "Y offset",
      "Z offset",
      "Max X",
      "Min X",
      "Max Y",
      "Min Y",
      "Max Z",
      "Min Z"
    ),
    Format = c(
      "char[4]", "unsigned short", "unsigned short", "unsigned long",
      "unsigned short", "unsigned short", "unsigned char[8]",
      "unsigned char", "unsigned char", "char[32]", "char[32]",
      "unsigned short", "unsigned short", "unsigned short",
      "unsigned long", "unsigned long", "unsigned char",
      "unsigned short", "unsigned long", "unsigned long[5]",
      "double", "double", "double", "double", "double", "double",
      "double", "double", "double", "double", "double", "double"
    ),
    Size = c(
      "4 bytes", "2 bytes", "2 bytes", "4 bytes", "2 bytes",
      "2 bytes", "8 bytes", "1 byte", "1 byte", "32 bytes",
      "32 bytes", "2 bytes", "2 bytes", "2 bytes", "4 bytes",
      "4 bytes", "1 byte", "2 bytes", "4 bytes", "20 bytes",
      "8 bytes", "8 bytes", "8 bytes", "8 bytes", "8 bytes",
      "8 bytes", "8 bytes", "8 bytes", "8 bytes", "8 bytes",
      "8 bytes", "8 bytes"
    ),
    Required = rep("*", 32)
  ), row.names = 2:33, class = "data.frame")

  hd$what <- ""
  hd$what[grep("unsigned", hd$Format)] <- "integer"
  hd$what[grep("char", hd$Format)] <- "raw"
  hd$what[grep("short", hd$Format)] <- "integer"
  hd$what[grep("long", hd$Format)] <- "integer"
  hd$what[grep("double", hd$Format)] <- "numeric"

  hd$signed <- TRUE
  hd$signed[grep("unsigned", hd$Format)] <- FALSE

  hd$n <- as.numeric(gsub("[[:alpha:][:punct:]]", "", hd$Format))
  hd$n[is.na(hd$n)] <- 1

  hd$Hsize <- as.numeric(gsub("[[:alpha:]]", "", hd$Size))
  hd$Rsize <- hd$Hsize / hd$n

  hd$Rsize[hd$what == "raw"] <- 1
  hd$n[hd$what == "raw"] <- hd$Hsize[hd$what == "raw"]

  hd
}

#' Write unsigned 1-byte integer
#' @keywords internal
.write_u1 <- function(value, con) {
  writeBin(as.integer(value), con, size = 1, endian = "little")
}

#' Write unsigned 2-byte integer
#' @keywords internal
.write_u2 <- function(value, con) {
  writeBin(as.integer(value), con, size = 2, endian = "little")
}

#' Write unsigned 4-byte integer
#' @keywords internal
.write_u4 <- function(value, con) {
  writeBin(as.integer(value), con, size = 4, endian = "little")
}

#' Write signed 1-byte integer
#' @keywords internal
.write_i1 <- function(value, con) {
  writeBin(as.integer(value), con, size = 1, endian = "little")
}

#' Write signed 4-byte integer
#' @keywords internal
.write_i4 <- function(value, con) {
  writeBin(as.integer(value), con, size = 4, endian = "little")
}

#' Write 8-byte floating-point value
#' @keywords internal
.write_f8 <- function(value, con) {
  writeBin(as.numeric(value), con, size = 8, endian = "little")
}

#' Write fixed-length raw vector
#' @keywords internal
.write_raw_fixed <- function(value, con, n) {
  value <- as.raw(value)
  if (length(value) < n) value <- c(value, raw(n - length(value)))
  if (length(value) > n) value <- value[seq_len(n)]
  writeBin(value, con, size = 1, endian = "little")
}

#' Write fixed-length character string
#' @keywords internal
.write_char_fixed <- function(value, con, n) {
  value <- charToRaw(sprintf(paste0("%-", n, "s"), value))
  value <- value[seq_len(n)]
  writeBin(value, con, size = 1, endian = "little")
}

#' Write LiDAR Data to LAS Format
#'
#' @description
#' Writes a LiDAR point cloud stored as a data.table to a LAS file.
#'
#' @details
#' This function writes an ASPRS LAS version 1.2 file using Point Data Record
#' Format 0. The input must be a data.table containing at least the columns
#' X, Y, and Z.
#'
#' @param x A data.table with columns X, Y, and Z. Optional columns are
#' Intensity, ReturnNumber, NumberOfReturns, Classification, ScanAngleRank,
#' UserData, and PointSourceID.
#' @param LASfile Character. Output LAS file path.
#' @param scale Numeric vector of length 3. Scale factors for X, Y, and Z.
#' Default is c(0.001, 0.001, 0.001).
#'
#' @return Invisibly returns the output LAS file path.
#'
#' @references
#' American Society for Photogrammetry and Remote Sensing, ASPRS.
#' LAS Specification, Version 1.2.
#'
#' @importFrom data.table is.data.table
#'
#' @keywords internal
writeLAS <- function(x, LASfile, scale = c(0.001, 0.001, 0.001)) {

  if (!data.table::is.data.table(x)) {
    stop("x must be a data.table.")
  }

  if (!is.character(LASfile) || length(LASfile) != 1) {
    stop("LASfile must be a single output file path.")
  }

  if (!all(c("X", "Y", "Z") %in% names(x))) {
    stop("x must contain the columns X, Y, and Z.")
  }

  if (nrow(x) == 0) {
    stop("x has no points.")
  }

  if (!all(sapply(x[, .SD, .SDcols = c("X", "Y", "Z")], is.numeric))) {
    stop("X, Y, and Z must be numeric.")
  }

  if (!is.numeric(scale) || length(scale) != 3) {
    stop("scale must be a numeric vector of length 3.")
  }

  n <- nrow(x)

  get_col <- function(name, default) {
    if (name %in% names(x)) x[[name]] else rep(default, n)
  }

  Intensity <- as.integer(get_col("Intensity", 0))
  ReturnNumber <- as.integer(get_col("ReturnNumber", 1))
  NumberOfReturns <- as.integer(get_col("NumberOfReturns", 1))
  Classification <- as.integer(get_col("Classification", 1))
  ScanAngleRank <- as.integer(get_col("ScanAngleRank", 0))
  UserData <- as.integer(get_col("UserData", 0))
  PointSourceID <- as.integer(get_col("PointSourceID", 0))

  Intensity <- pmin(pmax(Intensity, 0), 65535)
  ReturnNumber <- pmin(pmax(ReturnNumber, 1), 5)
  NumberOfReturns <- pmin(pmax(NumberOfReturns, 1), 5)
  Classification <- pmin(pmax(Classification, 0), 255)
  ScanAngleRank <- pmin(pmax(ScanAngleRank, -128), 127)
  UserData <- pmin(pmax(UserData, 0), 255)
  PointSourceID <- pmin(pmax(PointSourceID, 0), 65535)

  offset <- c(min(x$X), min(x$Y), min(x$Z))

  Xi <- as.integer(round((x$X - offset[1]) / scale[1]))
  Yi <- as.integer(round((x$Y - offset[2]) / scale[2]))
  Zi <- as.integer(round((x$Z - offset[3]) / scale[3]))

  creation_date <- as.POSIXlt(Sys.Date())

  con <- file(LASfile, open = "wb")
  on.exit(close(con))

  writeBin(charToRaw("LASF"), con, size = 1)
  .write_u2(0, con)
  .write_u2(0, con)
  .write_u4(0, con)
  .write_u2(0, con)
  .write_u2(0, con)
  .write_raw_fixed(raw(8), con, 8)
  .write_u1(1, con)
  .write_u1(2, con)
  .write_char_fixed("ICESat2VegR", con, 32)
  .write_char_fixed("writeLAS", con, 32)
  .write_u2(creation_date$yday + 1L, con)
  .write_u2(creation_date$year + 1900L, con)
  .write_u2(227L, con)
  .write_u4(227L, con)
  .write_u4(0L, con)
  .write_u1(0L, con)
  .write_u2(20L, con)
  .write_u4(n, con)

  returns_count <- tabulate(ReturnNumber, nbins = 5)

  for (j in seq_len(5)) {
    .write_u4(returns_count[j], con)
  }

  .write_f8(scale[1], con)
  .write_f8(scale[2], con)
  .write_f8(scale[3], con)

  .write_f8(offset[1], con)
  .write_f8(offset[2], con)
  .write_f8(offset[3], con)

  .write_f8(max(x$X), con)
  .write_f8(min(x$X), con)
  .write_f8(max(x$Y), con)
  .write_f8(min(x$Y), con)
  .write_f8(max(x$Z), con)
  .write_f8(min(x$Z), con)

  for (i in seq_len(n)) {

    .write_i4(Xi[i], con)
    .write_i4(Yi[i], con)
    .write_i4(Zi[i], con)

    .write_u2(Intensity[i], con)

    return_byte <- bitwOr(
      bitwAnd(ReturnNumber[i], 7),
      bitwShiftL(bitwAnd(NumberOfReturns[i], 7), 3)
    )

    .write_u1(return_byte, con)
    .write_u1(Classification[i], con)
    .write_i1(ScanAngleRank[i], con)
    .write_u1(UserData[i], con)
    .write_u2(PointSourceID[i], con)
  }

  invisible(LASfile)
}

#' Convert ICESat-2 ATL03 and ATL08 data.table objects to LAS files
#'
#' @description
#' Converts ICESat-2 ATL03 or ATL08 point data stored as a data.table with
#' longitude, latitude, and height values into one or more LAS files.
#'
#' @details
#' The input data must contain at least three columns representing longitude,
#' latitude, and height. The function renames the first three columns to X, Y,
#' and Z internally. The original X and Y values are assumed to be longitude
#' and latitude in EPSG:4326.
#'
#' Because LAS files should store projected coordinates, the function identifies
#' the appropriate UTM zone for each point, splits the dataset by UTM zone,
#' projects each subset to the corresponding UTM coordinate system, and writes
#' one LAS file per UTM zone.
#'
#'
#' @param dt A data.table with at least three columns: longitude, latitude,
#' and height. The first three columns are interpreted as longitude, latitude,
#' and height.
#' @param output Character. Output LAS file path. If the data span multiple UTM
#' zones, the EPSG code is appended to each output file name.
#'
#' @return Invisibly returns a character vector with the written LAS file paths.
#'
#' @examples
#' \dontrun{
#' library(data.table)
#'
#' dt <- data.table(
#'   lon = runif(1000, -52.5, -52.4),
#'   lat = runif(1000, -15.9, -15.8),
#'   h = runif(1000, 0, 35)
#' )
#'
#' dt_to_las(dt, "icesat2_output.las")
#' }
#'
#' @importFrom data.table as.data.table copy setnames
#' @importFrom terra project
#'
#' @include utmTools.R
#'
#' @export
dt_to_las <- function(dt, output) {

  if (!requireNamespace("data.table", quietly = TRUE)) {
    stop("Package 'data.table' is required.")
  }

  if (!requireNamespace("terra", quietly = TRUE)) {
    stop("Package 'terra' is required.")
  }

  if (!data.table::is.data.table(dt)) {
    dt <- data.table::as.data.table(dt)
  } else {
    dt <- data.table::copy(dt)
  }

  if (ncol(dt) < 3) {
    stop("dt must contain at least three columns: longitude, latitude, and height.")
  }

  if (!is.character(output) || length(output) != 1) {
    stop("output must be a single output file path.")
  }

  if (!grepl("\\.las$", output, ignore.case = TRUE)) {
    stop("output must end with '.las'. This function writes LAS files, not LAZ files.")
  }

  data.table::setnames(dt, names(dt)[1:3], c("X", "Y", "Z"))

  if (!all(vapply(dt[, .SD, .SDcols = c("X", "Y", "Z")], is.numeric, logical(1)))) {
    stop("The first three columns of dt must be numeric.")
  }

  X <- Y <- Z <- mask <- epsg <- NULL

  maskZones <- latLongToUtmMask(dt$Y, dt$X)

  message("====================================================")
  message(sprintf(
    "The provided data will be split into %s UTM zone(s)",
    nrow(maskZones)
  ))
  message("====================================================")

  written_files <- character(nrow(maskZones))

  for (ii in seq_len(nrow(maskZones))) {

    dtLocal <- data.table::copy(dt[maskZones[ii, mask][[1]]])
    epsgCode <- maskZones[ii, epsg]
    epsgText <- sprintf("epsg:%s", epsgCode)

    coordinates <- terra::project(
      as.matrix(dtLocal[, list(X, Y)]),
      from = "epsg:4326",
      to = epsgText
    )

    dtLocal[, X := coordinates[, 1]]
    dtLocal[, Y := coordinates[, 2]]

    localOutput <- gsub(
      "(\\.las)$",
      sprintf("_%s\\1", epsgCode),
      output,
      ignore.case = TRUE
    )

    message(sprintf("EPSG: %s, saved as %s", epsgCode, localOutput))

    writeLAS(
      x = dtLocal,
      LASfile = localOutput,
      scale = c(0.01, 0.01, 0.01)
    )

    written_files[ii] <- localOutput
  }

  invisible(written_files)
}
