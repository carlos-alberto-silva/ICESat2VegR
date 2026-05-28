#' Given a stack image raster from GEE
#' retrieve the point geometry with values for the images
#'
#' @param stack A single image or a vector/list of images from Earth Engine.
#' @param geom A geometry from [`terra::SpatVector-class`] read with [`terra::vect`].
#' @param scale The scale in meters for the extraction (image resolution).
#' @param chunk_size If the number of observations is greater than 1000,
#' it is recommended to chunk the results for not running out of memory within GEE server,
#' default is chunk by 1000.
#'
#' @return A [data.table::data.table-class] with the properties extracted
#' from the ee images.
#'
#' @export
seg_ancillary_extract <- function(stack, geom, scale = 10, chunk_size = 1000) {
  n <- nrow(geom)
  dts <- list()
  k <- 0L

  for (ii in seq.int(1L, n, by = chunk_size)) {
    tail <- min(ii + chunk_size - 1L, n)
    message(sprintf("Processing %d-%d of %d", ii, tail, n))

    sampled <- extract(stack, geom[ii:tail], scale)

    sz <- tryCatch(sampled$size()$getInfo(), error = function(e) 0L)
    if (is.null(sz) || sz == 0L) next

    k <- k + 1L
    dts[[k]] <- ee_to_dt(sampled)
  }

  if (length(dts) == 0L) return(data.table::data.table())
  data.table::rbindlist(dts, use.names = TRUE, fill = TRUE)
}


#' Given a geometry with point samples and images from Earth Engine
#' retrieve the point geometry with values for the images
#'
#' @param stack A single image or a vector/list of images from Earth Engine.
#' @param geom A geometry from [`terra::SpatVector-class`] read with [`terra::vect()`].
#' @param scale The scale in meters for the extraction (image resolution).
#'
#' @return An ee.FeatureCollection with the properties extracted from the
#' stack of images from ee.
#'
#' @export
extract <- function(stack, geom, scale) {
  img <- .compose_ee_image(stack)    # ee.Image
  fc  <- .points_to_ee_fc(geom)      # ee.FeatureCollection of points

  img$sampleRegions(
    collection = fc,
    scale      = as.integer(scale),
    tileScale  = 16
  )
}
# extract <- function(stack, geom, scale) {
#   final <- c(stack)
#   tempjson <- tempfile(fileext = ".geojson")
#   terra::writeVector(geom, tempjson, filetype = "geojson")
#   parsed <- jsonlite::parse_json(readLines(tempjson))
#   geojson <- ee$FeatureCollection(parsed)

#   sampled <- final$sampleRegions(
#     collection = geojson,
#     scale = scale,
#     tileScale = 16
#   )

#   return(sampled)
# }



# Compose a single ee.Image from common inputs
.compose_ee_image <- function(stack) {
  ee <- reticulate::import("ee", delay_load = FALSE)

  # Already a Python EE object?
  if (inherits(stack, "python.builtin.object")) {
    cls <- try(stack$`__class__`$`__name__`, silent = TRUE)
    if (!inherits(cls, "try-error") && grepl("ImageCollection", cls, ignore.case = TRUE)) {
      return(stack$median())  # reduce collection (change reducer if needed)
    }
    return(stack)  # assume ee.Image
  }

  # List of ee.Image objects
  if (is.list(stack)) {
    imgs <- Filter(function(z) inherits(z, "python.builtin.object"), stack)
    if (!length(imgs)) stop("No ee.Image objects supplied in 'stack'.")
    img <- imgs[[1L]]
    if (length(imgs) > 1L) {
      for (j in 2:length(imgs)) img <- img$addBands(imgs[[j]])
    }
    return(img)
  }

  stop("Unsupported 'stack' type. Provide an ee.Image, an ImageCollection, or a list of ee.Image objects.")
}

# Turn sf / SpatVector points into ee.FeatureCollection
.points_to_ee_fc <- function(geom) {
  ee <- reticulate::import("ee", delay_load = FALSE)

  if (inherits(geom, "sf")) {
    tf <- tempfile(fileext = ".geojson")
    suppressMessages(sf::st_write(geom, tf, quiet = TRUE))
    return(ee$FeatureCollection(jsonlite::read_json(tf, simplifyVector = FALSE)))
  }
  if (inherits(geom, "SpatVector")) {
    tf <- tempfile(fileext = ".geojson")
    terra::writeVector(geom, tf, filetype = "geojson")
    return(ee$FeatureCollection(jsonlite::read_json(tf, simplifyVector = FALSE)))
  }
  stop("geom must be sf or terra::SpatVector (POINT geometries) for sampling.")
}
