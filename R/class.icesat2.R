#' @include class_tools.R
.datatable.aware <- TRUE

# Base class for icesat H5 files
setClass(
  Class = "icesat2.granules_cloud",
  slots = list(granules = "ANY")
)

setRefClass("icesat2.granule_cloud")

#' Subset Granules
#'
#' This method subsets the `granules` slot of an
#' `icesat2.granules_cloud` object.
#'
#' @param x An object of class `icesat2.granules_cloud`.
#' @param i Subset index.
#' @param j Unused, just to match the generic signature.
#' @param drop Unused, just to match the generic signature.
#' @param ... Additional arguments (not used).
#' @return An object of class `icesat2.granule_cloud`.
#' @export
setMethod(
  "[",
  signature = c("icesat2.granules_cloud"),
  definition = function(x, i = NULL, ...) {
    tryCatch(
      {
        res <- x@granules[i - 1]
      },
      error = function(e) {
        res <<- x@granules[i]
      },
      finally = {
        prepend_class(res, "icesat2.granule_cloud")
        return(res)
      }
    )
  }
)


#' Extract Granules
#'
#' This method extracts a single element from the `granules`
#' slot of an `icesat2.granules_cloud` object.
#'
#' @param x An object of class `icesat2.granules_cloud`.
#' @param i Extraction index.
#' @param j Unused, just to match the generic signature.
#' @param ... Additional arguments (not used).
#'
#' @return An object of class `icesat2.granule_cloud`.
#'
#' @examples
#' \dontrun{
#' granule <- new("icesat2.granules_cloud")
#' extracted_granule <- granule[[1]]
#' }
#' @export
setMethod(
  "[[",
  signature = c("icesat2.granules_cloud", "ANY", "ANY"),
  definition = function(x, i = NULL, ...) {
    x[i]
  }
)


#' @importFrom hdf5r H5File
setRefClass("icesat2.hdf5r")

#' Base class for all ICESat2VegR package's H5 files for generic functions
#' that can be run on any H5
#' @export
setClass(
  Class = "icesat2.h5",
  slots = list(h5 = "ANY")
)
icesat2.h5 <- new("icesat2.h5")

#' Class for ICESat-2 ATL03
#'
#' @slot h5 Object of class [`H5File`][hdf5r::H5File-class] from
#' `hdf5r` package containing the ICESat-2 Global Geolocated
#' Photon Data (ATL03)
#'
#' @seealso [`H5File`][hdf5r::H5File-class] in the `hdf5r` package and
#' \url{https://icesat-2.gsfc.nasa.gov/sites/default/files/page_files/ICESat2_ATL03_ATBD_r006.pdf}
#'
#' @import methods
#' @export
icesat2.atl03_h5 <- setClass(
  Class = "icesat2.atl03_h5",
  contains = "icesat2.h5"
)

#' Class for ICESat-2 ATL08
#'
#' @slot h5 Object of class [`H5File`][hdf5r::H5File-class] from
#' `hdf5r` package containing the
#' ICESat-2 Land and Vegetation Along-Track Products (ATL08)

#'
#' @seealso [`H5File`][hdf5r::H5File-class] in the `hdf5r` package and
#' \url{https://icesat-2.gsfc.nasa.gov/sites/default/files/page_files/ICESat2_ATL08_ATBD_r006.pdf}
#'
#' @import methods
#' @export
icesat2.atl08_h5 <- setClass(
  Class = "icesat2.atl08_h5",
  contains = "icesat2.h5"
)

setRefClass("icesat2.h5_cloud")

#' Class for ATL08 attributes
#'
#' @seealso [`data.table`][data.table::data.table-class] in the `data.table` package and
#' \url{https://icesat-2.gsfc.nasa.gov/sites/default/files/page_files/ICESat2_ATL08_ATBD_r006.pdf}
#'
#' @import methods
#' @export
setRefClass("icesat2.atl08_dt")

#' Class for ATL03 attributes
#'
#' @seealso [`data.table`][data.table::data.table-class] in the `data.table` package and
#' \url{https://icesat-2.gsfc.nasa.gov/sites/default/files/page_files/ICESat2_ATL03_ATBD_r006.pdf}
#'
#' @import methods
#' @importClassesFrom data.table data.table
#' @export
setRefClass("icesat2.atl03_dt")

#' Class for ATL03 segment attributes
#'
#' @seealso [`data.table`][data.table::data.table-class] in the `data.table` package and
#' \url{https://icesat-2.gsfc.nasa.gov/sites/default/files/page_files/ICESat2_ATL03_ATBD_r006.pdf}
#'
#' @import methods
#' @importClassesFrom data.table data.table
#' @export
setRefClass("icesat2.atl03_seg_dt")

#' Class for joined ATL03 and ATL08 attributes
#'
#' @seealso [`data.table`][data.table::data.table-class] in the `data.table` package and
#' \url{https://icesat-2.gsfc.nasa.gov/sites/default/files/page_files/ICESat2_ATL03_ATBD_r006.pdf}
#'
#' @import methods
#' @export
setRefClass("icesat2.atl03atl08_dt")


h5closeall <- function(con, ...) {
  con$close_all()
}


#' Safely closes the ICESat2.h5 base classes
#'
#' @description
#' Closing files will avoid locking HDF5 ATL03 files.
#'
#' @param con An object of class `ICESat2.h5`
#' @param ... Inherited from base
#'
#' @rdname close
#' @export
setMethod("close", signature = c("icesat2.h5"), h5closeall)


#' Plot photons from ATL03 and ATL08 joined products
#'
#' @description This function plots ATL03 and ATL08 joined data along track
#'
#' @param x An object of class [`ICESat2VegR::icesat2.atl03atl08_dt-class`]
#' @param y should be missing in this case, defaults to "ph_h"
#' @param beam Character vector indicating only one beam to process
#' ("gt1l", "gt1r", "gt2l", "gt2r", "gt3l", "gt3r").
#' Default is "gt1r"
#' @param ... will be passed to the main plot, such as:
#'
#'  - `colors`: A vector containing colors for plotting noise, terrain,
#' vegetation and top canopy photons (e.g.
#' c("gray", "#bd8421", "forestgreen", "green"))
#'  - `legend`: the position of the legend. 'bottomleft',
#' 'bottomright', 'topleft', 'topright' or FALSE to omit
#'
#' @return No return value
#'
#' @examples
#' # Specifying the path to ATL03 file
#' atl03_path <- system.file("extdata",
#'   "atl03_clip.h5",
#'   package = "ICESat2VegR"
#' )
#'
#' # Specifying the path to ATL08 file
#' atl08_path <- system.file("extdata",
#'   "atl08_clip.h5",
#'   package = "ICESat2VegR"
#' )
#'
#' # Reading ATL03 data (h5 file)
#' atl03_h5 <- ATL03_read(atl03_path = atl03_path)
#'
#' # Reading ATL08 data (h5 file)
#' atl08_h5 <- ATL08_read(atl08_path = atl08_path)
#'
#' # Extracting ATL03 and ATL08 photons and heights
#' atl03_atl08_dt <- ATL03_ATL08_photons_attributes_dt_join(atl03_h5, atl08_h5)
#'
#' plot(
#'   atl03_atl08_dt,
#'   "ph_h",
#'   pch = 16, cex = 0.5,
#'   beam = "gt1r",
#'   colors = c("gray", "#bd8421", "forestgreen", "green")
#' )
#'
#' close(atl03_h5)
#' close(atl08_h5)
#' @rdname plot
#' @export
setMethod(
  f = "plot",
  signature("icesat2.atl03atl08_dt", "missing"),
  definition = function(x, ...) {
    plot(x, y = "ph_h", ...)
  }
)

#' @import grDevices
#' @keywords internal
plot_atl03_atl08_join_dt <- function(x, y, beam, col, xlim = NULL, ylim = NULL, colors = c("gray", "#bd8421", "forestgreen", "green"), legend = "topleft", ...) {
  .SD <- data.table::.SD
  if (is.null(beam)) {
    the_beam <- unique(x$beam)[1]
  } else {
    the_beam <- beam
  }
  xdt <- x[beam %in% the_beam, .SD, .SDcols = c("dist_ph_along", y, "classed_pc_flag")][]

  mask <- rep(TRUE, nrow(xdt))
  if (!is.null(xlim)) {
    mask <- mask & xdt$dist_ph_along >= xlim[1] & xdt$dist_ph_along <= xlim[2]
  }
  if (!is.null(ylim)) {
    mask <- mask & xdt[, 2] >= ylim[1] & xdt[, 2] <= ylim[2]
  }

  mask[!stats::complete.cases(mask)] <- FALSE
  mask <- (seq_along(xdt$dist_ph_along))[mask]
  dt_masked <- xdt[mask]

  if (is.null(xlim)) {
    xlim <- range(dt_masked$dist_ph_along)
  }
  if (is.null(ylim)) {
    ylim <- range(dt_masked[[y]])
  }


  old_palette <- grDevices::palette()
  on.exit(grDevices::palette(old_palette))
  grDevices::palette(colors)

  suppressWarnings({
    plot(
      x = dt_masked$dist_ph_along,
      y = dt_masked[, get(y)],
      col = dt_masked$classed_pc_flag + 1,
      xlab = "Distance along-track (m)",
      ylab = paste(y, "(m)"),
      xlim = xlim,
      ylim = ylim,
      ...
    )
    graphics::legend(legend,
      legend = c(
        "ATL03 unclassified",
        "ATL03 Terrain",
        "ATL03 Vegetation",
        "ATL03 Top canopy"
      ), pch = 16, col = colors, bty = "n"
    )
  })
}

#' Plot ATL03_ATL08 attributes
#'
#' @description This function plots ATL08 attributes along track
#'
#' @param x An object of class [icesat2.atl03atl08_dt-class]
#' @param y The attribute name for y axis
#' @param beam Character vector indicating only one beam to process
#' ("gt1l", "gt1r", "gt2l", "gt2r", "gt3l", "gt3r").
#' Default is "gt1r"
#' @param xlim The x limits to use for the plot
#' @param ylim the y limits to use for the plot
#' @param ... will be passed to the main plot, such as:
#'
#'  - `colors`: A vector containing colors for plotting noise, terrain,
#' vegetation and top canopy photons (e.g. c("gray", "#bd8421",
#' "forestgreen", "green")
#'  - `legend`: the position of the legend. 'bottomleft', 'bottomright',
#' 'topleft', 'topright' or FALSE to omit
#'
#' @return No return value
#'
#' @examples
#' # Specifying the path to ATL03 file
#' atl03_path <- system.file("extdata",
#'   "atl03_clip.h5",
#'   package = "ICESat2VegR"
#' )
#'
#' # Specifying the path to ATL08 file
#' atl08_path <- system.file("extdata",
#'   "atl08_clip.h5",
#'   package = "ICESat2VegR"
#' )
#'
#' # Reading ATL03 data (h5 file)
#' atl03_h5 <- ATL03_read(atl03_path = atl03_path)
#'
#' # Reading ATL08 data (h5 file)
#' atl08_h5 <- ATL08_read(atl08_path = atl08_path)
#'
#' # Extracting ATL03 and ATL08 photons and heights
#' atl03_atl08_dt <- ATL03_ATL08_photons_attributes_dt_join(atl03_h5, atl08_h5)
#'
#' plot(
#'   atl03_atl08_dt,
#'   "h_ph",
#'   colors = c("gray", "#bd8421", "forestgreen", "green"),
#'   pch = 16, cex = 0.5
#' )
#'
#' close(atl03_h5)
#' close(atl08_h5)
#' @rdname plot
#' @export
setMethod(
  "plot",
  signature("icesat2.atl03atl08_dt", "character"),
  function(x, y, beam = NULL,
           col = c("gray", "#bd8421", "forestgreen", "green"), ...) {
    if (!is(x, "icesat2.atl03atl08_dt")) {
      print(
        "Invalid input file.
It should be an object of class 'icesat2.atl03atl08_dt'"
      )
    } else {
      plot_atl03_atl08_join_dt(x, y, beam, col, ...)
    }
  }
)

#' Plot ATL08 attributes
#'
#' @description This function plots ATL08 attributes along track
#'
#' @param x An object of class [`ICESat2VegR::icesat2.atl08_dt-class`]
#' @param y The attribute name for y axis
#' @param beam Character vector indicating only one beam to process
#' ("gt1l", "gt1r", "gt2l", "gt2r", "gt3l", "gt3r").
#' Default is "gt1r"
#' @param xlim The x limits to use for the plot
#' @param ylim the y limits to use for the plot
#' @param ... will be passed to the main plot, such as:
#'
#'  - `colors`: A vector containing colors for plotting noise, terrain,
#' vegetation and top canopy photons (e.g. c("gray", "#bd8421",
#' "forestgreen", "green")
#'  - `legend`: the position of the legend. 'bottomleft', 'bottomright',
#' 'topleft', 'topright' or FALSE to omit
#'
#' @return No return value
#'
#' @examples
#' # Specifying the path to ATL08 file
#' atl08_path <- system.file("extdata",
#'   "atl08_clip.h5",
#'   package = "ICESat2VegR"
#' )
#'
#' # Reading ATL08 data (h5 file)
#' atl08_h5 <- ATL08_read(atl08_path = atl08_path)
#'
#' # Extracting ATL08 attributes
#' atl08_seg_dt <- ATL08_seg_attributes_dt(atl08_h5 = atl08_h5)
#'
#' plot(
#'   atl08_seg_dt,
#'   beam = "gt1r",
#'   col = "gray"
#' )
#'
#' close(atl08_h5)
#' @rdname plot
#' @export
setMethod(
  "plot",
  signature("icesat2.atl08_dt", "character"),
  function(x, y, beam = NULL,
           col = c("gray", "#bd8421", "forestgreen", "green"),
           xlim = NULL,
           ylim = NULL, ...) {
    if (!is(x, "icesat2.atl08_dt")) {
      print(
        "Invalid input file. It should be an object of class 'icesat2.atl08_dt'"
      )
    } else {
      plot_atl03_atl08_join_dt(x, y, beam, col, xlim, ylim, ...)
    }
  }
)


#' Plot ATL03 photons
#'
#' @description This function plots ATL03 photons along track
#'
#' @param x An object of class [`ICESat2VegR::icesat2.atl03_dt-class`]
#' @param beam Character vector indicating only one beam to process
#' ("gt1l", "gt1r", "gt2l", "gt2r", "gt3l", "gt3r").
#' Default is "gt1r"
#' @param col  Color for plotting the photons. Default is "gray"
#' @param ... will be passed to the main plot
#'
#' @return No return value
#'
#' @examples
#' # Specifying the path to ATL03 file
#' atl03_path <- system.file("extdata",
#'   "atl03_clip.h5",
#'   package = "ICESat2VegR"
#' )
#'
#' # Reading ATL03 data (h5 file)
#' atl03_h5 <- ATL03_read(atl03_path = atl03_path)
#'
#' # Extracting atl03 and atl03 photons and heights
#' atl03_photons_dt <- ATL03_seg_metadata_dt(
#'   atl03_h5 = atl03_h5,
#'   attributes = c(
#'     "reference_photon_lon",
#'     "reference_photon_lat",
#'     "segment_dist_x",
#'     "h_ph"
#'   )
#' )
#'
#' plot(
#'   atl03_photons_dt,
#'   "h_ph",
#'   col = "gray",
#'   pch = 16,
#'   cex = 0.5
#' )
#'
#' close(atl03_h5)
#' @export
#' @rdname plot
setMethod(
  f = "plot",
  signature("icesat2.atl03_dt", y = "character"),
  definition = function(x, y, col = "gray", ...) {
    beam <- NA

    x <- x[x$beam == beam, ]
    suppressWarnings({
      plot(
        x = x$segment_dist_x,
        y = x$h_ph,
        col = col,
        xlab = "Distance along-track (m)",
        ylab = "Elevation (m)",
        ...
      )
      graphics::legend(
        "topleft",
        legend = c("Noise", "Terrain", "Vegetation", "Top canopy"),
        pch = 16,
        col = col,
        bty = "n"
      )
    })
  }
)

genericICESatC <- function(classname, x, ...) {
  function(x, ...) {
    dt_list <- list(..., x)
    dt <- data.table::rbindlist(dt_list)
    prepend_class(dt, classname)
    dt
  }
}

#' @export
"c.icesat2.atl03_seg_dt" <- genericICESatC("icesat2.atl03_seg_dt")
#' @export
"c.icesat2.atl08_dt" <- genericICESatC("icesat2.atl08_dt")


#' Row-bind a list of objects while preserving class
#'
#' @description
#' Convenience wrapper around [`data.table::rbindlist()`] that
#' **row-binds a list of homogeneous objects and then restores their class**.
#' This is especially useful for custom S3/S4-like classes used in ICESat2VegR
#' (e.g. `icesat2.atl03_dt`, `icesat2.atl08_dt`), where a plain
#' `data.table::rbindlist()` call would drop the original class attribute.
#'
#' @param l A list whose elements are all of the *same* class, typically
#'   `data.table`, `data.frame`, or a custom subclass (e.g. `icesat2.*_dt`).
#'   All elements must share identical `class()`; otherwise an error is raised.
#' @param ... Additional arguments passed directly to
#'   [`data.table::rbindlist()`], such as `use.names` or `fill`.
#'
#' @return
#' A single `data.table` (the result of `data.table::rbindlist()`), with its
#' `class` attribute reset to match the class of the input elements.
#'
#' @export
rbindlist2 <- function(l, ...) {
  classes <- unique(lapply(l, class))
  n_classes <- length(classes)

  stopifnot(
    "All elements should belong to the same class(es)!" = n_classes == 1
  )

  # Row-bind with data.table::rbindlist, honoring any extra arguments
  result <- data.table::rbindlist(l, ...)

  # Restore original class (e.g. icesat2.atl03_dt, icesat2.atl08_dt, etc.)
  attr(result, "class") <- classes[[1]]

  result
}
