#' @importFrom hdf5r H5File
setRefClass("H5File")
requireNamespace("data.table")

#' Class for ICESat-2 ATL08
#'
#' @slot h5 Object of class [`H5File`][hdf5r::H5File-class] from `hdf5r` package containing the
#' ICESat-2 Land and Vegetation Along-Track Products (ATL08)

#'
#' @seealso [`H5File`][hdf5r::H5File-class] in the `hdf5r` package and
#' \url{https://icesat-2.gsfc.nasa.gov/sites/default/files/page_files/ICESat2_ATL08_ATBD_r006.pdf}
#'
#' @import methods
#' @export
icesat2.atl08_h5 <- setClass(
  Class = "icesat2.atl08_h5",
  slots = list(h5 = "H5File")
)

#' @importFrom hdf5r H5File
setRefClass("H5File")
requireNamespace("data.table")

#' Class for ICESat-2 ATL03
#'
#' @slot h5 Object of class [`H5File`][hdf5r::H5File-class] from `hdf5r` package containing the
#' ICESat-2 Global Geolocated Photon Data (ATL03)

#'
#' @seealso [`H5File`][hdf5r::H5File-class] in the `hdf5r` package and
#' \url{https://icesat-2.gsfc.nasa.gov/sites/default/files/page_files/ICESat2_ATL03_ATBD_r006.pdf}
#'
#' @import methods
#' @export
icesat2.atl03_h5 <- setClass(
  Class = "icesat2.atl03_h5",
  slots = list(h5 = "H5File")
)


#' @importFrom data.table data.table
setRefClass("data.table")
requireNamespace("data.table")

#' Class for ATL08 attributes
#'
#' @slot data.table Object of class [`data.table`][data.table::data.table-class] from `data.table` package containing the
#' Extracted ICESat-2 ATL08 attributes

#'
#' @seealso [`data.table`][data.table::data.table-class] in the `data.table` package and
#' \url{https://icesat-2.gsfc.nasa.gov/sites/default/files/page_files/ICESat2_ATL03_ATBD_r006.pdf}
#'
#' @import methods
#' @export
icesat2.atl08_dt <- setClass(
  Class = "icesat2.atl08_dt",
  slots = list(dt = "data.table")
)

#' @importFrom data.table data.table
setRefClass("data.table")
requireNamespace("data.table")

#' Class for joined ATL03 and ATL08 attributes
#'
#' @slot data.table Object of class [`data.table`][data.table::data.table-class] from `data.table` package containing the
#' Joined ICESat-2 ATL03 ATL08 attributes

#'
#' @seealso [`data.table`][data.table::data.table-class] in the `data.table` package and
#' \url{https://icesat-2.gsfc.nasa.gov/sites/default/files/page_files/ICESat2_ATL03_ATBD_r006.pdf}
#'
#' @import methods
#' @export
icesat2.atl03atl08_dt <- setClass(
  Class = "icesat2.atl03atl08_dt",
  slots = list(dt = "data.table")
)


h5closeall <- function(con, ...) {
  try(con@h5$close_all(), silent = TRUE)
}


#' Safely closes the [`rICESat2Veg::icesat2.atl03_h5-class`]
#'
#' @description
#' Closing files will avoid locking HDF5 ATL03 files.
#'
#' @param con An object of class `icesat2.atl03_h5`
#' @param ... Inherited from base
#'
#' @export
#' @rdname close
#' @method close icesat2.atl03_h5
setMethod("close", signature = c("icesat2.atl03_h5"), h5closeall)


#' Safely closes the [`rICESat2Veg::icesat2.atl08_h5-class`]
#'
#' @description
#' Closing files will avoid locking HDF5 ATL08 files.
#'
#' @param con An object of class `icesat2.atl08_h5`
#' @param ... Inherited from base
#' @method close icesat2.atl08_h5
#' @rdname close
setMethod("close", signature = c("icesat2.atl08_h5"), h5closeall)


#' Plot photons from ATL03 and ATL08 joined products
#'
#' @description This function plots photons along track
#'
#' @param x An object of class [`rICESat2Veg::icesat2.atl03atl08_dt-class`]
#' @param attribute photon attribute to be plot (ph_h or h_ph)
#' @param colors A vector containing colors for plotting noise, terrain, vegetation and top canopy photons
#' (e.g. c("gray", "#bd8421", "forestgreen", "green")
#' @param ... will be passed to the main plot
#'
#' @return No return value
#'
#' @examples
#'# Specifying the path to ATL03 and ATL08 file (zip file)
#'outdir = tempdir()
#'atl03_zip <- system.file("extdata",
#'                   "ATL03_20220401221822_01501506_005_01.zip",
#'                   package="rICESat2Veg")
#'
#'atl08_zip <- system.file("extdata",
#'                   "ATL08_20220401221822_01501506_005_01.zip",
#'                   package="rICESat2Veg")
#'
#'# Unzipping ATL03 file
#'atl03_path <- unzip(atl03_zip,exdir = outdir)
#'
#'# Unzipping ATL08 file
#'atl08_path <- unzip(atl08_zip,exdir = outdir)
#'
#'# Reading ATL03 data (h5 file)
#atl03_h5<-ATL03_read(atl03_path=atl03_path)
#'
#'# Reading ATL08 data (h5 file)
#atl08_h5<-ATL08_read(atl08_path=atl08_path)
#'
#'# Extracting ATL03 and ATL08 photons and heights
#'atl03_atl08_dt<-ATL03_ATL08_join_dt(atl03_h5,atl08_h5)
#'
#'rICESat2Veg::plot(atl03_atl08_dt=atl03_atl08_dt,attribute="ph_h",
#'                 colors = c("gray", "#bd8421", "forestgreen", "green"),
#'                 pch = 16, cex = 0.5)
#'
#'close(atl03_h5)
#'close(atl08_h5)
#' @export
#' @method plot icesat2.atl03atl08_dt
#' @rdname plot
setMethod(
  f = "plot",
  signature("icesat2.atl03atl08_dt", y = "missing", ),
  definition = function(x, attribute, colors, ...) {
    if (!is(x, "icesat2.atl03atl08_dt")) {
      print("Invalid input file. It should be an object of class 'icesat2.atl03atl08_dt' ")
    } else {


      #colors <- c("gray", "#bd8421", "forestgreen", "green")
      colorMap <- colors[x$classed_pc_flag + 1]
      x<-x$dist_along
      y<-x[,attribute]

        suppressWarnings({
          plot(x=x, y=y, col = colorMap,...)
          legend("topleft", legend=c("Noise","Terrain", "Vegetation", "Top canopy"), pch=16, col=colors, bty="n")
        })

    }
  }
)
