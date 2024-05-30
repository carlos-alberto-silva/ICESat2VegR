setRefClass("icesat2.predict_h5")

#' S4 generic method for predicting using H5 model
#'
#' @param model The trained model object
#' @param seg_dt The input segmentation data
#' @param output The output file path
#'
#' @return An [`icesat2.predict_h5-class`], which is an
#' h5 file with latitude, longitude and predicted datasets.
#'
#' @examples
#' atl03_path <- system.file(
#'   "extdata",
#'   "atl03_clip.h5",
#'   package = "ICESat2VegR"
#' )
#'
#' atl03_h5 <- ATL03_read(atl03_path = atl03_path)
#' atl03_seg_dt <- ATL03_seg_attributes_dt(atl03_h5)
#'
#' linear_model <- stats::lm(h_ph ~ segment_ph_cnt, data = atl03_seg_dt)
#' output_h5 <- tempfile(fileext = ".h5")
#' predicted_h5 <- predict_h5(linear_model, atl03_seg_dt, output_h5)
#'
#' # List datasets
#' predicted_h5$ls()$name
#'
#' # See predicted values
#' head(predicted_h5[["predicted"]][])
#'
#' # Close the file
#' close(predicted_h5)
#'
#' @export
setGeneric("predict_h5", function(model, seg_dt, output) {
  standardGeneric("predict_h5")
})
#' Method for predicting using atl03 segments data
#'
#' @param model The trained model object
#' @param seg_dt The input segmentation data
#' @param output The output file path
#'
#' @return An [`icesat2.predict_h5-class`], which is an
#' h5 file with latitude, longitude and predicted datasets.
#'
#' @examples
#' atl03_path <- system.file(
#'   "extdata",
#'   "atl03_clip.h5",
#'   package = "ICESat2VegR"
#' )
#'
#' atl03_h5 <- ATL03_read(atl03_path = atl03_path)
#' atl03_seg_dt <- ATL03_seg_attributes_dt(atl03_h5)
#'
#' linear_model <- stats::lm(h_ph ~ segment_ph_cnt, data = atl03_seg_dt)
#' output_h5 <- tempfile(fileext = ".h5")
#' predicted_h5 <- predict_h5(linear_model, atl03_seg_dt, output_h5)
#'
#' # List datasets
#' predicted_h5$ls()$name
#'
#' # See predicted values
#' head(predicted_h5[["predicted"]][])
#'
#' # Close the file
#' close(predicted_h5)
#' @export
setMethod(
  "predict_h5",
  signature(model = "ANY", seg_dt = "icesat2.atl03_seg_dt", output = "character"),
  function(model, seg_dt, output) {
    h5 <- hdf5r::H5File$new(output, "w")

    h5[["prediction"]] <- predict(model, seg_dt)
    h5[["longitude"]] <- seg_dt[["reference_photon_lon"]]
    h5[["latitude"]] <- seg_dt[["reference_photon_lat"]]

    h5$close()
    h5 <- hdf5r::H5File$new(output, "r")
    prepend_class(h5, "icesat2.predict_h5")

    return(h5)
  }
)


#' Method for predicting using atl08 segments data
#'
#' @param model The trained model object
#' @param seg_dt The input segmentation data
#' @param output The output file path
#'
#' @return An [`icesat2.predict_h5-class`], which is an
#' h5 file with latitude, longitude and predicted datasets.
#'
#' @examples
#' atl08_path <- system.file(
#'   "extdata",
#'   "atl08_clip.h5",
#'   package = "ICESat2VegR"
#' )
#' atl08_h5 <- ATL08_read(atl08_path = atl08_path)
#' atl08_dt <- ATL08_seg_attributes_dt(atl08_h5)
#' linear_model <- stats::lm(h_canopy ~ canopy_openness, data = atl08_dt)
#' output_h5 <- tempfile(fileext = ".h5")
#' predicted_h5 <- predict_h5(linear_model, atl08_dt, output_h5)
#'
#' # List datasets
#' predicted_h5$ls()$name
#'
#' # See predicted values
#' head(predicted_h5[["predicted"]][])
#'
#' # Close the file
#' close(predicted_h5)
#'
#' @export
setMethod(
  "predict_h5",
  signature(model = "ANY", seg_dt = "icesat2.atl08_dt", output = "character"),
  function(model, seg_dt, output) {
    h5 <- hdf5r::H5File$new(output, "w")

    h5[["prediction"]] <- predict(model, seg_dt)
    h5[["longitude"]] <- seg_dt[["longitude"]]
    h5[["latitude"]] <- seg_dt[["latitude"]]

    h5$close_all()
    h5 <- hdf5r::H5File$new(output, "r")
    prepend_class(h5, "icesat2.predict_h5")

    return(h5)
  }
)

#' @export
"length.icesat2.predict_h5" <- function(x) {
  exp(sum(log(x[["latitude"]]$dims)))
}

#' @export
setMethod(
  "close",
  signature = "icesat2.predict_h5",
  definition = function(con, ...) {
    con$close_all()
  }
)
