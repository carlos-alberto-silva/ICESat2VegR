setRefClass("icesat2.predict_h5")

internal_predict_h5 <- new.env()
internal_predict_h5$path <- ""
internal_predict_h5$h5 <- NULL

#' Model prediction over data.tables using HDF5 file as output
#'
#' Model prediction over a data.table from ATL03 or ATL08 data
#' containing geolocation data.
#' It can both append results to an existing HDF5 file or
#' create a new file, allowing to incrementally add predictions
#' to the file to avoid memory issues.
#'
#' @param model The trained model object
#' @param dt The input data.table to run the model
#' @param output The output HDF5 file path
#'
#' @return An [`icesat2.predict_h5-class`], which is an
#' h5 file with latitude, longitude and prediction datasets.
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
#' head(predicted_h5[["prediction"]][])
#'
#' # Close the file
#' close(predicted_h5)
#'
#' @export
setGeneric("predict_h5", function(model, dt, output) {
  standardGeneric("predict_h5")
})

#' S4 method for predicting using HDF5 file as output
#'
#' This method is used to predict using a trained model and
#' save the results in HDF5 file.
#'
#' @param model The trained model object
#' @param dt The input data.table to run the model
#' @param output The output file path
#'
#' @return An [`icesat2.predict_h5-class`], which is an
#' h5 file with latitude, longitude and prediction
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
#' head(predicted_h5[["prediction"]][])
#'
#' # Close the file
#' close(predicted_h5)
#' @export
setMethod(
  "predict_h5",
  signature(model = "ANY", dt = "icesat2.atl03_seg_dt", output = "character"),
  function(model, dt, output) {
    h5 <- hdf5r::H5File$new(output, "w")

    h5[["prediction"]] <- predict(model, dt)
    h5[["longitude"]] <- dt[["reference_photon_lon"]]
    h5[["latitude"]] <- dt[["reference_photon_lat"]]

    h5$close()
    h5 <- hdf5r::H5File$new(output, "r")
    prepend_class(h5, "icesat2.predict_h5")

    return(h5)
  }
)


#' S4 method for predicting using HDF5 file as output
#'
#' This method is used to predict using a trained model and
#' save the results in HDF5 file.
#'
#' @param model The trained model object
#' @param dt The input data.table to run the model
#' @param output The output file path
#'
#' @details
#' This method is used to predict using a trained model and
#' save the results in an HDF5 file.
#'
#' @return An [`icesat2.predict_h5-class`], which is an
#' h5 file with latitude, longitude and prediction datasets.
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
#' head(predicted_h5[["prediction"]][])
#'
#' # Close the file
#' close(predicted_h5)
#'
#' @export
setMethod(
  "predict_h5",
  signature(model = "ANY", dt = "icesat2.atl08_dt", output = "character"),
  function(model, dt, output) {
    if (internal_predict_h5$path == output) {
      h5 <- internal_predict_h5$h5$close_all()
    }

    if (!file.exists(output)) {
      h5 <- hdf5r::H5File$new(output, "w")
    } else {
      h5 <- hdf5r::H5File$new(output, "r+")
    }
    internal_predict_h5$path <- output

    if (h5$exists("prediction")) {
      append_range <- seq(h5[["prediction"]]$dims + 1, length.out = nrow(dt))
      h5[["prediction"]][append_range] <- predict(model, dt)
      h5[["longitude"]][append_range] <- dt[["longitude"]]
      h5[["latitude"]][append_range] <- dt[["latitude"]]
    } else {
      h5[["prediction"]] <- predict(model, dt)
      h5[["longitude"]] <- dt[["longitude"]]
      h5[["latitude"]] <- dt[["latitude"]]
    }

    h5$close_all()
    h5 <- hdf5r::H5File$new(output, "r")
    prepend_class(h5, "icesat2.predict_h5")
    internal_predict_h5$h5 <- h5

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
