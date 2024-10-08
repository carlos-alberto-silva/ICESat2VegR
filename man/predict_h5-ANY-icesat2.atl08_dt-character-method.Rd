% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/predict_h5.R
\name{predict_h5,ANY,icesat2.atl08_dt,character-method}
\alias{predict_h5,ANY,icesat2.atl08_dt,character-method}
\title{S4 method for predicting using HDF5 file as output}
\usage{
\S4method{predict_h5}{ANY,icesat2.atl08_dt,character}(model, dt, output)
}
\arguments{
\item{model}{The trained model object}

\item{dt}{The input data.table to run the model}

\item{output}{The output file path}
}
\value{
An \code{\linkS4class{icesat2.predict_h5}}, which is an
h5 file with latitude, longitude and prediction datasets.
}
\description{
This method is used to predict using a trained model and
save the results in HDF5 file.
}
\details{
This method is used to predict using a trained model and
save the results in an HDF5 file.
}
\examples{
atl08_path <- system.file(
  "extdata",
  "atl08_clip.h5",
  package = "ICESat2VegR"
)
atl08_h5 <- ATL08_read(atl08_path = atl08_path)
atl08_dt <- ATL08_seg_attributes_dt(atl08_h5)
linear_model <- stats::lm(h_canopy ~ canopy_openness, data = atl08_dt)
output_h5 <- tempfile(fileext = ".h5")
predicted_h5 <- predict_h5(linear_model, atl08_dt, output_h5)

# List datasets
predicted_h5$ls()$name

# See predicted values
head(predicted_h5[["prediction"]][])

# Close the file
close(predicted_h5)

}
