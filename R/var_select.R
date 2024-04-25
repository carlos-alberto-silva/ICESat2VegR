# #' Uses Google Earth Engine for selecting the variables
# #' based on a forward feature selection wrapper
# #'
# #' @param x [`data.frame-class`]. The vectors of features available for predicion
# #'
# #' @return A [`data.frame-class`] containing the selected properties and the mean
# #' RMSEs from the bootstrapping runs for each variable selected.
# #'
# #' @include model_tools.R
# #' @export
# var_select <- function(x, y, method = "forward", nboots = 10, nTrees = 100, train_split = 0.7, delta = 0.01) {
#   ee_bandNames <- colnames(x)
#   n <- nrow(x)
#   train_size <- as.integer(round(n * train_split))
#   validation_size <- n - train_size
#   x <- build_fc(x, y)
#   y <- names(y)
# 
#   selected_properties <- ee$List(list())
#   rmses <- ee$List(list())
# 
#   minRmse <- ee_number(1e10)
#   lastRmse <- ee_number(1e11)
#   while (((1 - minRmse / lastRmse) >= delta)$getInfo() == 1) {
#     result <- list(
#       property = ee$List(list()),
#       rmse = ee$List(list())
#     )
# 
#     for (band in ee_bandNames) {
#       current <- selected_properties$add(band)
# 
#       rmseList <- ee$List(list())
#       message(gettextf("Testing %s", list(current$getInfo())), appendLF = TRUE)
#       for (i in 1:nboots) {
#         message(gettextf("\rBoot %d/%d", i, nboots), appendLF = FALSE)
#         x <- x$randomColumn(seed = floor(runif(1) * INT_MAX))
#         train_sample <- x$limit(train_size, "random")$select(current$add(y))
#         validation_sample <- x$limit(validation_size, "random", ascending = FALSE)$select(current$add(y))
# 
# 
#         randomForestClassifier <- randomForestRegression(train_sample, property = y, train_properties = current, nTrees = nTrees)
#         classification <- validation_sample$classify(randomForestClassifier)
#         classification2 <- classification$map(function(f) {
#           f[["sqerror"]] <- (f[["h_canopy"]] - f[["classification"]])^2
#           return(f)
#         })
#         rmse <- sqrt(classification2$aggregate_mean("sqerror"))
#         rmseList[[]] <- rmse
#       }
#       message(appendLF = TRUE)
# 
#       meanRmse <- rmseList$reduce(ee$Reducer$mean())
#       result$rmse <- result$rmse$add(meanRmse)
#       result$property <- result$property$add(band)
#     }
#     lastRmse <- minRmse
#     minRmse <- result$rmse$reduce(ee$Reducer$min())
#     minIdx <- result$rmse$indexOf(minRmse)
#     minRmse <- result$rmse$getNumber(minIdx)
#     # ee$Dictionary(result)$getInfo()
#     rmses <- rmses$add(minRmse)
# 
#     # rmses$getInfo()
# 
#     bandAdd <- result$property$getString(minIdx)
#     selected_properties <- selected_properties$add(bandAdd)
#     ee_bandNames <- ee_bandNames[ee_bandNames != bandAdd$getInfo()]
#   }
# 
#   final_result <- ee$Dictionary(list(
#     properties = selected_properties,
#     rmse = rmses
#   ))
# 
#   df <- data.frame(final_result$getInfo())
#   return(df)
# }