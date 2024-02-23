#' Given a geometry with point samples and images from Earth Engine
#' retrieve the point geometry with values for the images
#'
#' @param images A single image or a vector/list of images from Earth Engine.
#' @param geom A geometry from [`terra::SpatVector-class`] read with [`terra::vect`].
#' @param scale The scale in meters for the extraction (image resolution).
#'
#' @export
extract <- function(images, geom, scale) {
    final <- c(images)
    tempjson <- tempfile(fileext = ".geojson")
    terra::writeVector(geom, tempjson, filetype = "geojson")
    parsed <- jsonlite::parse_json(readLines(tempjson))
    geojson <- ee$FeatureCollection(parsed)

    sampled <- final$sampleRegions(
        collection = geojson,
        scale = scale
    )

    return(sampled)
}

#' @export
ee_to_df <- function(sampled) {
    ids <- sampled$aggregate_array("system:index")$getInfo()
    ids <- as.integer(gsub("_0", "", ids)) + 1
    df <- data.frame(ids = ids)
    df[["ids"]] <- ids

    column_df <- sampled$first()$propertyNames()$getInfo()[-1]
    nested_list <- sampled$reduceColumns(ee$Reducer$toList(length(column_df)), column_df)$values()$get(0)

    pymain <- reticulate::import_main()
    pymain$nested_list <- nested_list
    reticulate::py_run_string("import numpy as np")
    reticulate::py_run_string("result = np.array(nested_list.getInfo())")


    df2 <- data.frame(pymain$result)
    names(df2) <- column_df

    return(cbind(df, df2))
}

#' @export
randomForestRegression <- function(
    featurecollection,
    property_name,
    nTrees = 500,
    mtry = NULL,
    nodesize = 5,
    sampsize = .632,
    seed = floor((runif(1) * 2147483647 * 2) - 2147483647)) {
    if (is.null(mtry)) {
        mtry <- max(1, floor(length(featurecollection$first()$propertyNames()$getInfo()) / 3))
        classifier <- ee$Classifier$smileRandomForest(
            numberOfTrees = nTrees,
            variablesPerSplit = mtry,
            minLeafPopulation = nodesize,
            bagFraction = sampsize,
            seed = seed
        )$
            train(
            features = featurecollection,
            classProperty = "h_canopy"
        )$
            setOutputMode("REGRESSION")

        return(classifier)
    }

    ee$Classifier$smileRandomForest(
        numberOfTrees = nTrees,
        variablesPerSplit = mtry,
        minLeafPopulation = nodesize,
        bagFraction = sampsize,
        seed = seed
    )$
        train(
        features = featurecollection,
        classProperty = "h_canopy"
    )$
        setOutputMode("REGRESSION")
}


setRefClass("ee.Classifier")
#' @export
"predict.ee.Classifier" <- function(x, data, ...) {
    predicted <- data$classify(x)
    classification <- predicted$aggregate_array("classification")$getInfo()
    return(classification)
}
