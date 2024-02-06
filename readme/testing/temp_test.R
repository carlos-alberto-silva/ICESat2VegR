# source("R/clipTools.R")
# source("R/ATL08_h5_clipBox.R")

# library(ICESat2VegR)

# atl03_path <- "../inst/exdata/ATL03_20220401221822_01501506_005_01.h5"

# atl03 <- ATL03_read(atl03_path)
# vect <- terra::vect("inst/exdata/polygons.shp")
# output <- file.path("C:/Users/caiohamamura/Desktop/saida", "output.h5")
# polygon_id <- "id"

# bbox <- terra::ext(vect)

# out_h5 <- ATL03_h5_clipBox(atl03, output, bbox)
# close(out_h5)
# file.remove(out_h5)
# out_h5 <- clip(atl03, output, bbox)
# close(out_h5)
# file.remove(output)


# out_h5 <- clip(atl03, output, vect, "id")
# lapply(out_h5, close)





library(ICESat2VegR)

atl08_path <- "../inst/exdata/ATL08_20220401221822_01501506_005_01.h5"

atl08 <- ATL08_read(atl08_path)
vect <- terra::vect("inst/exdata/polygons.shp")
output <- file.path("C:/Users/caiohamamura/Desktop/saida", "output2.h5")
polygon_id <- "id"


landSegmentsMask_fn = landsegmentsMask_bbox
bbox <- terra::ext(vect)
clip_obj <- bbox
# out_h5 <- ATL08_h5_clipBox(atl08, output, bbox)
# close(out_h5)
# file.remove(output)


out_h5 <- ATL08_h5_clipGeometry(atl08, output, vect, "id")

lapply(out_h5, close)
