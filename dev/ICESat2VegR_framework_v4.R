#=============================================================================
# ICESat2VegR - Manuscript Reproducibility Pipeline (Fully Commented, No Roxygen)
#=============================================================================
# This script sets up the environment, discovers/downloads ICESat-2 granules,
# aggregates 10 m segments, builds predictors from Earth Engine, fits a model,
# and visualizes predictions. Every line is commented for clarity.
#=============================================================================


#=============================================================================
# 1) Repositories (r-universe + CRAN)
#=============================================================================
repos <- c("https://caiohamamura.r-universe.dev", "https://cloud.r-project.org")    # Repos for ICESat2VegR (r-universe) and CRAN


#=============================================================================
# 2) Required R packages (install if missing)
#=============================================================================
need_pkgs <- c(                                                                       # Packages required by the workflow
  "ICESat2VegR",  # ICESat-2 vegetation processing
  "reticulate",   # Python <-> R interface
  "leaflet",      # interactive maps
  "sf",           # spatial vectors
  "terra",        # rasters & vectors
  "data.table",   # fast tables
  "dplyr"         # tidy helpers
)
missing <- need_pkgs[!need_pkgs %in% rownames(installed.packages())]                 # Identify missing packages
if (length(missing)) {                                                                # If any are missing
  message("Installing missing R packages: ", paste(missing, collapse = ", "))         # Inform the user
  install.packages(missing, repos = repos, dependencies = TRUE)                       # Install from repos with dependencies
}


#=============================================================================
# 3) Load libraries quietly
#=============================================================================
suppressPackageStartupMessages({                                                       # Suppress startup messages for clean logs
  library(reticulate)                                                                  # Interface to Python
  library(ICESat2VegR)                                                                 # ICESat-2 vegetation tools
  library(leaflet)                                                                     # Interactive maps
  library(sf)                                                                          # Simple Features for vector data
  library(terra)                                                                       # Raster + vector geospatial ops
  library(data.table)                                                                  # Fast data tables
  library(dplyr)                                                                       # Tidy verbs
})
cat("Loading libraries... done.\n")                                                    # Confirm loading


#=============================================================================
# 4) Configure Python environment (optional interactive helper)
#=============================================================================
ICESat2VegR_configure()


#=============================================================================
# 5) Verify configuration status
#=============================================================================

# Helper safely() runs an expression and returns a default value if an error occurs
safely <- function(expr, default = NA) tryCatch(expr, error = function(e) default)

# Check if 'reticulate' is installed
have_reticulate <- requireNamespace("reticulate", quietly = TRUE)

# Attempt to discover which Python executable is being used
py_path <- if (have_reticulate) safely(reticulate::py_discover_config()$python, NA) else NA

# Verify if Python is properly initialized and available
py_ready <- have_reticulate && !inherits(try(reticulate::py_config(), silent = TRUE), "try-error")

# Build a list summarizing the environment and module status
status <- list(
  ICESat2VegR_loaded = "package:ICESat2VegR" %in% search(),  # Is package loaded in current session?
  python_used        = py_path,                              # Path to active Python binary
  h5py               = if (py_ready) reticulate::py_module_available("h5py")        else FALSE,
  earthaccess        = if (py_ready) reticulate::py_module_available("earthaccess") else FALSE,
  ee                 = if (py_ready) reticulate::py_module_available("ee")          else FALSE
)

# Print the collected environment information to the console
print(status)


#=============================================================================
# 6) Authenticate & initialize Earth Engine
#=============================================================================
ee_initialize(project = "ee-carlossilvaengflorestal")                                   # Initialize Earth Engine under a specific GCP project

#=============================================================================
# 7) Helper functions
#=============================================================================

Generate_EmbeddingAndAncillaryEEstack <- function(                                     # Create an EE stack (AlphaEarth + terrain + optional lon/lat)
  geom,                                                                                # AOI geometry (sf/sfc, SpatVector, or EE geometry)
  startYY, endYY,                                                                      # Inclusive year range
  mask_outside = TRUE,                                                                 # Mask pixels outside AOI
  terrain_scale = 30,                                                                  # Reproject terrain to ~30 m
  multiply_slope_aspect_by10 = TRUE,                                                   # Keep slope/aspect scaling convention
  add_lonlat = TRUE                                                                     # Add lon/lat bands
) {
  ee <- reticulate::import("ee", delay_load = FALSE)                                    # Import Python Earth Engine
  .ee_ping(ee)                                                            # Ping EE session (internal helper)

  startYY <- as.integer(startYY); endYY <- as.integer(endYY)                            # Coerce years to integers
  if (!is.finite(startYY) || !is.finite(endYY) || endYY < startYY)                      # Validate year range
    stop("Invalid years: startYY/endYY")                                                # Abort if invalid

  ee_geom <- ICESat2VegR:::.as_ee_geom(geom)                                            # Convert R geometry to EE geometry

  emb_ic <- ee$ImageCollection("GOOGLE/SATELLITE_EMBEDDING/V1/ANNUAL")$                # AlphaEarth annual embeddings collection
    filterDate(sprintf("%04d-01-01", startYY), sprintf("%04d-12-31", endYY))$           # Filter by year range
    filterBounds(ee_geom)                                                               # Filter by AOI
  embedding <- emb_ic$median()$clip(ee_geom)                                            # Median composite and clip to AOI

  elevation <- tryCatch({                                                               # Try 3DEP 1 m elevation
    ee$ImageCollection("USGS/3DEP/1m")$mosaic()$select("elevation")$clip(ee_geom)
  }, error = function(e) NULL)
  if (is.null(elevation)) {                                                             # Fallback to 3DEP 10 m
    elevation <- tryCatch({
      ee$ImageCollection("USGS/3DEP/10m")$mosaic()$select("elevation")$clip(ee_geom)
    }, error = function(e) NULL)
  }
  if (is.null(elevation)) {                                                             # Fallback to NASADEM
    elevation <- ee$Image("NASA/NASADEM_HGT/001")$select("elevation")$clip(ee_geom)
  }

  elev_proj <- elevation$reproject(crs = "EPSG:4326", scale = terrain_scale)            # Reproject elevation to target scale
  slope  <- ee$Terrain$slope(elev_proj)$rename("slope")                                  # Compute slope
  aspect <- ee$Terrain$aspect(elev_proj)$rename("aspect")                                # Compute aspect
  if (isTRUE(multiply_slope_aspect_by10)) {                                              # Optionally scale values by 10
    slope  <- slope$multiply(10)
    aspect <- aspect$multiply(10)
  }
  slope  <- slope$clip(ee_geom)                                                          # Clip slope to AOI
  aspect <- aspect$clip(ee_geom)                                                         # Clip aspect to AOI

  if (isTRUE(add_lonlat)) {                                                              # Optionally add lon/lat bands
    lonlat <- ee$Image$pixelLonLat()$clip(ee_geom)
    lon    <- lonlat$select("longitude")$rename("lon")
    lat    <- lonlat$select("latitude")$rename("lat")
  }

  stack <- embedding$                                                                    # Start with embedding
    addBands(elevation$rename("elevation"))$                                             # Add elevation
    addBands(slope)$                                                                      # Add slope
    addBands(aspect)                                                                      # Add aspect
  if (isTRUE(add_lonlat)) stack <- stack$addBands(lon)$addBands(lat)                     # Add lon/lat if requested

  if (isTRUE(mask_outside)) {                                                            # Mask outside AOI
    mask <- ee$Image$constant(1)$clip(ee_geom)
    stack <- stack$updateMask(mask)
  }

  stack                                                                                  # Return the final EE image stack
}


sample_granules_by_year <- function(urls, n_per_year = 5, seed = NULL) {                # Sample up to n_per_year ATL URLs per year
  if (is.matrix(urls) || is.data.frame(urls)) urls <- as.character(urls[, 1])           # Normalize to character vector
  stopifnot(is.character(urls), length(urls) > 0)                                       # Basic checks

  get_year <- function(u) {                                                              # Helper to parse year from URL or filename
    m <- regexpr("/(20\\d{2})/", u, perl = TRUE)                                        # Prefer /YYYY/ in path
    if (m[1] > 0) return(substr(u, m[1] + 1, m[1] + attr(m, "match.length") - 2))       # Extract year from path
    m2 <- regexpr("ATL0[38]_(20\\d{2})\\d{10}", u, perl = TRUE)                         # Fallback: ATL0X_YYYYMMDD...
    if (m2[1] > 0) return(substr(u, m2[1] + 6, u, m2[1] + 9))                           # Extract year from filename
    NA_character_                                                                        # Return NA if not found
  }

  years <- vapply(urls, get_year, character(1))                                         # Vectorize year extraction
  keep  <- !is.na(years)                                                                 # Keep only URLs with a year
  if (!all(keep)) warning("Dropped ", sum(!keep), " URL(s) with no detectable year.")    # Warn if any dropped

  urls  <- urls[keep]                                                                    # Filter URLs
  years <- years[keep]                                                                   # Filter corresponding years

  if (!is.null(seed)) set.seed(seed)                                                     # Set RNG seed if provided
  split_urls <- split(urls, years)                                                       # Split by year

  sampled <- lapply(names(split_urls), function(y) {                                     # For each year
    u <- sort(unique(split_urls[[y]]))                                                   # Unique sorted URLs
    if (length(u) <= n_per_year) {                                                       # If few URLs, keep all
      data.frame(year = as.integer(y), url = u, stringsAsFactors = FALSE)               # Return as data.frame
    } else {                                                                             # Otherwise sample
      pick <- sample(u, n_per_year)
      data.frame(year = as.integer(y), url = pick, stringsAsFactors = FALSE)
    }
  })

  out <- do.call(rbind, sampled)                                                         # Row-bind per-year data.frames
  rownames(out) <- NULL                                                                  # Clean row names
  out[order(out$year, out$url), ]                                                        # Order by year, then URL
}


safe_write_geojson <- function(dt, path, xcol = "lon", ycol = "lat",                    # Write GeoJSON from table with lon/lat
                               crs = "EPSG:4326", overwrite = TRUE) {
  has_sf <- requireNamespace("sf", quietly = TRUE)                                       # Check if sf is available
  try({                                                                                  # First attempt: use terra via ICESat2VegR::to_vect
    sv <- ICESat2VegR::to_vect(dt, xcol = xcol, ycol = ycol, crs = crs)                  # Table -> SpatVector
    terra::writeVector(sv, path, overwrite = overwrite)                                  # Write GeoJSON
    return(invisible(TRUE))                                                              # Success
  }, silent = TRUE)
  if (has_sf) {                                                                          # Fallback: use sf
    sf_pts <- sf::st_as_sf(as.data.frame(dt), coords = c(xcol, ycol), crs = 4326, remove = FALSE)  # Create sf points
    sf::st_write(sf_pts, path, delete_dsn = TRUE, quiet = TRUE)                          # Write GeoJSON
    return(invisible(TRUE))                                                              # Success
  }
  stop("Failed to write GeoJSON with terra, and sf is not available.")                   # Give up if both fail
}


get_timestamp <- function(file_path) {                                                   # Extract 14-digit timestamp from ATL filename
  sub(".*_(\\d{14})_.*", "\\1", basename(file_path))                                     # Regex extract YYYYMMDDHHMMSS
}


#=============================================================================
# 8) AOI & ICESat-2 Granule Discovery/Download
#=============================================================================

boundary <- sf::st_read("C:\\ICESat2VegR\\01_data\\01_shp\\hurricanemichael_4326.geojson", quiet = TRUE)  # Load AOI GeoJSON

box <- sf::st_bbox(boundary)                                                             # Compute AOI bounding box
box_extent <- sf::st_as_sfc(box)                                                         # Convert bbox to geometry
plot(box_extent, col = NA, border = "blue", main = "AOI and Bounding Box")               # Plot bbox in blue
plot(boundary$geometry, add = TRUE, border = "red", lty = 2, lwd = 2)                    # Overlay AOI boundary in red dashed

lower_left_lon  <- box["xmin"]; lower_left_lat  <- box["ymin"]                           # Lower-left corner lon/lat
upper_right_lon <- box["xmax"]; upper_right_lat <- box["ymax"]                           # Upper-right corner lon/lat
daterange <- c("2018-01-01", "2024-12-31")                                               # Date range for granule search

atl03_granules_find <- ICESat2VegR::ATLAS_dataFinder(                                    # Search ATL03 granules intersecting AOI and dates
  short_name = "ATL03",
  lower_left_lon, lower_left_lat, upper_right_lon, upper_right_lat,
  version = "007", daterange = daterange, persist = TRUE, cloud_computing = FALSE
)
atl08_granules_find <- ICESat2VegR::ATLAS_dataFinder(                                    # Search ATL08 granules intersecting AOI and dates
  short_name = "ATL08",
  lower_left_lon, lower_left_lat, upper_right_lon, upper_right_lat,
  version = "007", daterange = daterange, persist = TRUE, cloud_computing = FALSE
)

atl03_granules_find2 <- sample_granules_by_year(atl03_granules_find, n_per_year = 5, seed = 42)  # Sample up to 5 ATL03 per year
print(head(atl03_granules_find2)); print(table(atl03_granules_find2$year))              # Inspect sampled ATL03 distribution

atl08_granules_find2 <- atl03_granules_find2                                             # Mirror ATL03 selection to ATL08
atl08_granules_find2$url <- gsub("ATL03", "ATL08", atl03_granules_find2$url)             # Replace product code to create ATL08 URLs
print(head(atl08_granules_find2)); print(table(atl08_granules_find2$year))              # Inspect ATL08 URLs and distribution

outdir <- file.path("C:\\ICESat2VegR\\03_data\\02_h5")                                   # Output folder for HDF5 downloads
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)                               # Ensure folder exists
ICESat2VegR::ATLAS_dataDownload(atl03_granules_find2$url, outdir)                        # Download ATL03 granules
ICESat2VegR::ATLAS_dataDownload(atl08_granules_find2$url, outdir)                        # Download ATL08 granules


#=============================================================================
# 9) Aggregation of ATL03/ATL08 into 10 m segments
#=============================================================================

outdir <- file.path("C:\\ICESat2VegR\\01_data\\02_h5")                                   # Input/output folder for HDF5/GeoJSON
atl03_files <- list.files(outdir, pattern = "ATL03.*h5", full.names = TRUE)              # List ATL03 HDF5 files
atl08_files <- list.files(outdir, pattern = "ATL08.*h5", full.names = TRUE)              # List ATL08 HDF5 files

timestamps_atl08 <- sapply(atl08_files, get_timestamp)                                   # Extract timestamps from ATL08 names
timestamps_atl03 <- sapply(atl03_files, get_timestamp)                                   # Extract timestamps from ATL03 names
timestamps_common <- intersect(timestamps_atl08, timestamps_atl03)                       # Intersect timestamps to get paired granules
message("Number of common timestamps: ", length(timestamps_common))                       # Report how many pairs are available

if (!inherits(boundary, "SpatVector")) {                                                 # Ensure AOI is a SpatVector
  if (inherits(boundary, c("sf", "sfc"))) {                                              # If sf, convert to SpatVector
    boundary <- terra::vect(boundary)
  } else {                                                                               # Otherwise error
    stop("`boundary` must be terra::SpatVector or sf/sfc.")
  }
}

seg_10m_metrics <- data.table::data.table()                                              # Accumulator for per-segment metrics

for (ts in timestamps_common) {                                                          # Loop over each common timestamp
  year_ts <- substring(as.character(ts), 1, 4)                                           # Extract year from timestamp
  atl03 <- atl08 <- atl03_clip <- atl08_clip <- NULL                                     # Placeholders for HDF5 handles

  tryCatch({                                                                              # Robust processing per pair
    message("Processing granule: ", ts)                                                  # Log progress

    atl03_file <- atl03_files[grepl(ts, atl03_files)][1L]                                # Find ATL03 file matching timestamp
    atl08_file <- atl08_files[grepl(ts, atl08_files)][1L]                                # Find ATL08 file matching timestamp
    if (!length(atl03_file) || !length(atl08_file)) stop("Matching ATL03/ATL08 not found for ", ts)  # Sanity check

    atl03 <- ICESat2VegR::ATL03_read(atl03_file)                                         # Open ATL03 HDF5
    atl08 <- ICESat2VegR::ATL08_read(atl08_file)                                         # Open ATL08 HDF5

    bb <- as.numeric(sf::st_bbox(terra::vect(boundary)))                                 # AOI bbox as numeric
    clip_region <- terra::ext(bb["xmin"], bb["xmax"], bb["ymin"], bb["ymax"])            # Terra extent from bbox values

    output03 <- tempfile(fileext = ".h5")                                                # Temp file for clipped ATL03
    output08 <- tempfile(fileext = ".h5")                                                # Temp file for clipped ATL08
    atl03_clip <- ICESat2VegR::ATL03_h5_clipBox(atl03, output03, clip_region)            # Clip ATL03 by bbox
    atl08_clip <- ICESat2VegR::ATL08_h5_clipBox(atl08, output08, clip_region)            # Clip ATL08 by bbox

    dt <- ICESat2VegR::ATL03_ATL08_photons_attributes_dt_join(atl03_clip, atl08_clip)    # Join photons with ATL08 attributes
    if (is.null(dt) || nrow(dt) == 0L) stop("Null or empty join result")                 # Stop if nothing to process

    seg10 <- ICESat2VegR::ATL03_ATL08_segment_create(                                    # Create 10 m segments
      dt, 10, centroid = "mean", output = NA, overwrite = FALSE
    )

    stats10 <- ICESat2VegR::ATL03_ATL08_compute_seg_attributes_dt_segStat(               # Compute per-segment metrics
      seg10,
      list(
        h_canopy_gt0   = quantile(ph_h[ph_h > 0 & classed_pc_flag %in% c(2,3)], 0.98),   # 98th percentile canopy height for canopy photons
        n_ground       = sum(classed_pc_flag == 1),                                      # Count ground photons
        n_mid_canopy   = sum(classed_pc_flag == 2),                                      # Count mid-canopy photons
        n_top_canopy   = sum(classed_pc_flag == 3),                                      # Count top-canopy photons
        n_canopy_total = sum(classed_pc_flag >= 2),                                      # Count all canopy photons
        mean_solar     = mean(solar_elevation),                                          # Mean solar elevation
        night_flag2    = as.integer(mean(night_flag) > 0.5)                              # Night flag if majority night
      ),
      ph_class = c(2, 3)                                                                 # Use canopy photons for stats
    )

    stats10_clip <- ICESat2VegR::ATL03_ATL08_seg_attributes_dt_clipGeometry(             # Clip segments to AOI polygon
      stats10, boundary
    )
    if (nrow(stats10_clip) == 0L) {                                                      # Skip if empty after clip
      message("⚠ No valid data after AOI clipping.")
      next
    }

    stats10_clip[, `:=`(lon = longitude, lat = latitude)]                                # Add lon/lat columns used for exports
    stats10_clip <- stats10_clip[h_canopy_gt0 <= 50 & is.finite(lon) & is.finite(lat)]   # QC: realistic height + finite coords
    if (nrow(stats10_clip) == 0L) {                                                      # Skip if empty after QC
      message("⚠ No finite coordinates after filters.")
      next
    }

    stats10_clip$year <- year_ts                                                         # Attach year to segment records
    seg_10m_metrics <- data.table::rbindlist(                                            # Accumulate into master table
      list(seg_10m_metrics, stats10_clip), fill = TRUE
    )

    outdir_sub <- file.path(outdir, "ATL_joined_10")                                     # Per-granule GeoJSON output folder
    dir.create(outdir_sub, showWarnings = FALSE, recursive = TRUE)                       # Ensure folder exists
    f1 <- file.path(outdir_sub, paste0("ATL03_ATL08_", ts, "_10.geojson"))               # Filename for this granule
    safe_write_geojson(stats10_clip, f1, xcol = "lon", ycol = "lat",                     # Write GeoJSON for this granule
                       crs = "EPSG:4326", overwrite = TRUE)

  }, error = function(e) {                                                               # On error during processing
    message("❌ Skipping ", ts, " due to error: ", e$message)                            # Log and skip
  }, finally = {                                                                         # Always attempt to close files
    try({ if (!is.null(atl03)) close(atl03) }, silent = TRUE)                            # Close ATL03
    try({ if (!is.null(atl08)) close(atl08) }, silent = TRUE)                            # Close ATL08
    try({ if (!is.null(atl03_clip)) close(atl03_clip) }, silent = TRUE)                  # Close clipped ATL03
    try({ if (!is.null(atl08_clip)) close(atl08_clip) }, silent = TRUE)                  # Close clipped ATL08
  })
}

all_geojson <- file.path(outdir, "ATL03_ATL08_all_10.geojson")                           # Path for merged GeoJSON
safe_write_geojson(seg_10m_metrics, all_geojson, xcol = "lon", ycol = "lat",             # Write merged GeoJSON
                   crs = "EPSG:4326", overwrite = TRUE)


#=============================================================================
# 10) Sampling & EE-ready geometry
#=============================================================================

all_geojson <- file.path("C:\\ICESat2VegR\\01_data\\02_h5", "ATL03_ATL08_all_10.geojson") # Path to merged segments file
seg_10m_metrics <- sf::st_read(all_geojson, quiet = TRUE)                                  # Read merged segments

head(seg_10m_metrics)

seg_10m_metrics_sub0 <- seg_10m_metrics[                                                   # QC: enough canopy photons & reasonable height
  seg_10m_metrics$n_canopy_total > 5 & seg_10m_metrics$h_canopy_gt0 < 40, ]

seg_10m_metrics_sub <- sf::st_drop_geometry(seg_10m_metrics_sub0)

head(seg_10m_metrics_sub)


class(seg_10m_metrics_sub) <- c("icesat2.atl03_atl08_seg_dt", "data.table") # Preserve class expected by ICESat2VegR

#require(ICESat2VegR)
s <- ICESat2VegR::sample(                                                                 # Spatially spaced sample of segments
  seg_10m_metrics_sub,
  method = ICESat2VegR::spacedSampling(size = 5000, radius = 0.001)
)

plot(s)                                                        # Quick scatter check of sampled points

#s<-seg_10m_metrics_sub[sample(1:nrow(seg_10m_metrics_sub0),1000),]
head(s)
s_vect <- ICESat2VegR::to_vect(s)                                                         # Convert sampled table to SpatVector

plot(s_vect)

outdir_acc <- file.path("C:/ICESat2VegR/03_Results/Accuracy")                             # Output folder for accuracy results
if (!dir.exists(outdir_acc)) dir.create(outdir_acc, recursive = TRUE)                     # Ensure it exists

boundary <- sf::st_read("C:\\ICESat2VegR\\01_data\\01_shp\\hurricanemichael_4326.geojson", quiet = TRUE)  # Load AOI GeoJSON
boundary_ee <- vect_as_ee(boundary)                                   # Convert to EE geometry

# boundary_sf <- sf::st_as_sf(boundary)                                        # Ensure AOI is sf and valid
# boundary_valid <- sf::st_make_valid(boundary_sf)                                          # Fix invalid geometries
# boundary_union <- sf::st_union(boundary_valid)                                            # Dissolve to a single geometry
# boundary_simple <- sf::st_simplify(boundary_union, dTolerance = 0.001)                    # Simplify for EE efficiency
# boundary_ee <- vect_as_ee(boundary_simple)                                   # Convert to EE geometry
geom <- boundary_ee                                                                        # Alias used later
aoi_ee <- geom                                                                             # AOI EE geometry used for mapping/prediction


#=============================================================================
# 11) Build predictors (AlphaEarth + terrain) & extract samples
#=============================================================================

yrs <- sort(unique(s_vect$year))                                                           # Unique years in sampled points
sampling <- vector("list", length(yrs))                                                    # List to collect per-year extracts
names(sampling) <- yrs                                                                     # Name list by year for clarity

for (k in seq_along(yrs)) {                                                                # Loop over years
  yr <- yrs[k]                                                                             # Year value
  sf_part <- s_vect[s_vect$year == yr, ]                                                   # Subset points for this year
  #plot(sf_part)                                                                            # Optional quick plot for QC

  cat(sprintf("Chunk %d/%d - year=%s - n=%d\n",                                            # Progress message
              k, length(yrs), as.character(yr), nrow(sf_part)))

  combined <- Generate_EmbeddingAndAncillaryEEstack(boundary, yr, yr)               # Build EE stack for this year

  extracted <- ICESat2VegR::seg_ancillary_extract(                                  # Extract predictor values at points
    combined, geom = sf_part, scale = 10
  )

  sampling[[k]] <- cbind(year = yr, extracted)                                             # Tag with year and store
}

sampling_df <- dplyr::bind_rows(sampling)                                                  # Bind all yearly extracts
head(sampling_df)                                                                           # Peek at the dataset


#=============================================================================
# 12) Canopy height modeling (feature selection + RF)
#=============================================================================

x <- sampling_df %>%                                                                       # Build feature matrix
  dplyr::select(dplyr::starts_with("A"), elevation, slope, aspect, lon, lat) %>%           # Use AlphaEarth + terrain + lon/lat
  data.frame()                                                                              # Coerce to base data.frame
y <- sampling_df %>%                                                                       # Response vector
  dplyr::select(h_canopy_gt0) %>%                                                           # Canopy height metric
  data.frame()                                                                              # Coerce to base data.frame

nrow(x)
nrow(y)
head(y)

res_f <- varSel(
     xdata     = x,
     ydata     = y$h_canopy_gt0,
     imp.scale = "mir",
     ntree     = 100,
     r         = c(0.25, 0.50, 0.75),
     min.imp =0.2
   )

res_f$selvars
res_f$importance
print(res_f$importance)


#'
## Plot selected-variable importance (most important at top)
plot(res_f, which = "importance")

x <- sampling_df %>% dplyr::select(res_f$selvars)                                          # Restrict features to selected subset
y <- sampling_df$h_canopy_gt0                                                               # Response as numeric vector

fit_k_fold <- fit_model(                                                       # 5-fold CV random forest
  x, y,
  rf_args = list(ntree = 100, mtry = NULL, nodesize = 5, sampsize = NULL),
  test = list(method = "k-fold", k = 5, seed = 42)
)
print(fit_k_fold$stats_train)                                                               # Print cross-validated training stats


#=============================================================================
# 13) Map prediction for 2018 and visualize
#=============================================================================


suppressPackageStartupMessages({                                                       # Suppress startup messages for clean logs
  library(reticulate)                                                                  # Interface to Python
  library(ICESat2VegR)                                                                 # ICESat-2 vegetation tools
  library(leaflet)                                                                     # Interactive maps
  library(sf)                                                                          # Simple Features for vector data
  library(terra)                                                                       # Raster + vector geospatial ops
  library(data.table)                                                                  # Fast data tables
  library(dplyr)                                                                       # Tidy verbs
})
cat("Loading libraries... done.\n")                                                    # Confirm loading

boundary <- sf::st_read("C:\\ICESat2VegR\\01_data\\01_shp\\hurricanemichael_4326.geojson", quiet = TRUE)  # Load AOI GeoJSON

ICESat2VegR_configure()
ee_initialize(project = "ee-carlossilvaengflorestal")                                   # Initialize Earth Engine under a specific GCP project

Generate_EmbeddingAndAncillaryEEstack <- function(                                     # Create an EE stack (AlphaEarth + terrain + optional lon/lat)
  geom,                                                                                # AOI geometry (sf/sfc, SpatVector, or EE geometry)
  startYY, endYY,                                                                      # Inclusive year range
  mask_outside = TRUE,                                                                 # Mask pixels outside AOI
  terrain_scale = 30,                                                                  # Reproject terrain to ~30 m
  multiply_slope_aspect_by10 = TRUE,                                                   # Keep slope/aspect scaling convention
  add_lonlat = TRUE                                                                     # Add lon/lat bands
) {
  ee <- reticulate::import("ee", delay_load = FALSE)                                    # Import Python Earth Engine
  .ee_ping(ee)                                                            # Ping EE session (internal helper)

  startYY <- as.integer(startYY); endYY <- as.integer(endYY)                            # Coerce years to integers
  if (!is.finite(startYY) || !is.finite(endYY) || endYY < startYY)                      # Validate year range
    stop("Invalid years: startYY/endYY")                                                # Abort if invalid

  ee_geom <- ICESat2VegR:::.as_ee_geom(geom)                                            # Convert R geometry to EE geometry

  emb_ic <- ee$ImageCollection("GOOGLE/SATELLITE_EMBEDDING/V1/ANNUAL")$                # AlphaEarth annual embeddings collection
    filterDate(sprintf("%04d-01-01", startYY), sprintf("%04d-12-31", endYY))$           # Filter by year range
    filterBounds(ee_geom)                                                               # Filter by AOI
  embedding <- emb_ic$median()$clip(ee_geom)                                            # Median composite and clip to AOI

  elevation <- tryCatch({                                                               # Try 3DEP 1 m elevation
    ee$ImageCollection("USGS/3DEP/1m")$mosaic()$select("elevation")$clip(ee_geom)
  }, error = function(e) NULL)
  if (is.null(elevation)) {                                                             # Fallback to 3DEP 10 m
    elevation <- tryCatch({
      ee$ImageCollection("USGS/3DEP/10m")$mosaic()$select("elevation")$clip(ee_geom)
    }, error = function(e) NULL)
  }
  if (is.null(elevation)) {                                                             # Fallback to NASADEM
    elevation <- ee$Image("NASA/NASADEM_HGT/001")$select("elevation")$clip(ee_geom)
  }

  elev_proj <- elevation$reproject(crs = "EPSG:4326", scale = terrain_scale)            # Reproject elevation to target scale
  slope  <- ee$Terrain$slope(elev_proj)$rename("slope")                                  # Compute slope
  aspect <- ee$Terrain$aspect(elev_proj)$rename("aspect")                                # Compute aspect
  if (isTRUE(multiply_slope_aspect_by10)) {                                              # Optionally scale values by 10
    slope  <- slope$multiply(10)
    aspect <- aspect$multiply(10)
  }
  slope  <- slope$clip(ee_geom)                                                          # Clip slope to AOI
  aspect <- aspect$clip(ee_geom)                                                         # Clip aspect to AOI

  if (isTRUE(add_lonlat)) {                                                              # Optionally add lon/lat bands
    lonlat <- ee$Image$pixelLonLat()$clip(ee_geom)
    lon    <- lonlat$select("longitude")$rename("lon")
    lat    <- lonlat$select("latitude")$rename("lat")
  }

  stack <- embedding$                                                                    # Start with embedding
    addBands(elevation$rename("elevation"))$                                             # Add elevation
    addBands(slope)$                                                                      # Add slope
    addBands(aspect)                                                                      # Add aspect
  if (isTRUE(add_lonlat)) stack <- stack$addBands(lon)$addBands(lat)                     # Add lon/lat if requested

  if (isTRUE(mask_outside)) {                                                            # Mask outside AOI
    mask <- ee$Image$constant(1)$clip(ee_geom)
    stack <- stack$updateMask(mask)
  }

  stack                                                                                  # Return the final EE image stack
}

writeRDS

saveRDS(fit_k_fold$model, file = "C:/Users/c.silva/Documents/ICESat2VegR/dev/my_model.rds")
getwd()
rfmodel <- readRDS("C:/Users/c.silva/Documents/ICESat2VegR/dev/my_model.rds")



clip_region <- terra::ext(-83.2, -83.18, 32.12, 32.16)                                     # Example extent (min/max lon/lat)
clip_region <- terra::as.polygons(clip_region, crs = "EPSG:4326")                          # Convert extent to polygon for display (optional)

stack_2018 <- Generate_EmbeddingAndAncillaryEEstack(boundary, 2018, 2018)                  # Build predictors for 2018

min_hcanopy <- 0                                                                           # Legend min for map
max_hcanopy <- 20                                                                          # Legend max for map

#s_vect <- ICESat2VegR::to_vect(s)                                                          # Convert sampled points for mapping
#head(s)

forest_height_palette <- c("#ffffff", "#99cc99", "#006600", "#004d00")                     # Simple palette white->green->dark
# palette_colors <- leaflet::colorNumeric(forest_height_palette,                             # Colors for point overlay (plet map)
#                                         range(s$h_canopy_gt0))(s[order(s$h_canopy_gt0), s$h_canopy_gt0])
d <- sf::st_drop_geometry(s)
pal <- leaflet::colorNumeric(forest_height_palette, domain = d$h_canopy_gt0, na.color = "transparent")
palette_colors <- pal(d$h_canopy_gt0)


boundary_ee <- vect_as_ee(boundary)                                   # Convert to EE geometry
geom <- boundary_ee                                                                        # Alias used later
aoi_ee <- geom                                                                             # AOI EE geometry used for mapping/prediction


pred_2018 <- map_create(                                                      # Predict with trained model on EE stack
  model   = rfmodel,
  stack   = stack_2018,
  aoi     = aoi_ee,
  reducer = "mosaic",
  mode    = "auto",
  to_float = TRUE
)


devtools::load_all(getwd())
# or
library(ICESat2VegR)


pal_fun  <- leaflet::colorNumeric(forest_height_palette,                                   # Legend palette for leaflet map
                                  domain = c(min_hcanopy, max_hcanopy))
pal_fun2 <- leaflet::colorNumeric(c("#00441B","#1B7837","#A6DBA0","#E7E1EF","#762A83"),    # Alternative 5-color palette
                                  domain = c(min_hcanopy, max_hcanopy))

centroid <- mean(terra::ext(-83.2, -83.18, 32.12, 32.16))                                  # Quick centroid for initial map view

modelled_map <- leaflet::leaflet() |>                                                       # Start a leaflet map
  ICESat2VegR::addEEImage(                                                                  # Add EE tile from prediction
    pred_2018,
    bands   = "prediction_layer",
    group   = "map",
    min     = min_hcanopy,
    max     = max_hcanopy,
    palette = forest_height_palette
  ) |>
  leaflet::addLegend(                                                                       # Add legend
    pal      = pal_fun,
    values   = c(min_hcanopy, max_hcanopy),
    opacity  = 1,
    title    = "h_canopy",
    position = "bottomleft"
  ) |>
  leaflet::setView(lng = centroid[1], lat = centroid[2], zoom = 12) |>                      # Center map
  leaflet::addLayersControl(                                                                # Add layer control
    overlayGroups = c("map"),
    options = leaflet::layersControlOptions(collapsed = FALSE)
  )

modelled_map                                                                                # Render the map

out_tif <- ICESat2VegR::map_download(                                                       # Export predicted raster to Drive and local temp path
  ee_image = pred_2018, method = "drive",
  region = aoi_ee, scale = 10,
  file_name_prefix = "prediction_2018",
  dsn = file.path(tempdir(),"prediction_2018.tif"),
  drive_folder = "EE_Exports", monitor = TRUE
)
terra::plot(out_tif)                                                                         # Plot local copy (if returned)

modelled_map2 <- leaflet() |>
  ICESat2VegR::addEEImage(                                                                   # Add the EE prediction
    pred_2018,
    bands   = "prediction_layer",
    group   = "map",
    min     = min_hcanopy,
    max     = max_hcanopy,
    palette = colorRampPalette(
      c("#00441B","#1B7837","#A6DBA0","#E7E1EF","#762A83")
    )(128)
  ) |>
  leaflet::addLegend(                                                                        # Add legend
    pal      = pal_fun2,
    values   = c(min_hcanopy, max_hcanopy),
    opacity  = 1,
    title    = "h_canopy",
    position = "bottomleft"
  ) |>
  leaflet::setView(lng = centroid[1], lat = centroid[2], zoom = 12) |>                       # Center map
  leaflet::addLayersControl(                                                                 # Layers control
    overlayGroups = c("map"),
    options = leaflet::layersControlOptions(collapsed = FALSE)
  )

modelled_map2                                                                                # Render the composite map


m <- map_view(list(
  list(type="ee_image", x=pred_2018, bands="prediction_layer", aoi=aoi_ee,
       palette   = colorRampPalette(
         c("#00441B","#1B7837","#A6DBA0","#E7E1EF","#762A83")
       )(128),
       group="Height (m)", min_value=0, max_value=20,
       legend=list(title="Height (m)", auto="quantile")),
  list(type="vector", vect=boundary, group="Site Boundary", color_field="site")
))
m

# # Download (Drive)
out <- map_download(ee_image = pred_2018, method = "drive",
                    region = aoi_ee, scale = 10,
                    file_name_prefix = "prediction_2018",
                    dsn = file.path(tempdir(),"prediction_2018.tif"),
                    drive_folder = "EE_Exports", monitor = TRUE)
terra::plot(out)


atl03_path <- system.file("extdata",
  "atl03_clip.h5",
  package = "ICESat2VegR"
)

atl03_path<-"C:\\Users\\c.silva\\Documents\\ICESat2VegR\\inst\\exdata\\atl03_clip.h5"
# Reading ATL03 data (h5 file)
atl03_h5 <- ATL03_read(atl03_path = atl03_path)

# Extracting ATL03 segment attributes
atl03_segment_dt <- ATL03_seg_metadata_dt(atl03_h5 = atl03_h5)
atl03_dt <- ATL03_photons_attributes_dt(atl03_h5 = atl03_h5)

head(atl03_segment_dt)
head(atl03_dt)

close(atl03_h5)
