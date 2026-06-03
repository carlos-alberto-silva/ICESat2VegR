test_that("README local workflow runs end-to-end with bundled fixtures", {
  skip_if_not_installed("sf")
  skip_if_not_installed("terra")

  outdir <- withr::local_tempdir()
  atl03_path <- system.file("extdata", "atl03_clip.h5", package = "ICESat2VegR")
  atl08_path <- system.file("extdata", "atl08_clip.h5", package = "ICESat2VegR")
  aoi_path <- system.file("extdata", "aoi_4326.geojson", package = "ICESat2VegR")
  skip_if(
    atl03_path == "" || atl08_path == "" || aoi_path == "",
    "Bundled README fixtures are unavailable."
  )

  file.copy(atl03_path, file.path(outdir, "ATL03_20250801000000_fixture.h5"))
  file.copy(atl08_path, file.path(outdir, "ATL08_20250801000000_fixture.h5"))

  boundary <- sf::st_read(aoi_path, quiet = TRUE)
  box <- sf::st_bbox(boundary)
  lower_left_lon <- box["xmin"]
  lower_left_lat <- box["ymin"]
  upper_right_lon <- box["xmax"]
  upper_right_lat <- box["ymax"]
  daterange <- c("2025-08-01", "2025-08-31")

  expect_true(lower_left_lon < upper_right_lon)
  expect_true(lower_left_lat < upper_right_lat)
  expect_equal(daterange, c("2025-08-01", "2025-08-31"))

  atl03_files <- list.files(outdir, pattern = "ATL03.*h5", full.names = TRUE)
  atl08_files <- list.files(outdir, pattern = "ATL08.*h5", full.names = TRUE)
  get_timestamp <- function(f) sub(".*_(\\d{14})_.*", "\\1", basename(f))
  timestamps_common <- intersect(
    vapply(atl08_files, get_timestamp, character(1)),
    vapply(atl03_files, get_timestamp, character(1))
  )

  expect_length(timestamps_common, 1)

  ts <- timestamps_common[1]
  atl03 <- ATL03_read(atl03_files[grepl(ts, atl03_files)])
  atl08 <- ATL08_read(atl08_files[grepl(ts, atl08_files)])
  on.exit({ close(atl03); close(atl08) }, add = TRUE)

  joined <- ATL03_ATL08_photons_attributes_dt_join(atl03, atl08)
  seg20 <- ATL03_ATL08_segment_create(
    joined,
    20,
    centroid = "mean",
    output = NA,
    overwrite = FALSE
  )
  stats20 <- ATL03_ATL08_compute_seg_attributes_dt_segStat(
    seg20,
    list(
      rh98           = stats::quantile(
        ph_h[ph_h > 0 & classed_pc_flag %in% c(2, 3)],
        0.98,
        na.rm = TRUE
      ),
      n_ground       = sum(classed_pc_flag == 1),
      n_top_canopy   = sum(classed_pc_flag == 3),
      n_canopy_total = sum(classed_pc_flag >= 2),
      mean_solar     = mean(solar_elevation, na.rm = TRUE),
      night_flag2    = as.integer(mean(night_flag, na.rm = TRUE) > 0.5)
    ),
    ph_class = c(2, 3)
  )
  stats20_vect <- to_vect(stats20)

  expect_s4_class(stats20_vect, "SpatVector")
  expect_true(all(c("rh98", "n_canopy_total", "mean_solar") %in% names(stats20)))
  expect_gt(nrow(stats20), 0)

  stats20_clip <- suppressWarnings(
    ATL03_ATL08_seg_attributes_dt_clipGeometry(stats20, boundary)
  )
  if (nrow(stats20_clip) == 0 || !all(c("longitude", "latitude") %in% names(stats20_clip))) {
    stats20_clip <- stats20
  }
  stats20_clip$year <- as.integer(substr(ts, 1, 4))
  stats20_clip <- stats20_clip[stats20_clip$rh98 <= 50 | is.na(stats20_clip$rh98)]
  stats20_clip_sf <- sf::st_as_sf(
    as.data.frame(stats20_clip),
    coords = c("longitude", "latitude"),
    crs = 4326,
    remove = FALSE
  )
  geojson_path <- file.path(outdir, paste0("ATL03_ATL08_", ts, "_20.geojson"))
  sf::st_write(stats20_clip_sf, geojson_path, delete_dsn = TRUE, quiet = TRUE)

  expect_true(file.exists(geojson_path))
  expect_gt(file.info(geojson_path)$size, 0)
})

test_that("README modelling phase works with bundled example segments and synthetic predictors", {
  skip_if_not_installed("randomForest")

  seg_path <- system.file(
    "extdata",
    "ATL03_ATL08_example_segments.geojson",
    package = "ICESat2VegR"
  )
  skip_if(seg_path == "", "Bundled example segment fixture is unavailable.")

  data_raw <- sf::read_sf(seg_path, quiet = TRUE)
  df_dt <- data.table::as.data.table(sf::st_drop_geometry(data_raw))
  class(df_dt) <- c("icesat2.atl03_atl08_seg_dt", "data.table", "data.frame")

  set.seed(1)
  df_sampled <- ICESat2VegR::sample(
    df_dt,
    method = randomSampling(min(50, nrow(df_dt)))
  )
  sampling_df <- data.frame(
    rh98 = df_sampled$rh98,
    A00 = df_sampled$longitude,
    A20 = df_sampled$latitude,
    A40 = df_sampled$n_canopy_total,
    slope = seq_len(nrow(df_sampled)) / nrow(df_sampled),
    aspect = rev(seq_len(nrow(df_sampled))) / nrow(df_sampled),
    elevation = df_sampled$rh98 + 100
  )
  sampling_df <- stats::na.omit(sampling_df)
  skip_if(nrow(sampling_df) < 10, "Not enough complete rows for modelling E2E.")

  x <- sampling_df[, c("A00", "A20", "A40", "slope", "aspect", "elevation")]
  y <- sampling_df$rh98

  fit_rf <- fit_model(
    x = x,
    y = y,
    rf_args = list(ntree = 20),
    test = list(method = "split", test_size = 0.3, seed = 42),
    verbose = FALSE,
    list_test_models = FALSE
  )
  pred_test <- stats::predict(fit_rf$model, newdata = x[fit_rf$test_index, , drop = FALSE])
  results_test <- fit_metrics(y[fit_rf$test_index], as.numeric(pred_test))

  expect_s3_class(fit_rf$model, "randomForest")
  expect_true(all(c("rmse", "mae", "bias", "adj_r2") %in% results_test$stat))
  expect_equal(nrow(results_test), 8)
})
