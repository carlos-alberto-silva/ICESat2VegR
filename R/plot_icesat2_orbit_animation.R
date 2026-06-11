#' Plot ICESat-2 Orbital Tracks as a 3D Globe Animation
#'
#' @description
#' Creates an interactive HTML animation of ICESat-2 orbital tracks from one or
#' more KML files, or from a \code{terra::SpatVector} returned by
#' \code{rgt_extract()}.
#'
#' @param kml_files Character vector containing paths to one or more KML files.
#' If \code{NULL} and \code{rgt = NULL}, example KML files distributed with the
#' package are used.
#' @param rgt A \code{terra::SpatVector} returned by \code{rgt_extract()}.
#' If provided, this object is used instead of \code{kml_files}.
#' @param output_dir Character. Directory where the HTML animation and required
#' assets will be written. Default is \code{tempdir()}.
#' @param output_file Character. Name of the output HTML file.
#' @param track_speed Numeric. Initial orbital track animation speed.
#' Default is \code{2}.
#' @param earth_rotation_speed Numeric. Initial Earth rotation speed.
#' Default is \code{2}.
#' @param launch Logical. If \code{TRUE}, launches a local web server and opens
#' the animation in the default browser.
#'
#' @return
#' Invisibly returns the path to the generated HTML file.
#'
#' @details
#' The function generates an interactive Three.js 3D visualization of ICESat-2
#' orbital tracks on a rotating Earth. Tracks can be supplied either as KML files
#' or directly as the \code{terra::SpatVector} output from \code{rgt_extract()}.
#'
#' @references
#' NASA ICESat-2 Reference Ground Tracks:
#' \url{https://icesat-2.gsfc.nasa.gov/science/specs}
#'
#' @seealso
#' \code{\link{rgt_extract}},
#' \url{https://icesat-2.gsfc.nasa.gov}
#'
#' @examples
#' \dontrun{
#' plot_icesat2_orbit_animation()
#'
#' rgt_line <- rgt_extract(h5 = atl03_h5, line = TRUE)
#'
#' plot_icesat2_orbit_animation(
#'   rgt = rgt_line,
#'   launch = TRUE
#' )
#'
#' plot_icesat2_orbit_animation(
#'   kml_files = system.file("extdata", "IS2_RGT_0001_cycle33_12-Sep-2026.kml", package="ICESat2VegR"),
#'   output_dir = tempdir(),
#'   launch = TRUE
#' )
#' }
#'
#' @importFrom sf st_as_sf st_cast st_coordinates st_drop_geometry st_layers st_read
#' @importFrom dplyr mutate select arrange desc row_number bind_rows group_split
#' @importFrom purrr map_dfr
#' @importFrom stringr str_extract
#' @importFrom jsonlite toJSON
#' @export
plot_icesat2_orbit_animation <- function(
    kml_files = NULL,
    rgt = NULL,
    output_dir = tempdir(),
    output_file = "ICESat2_orbit_animation.html",
    track_speed = 2,
    earth_rotation_speed = 2,
    launch = TRUE
) {

  requireNamespace("sf")
  requireNamespace("dplyr")
  requireNamespace("purrr")
  requireNamespace("stringr")
  requireNamespace("jsonlite")
  requireNamespace("servr")
  requireNamespace("terra")

  earth_texture <- "Stylized_World_Topo_5400x2700.jpeg"

  extdata_dir <- system.file(
    "extdata",
    package = "ICESat2VegR"
  )

  if (extdata_dir == "") {
    stop("Could not find the package extdata folder.")
  }

  earth_texture_path <- file.path(
    extdata_dir,
    earth_texture
  )

  if (!file.exists(earth_texture_path)) {
    stop("Earth texture not found in inst/extdata: ", earth_texture)
  }

  if (!dir.exists(output_dir)) {
    dir.create(
      output_dir,
      recursive = TRUE
    )
  }

  file.copy(
    earth_texture_path,
    file.path(output_dir, earth_texture),
    overwrite = TRUE
  )

  all_tracks <- list()

  if (!is.null(rgt)) {

    if (!inherits(rgt, "SpatVector")) {
      stop("rgt must be a terra::SpatVector returned by rgt_extract().")
    }

    rgt_sf <- sf::st_as_sf(rgt)
    rgt_sf <- sf::st_cast(rgt_sf, "POINT")

    coords <- sf::st_coordinates(rgt_sf)

    track_df <- data.frame(
      track_id = "rgt_extract",
      track_order = 1L,
      step = seq_len(nrow(coords)),
      time_label = NA_character_,
      lon = coords[, 1],
      lat = coords[, 2]
    )

  } else {

    if (is.null(kml_files)) {
      kml_files <- list.files(
        extdata_dir,
        pattern = "\\.kml$",
        full.names = TRUE
      )
    }

    if (length(kml_files) == 0) {
      stop(
        "No KML files provided, no rgt object supplied, ",
        "and no example KML files found in inst/extdata."
      )
    }

    if (any(!file.exists(kml_files))) {
      stop("Some KML files do not exist.")
    }

    kml_files <- sort(kml_files)

    for (i in seq_along(kml_files)) {

      kf <- kml_files[i]
      layers <- sf::st_layers(kf)$name

      track_pts <- purrr::map_dfr(
        layers,
        function(layer_name) {
          x <- sf::st_read(
            kf,
            layer = layer_name,
            quiet = TRUE
          )

          x$layer_name <- layer_name
          x
        }
      )

      track_pts <- sf::st_cast(
        track_pts,
        "POINT"
      )

      coords <- sf::st_coordinates(track_pts)

      one_track <- track_pts |>
        dplyr::mutate(
          lon = coords[, 1],
          lat = coords[, 2],
          track_id = tools::file_path_sans_ext(
            basename(kf)
          ),
          track_order = i,
          step = dplyr::row_number(),
          time_label = stringr::str_extract(
            layer_name,
            "\\d{2}:\\d{2}:\\d{2}"
          )
        ) |>
        sf::st_drop_geometry() |>
        dplyr::select(
          track_id,
          track_order,
          step,
          time_label,
          lon,
          lat
        )

      if (nrow(one_track) > 2) {

        first_pt <- one_track[1, ]
        last_pt <- one_track[nrow(one_track), ]

        close_dist <- sqrt(
          (first_pt$lon - last_pt$lon)^2 +
            (first_pt$lat - last_pt$lat)^2
        )

        if (close_dist < 0.01) {
          one_track <- one_track[-nrow(one_track), ]
        }
      }

      one_track <- one_track |>
        dplyr::mutate(
          step = dplyr::row_number()
        )

      all_tracks[[i]] <- one_track
    }

    track_df <- dplyr::bind_rows(all_tracks)
  }

  distance_deg <- function(lon1, lat1, lon2, lat2) {
    sqrt((lon1 - lon2)^2 + (lat1 - lat2)^2)
  }

  track_list <- track_df |>
    dplyr::arrange(track_order, step) |>
    dplyr::group_split(track_order)

  ordered_tracks <- list()
  ordered_tracks[[1]] <- track_list[[1]]

  if (length(track_list) > 1) {

    for (i in 2:length(track_list)) {

      previous_track <- ordered_tracks[[i - 1]]
      current_track <- track_list[[i]]

      previous_end <- previous_track[nrow(previous_track), ]
      current_start <- current_track[1, ]
      current_end <- current_track[nrow(current_track), ]

      dist_to_start <- distance_deg(
        previous_end$lon,
        previous_end$lat,
        current_start$lon,
        current_start$lat
      )

      dist_to_end <- distance_deg(
        previous_end$lon,
        previous_end$lat,
        current_end$lon,
        current_end$lat
      )

      if (dist_to_end < dist_to_start) {
        current_track <- current_track |>
          dplyr::arrange(dplyr::desc(step)) |>
          dplyr::mutate(
            step = dplyr::row_number()
          )
      }

      ordered_tracks[[i]] <- current_track
    }
  }

  track_df <- dplyr::bind_rows(ordered_tracks) |>
    dplyr::mutate(
      global_step = dplyr::row_number()
    )

  track_json <- jsonlite::toJSON(
    track_df,
    dataframe = "rows",
    auto_unbox = TRUE
  )

  html_code <- paste0('
<!DOCTYPE html>
<html>
<head>
<meta charset="utf-8">
<title>ICESat-2 Orbit Animation</title>
</head>
<body>
<script>
const trackData = ', track_json, ';
</script>
</body>
</html>
')

  out_file <- file.path(
    output_dir,
    output_file
  )

  writeLines(
    html_code,
    out_file
  )

  if (launch) {
    try(
      servr::daemon_stop(),
      silent = TRUE
    )

    servr::httd(
      dir = output_dir,
      port = 6733
    )

    utils::browseURL(
      paste0(
        "http://127.0.0.1:6733/",
        output_file
      )
    )
  }

  invisible(out_file)
}
