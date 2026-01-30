# ============================================================
# All-in-one: EE + Leaflet + Exports
# ============================================================

# ---------- tiny utilities ----------
#' Ensure a required namespace is available
#' @param pkg Package name (character).
#' @keywords internal
.must_have <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE))
    stop("Package '", pkg, "' is required.")
}

#' Null-coalescing operator
#' @param a,b Values where `a` is returned unless it is NULL, otherwise `b`.
#' @keywords internal
`%||%` <- function(a, b) if (!is.null(a)) a else b

#' Random identifier
#' @param prefix Character prefix for the id.
#' @return Character scalar.
#' @keywords internal
.rand_id <- function(prefix = "id") paste0(prefix, "_", paste(sample(c(letters, LETTERS, 0:9), 8, TRUE), collapse = ""))

# ============================================================
# AOI helpers
# ============================================================

#' Convert a terra extent or spatial object to an Earth Engine Rectangle
#'
#' @description
#' Converts a `terra::ext`, `terra::SpatVector`, or `terra::SpatRaster`
#' into a Google Earth Engine geometry of type `ee$Geometry$Rectangle`.
#' The resulting rectangle is always expressed in geographic coordinates
#' (EPSG:4326), regardless of the input object's projection.
#'
#' This function is mainly used internally to standardize spatial inputs
#' before Earth Engine requests (e.g., filtering, clipping, or exporting data).
#'
#' @param x A spatial object of class `terra::ext`, `terra::SpatVector`, or
#'   `terra::SpatRaster`. Its extent is extracted and converted into an
#'   Earth Engine bounding box.
#'
#' @return
#' A Python object representing `ee$Geometry$Rectangle`, suitable for use
#' with Earth Engine Python API functions accessed via `reticulate`.
#'
#' @examples
#' \dontrun{
#' library(terra)
#'
#' # Create an extent over Gainesville, FL
#' bb <- ext(c(-82.4, -82.2, 29.6, 29.8))
#'
#' # Convert to EE geometry
#' aoi <- ext_to_ee(bb)
#' }
#'
#' @export
ext_to_ee <- function(x) {
  .must_have("terra")
  bb <- terra::ext(x)
  ee$Geometry$Rectangle(c(bb$xmin, bb$ymin, bb$xmax, bb$ymax), NULL, FALSE)
}


#' Convert an EE Rectangle to an sf polygon (EPSG:4326)
#'
#' @param aoi An `ee$Geometry$Rectangle`.
#'
#' @return An `sf` polygon with CRS EPSG:4326.
#' @export
ee_rect_to_sf <- function(aoi) {
  .must_have("sf")
  coords <- aoi$coordinates()$getInfo()[[1]]
  mat <- do.call(rbind, lapply(coords, function(p) c(p[[1]], p[[2]])))
  pol <- sf::st_polygon(list(mat))
  sf::st_sf(geometry = sf::st_sfc(pol, crs = 4326))
}

#' Fit a leaflet map to an EE Rectangle AOI
#'
#' @param m A `leaflet` map widget.
#' @param aoi An `ee$Geometry$Rectangle`.
#'
#' @return The input `leaflet` map with bounds set to the AOI.
#' @keywords internal
fit_map_to_aoi <- function(m, aoi) {
  .must_have("leaflet")
  bb <- aoi$coordinates()$getInfo()[[1]]
  xmin <- bb[[1]][[1]]; ymin <- bb[[1]][[2]]
  xmax <- bb[[3]][[1]]; ymax <- bb[[3]][[2]]
  leaflet::fitBounds(m, lng1 = xmin, lat1 = ymin, lng2 = xmax, lat2 = ymax)
}

# ============================================================
# Model → EE prediction image ("prediction_layer")
# ============================================================

#' Create a Prediction Map in Google Earth Engine Using a Fitted Random Forest Model
#'
#' @description
#' Applies a fitted Random Forest model (from R's \code{randomForest} package)
#' to a Google Earth Engine (\strong{EE}) image or image collection and returns
#' an EE image containing predicted values.
#'
#' The model is converted to an Earth Engine estimator using
#' \code{\link{build_ee_forest}}, which internally serializes the R forest
#' through the ICESat2VegR C++ module (\code{icesat2_module}) and constructs an
#' EE \code{Classifier} via \code{ee$Classifier$decisionTreeEnsemble()}.
#'
#' @param model A fitted \code{randomForest::randomForest} model, or a wrapper
#'   object containing an element named \code{model}.
#' @param stack An \code{ee$Image} or \code{ee$ImageCollection} providing the
#'   predictor variables (bands) matching those used to train the model.
#' @param aoi Optional EE geometry. If provided, the output map is clipped to it.
#' @param reducer Aggregation method when \code{stack} is an image collection.
#'   One of: \code{"mosaic"} (default) or \code{"median"}.
#' @param mode Controls how the estimator is applied. Currently supported:
#'   \itemize{
#'     \item \code{"auto"}       — (default) applies the classifier with
#'       \code{img$classify()}.
#'     \item \code{"classifier"} — explicitly applies \code{img$classify()}.
#'   }
#'   The \code{"regressor"} option is not supported, as the current
#'   implementation only builds an Earth Engine classifier.
#' @param to_float Logical; if \code{TRUE} (default) converts the output band to
#'   32-bit float.
#'
#' @details
#' This function works for both:
#' \itemize{
#'   \item EE images: predictors already merged into one raster.
#'   \item EE image collections: predictions computed for each image and reduced
#'     using mosaic or median.
#' }
#'
#' The output band is always named \code{"prediction_layer"}.
#'
#' Properties from a zero-pixel placeholder image are copied to preserve band
#' metadata.
#'
#' @return An Earth Engine \code{ee$Image} containing model predictions.
#'
#' @examples
#' \dontrun{
#'   library(ICESat2VegR)
#'   library(randomForest)
#'
#'   # Fit a model in R
#'   data(airquality)
#'   air <- na.omit(airquality)
#'
#'   rf_model <- randomForest(
#'     x = air[, 2:6],
#'     y = air[, 1],
#'     ntree = 100,
#'     importance = TRUE
#'   )
#'
#'   # Suppose 'stack_ee' is an Earth Engine image with matching bands
#'   pred <- map_create(
#'     model = rf_model,
#'     stack = stack_ee,
#'     aoi   = NULL,
#'     mode  = "auto"
#'   )
#'
#'   # Now pred is: ee$Image with one band ("prediction_layer")
#' }
#'
#' @seealso
#'   \code{\link{build_ee_forest}}
#'   \code{randomForest::randomForest}
#'   Earth Engine API reference: <https://developers.google.com/earth-engine>
#'
#' @export
map_create <- function(model,
                       stack,
                       aoi      = NULL,
                       reducer  = c("mosaic", "median"),
                       mode     = c("auto", "classifier"),
                       to_float = TRUE) {

  reducer <- match.arg(reducer)
  mode    <- match.arg(mode)

  # Unwrap nested model if needed (e.g., k-fold or caret objects)
  if (!is.null(model[["model"]])) {
    model <- model[["model"]]
  }

  # --------------------------------------------------------------------------
  # 1) Convert R randomForest to EE classifier
  # --------------------------------------------------------------------------
  est <- build_ee_forest(model)

  is_py <- inherits(est, "python.builtin.object")
  classes <- if (is_py) paste(class(est), collapse = " ") else ""

  if (!is_py || !grepl("ee.classifier.Classifier", classes, ignore.case = TRUE)) {
    stop(
      "build_ee_forest() must return an ee.Classifier object.\n",
      "Got object of class: ", classes
    )
  }

  # We currently always classify; mode is kept for API compatibility
  apply_model <- function(img) {
    img$classify(est)
  }

  has <- function(o, a) reticulate::py_has_attr(o, a)

  # --------------------------------------------------------------------------
  # 2) Apply to Image or ImageCollection
  # --------------------------------------------------------------------------
  if (has(stack, "visualize") && has(stack, "select")) {
    # ee$Image
    out <- apply_model(stack)

  } else if (has(stack, "map") && has(stack, "mosaic")) {
    # ee$ImageCollection
    ic  <- stack$map(apply_model)
    out <- if (reducer == "mosaic") ic$mosaic() else ic$median()

  } else {
    stop("`stack` must be an ee$Image or ee$ImageCollection.")
  }

  # --------------------------------------------------------------------------
  # 3) Optional AOI, naming, dtype, metadata
  # --------------------------------------------------------------------------
  if (!is.null(aoi)) {
    out <- out$clip(aoi)
  }

  out <- out$rename("prediction_layer")

  if (isTRUE(to_float)) {
    out <- out$toFloat()
  }

  # Copy basic properties (band name, coordinate system, etc.)
  out <- ee$Image(out$copyProperties(ee$Image(0)))

  ee$Image(out)
}

# ============================================================
# Leaflet painters and composer
# ============================================================

#' Register an internal tile-layer prefix per group
#' @param m Leaflet map.
#' @param group Group name.
#' @param prefix Prefix to register.
#' @return Map with attribute updated.
#' @keywords internal
.register_group_prefix <- function(m, group, prefix) {
  reg <- attr(m, "ee_layer_prefixes")
  if (is.null(reg)) reg <- list()
  reg[[group]] <- unique(c(reg[[group]], prefix))
  attr(m, "ee_layer_prefixes") <- reg
  m
}

#' Compute percentiles for a band within an AOI (server-side reduce)
#' @param img ee Image.
#' @param band Band name.
#' @param aoi EE geometry.
#' @param probs Numeric percentiles (0–100).
#' @param scale Optional scale.
#' @param maxPixels Max pixels for reduceRegion.
#' @return Numeric vector of percentiles.
#' @keywords internal
.ee_percentiles <- function(img, band, aoi, probs = c(2, 98), scale = NULL, maxPixels = 1e8) {
  red <- ee$Reducer$percentile(as.list(as.numeric(probs)))
  dct <- img$select(band)$reduceRegion(reducer = red, geometry = aoi, scale = scale, maxPixels = maxPixels)$getInfo()
  as.numeric(unlist(dct, use.names = TRUE))
}

#' Add an EE raster layer to a leaflet map (single tile source)
#' @keywords internal
.add_ee_image_layer <- function(map, x, bands,
                                aoi = NULL, group = NULL,
                                min_value = 0, max_value = 1,
                                palette = c("#00441B","#1B7837","#A6DBA0","#E7E1EF","#762A83"),
                                is_class = FALSE, scale_to_int = TRUE, int_factor = 1000L) {
  .must_have("leaflet"); .must_have("reticulate")
  has <- function(o,a) reticulate::py_has_attr(o,a)
  if (!(has(x,"visualize") && has(x,"select"))) {
    if (has(x,"median")) x <- x$median() else if (has(x,"mosaic")) x <- x$mosaic() else stop("x must be ee.Image or reducible IC")
  }
  avail <- try(x$bandNames()$getInfo(), silent = TRUE)
  if (missing(bands) || is.null(bands)) {
    bands <- if (!inherits(avail, "try-error")) avail else stop("Provide `bands` (1 or 3).")
  } else if (is.numeric(bands)) {
    if (inherits(avail, "try-error")) stop("Numeric `bands` require bandNames().")
    bands <- as.character(avail[bands])
  } else bands <- as.character(bands)

  if (!(length(bands) %in% c(1L,3L))) stop("`bands` must have length 1 or 3.")
  is_rgb <- length(bands) == 3L
  if (is_rgb) is_class <- FALSE
  if (!is_rgb && length(palette) > 256L) palette <- palette[seq_len(256L)]

  img <- x$select(bands)
  if (!is.null(aoi)) img <- img$clip(aoi)
  if (is_rgb) img <- img$toFloat() else if (is_class) img <- img$toUint8() else if (isTRUE(scale_to_int)) img <- img$multiply(as.numeric(int_factor))$round()$toInt16() else img <- img$toFloat()
  img <- ee$Image(img$copyProperties(ee$Image(0)))

  if (!is_rgb && !is_class && isTRUE(scale_to_int)) { vmin <- as.numeric(round(min_value * int_factor)); vmax <- as.numeric(round(max_value * int_factor)) } else { vmin <- min_value; vmax <- max_value }
  vis <- if (is_rgb) img$visualize(min = vmin, max = vmax) else img$visualize(min = vmin, max = vmax, palette = palette)

  prefix <- paste0(group %||% "layer", "_", .rand_id(""))
  url <- vis$getMapId()$tile_fetcher$url_format
  map <- leaflet::addTiles(map, urlTemplate = url, group = group, layerId = paste0(prefix, "_0"), options = leaflet::tileOptions(opacity = 1))
  .register_group_prefix(map, group %||% "layer", prefix)
}

#' Add a tiled EE raster layer to a leaflet map (large AOIs)
#' @keywords internal
.add_ee_image_tiled <- function(map, x, bands,
                                aoi, group = NULL,
                                min_value = 0, max_value = 1,
                                palette = c("#00441B","#1B7837","#A6DBA0","#E7E1EF","#762A83"),
                                is_class = FALSE, scale_to_int = TRUE, int_factor = 1000L,
                                nx = NULL, ny = NULL) {
  .must_have("leaflet"); .must_have("reticulate")
  rect_coords <- function(rect) rect$coordinates()$getInfo()[[1]]
  mk_grid <- function(rect, nx, ny) {
    c0 <- rect_coords(rect)
    xmin <- c0[[1]][[1]]; ymin <- c0[[1]][[2]]
    xmax <- c0[[3]][[1]]; ymax <- c0[[3]][[2]]
    xs <- seq(xmin, xmax, length.out = nx + 1L)
    ys <- seq(ymin, ymax, length.out = ny + 1L)
    out <- vector("list", nx*ny); k <- 0L
    for (i in seq_len(nx)) for (j in seq_len(ny)) { k <- k + 1L; out[[k]] <- ee$Geometry$Rectangle(c(xs[i], ys[j], xs[i+1L], ys[j+1L]), NULL, FALSE) }
    out
  }
  lean_tile <- function(img, bands, aoi_one) {
    out <- img$select(bands); if (!is.null(aoi_one)) out <- out$clip(aoi_one)
    if (length(bands) == 3L) out <- out$toFloat() else if (is_class) out <- out$toUint8() else if (isTRUE(scale_to_int)) out <- out$multiply(as.numeric(int_factor))$round()$toInt16() else out <- out$toFloat()
    ee$Image(out$copyProperties(ee$Image(0)))
  }
  add_one <- function(map, img, prefix, idx) {
    if (length(bands)!=3L && !is_class && isTRUE(scale_to_int)) { vmin <- as.numeric(round(min_value * int_factor)); vmax <- as.numeric(round(max_value * int_factor)) } else { vmin <- min_value; vmax <- max_value }
    vis <- if (length(bands)==3L) img$visualize(min = vmin, max = vmax) else img$visualize(min = vmin, max = vmax, palette = palette)
    url <- vis$getMapId()$tile_fetcher$url_format
    leaflet::addTiles(map, urlTemplate = url, group = group, layerId = paste0(prefix, "_", idx), options = leaflet::tileOptions(opacity = 1))
  }
  grids <- if (!is.null(nx) && !is.null(ny)) list(c(nx,ny)) else lapply(2:12, \(k) c(k,k))
  prefix <- paste0(group %||% "layer", "_", .rand_id(""))
  out <- map; last_err <- NULL
  for (g in grids) {
    tiles <- mk_grid(aoi, g[1], g[2]); out <- map; failed <- FALSE; idx <- 0L
    for (r in tiles) {
      img_tile <- lean_tile(x, bands, r); idx <- idx + 1L
      res <- tryCatch(add_one(out, img_tile, prefix, idx), error = function(e) e)
      if (inherits(res, "error")) { last_err <- res; failed <- TRUE; break } else out <- res
    }
    if (!failed) { out <- .register_group_prefix(out, group %||% "layer", prefix); return(out) }
  }
  stop(last_err)
}

#' Add vector overlays (sf/terra) to a leaflet map with optional categorical legend
#'
#' @param map Leaflet map widget.
#' @param vect `sf` object or `terra::SpatVector`.
#' @param group Overlay group name.
#' @param border_color Border color for features.
#' @param border_weight Border width (pixels).
#' @param fill Logical; fill polygons/markers.
#' @param fill_color Fill color.
#' @param fill_opacity Fill opacity (0–1).
#' @param color_field Optional column used to color features by category.
#' @param palette Optional palette for categories; defaults to `rainbow(n)`.
#' @param legend_title Optional legend title.
#'
#' @return Leaflet map.
#' @keywords internal
.add_vector_overlay <- function(map, vect,
                                group = "overlay",
                                border_color = "#FF3B3B", border_weight = 2,
                                fill = FALSE, fill_color = "#FF3B3B", fill_opacity = 0.2,
                                color_field = NULL, palette = NULL, legend_title = NULL) {
  .must_have("leaflet"); .must_have("sf")
  if (inherits(vect, "SpatVector")) { .must_have("terra"); vect <- sf::st_as_sf(vect) }
  stopifnot(inherits(vect, "sf"))
  if (is.na(sf::st_crs(vect))) sf::st_crs(vect) <- 4326 else if (sf::st_crs(vect)$epsg != 4326) vect <- sf::st_transform(vect, 4326)

  has_cats <- !is.null(color_field) && (color_field %in% names(vect))
  if (has_cats) {
    vals <- as.factor(vect[[color_field]])
    if (is.null(palette)) palette <- grDevices::rainbow(length(levels(vals)))
    pal <- leaflet::colorFactor(palette, domain = levels(vals), ordered = FALSE)

    if (any(sf::st_geometry_type(vect) %in% c("POINT","MULTIPOINT"))) {
      map <- leaflet::addCircleMarkers(map, data = vect,
                                       lng = ~sf::st_coordinates(geometry)[,1], lat = ~sf::st_coordinates(geometry)[,2],
                                       color = ~pal(get(color_field)), weight = border_weight,
                                       fillColor = ~pal(get(color_field)), fillOpacity = fill_opacity, group = group, radius = 5)
    } else if (any(sf::st_geometry_type(vect) %in% c("LINESTRING","MULTILINESTRING"))) {
      map <- leaflet::addPolylines(map, data = vect, color = ~pal(get(color_field)), weight = border_weight, group = group)
    } else {
      map <- leaflet::addPolygons(map, data = vect, color = ~pal(get(color_field)), weight = border_weight,
                                  fill = fill, fillColor = ~pal(get(color_field)), fillOpacity = fill_opacity, group = group)
    }
    if (!is.null(legend_title)) {
      map <- leaflet::addLegend(map, pal = pal, values = levels(vals), title = legend_title, opacity = 1, group = group)
    }
    return(map)
  }

  # single-style
  if (any(sf::st_geometry_type(vect) %in% c("POINT","MULTIPOINT"))) {
    map <- leaflet::addCircleMarkers(map, data = vect,
                                     lng = ~sf::st_coordinates(geometry)[,1], lat = ~sf::st_coordinates(geometry)[,2],
                                     color = border_color, weight = border_weight, fillColor = fill_color, fillOpacity = fill_opacity, group = group, radius = 5)
  } else if (any(sf::st_geometry_type(vect) %in% c("LINESTRING","MULTILINESTRING"))) {
    map <- leaflet::addPolylines(map, data = vect, color = border_color, weight = border_weight, group = group)
  } else {
    map <- leaflet::addPolygons(map, data = vect, color = border_color, weight = border_weight,
                                fill = fill, fillColor = fill_color, fillOpacity = fill_opacity, group = group)
  }
  map
}

#' Attach per-group opacity sliders (client-side)
#'
#' @param m Leaflet map.
#' @param position Control position ("topright","topleft","bottomleft","bottomright").
#'
#' @return Leaflet map with interactive opacity controls.
#' @keywords internal
.attach_opacity_control <- function(m, position = "topright") {
  .must_have("htmlwidgets")
  reg <- attr(m, "ee_layer_prefixes")
  if (is.null(reg) || !length(reg)) return(m)
  groups <- names(reg)
  slider_divs <- paste0(
    vapply(groups, function(g) {
      id <- paste0("op_", g)
      sprintf('<div style="margin:6px 8px;">
                 <label style="display:block;margin-bottom:2px;font:12px/1.2 Arial;">%s opacity</label>
                 <input type="range" id="%s" min="0" max="1" step="0.01" value="1" style="width:140px;">
               </div>', g, id)
    }, character(1L)), collapse = ""
  )
  html <- sprintf('<div class="leaflet-control" style="background:rgba(255,255,255,0.9);padding:6px;border-radius:4px;box-shadow:0 1px 5px rgba(0,0,0,0.4);">%s</div>', slider_divs)

  m <- htmlwidgets::onRender(m, sprintf("
    function(el, x) {
      var map = this;
      var ctl = L.control({position: '%s'});
      ctl.onAdd = function() {
        var div = L.DomUtil.create('div');
        div.innerHTML = `%s`;
        L.DomEvent.disableClickPropagation(div);
        return div;
      };
      ctl.addTo(map);
      var groupPrefixes = %s;
      function setOpacityForGroup(group, val) {
        var prefixes = groupPrefixes[group] || [];
        map.eachLayer(function(l) {
          var id = (l && l.options) ? l.options.layerId : null;
          if (id && l.setOpacity) {
            for (var i=0;i<prefixes.length;i++){
              if (id.indexOf(prefixes[i]) === 0) l.setOpacity(val);
            }
          }
        });
      }
      Object.keys(groupPrefixes).forEach(function(g) {
        var slider = document.getElementById('op_' + g);
        if (slider) slider.addEventListener('input', function(e) { setOpacityForGroup(g, parseFloat(e.target.value)); });
      });
    }",
                                        position, gsub("`", "\\\\`", html), jsonlite::toJSON(reg, auto_unbox = TRUE)))
  m
}

#' Compose a leaflet map from multiple layers (EE rasters and/or vectors)
#'
#' @description High-level map composer capable of adding Earth Engine raster
#' tiles (optionally tiled by AOI grid) and vector overlays (`sf`/`SpatVector`),
#' with layer control, legends, and per-group opacity sliders.
#'
#' @param layers A list of layer specs. For rasters: `list(type="ee_image", x, bands,
#'   aoi=NULL, group=NULL, min_value=0, max_value=1, palette=NULL, is_class=FALSE,
#'   scale_to_int=TRUE, int_factor=1000L, tile=FALSE, nx=NULL, ny=NULL,
#'   legend=list(title=NULL, position=\"bottomright\", opacity=1, auto=NULL, probs=c(2,98), scale=NULL))`.
#'   For vectors: `list(type="vector", vect, group=NULL, border_color, border_weight,
#'   fill, fill_color, fill_opacity, color_field, palette, legend_title)`.
#' @param base_tiles One of `"OSM"`, `"Carto.Light"`, `"Carto.Dark"`.
#' @param add_layers_control Logical; add layers control panel.
#' @param add_opacity_controls Logical; add per-group opacity sliders.
#' @param fit_to Optional EE Rectangle to fit the initial view.
#'
#' @return A `leaflet` htmlwidget.
#' @export
#' @examples
#' \dontrun{
#' m <- map_view(list(
#'   list(type="ee_image", x=img, bands="prediction_layer", aoi=aoi,
#'        group="Prediction", min_value=0, max_value=30,
#'        legend=list(title="Height (m)", auto="quantile")),
#'   list(type="vector", vect=polys_sf, group="Sites", color_field="site")
#' ))
#' }
map_view <- function(layers,
                     base_tiles = c("OSM", "Carto.Light", "Carto.Dark"),
                     add_layers_control = TRUE,
                     add_opacity_controls = TRUE,
                     fit_to = NULL) {
  .must_have("leaflet")
  m <- leaflet::leaflet()
  base_tiles <- match.arg(base_tiles)
  if (base_tiles == "OSM") m <- leaflet::addTiles(m) else
    if (base_tiles == "Carto.Light") m <- leaflet::addProviderTiles(m, "CartoDB.Positron") else
      m <- leaflet::addProviderTiles(m, "CartoDB.DarkMatter")

  overlay_groups <- character(0)
  for (L in layers) {
    stopifnot(is.list(L), !is.null(L$type))
    typ <- L$type
    grp <- L$group %||% paste0("layer_", length(overlay_groups)+1)
    overlay_groups <- unique(c(overlay_groups, grp))

    if (identical(typ, "ee_image")) {
      x <- L$x; bands <- L$bands; aoi <- L$aoi %||% NULL
      tile <- isTRUE(L$tile); nx <- L$nx %||% NULL; ny <- L$ny %||% NULL
      minv <- L$min_value %||% 0; maxv <- L$max_value %||% 1
      pal  <- L$palette %||% c("#00441B","#1B7837","#A6DBA0","#E7E1EF","#762A83")
      is_cls <- isTRUE(L$is_class)
      s2i <- if (is.null(L$scale_to_int)) TRUE else isTRUE(L$scale_to_int); if (is_cls) s2i <- FALSE
      intf <- as.integer(L$int_factor %||% 1000L)

      if (isTRUE(tile) && !is.null(aoi)) {
        m <- .add_ee_image_tiled(m, x, bands, aoi, grp, minv, maxv, pal, is_cls, s2i, intf, nx, ny)
      } else {
        m <- .add_ee_image_layer(m, x, bands, aoi, grp, minv, maxv, pal, is_cls, s2i, intf)
      }

      # continuous legend
      if (!is_cls && length(bands) == 1L && !is.null(L$legend)) {
        lg <- L$legend; pos <- lg$position %||% "bottomright"; op <- lg$opacity %||% 1; ttl <- lg$title %||% ""
        if (!is.null(lg$auto) && lg$auto == "quantile" && !is.null(aoi)) {
          qs <- try(.ee_percentiles(x, bands, aoi, probs = lg$probs %||% c(2,98), scale = lg$scale %||% NULL), silent = TRUE)
          if (!inherits(qs, "try-error") && length(qs) >= 2) {
            mn <- min(qs, na.rm = TRUE); mx <- max(qs, na.rm = TRUE)
            pal_fn <- leaflet::colorNumeric(pal, domain = c(mn, mx))
            m <- leaflet::addLegend(m, pal = pal_fn, values = c(mn, mx), title = ttl, position = pos, opacity = op, group = grp)
          }
        } else {
          pal_fn <- leaflet::colorNumeric(pal, domain = c(minv, maxv))
          m <- leaflet::addLegend(m, pal = pal_fn, values = c(minv, maxv), title = ttl, position = pos, opacity = op, group = grp)
        }
      }

    } else if (identical(typ, "vector")) {
      m <- .add_vector_overlay(
        map = m, vect = L$vect, group = grp,
        border_color = L$border_color %||% "#FF3B3B",
        border_weight= L$border_weight %||% 2,
        fill         = L$fill %||% FALSE,
        fill_color   = L$fill_color %||% "#FF3B3B",
        fill_opacity = L$fill_opacity %||% 0.2,
        color_field  = L$color_field %||% NULL,
        palette      = L$palette %||% NULL,
        legend_title = L$legend_title %||% NULL
      )
    } else stop("Unknown layer type: ", typ)
  }

  if (isTRUE(add_layers_control) && length(overlay_groups)) {
    m <- leaflet::addLayersControl(m, overlayGroups = overlay_groups, options = leaflet::layersControlOptions(collapsed = FALSE))
  }
  if (!is.null(fit_to)) m <- fit_map_to_aoi(m, fit_to)
  if (isTRUE(add_opacity_controls)) m <- .attach_opacity_control(m, "topright")
  m
}

# ============================================================
# EE Export helpers: create tasks
# ============================================================

#' Create an unstarted Drive export task for an EE image
#'
#' @description Thin wrapper over `ee$batch$Export$image$toDrive()` that supports
#' optional time-based prefixes. Returns an unstarted task.
#'
#' @param image An `ee$Image`.
#' @param description Task description.
#' @param folder Drive folder name.
#' @param fileNamePrefix File name prefix; if `timePrefix=TRUE`, the timestamp is appended.
#' @param timePrefix Logical; append a timestamp to `fileNamePrefix`.
#' @param dimensions,region,scale,crs,crsTransform,maxPixels,shardSize,fileDimensions,
#'   skipEmptyTiles,fileFormat,formatOptions Passed to EE export.
#'
#' @return An unstarted EE `Task` (Python object).
#' @keywords internal
ee_image_to_drive <- function(image,
                              description = "myExportImageTask",
                              folder = "EE_Exports",
                              fileNamePrefix = NULL,
                              timePrefix = FALSE,
                              dimensions = NULL,
                              region = NULL,
                              scale = NULL,
                              crs = NULL,
                              crsTransform = NULL,
                              maxPixels = 1e10,
                              shardSize = NULL,
                              fileDimensions = NULL,
                              skipEmptyTiles = NULL,
                              fileFormat = "GeoTIFF",
                              formatOptions = NULL) {
  if (isTRUE(timePrefix)) {
    timePrefix_chr <- gsub("\\s","_", format(Sys.time(), "%Y_%m_%d_%H_%M_%S"))
    fileNamePrefix <- if (is.null(fileNamePrefix)) sprintf("%s_%s", description, timePrefix_chr) else sprintf("%s_%s", fileNamePrefix, timePrefix_chr)
  }
  ee$batch$Export$image$toDrive(
    image = image,
    description = description,
    folder = folder,
    fileNamePrefix = fileNamePrefix,
    dimensions = dimensions,
    region = region,
    scale = scale,
    crs = crs,
    crsTransform = crsTransform,
    maxPixels = maxPixels,
    shardSize = shardSize,
    fileDimensions = fileDimensions,
    skipEmptyTiles = skipEmptyTiles,
    fileFormat = fileFormat,
    formatOptions = formatOptions
  )
}

#' Create an unstarted Cloud Storage export task for an EE image
#'
#' @description Thin wrapper over `ee$batch$Export$image$toCloudStorage()`.
#'
#' @inheritParams ee_image_to_drive
#' @param bucket Cloud Storage bucket name (required).
#'
#' @return An unstarted EE `Task` (Python object).
#' @keywords internal
ee_image_to_gcs <- function(image,
                            description = "myExportImageTask",
                            bucket = NULL,
                            fileNamePrefix = NULL,
                            timePrefix = FALSE,
                            dimensions = NULL,
                            region = NULL,
                            scale = NULL,
                            crs = NULL,
                            crsTransform = NULL,
                            maxPixels = 1e10,
                            shardSize = NULL,
                            fileDimensions = NULL,
                            skipEmptyTiles = NULL,
                            fileFormat = "GeoTIFF",
                            formatOptions = NULL) {
  if (is.null(bucket)) stop("Cloud Storage bucket was not defined")
  if (isTRUE(timePrefix)) {
    timePrefix_chr <- gsub("\\s","_", format(Sys.time(), "%Y_%m_%d_%H_%M_%S"))
    fileNamePrefix <- if (is.null(fileNamePrefix)) sprintf("%s_%s", description, timePrefix_chr) else sprintf("%s_%s", fileNamePrefix, timePrefix_chr)
  }
  ee$batch$Export$image$toCloudStorage(
    image = image,
    description = description,
    bucket = bucket,
    fileNamePrefix = fileNamePrefix,
    dimensions = dimensions,
    region = region,
    scale = scale,
    crs = crs,
    crsTransform = crsTransform,
    maxPixels = maxPixels,
    shardSize = shardSize,
    fileDimensions = fileDimensions,
    skipEmptyTiles = skipEmptyTiles,
    fileFormat = fileFormat,
    formatOptions = formatOptions
  )
}

#' Create an unstarted Asset export task for an EE image
#'
#' @description Thin wrapper over `ee$batch$Export$image$toAsset()`. Optionally
#' deletes an existing asset when `overwrite=TRUE` (requires appropriate perms).
#'
#' @param assetId Destination asset ID (e.g., `"users/me/my_asset"`).
#' @param overwrite Logical; if TRUE, attempt to delete existing `assetId` first.
#' @inheritParams ee_image_to_drive
#'
#' @return An unstarted EE `Task` (Python object).
#' @keywords internal
ee_image_to_asset <- function(image,
                              description = "myExportImageTask",
                              assetId = NULL,
                              overwrite = FALSE,
                              pyramidingPolicy = NULL,
                              dimensions = NULL,
                              region = NULL,
                              scale = NULL,
                              crs = NULL,
                              crsTransform = NULL,
                              maxPixels = 1e10) {
  if (isTRUE(overwrite) && !is.null(assetId)) {
    # best-effort delete (requires ee asset perms)
    try({ ee$data$deleteAsset(assetId) }, silent = TRUE)
  }
  ee$batch$Export$image$toAsset(
    image = image,
    description = description,
    assetId = assetId,
    pyramidingPolicy = pyramidingPolicy,
    dimensions = dimensions,
    region = region,
    scale = scale,
    crs = crs,
    crsTransform = crsTransform,
    maxPixels = maxPixels
  )
}

# ============================================================
# Task control: start safely, monitor, status
# ============================================================

#' Convert an EE Task status to an R list
#'
#' @param task EE `Task` object.
#' @return A named list with task fields (e.g., `state`, `description`, `destination_uris`).
#' @keywords internal
.ee_status_to_list <- function(task) {
  ee <- reticulate::import("ee", delay_load = TRUE)
  st  <- ee$batch$Task$status(task)
  if (inherits(st, "python.builtin.object")) reticulate::py_to_r(st) else st
}

#' Safe list extractor
#' @keywords internal
.lst_get <- function(x, key, default = NULL) if (is.null(x) || is.null(x[[key]])) default else x[[key]]

#' Start an EE Task if not already running
#'
#' @param task EE `Task` object (unstarted or already active).
#'
#' @return Invisibly returns `task`. Messages current state.
#' @keywords internal
ee_task_start_safe <- function(task) {
  st <- .ee_status_to_list(task)
  state <- .lst_get(st, "state", NA_character_)
  if (!is.na(state) && state %in% c("READY","RUNNING","COMPLETED")) {
    message("▶ Task already ", state, ": ", .lst_get(st, "description", "<no description>"))
    return(invisible(task))
  }
  task$start()
  message("▶ Task started: ", .lst_get(st, "description", ""))
  invisible(task)
}

#' Lightweight task monitor (poll-only)
#'
#' @description Polls the task state at fixed intervals until the task leaves
#' `READY`/`RUNNING`, or until `max_attempts` is reached.
#'
#' @param task EE `Task` to monitor.
#' @param task_time Seconds between polls.
#' @param quiet Logical; if `FALSE`, prints state.
#' @param max_attempts Maximum number of polls (use `Inf` to wait indefinitely).
#'
#' @return The final task status as a list.
#' @keywords internal
ee_monitoring <- function(task, task_time = 5, quiet = FALSE, max_attempts = Inf) {
  if (missing(task)) stop("Provide a task (ee$batch$Task).")
  attempts <- 0L
  repeat {
    st <- .ee_status_to_list(task)
    state <- .lst_get(st, "state", NA_character_)
    if (!quiet && !is.na(state)) cat("State:", state, "\n")
    if (!is.na(state) && !(state %in% c("READY","RUNNING"))) break
    attempts <- attempts + 1L
    if (attempts >= max_attempts) break
    Sys.sleep(task_time)
  }
  invisible(.ee_status_to_list(task))
}

#' One-shot task status check
#'
#' @param task EE `Task`.
#' @param quiet Logical; if `TRUE`, suppress printing.
#'
#' @return The task status (list).
#' @keywords internal
ee_check_task_status <- function(task, quiet = TRUE) {
  ee_monitoring(task, task_time = 1, quiet = quiet, max_attempts = 1)
}

# ============================================================
# Drive/GCS fetchers (no rgee)
# ============================================================

#' Extract a sortable timestamp from a googledrive resource row
#' @keywords internal
.gd_get_time <- function(drive_resource_row) {
  if (is.null(drive_resource_row)) return(NA_real_)
  mt <- drive_resource_row$modifiedTime
  ct <- drive_resource_row$createdTime
  ts <- if (!is.null(mt)) mt else ct
  if (is.null(ts)) return(NA_real_)
  as.numeric(as.POSIXct(ts, format = "%Y-%m-%dT%H:%M:%OSZ", tz = "UTC"))
}

#' Order googledrive results by recency
#' @keywords internal
.gd_order_by_time <- function(drib) {
  if (!"drive_resource" %in% names(drib)) return(order(drib$name))
  times_num <- vapply(drib$drive_resource, .gd_get_time, numeric(1))
  if (all(is.na(times_num))) order(drib$name) else order(times_num, decreasing = TRUE, na.last = TRUE)
}

#' Find, wait, and download the most recent Drive file with a given prefix
#'
#' @description After a Drive export task is `COMPLETED`, this utility queries
#' Google Drive for files whose names contain `file_name_prefix`, waits until
#' results are indexed, and downloads the most recent match to `dsn`.
#'
#' @param task EE `Task` (should be in `COMPLETED` state).
#' @param dsn Destination path on disk (GeoTIFF recommended).
#' @param file_name_prefix The export prefix used in the task.
#' @param overwrite Overwrite existing file at `dsn`.
#' @param poll_drive_secs Seconds between Drive searches.
#' @param max_wait_secs Maximum seconds to wait before giving up.
#' @param verbose Logical; print progress messages.
#'
#' @return If `terra` is available, a `SpatRaster`; otherwise the `dsn` path.
#' @keywords internal
ee_drive_fetch_completed <- function(task,
                                     dsn,
                                     file_name_prefix,
                                     overwrite = TRUE,
                                     poll_drive_secs = 10,
                                     max_wait_secs = Inf,
                                     verbose = TRUE) {
  if (!requireNamespace("googledrive", quietly = TRUE))
    stop("Install 'googledrive' and run googledrive::drive_auth().")
  if (missing(dsn) || !nzchar(dsn)) stop("Provide a valid output path in `dsn`.")
  out_dir <- dirname(dsn); if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
  if (file.exists(dsn) && !overwrite) stop("File exists: ", dsn)

  st <- .ee_status_to_list(task)
  state <- .lst_get(st, "state", NA_character_)
  if (!identical(state, "COMPLETED"))
    stop("Task is not COMPLETED (state=", state, ").")

  # Try to parse export folder id from destination_uris
  dest_uris <- .lst_get(st, "destination_uris", NULL)
  folder_id <- NULL
  if (!is.null(dest_uris) && length(dest_uris) >= 1) {
    m <- regexec("folders/([A-Za-z0-9_-]+)", dest_uris[[1]])
    hit <- regmatches(dest_uris[[1]], m)[[1]]
    if (length(hit) >= 2) folder_id <- hit[2]
  }
  if (missing(file_name_prefix) || !nzchar(file_name_prefix))
    stop("Provide the exact `file_name_prefix` used in the export task.")

  started <- Sys.time()
  repeat {
    found <- try({
      if (!is.null(folder_id)) {
        googledrive::drive_find(q = sprintf("name contains '%s' and '%s' in parents", file_name_prefix, folder_id))
      } else {
        googledrive::drive_find(q = sprintf("name contains '%s'", file_name_prefix))
      }
    }, silent = TRUE)
    if (!inherits(found, "try-error") && nrow(found) > 0) break
    if (verbose) message("… waiting for Drive to index files with prefix: ",
                         file_name_prefix,
                         if (!is.null(folder_id)) paste0(" in folder [", folder_id, "]") else " in My Drive")
    if (difftime(Sys.time(), started, units = "secs") > max_wait_secs)
      stop("Timed out waiting for Drive to index the exported file(s).")
    Sys.sleep(as.integer(poll_drive_secs))
  }

  ord <- .gd_order_by_time(found)
  found <- found[ord, , drop = FALSE]
  if (verbose) message("⬇️  Downloading: ", found$name[1], " → ", dsn)
  googledrive::drive_download(found[1, ], path = dsn, overwrite = overwrite)
  if (requireNamespace("terra", quietly = TRUE)) terra::rast(dsn) else dsn
}

#' Alias of `ee_drive_fetch_completed()`
#'
#' @inheritParams ee_drive_fetch_completed
#' @return See `ee_drive_fetch_completed()`.
#' @keywords internal
ee_drive_to_local <- function(task, dsn, file_name_prefix, overwrite = TRUE,
                              poll_drive_secs = 10, max_wait_secs = Inf, verbose = TRUE) {
  ee_drive_fetch_completed(task, dsn, file_name_prefix, overwrite, poll_drive_secs, max_wait_secs, verbose)
}

#' Find, wait, and download the most recent GCS object with a given prefix
#'
#' @description After a GCS export task is `COMPLETED`, this utility lists
#' objects in the target bucket with names starting with `file_name_prefix`,
#' waits for availability, and downloads the largest (by size) to `dsn`.
#'
#' @param bucket Cloud Storage bucket name. If omitted, the function attempts
#'   to infer it from the task `destination_uris`.
#' @inheritParams ee_drive_fetch_completed
#' @param poll_secs Seconds between bucket listings.
#'
#' @return If `terra` is available, a `SpatRaster`; otherwise the `dsn` path.
#' @keywords internal
ee_gcs_fetch_completed <- function(task,
                                   dsn,
                                   file_name_prefix,
                                   bucket = NULL,
                                   overwrite = TRUE,
                                   poll_secs = 10,
                                   max_wait_secs = Inf,
                                   verbose = TRUE) {
  if (!requireNamespace("googleCloudStorageR", quietly = TRUE))
    stop("Install 'googleCloudStorageR' and run googleCloudStorageR::gcs_auth().")
  if (missing(dsn) || !nzchar(dsn)) stop("Provide a valid output path in `dsn`.")
  if (file.exists(dsn) && !overwrite) stop("File exists: ", dsn)
  out_dir <- dirname(dsn); if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

  st <- .ee_status_to_list(task)
  state <- .lst_get(st, "state", NA_character_)
  if (!identical(state, "COMPLETED"))
    stop("Task is not COMPLETED (state=", state, ").")

  # infer bucket/prefix if possible
  if (is.null(bucket)) {
    dest_uris <- .lst_get(st, "destination_uris", NULL)
    if (!is.null(dest_uris) && length(dest_uris) >= 1) {
      m <- regexec("^gs://([^/]+)/(.+)$", dest_uris[[1]])
      hit <- regmatches(dest_uris[[1]], m)[[1]]
      if (length(hit) >= 3) {
        bucket <- hit[2]
        if (missing(file_name_prefix) || !nzchar(file_name_prefix)) file_name_prefix <- hit[3]
      }
    }
  }
  if (is.null(bucket)) stop("Could not infer GCS bucket; please provide `bucket`.")

  started <- Sys.time()
  repeat {
    lst <- try(googleCloudStorageR::gcs_list_objects(bucket = bucket, prefix = file_name_prefix), silent = TRUE)
    if (!inherits(lst, "try-error") && !is.null(lst) && nrow(lst) > 0) break
    if (verbose) message("… waiting for GCS to list objects with prefix: ", file_name_prefix, " in bucket ", bucket)
    if (difftime(Sys.time(), started, units = "secs") > max_wait_secs)
      stop("Timed out waiting for GCS to list the exported file(s).")
    Sys.sleep(as.integer(poll_secs))
  }

  cand <- subset(lst, grepl("\\.tif(f)?$", name, ignore.case = TRUE))
  if (!nrow(cand)) cand <- lst
  cand <- cand[order(cand$size, decreasing = TRUE), , drop = FALSE]
  if (verbose) message("⬇️  Downloading: gs://", bucket, "/", cand$name[1], " → ", dsn)
  googleCloudStorageR::gcs_get_object(object_name = cand$name[1], bucket = bucket, saveToDisk = dsn, overwrite = overwrite)
  if (requireNamespace("terra", quietly = TRUE)) terra::rast(dsn) else dsn
}

#' Alias of `ee_gcs_fetch_completed()`
#'
#' @inheritParams ee_gcs_fetch_completed
#' @return See `ee_gcs_fetch_completed()`.
#' @keywords internal
ee_gcs_to_local <- function(task, dsn, file_name_prefix, bucket = NULL,
                            overwrite = TRUE, poll_secs = 10, max_wait_secs = Inf, verbose = TRUE) {
  ee_gcs_fetch_completed(task, dsn, file_name_prefix, bucket, overwrite, poll_secs, max_wait_secs, verbose)
}

# ============================================================
# Unified wrapper: map_download()
# ============================================================

#' map_download: create task → start → (monitor) → download/return id
#'
#' @description One function to export an EE image to Drive, Cloud Storage, or
#' Asset; optionally monitor the task; and, for Drive/GCS, download the result
#' locally. Returns a local file (or `SpatRaster` if `terra` is available) for
#' Drive/GCS, or the `asset_id` (invisible) for Asset exports.
#'
#' @param ee_image An `ee$Image` to export.
#' @param method One of `"drive"`, `"gcs"`, or `"asset"`.
#' @param region EE geometry/feature collection defining the export region.
#' @param scale Numeric pixel size in meters.
#' @param file_name_prefix Export file prefix (used to search/download).
#' @param dsn Destination path on disk (Drive/GCS methods).
#' @param drive_folder Drive folder name for Drive exports.
#' @param gcs_bucket Cloud Storage bucket name for GCS exports.
#' @param asset_id Destination asset id for Asset exports.
#' @param monitor Logical; if `TRUE`, poll the task until completion.
#' @param task_time Seconds between polls when monitoring.
#' @param ... Additional arguments forwarded to the corresponding export helper
#' (`ee_image_to_drive()`, `ee_image_to_gcs()`, or `ee_image_to_asset()`), e.g.,
#' CRS, pyramiding policy, `maxPixels`, etc.
#'
#' @return For `"drive"`/`"gcs"`, a local path or a `SpatRaster` if `terra` is installed.
#' For `"asset"`, returns `invisible(asset_id)`.
#' @export
#' @examples
#' \dontrun{
#' out <- map_download(
#'   ee_image = img, method = "drive", region = aoi, scale = 10,
#'   file_name_prefix = "my_pred", dsn = "pred.tif",
#'   drive_folder = "EE_Exports", monitor = TRUE
#' )
#' }
map_download <- function(ee_image,
                         method = c("drive","gcs","asset"),
                         region,
                         scale = NULL,
                         file_name_prefix,
                         dsn = NULL,
                         drive_folder = NULL,
                         gcs_bucket = NULL,
                         asset_id = NULL,
                         monitor = TRUE,
                         task_time = 5,
                         ...) {
  method <- match.arg(method)
  if (missing(ee_image)) stop("`ee_image` is required.")
  if (missing(region) && method %in% c("drive","gcs","asset")) stop("`region` is required.")
  if (missing(file_name_prefix) || !nzchar(file_name_prefix)) stop("`file_name_prefix` is required.")
  if (method %in% c("drive","gcs") && (is.null(dsn) || !nzchar(dsn))) stop("For Drive/GCS, provide `dsn`.")
  if (region != NULL) {
    if(inherits(aoi_ee, "ee.featurecollection.FeatureCollection")) {
      region <- region$geometry()$dissolve()
    }
  }

  task <- switch(method,
                 drive = ee_image_to_drive(image = ee_image,
                                           folder = drive_folder %||% "EE_Exports",
                                           region = region, scale = scale,
                                           fileNamePrefix = file_name_prefix, ...),
                 gcs   = ee_image_to_gcs(image = ee_image,
                                         bucket = gcs_bucket,
                                         region = region, scale = scale,
                                         fileNamePrefix = file_name_prefix, ...),
                 asset = ee_image_to_asset(image = ee_image,
                                           assetId = asset_id,
                                           region  = region, scale = scale, ...)
  )

  ee_task_start_safe(task)
  if (isTRUE(monitor)) ee_monitoring(task, task_time = task_time, quiet = FALSE, max_attempts = Inf)

  if (method == "drive") {
    return(ee_drive_to_local(task, dsn, file_name_prefix, overwrite = TRUE,
                             poll_drive_secs = task_time, max_wait_secs = Inf, verbose = TRUE))
  } else if (method == "gcs") {
    return(ee_gcs_to_local(task, dsn, file_name_prefix, bucket = gcs_bucket,
                           overwrite = TRUE, poll_secs = task_time, max_wait_secs = Inf, verbose = TRUE))
  } else {
    invisible(asset_id)
  }
}

  
# =============================================================================
# Sample ATL granule URLs by year
# =============================================================================

#' Sample ICESat-2 ATL granule URLs by year
#'
#' @description
#' Given a vector (or first column of a matrix/data frame) of ICESat-2 ATL03/ATL08
#' granule URLs or file paths, this function:
#'
#' \itemize{
#'   \item Attempts to parse the acquisition year from each URL/path.
#'   \item Groups granules by year.
#'   \item Samples up to \code{n_per_year} unique URLs per year.
#' }
#'
#' This is useful for creating manageable subsets of ATL granules for testing,
#' model fitting, or cross-validation.
#'
#' @param urls Character vector, matrix, or data frame containing granule URLs
#'   or file paths. If a matrix or data frame is provided, the first column is
#'   used.
#' @param n_per_year Integer. Maximum number of URLs to sample per year
#'   (default \code{5}).
#' @param seed Optional integer seed passed to \code{set.seed()} to make the
#'   sampling reproducible. If \code{NULL} (default), the random state is
#'   left unchanged.
#'
#' @details
#' The year is extracted using the following heuristics:
#'
#' \enumerate{
#'   \item A path component of the form \code{/YYYY/} (e.g., \code{".../2020/..."}).
#'   \item A filename pattern of the form \code{"ATL0[38]_YYYYMMDD..."}.
#' }
#'
#' URLs for which no year can be detected are dropped with a warning.
#'
#' @return
#' A data frame with columns:
#' \itemize{
#'   \item \code{year}: integer acquisition year.
#'   \item \code{url}: the sampled URL or file path.
#' }
#' Rows are ordered by \code{year} and then \code{url}.
#'
#' @examples
#' \dontrun{
#'   urls <- c(
#'     "https://example.org/ATL03_20200101000000_001.h5",
#'     "https://example.org/ATL03_20200102000000_002.h5",
#'     "https://example.org/ATL03_20210101000000_003.h5"
#'   )
#'
#'   sample_df <- ee_sample_atl_granules_by_year(urls, n_per_year = 1, seed = 123)
#'   sample_df
#' }
#'
#' @export
sample_ATL_granules_by_year <- function(
    urls,
    n_per_year = 5,
    seed       = NULL
) {
  # Normalize to character vector
  if (is.matrix(urls) || is.data.frame(urls)) {
    urls <- as.character(urls[, 1])
  }

  stopifnot(is.character(urls), length(urls) > 0)

  # Helper to parse year from URL or filename
  get_year <- function(u) {
    # Prefer /YYYY/ in path
    m <- regexpr("/(20\\d{2})/", u, perl = TRUE)
    if (m[1] > 0) {
      yy <- substr(u, m[1] + 1, m[1] + attr(m, "match.length") - 2)
      return(yy)
    }

    # Fallback: ATL0X_YYYYMMDD...
    m2 <- regexpr("ATL0[38]_(20\\d{2})\\d{10}", u, perl = TRUE)
    if (m2[1] > 0) {
      # Extract the four-digit year starting at position m2 + 6
      yy <- substr(u, m2[1] + 6, m2[1] + 9)
      return(yy)
    }

    NA_character_
  }

  years <- vapply(urls, get_year, character(1))
  keep  <- !is.na(years)

  if (!all(keep)) {
    warning("Dropped ", sum(!keep), " URL(s) with no detectable year.")
  }

  urls  <- urls[keep]
  years <- years[keep]

  if (!is.null(seed)) {
    set.seed(seed)
  }

  split_urls <- split(urls, years)

  sampled <- lapply(names(split_urls), function(y) {
    u <- sort(unique(split_urls[[y]]))
    if (length(u) <= n_per_year) {
      data.frame(
        year = as.integer(y),
        url  = u,
        stringsAsFactors = FALSE
      )
    } else {
      pick <- sample(u, n_per_year)
      data.frame(
        year = as.integer(y),
        url  = pick,
        stringsAsFactors = FALSE
      )
    }
  })

  out <- do.call(rbind, sampled)
  rownames(out) <- NULL
  out[order(out$year, out$url), ]
}

# =============================================================================
# Write GeoJSON from a table with lon/lat
# =============================================================================

#' Safely write a GeoJSON file from a lon/lat table
#'
#' @description
#' Converts a data frame or similar tabular object with longitude/latitude
#' columns into a point layer and writes it to disk as a GeoJSON file.
#' The function first attempts to use \pkg{terra} via
#' \code{ICESat2VegR::to_vect()}, and falls back to \pkg{sf} if available.
#'
#' @param dt A data frame, \code{data.table}, or similar object containing at
#'   least two numeric columns for coordinates.
#' @param path Character. Path to the output GeoJSON file.
#' @param xcol,ycol Character. Names of the longitude and latitude columns in
#'   \code{dt}. Defaults are \code{"lon"} and \code{"lat"}.
#' @param crs Character. Coordinate reference system of the input coordinates.
#'   Default is \code{"EPSG:4326"} (longitude/latitude in WGS84).
#' @param overwrite Logical. If \code{TRUE} (default), allow overwriting an
#'   existing file.
#'
#' @return
#' Invisibly returns \code{TRUE} on success. An error is thrown if both the
#' \pkg{terra}-based and \pkg{sf}-based write attempts fail.
#'
#' @examples
#' \dontrun{
#'   pts <- data.frame(
#'     id  = 1:3,
#'     lon = c(-82.35, -82.34, -82.33),
#'     lat = c( 29.65,  29.66,  29.67)
#'   )
#'
#'   write_geojson_safe(pts, "points.geojson")
#' }
#'
#' @export
write_geojson <- function(
    dt,
    path,
    xcol      = "lon",
    ycol      = "lat",
    crs       = "EPSG:4326",
    overwrite = TRUE
) {
  has_sf <- requireNamespace("sf", quietly = TRUE)

  # First attempt: use terra via ICESat2VegR::to_vect
  ok <- try({
    sv <- ICESat2VegR::to_vect(dt, xcol = xcol, ycol = ycol, crs = crs)
    terra::writeVector(sv, path, overwrite = overwrite)
    TRUE
  }, silent = TRUE)

  if (isTRUE(ok)) {
    return(invisible(TRUE))
  }

  # Fallback: use sf if available
  if (has_sf) {
    df <- as.data.frame(dt)
    sf_pts <- sf::st_as_sf(df, coords = c(xcol, ycol), crs = 4326, remove = FALSE)
    sf::st_write(sf_pts, path, delete_dsn = TRUE, quiet = TRUE)
    return(invisible(TRUE))
  }

  stop("Failed to write GeoJSON with terra, and sf is not available.")
}

# =============================================================================
# Extract 14-digit timestamp from ATL filename
# =============================================================================

#' Extract the ICESat-2 timestamp from an ATL filename
#'
#' @description
#' Extracts a 14-digit timestamp (\code{YYYYMMDDHHMMSS}) from an ICESat-2 ATL
#' granule filename or path. Many ATL03/ATL08 filenames include a timestamp
#' segment such as \code{"ATL03_20200101000000_..."}; this helper retrieves
#' that substring.
#'
#' @param file_path Character vector of file paths or URLs. The timestamp is
#'   extracted from the basename of each path.
#'
#' @return
#' A character vector of the same length as \code{file_path}, where each
#' element is either the extracted 14-digit timestamp (\code{"YYYYMMDDHHMMSS"})
#' or, if no match is found, the original string (due to the use of \code{sub()}).
#'
#' @examples
#' \dontrun{
#'
#' # Example ATL03 granule from NSIDC Earthdata Cloud
#' url <- "https://data.nsidc.earthdatacloud.nasa.gov/
#' nsidc-cumulus-prod-protected/ATLAS/ATL03/006/2021/10/02/
#' ATL03_20211002001658_01461302_006_01.h5"
#'
#' # Extract timestamp from filename
#' ts <- atl_extract_timestamp(url)
#' ts
#' # [1] "20211002001658"
#'
#' }
#'
#' @keywords internal
atl_extract_timestamp <- function(file_path) {
  sub(".*_(\\d{14})_.*", "\\1", basename(file_path))
}
