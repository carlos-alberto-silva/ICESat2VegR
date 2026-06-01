ee_search_cache <- NULL

# From https://stackoverflow.com/a/51828654/2548351
unescape_html2 <- function(str) {
  html <- paste0("<x>", paste0(str, collapse = "#_|"), "</x>")
  parsed <- xml2::xml_text(xml2::read_html(html))
  strsplit(parsed, "#_|", fixed = TRUE)[[1]]
}

#' Search for Google Earth Engine datasets
#'
#' @description
#' Searches the Google Earth Engine public data catalog for datasets matching
#' one or more keywords. Results are filtered by matching the query terms
#' against dataset titles and/or descriptions, and are returned ordered by
#' relevance. The catalog is cached locally after the first call to avoid
#' repeated network requests.
#'
#' @param ... One or more character strings to search for within the dataset
#'   title and/or description. Multiple terms are combined using the
#'   \code{operator} argument.
#' @param title logical. Whether to search within dataset titles.
#'   Default is \code{TRUE}.
#' @param description logical. Whether to search within dataset descriptions.
#'   Default is \code{TRUE}.
#' @param operator character. How to combine multiple search terms.
#'   Use \code{"AND"} to require all terms to match, or \code{"OR"} to
#'   return datasets matching any term. Default is \code{"AND"}.
#' @param refresh_cache logical. Whether to refresh the local dataset catalog
#'   cache. Set to \code{TRUE} to force a fresh download of the catalog.
#'   Default is \code{FALSE}.
#'
#' @return A \code{\link[data.table]{data.table}} with three columns:
#'   \describe{
#'     \item{id}{Character. The dataset catalog ID used to load the dataset
#'       in GEE (e.g., \code{"NASA/NASADEM_HGT/001"}).}
#'     \item{title}{Character. The dataset title.}
#'     \item{description}{Character. A brief description of the dataset.}
#'   }
#'   Results are ordered by relevance to the search query.
#'
#' @details
#' The function queries the GEE public data catalog, which is cached locally
#' after the first call. Use \code{refresh_cache = TRUE} to force a fresh
#' download if the catalog may have been updated.
#'
#' The returned \code{id} column can be passed to
#' \code{\link{get_catalog_id}} to extract the GEE-formatted dataset ID
#' for use with \code{ee$Image()} or \code{ee$ImageCollection()}.
#'
#' @examples
#' \dontrun{
#'   Sys.setenv(EE_PROJECT = "your-ee-project-id")
#'   ICESat2VegR::ee_initialize()
#'
#'   # Search for NASA DEM datasets
#'   result <- search_datasets("nasa", "dem")
#'   result
#'
#'   # Search for HLS datasets
#'   hls <- search_datasets("hls")
#'   hls
#'
#'   # Search within title only
#'   landsat <- search_datasets("landsat",
#'     title       = TRUE,
#'     description = FALSE
#'   )
#'   landsat
#'
#'   # Search with OR operator
#'   dem_or_elevation <- search_datasets(
#'     "dem", "elevation",
#'     operator = "OR"
#'   )
#'   dem_or_elevation
#'
#'   # Get the catalog ID for use in GEE
#'   dem_id <- get_catalog_id(result)
#'   dem_id
#'
#'   # Load the dataset in GEE
#'   elevation <- ee$Image(dem_id)
#'   elevation
#' }
#'
#' @export
search_datasets <- function(..., title = TRUE, description = TRUE, operator = "and", refresh_cache = FALSE) {
  args <- list(...)
  stopifnot("Operator should be either 'AND' or 'OR'" = tolower(operator) %in% c("and", "or"))
  stopifnot(
    "At least one of `title` or `description` should be TRUE" =
      (inherits(title, "logical") && title[1]) ||
        (inherits(description, "logical") && description[1])
  )

  if (refresh_cache || is.null(ee_cache$search)) {
    "%>%" <- magrittr::"%>%"
    req <- httr2::request("https://developers.google.com/earth-engine/datasets/catalog") %>% httr2::req_perform()
    txt <- req %>% httr2::resp_body_string()

    ids <- sub(
      "ee-dataset[^a]+a href=\"/earth-engine/datasets/catalog/",
      "",
      regmatches(txt, gregexpr('ee-dataset[^a]+a href="(/earth-engine/datasets/catalog/[^"]+)', txt))[[1]]
    )
    the_title <- sub(
      'a href="/earth-engine/datasets/catalog/.*data-text="',
      "",
      regmatches(txt, gregexpr('a href="/earth-engine/datasets/catalog/[^"]+.*?data-text="([^"]+)', txt))[[1]]
    )
    the_title <- unescape_html2(the_title)
    the_description <- sub(
      'ee-dataset-description-snippet">',
      "",
      regmatches(txt, gregexpr('ee-dataset-description-snippet">[^<]+', txt))[[1]]
    )
    the_description <- unescape_html2(the_description)

    dt <- data.table::data.table(id = ids, title = the_title, description = the_description)
    ee_cache$search <- dt
  }
  dt <- ee_cache$search

  # Search params
  all_idx <- list()
  for (arg in args) {
    idx <- numeric()
    if (title) {
      idx <- grep(tolower(arg), tolower(dt$title))
    }
    if (description) {
      idx <- union(idx, grep(tolower(arg), tolower(dt$description)))
    }
    all_idx[[""]] <- list(idx = idx)
  }

  n <- i <- NA
  .N <- data.table::.N
  dt_idx <- data.table::rbindlist(all_idx)[
    , list(n = .N, i = min(.I)),
    by = c("idx")
  ][order(-n, i)]

  if (tolower(operator) == "and") {
    dt_idx <- dt_idx[n == length(args)]
  }
  dt[dt_idx$idx]
}

#' Retrieve the Google Earth Engine image catalog ID
#'
#' @description
#' Converts a dataset ID from the \code{\link{search_datasets}} result
#' into the GEE-formatted catalog ID required by \code{ee$Image()} or
#' \code{ee$ImageCollection()}. The function converts the underscore-separated
#' format returned by \code{\link{search_datasets}} into the slash-separated
#' format used by GEE (e.g., \code{"NASA_NASADEM_HGT_001"} becomes
#' \code{"NASA/NASADEM_HGT/001"}).
#'
#' @param id Either a \code{\link[data.table]{data.table}} returned by
#'   \code{\link{search_datasets}}, or a character string containing the
#'   dataset ID in the underscore-separated format.
#'
#' @return A character string containing the GEE-formatted catalog ID
#'   (slash-separated) ready for use with \code{ee$Image()} or
#'   \code{ee$ImageCollection()}.
#'
#' @details
#' This function is typically used together with \code{\link{search_datasets}}
#' to find and load GEE datasets:
#'
#' \enumerate{
#'   \item Use \code{\link{search_datasets}} to find the dataset.
#'   \item Pass the result to \code{get_catalog_id} to get the GEE ID.
#'   \item Use the ID with \code{ee$Image()} or \code{ee$ImageCollection()}.
#' }
#'
#' @examples
#' \dontrun{
#'   Sys.setenv(EE_PROJECT = "your-ee-project-id")
#'   ICESat2VegR::ee_initialize()
#'
#'   # Search and retrieve catalog ID for NASA DEM
#'   result <- search_datasets("nasa", "dem")
#'   dem_id <- get_catalog_id(result)
#'   dem_id
#'
#'   # Load the dataset in GEE
#'   elevation <- ee$Image(dem_id)
#'   elevation
#'
#'   # Search and retrieve catalog ID for HLS
#'   hls    <- search_datasets("hls")
#'   hls_id <- get_catalog_id(hls)
#'   hls_id
#'
#'   # Load as ImageCollection
#'   hls_ic <- ee$ImageCollection(hls_id)
#'   hls_ic
#' }
#'
#' @export
get_catalog_id <- function(id)
get_catalog_id <- function(id) {
  "%>%" <- magrittr::"%>%"
  if (inherits(id, "data.frame")) {
    id <- id[1]$id
  }
  id2 <- httr2::request(sprintf("https://developers.google.com/earth-engine/datasets/catalog/%s", id)) %>%
    httr2::req_perform() %>%
    httr2::resp_body_string()
  utils::strcapture('<code class="devsite-click-to-copy[^\\(]+\\("([^"]+)', id2, proto = list(character()))[[1]]
}
