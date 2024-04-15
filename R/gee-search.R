#' @export
ee_search_cache <- NULL

# From https://stackoverflow.com/a/51828654/2548351
unescape_html2 <- function(str) {
  html <- paste0("<x>", paste0(str, collapse = "#_|"), "</x>")
  parsed <- xml2::xml_text(xml2::read_html(html))
  strsplit(parsed, "#_|", fixed = TRUE)[[1]]
}

#' Search for Google Earth Engine datasets
#'
#' @param ... character arguments to search for within title and/or description
#' @param title logical. Whether should search within the title, default TRUE.
#' @param description logical. Whether should search within the description of the dataset, default TRUE.
#' @param operator character. Should be either "OR" or "AND" to tell if the search needs to include
#' all the queries "AND" or any of the queries "OR". Default "AND".
#' 
#' @return A `data.table` containing the id, title and description of the datasets that matched
#' the supplied query ordered by relevance.
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

#' Retrieve the Google Earth Engine image catalog id
#' 
#' @param id character. The id retrieved from the data.table resulting from [`search_datasets()`].
#' 
#' @return The catalog id to open within Google Earth Engine.
#' @export
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
