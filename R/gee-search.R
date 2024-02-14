cache <- NULL

# From https://stackoverflow.com/a/51828654/2548351
unescape_html2 <- function(str) {
    html <- paste0("<x>", paste0(str, collapse = "#_|"), "</x>")
    parsed <- xml2::xml_text(xml2::read_html(html))
    strsplit(parsed, "#_|", fixed = TRUE)[[1]]
}

#' @export
search_datasets <- function(..., refresh_cache = FALSE) {
    args <- list(...)

    if (refresh_cache || is.null(cache)) {
        "%>%" <- magrittr::"%>%"
        req <- httr2::request("https://developers.google.com/earth-engine/datasets/catalog") %>% httr2::req_perform()
        txt <- req %>% httr2::resp_body_string()
        ids <- sub("ee-dataset[^a]+a href=\"/earth-engine/datasets/catalog/", "", regmatches(txt, gregexpr('ee-dataset[^a]+a href="(/earth-engine/datasets/catalog/[^"]+)', txt))[[1]])
        title <- sub('a href="/earth-engine/datasets/catalog/.*data-text="', "", regmatches(txt, gregexpr('a href="/earth-engine/datasets/catalog/[^"]+.*?data-text="([^"]+)', txt))[[1]])
        title <- unescape_html2(title)
        description <- sub('ee-dataset-description-snippet">', "", regmatches(txt, gregexpr('ee-dataset-description-snippet">[^<]+', txt))[[1]])
        description <- unescape_html2(description)
        cache <<- data.table::data.table(id = ids, title = title, description = description)
    }
    dt <- cache
    for (arg in args) {
        dt <- dt[grep(tolower(arg), tolower(title))]
    }
    dt
}

#' @export
get_catalog_path <- function(id) {
    if (inherits(id, "data.frame")) {
        id <- id[1]$id
    }
    id2 <- httr2::request(sprintf("https://developers.google.com/earth-engine/datasets/catalog/%s", id)) %>%
        httr2::req_perform() %>%
        httr2::resp_body_string()
    utils::strcapture('<code class="devsite-click-to-copy[^\\(]+\\("([^"]+)', id2, proto = list(character()))[[1]]
}
