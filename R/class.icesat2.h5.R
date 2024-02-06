#' @import R6
#' @export
ICESat2.h5 <- R6Class("ICESat2.h5", list(
    h5 = NULL,
    getBeams = function() {
        groups <- self$h5$ls()
        grep("gt[1-3][lr]", groups, value = TRUE)
    },
    ls = function() {
    }
))
