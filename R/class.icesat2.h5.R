# ICESat2.h5 <- R6Class("ICESat2.h5", list(
#     h5 = NULL,
#     initialize = function(h5) {
#         self$h5 <- h5
#     },
#     ls = function() {
#         self$h5$ls()
#     },
#     getBeams = function() {
#         groups <- self$h5$ls()
#         grep("gt[1-3][lr]", groups, value = TRUE)
#     },
#     create_group = function(group) {
#         self$h5$create_group(group)
#     },
#     create_attr = function(attribute, data) {
#         self$h5$create_attr(attribute, data)
#     },
#     close_all = function() {
#         self$h5$close_all()
#     }
# ))

# "[[.ICESat2.h5" <- function(x, i) {
#     x$h5[[i]]
# }
