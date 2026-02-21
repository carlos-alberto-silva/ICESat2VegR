setwd("C:/Users/c.silva/Documents/ICESat2VegR")

# clean compiled bits in your tree
devtools::clean_dll()

# regenerate docs/namespace if you're using roxygen
Sys.setenv(ROXYGEN_RUN = "true")
devtools::document()
devtools::check(vignettes = FALSE)


# install (no vignettes, to be quick)
devtools::install(build_vignettes = FALSE, upgrade = "never")
require("ICESat2VegR")

#rm(list = c("fit_model"))

devtools::build(args = '--no-build-vignettes')

devtools::build_manual()
system("R CMD Rd2pdf . --output=coloringBookTools_manual.pdf")

# install.packages(c("tinytex","rmarkdown","Rdpack"))
# tinytext::install_tinytext()
#
# tools::texi2pdf(system.file("doc/manual.pdf", package = "stats"))

devtools::build(vignettes=FALSE,manual=TRUE)

require(ICESat2VegR)
atl08_path <- system.file("extdata",
                             "atl08_clip.h5",
                             package = "ICESat2VegR"
                           )
