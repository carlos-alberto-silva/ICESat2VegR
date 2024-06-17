# h/t to @jimhester and @yihui for this parse block:
# https://github.com/yihui/knitr/blob/dc5ead7bcfc0ebd2789fe99c527c7d91afb3de4a/Makefile#L1-L4
# Note the portability change as suggested in the manual:
# https://cran.r-project.org/doc/manuals/r-release/R-exts.html#Writing-portable-packages
PKGNAME = `sed -n "s/Package: *\([^ ]*\)/\1/p" DESCRIPTION`
PKGVERS = `sed -n "s/Version: *\([^ ]*\)/\1/p" DESCRIPTION`
input_vignettes = $(wildcard vignettes/*.Rmd)
output_vignettes = inst/doc/$(input_vignettes:.Rmd=.html)

all: check

build: install_deps
	R CMD build .

install_deps:
	Rscript \
	-e 'if (!requireNamespace("remotes")) install.packages("remotes")' \
	-e 'remotes::install_deps(dependencies = TRUE)'

install: build
	R CMD INSTALL $(PKGNAME)_$(PKGVERS).tar.gz

clean:
	@rm -rf $(PKGNAME)_$(PKGVERS).tar.gz $(PKGNAME).Rcheck
	./cleanup

vignettes: $(output_vignettes)

preprocess:
	autoreconf
	autoconf --output=configure.win configure.ac
	./cleanup

$(output_vignettes): $(input_vignettes)
	Rscript -e 'knitr::knit("vignettes/v08_GeeModelling.Rmd.orig", output = "vignettes/v08_GeeModelling.Rmd")'
	Rscript -e 'devtools::build_vignettes()'
	rm doc/*.Rmd
	rm doc/*.R
	rm -rf inst/doc
	mv doc inst/doc
