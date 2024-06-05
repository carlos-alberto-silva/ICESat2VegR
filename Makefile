# h/t to @jimhester and @yihui for this parse block:
# https://github.com/yihui/knitr/blob/dc5ead7bcfc0ebd2789fe99c527c7d91afb3de4a/Makefile#L1-L4
# Note the portability change as suggested in the manual:
# https://cran.r-project.org/doc/manuals/r-release/R-exts.html#Writing-portable-packages
PKGNAME = `sed -n "s/Package: *\([^ ]*\)/\1/p" DESCRIPTION`
PKGVERS = `sed -n "s/Version: *\([^ ]*\)/\1/p" DESCRIPTION`


all: check

build: install_deps
	R CMD build .

check1:
	Rscript -e "zz <- file('check_output.log', 'wt'); sink(zz, type = 'output'); sink(zz, type='message'); devtools::check(); sink(); close(zz)"

check2:
	Rscript -e "zz <- file('check_output2.log', 'wt'); sink(zz, type = 'output'); sink(zz, type='message'); devtools::check(); sink(); close(zz)"

install_deps:
	Rscript \
	-e 'if (!requireNamespace("remotes")) install.packages("remotes")' \
	-e 'remotes::install_deps(dependencies = TRUE)'

install: build
	R CMD INSTALL $(PKGNAME)_$(PKGVERS).tar.gz

clean:
	@rm -rf $(PKGNAME)_$(PKGVERS).tar.gz $(PKGNAME).Rcheck

vignettes: doc

doc: 
	Rscript -e 'devtools::build_vignettes()'
	rm doc/*.Rmd
	rm doc/*.R
	rm -rf inst/doc
	mv doc inst/doc
