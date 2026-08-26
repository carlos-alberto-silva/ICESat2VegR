## Resubmission

This is a resubmission. The previous submission (version 0.0.6) received the
following comments from the CRAN team, all of which have been addressed:

* Title shortened to under 65 characters.
* Replaced `:::` with `::` for the exported `.as_ee_geom()` helper referenced
  in documentation examples.
* `addEEImage()` is now exported; examples were removed from the two
  functions that are genuinely internal (`atl_extract_timestamp()`,
  `tryInitializeEarthEngine()`) rather than exporting them.
* `\dontrun{}` has been replaced with runnable examples or `\donttest{}`
  wherever the example only needs the bundled sample data (`inst/extdata`).
  `\dontrun{}` is kept only where an example genuinely requires Earth Engine
  authentication, NASA Earthdata Cloud credentials, or launches an
  interactive browser session -- none of which are available in an
  automated check environment. All examples (including `--run-donttest`)
  pass locally; the Earth Engine/Earthdata/Google Drive examples that remain
  `\dontrun{}` were additionally verified by hand against live services.
* Removed commented-out code from examples (`to_vect()` and others).
* Replaced `cat()` with `message()` for progress output in `fit_model.R`.
* `inst/scripts/upscaling_alphaearth_workflow.R` no longer calls
  `install.packages()`/`installed.packages()`, and now restores `par()`
  after changing it.
* Removed the unused `purrr` Import.

### Authors@R / copyright holders

Per the request to credit all authors, contributors, and copyright holders
of code included in the package, `Authors@R` now also lists:

* Sunil Arya and David Mount (`ctb`, `cph`), and the University of Maryland
  (`cph`) -- authors and copyright holder of the ANN (Approximate Nearest
  Neighbor) library bundled and compiled under `src/`.
* Chuck Gantz (`ctb`) -- author of the latitude/longitude to UTM conversion
  algorithm adapted in `R/utmTools.R` (already credited inline in that file).
* Cole Krehbiel (`ctb`) -- author of the NASA data download routine adapted
  in `R/ATLAS_dataDownload.R` (already credited via `@references` in that
  file).

## Test environments

* Windows 11 x64, R 4.5.2, local `R CMD check --as-cran`: 0 errors,
  0 warnings, 0 notes.

## Other notes

The words ATL, ICESat, and Geolocated are domain-specific terms related to
NASA ICESat-2 products.
