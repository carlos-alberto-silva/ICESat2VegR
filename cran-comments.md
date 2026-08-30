## Resubmission (hdf5r dependency decoupled)

This is a resubmission. The previous submission (version 0.0.1, submitted
2026-08-26) failed the incoming pretest on the Debian flavor with:

    Package required but not available: 'hdf5r'

This was not an issue in ICESat2VegR: `hdf5r` itself was (and, per
<https://cran.r-project.org/web/checks/check_results_hdf5r.html>, still is)
failing to install on CRAN's Debian check machines
(r-devel-linux-x86_64-debian-{gcc,clang}), because its C wrapper is
incompatible with the HDF5 2.1.0 build present on those machines
(`Wrapper_auto_H5FDfamily.c:23:22: error: implicit declaration of function
'H5FD_family_init'`). Uwe Ligges asked us to pursue a fix with the `hdf5r`
maintainer, but since that is a third-party timeline outside our control,
we instead removed the hard dependency:

* `hdf5r` has been moved from `Imports` to `Suggests`. It is only needed
  for reading/writing *local* HDF5 (`.h5`) files; the cloud-streaming read
  path (`h5py` via `reticulate`) does not use it at all.
* All internal calls to `hdf5r` are namespace-qualified (`hdf5r::...`); the
  blanket `import(hdf5r)` and the one remaining unqualified call have been
  removed from `NAMESPACE`/`R/`.
* Every entry point that actually opens a local `.h5` file
  (`ATL03_read()`/`ATL08_read()` via `ICESat2.h5_local$initialize()`, and
  the two `predict_h5()` methods) now calls a small internal
  `check_hdf5r()` guard first, which raises a clear, actionable error
  (`install.packages("hdf5r")`) if the package is missing, instead of a
  cryptic namespace-not-found error.
* Every `@examples` block that exercises the local-HDF5 path is wrapped in
  `if (requireNamespace("hdf5r", quietly = TRUE)) { ... }`, and the two
  `testthat` tests that read the bundled local fixtures gained
  `skip_if_not_installed("hdf5r")`.
* In practice this affects almost no real users: the HDF5 2.1.0
  incompatibility is specific to CRAN's own bleeding-edge Debian check
  machines. Windows/macOS CRAN binaries and ordinary Linux installs of
  `hdf5r` continue to work fine, so a normal `install.packages("hdf5r")`
  after installing ICESat2VegR restores full functionality. We plan to
  move `hdf5r` back to `Imports` once the upstream fix lands.

This resubmission was additionally verified against real data: real NASA
Earthdata downloads and cloud-streamed reads, and a full live Google Earth
Engine workflow (AlphaEarth embeddings -> Random Forest -> wall-to-wall
canopy-height prediction), not only the bundled example fixtures. That
process turned up and fixed two more small, unrelated bugs:

* `ee_check_task_status()` (and the shared `ee_monitoring()` helper) now
  validates that `task` is a live `ee.batch.Task` object, raising a clean
  R-level error instead of a raw Python `AttributeError` when given
  something else (e.g. a bare task-id string).
* `ICESat2.h5_local`'s exit finalizer used to call `close_all()`
  unconditionally, even for transient wrappers created by navigating into
  a sub-group/dataset. Since `close_all()` closes every open object in the
  whole file, garbage-collecting one of those short-lived wrappers at an
  unpredictable time could silently invalidate a *different*, still-live
  handle to the same file. Only the wrapper that actually opened the file
  now calls `close_all()`; sub-group/dataset wrappers only close their own
  identifier.

## Earlier resubmission (version 0.0.6 comments, already addressed)

The submission prior to the one above (version 0.0.6) received the
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

* Windows 11 x64, R 4.6.1, local `R CMD check --as-cran` (including
  `--run-donttest`): 0 errors, 0 warnings, 0 notes.
* Full test suite (`testthat`, 80 tests) passes with `hdf5r` installed.
* All exported functions additionally exercised against live services:
  real NASA Earthdata downloads/cloud reads and a full Google Earth Engine
  workflow, not only the bundled example fixtures.

## Other notes

The words ATL, ICESat, and Geolocated are domain-specific terms related to
NASA ICESat-2 products.
