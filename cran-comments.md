## Resubmission (ICESat2VegR 0.0.3)

Dear CRAN Team,

Thank you for running the checks and for reporting the example failures on
`ATL03_ATL08_photons_attributes_dt_join`. This submission fixes the cause of
those failures. The only functional change since 0.0.2 is the HDF5
handle-lifetime correction described below, plus two example/test updates;
there are no new features and no new dependencies.

## Cause of the `id is invalid` failures

Version 0.0.2 registered an R finalizer on the object wrapping an opened
HDF5 file. That finalizer called `hdf5r`'s `H5File$close_all()`, which closes
`every` open identifier belonging to that file, not only the identifier
owned by the wrapper. R runs finalizers at a moment chosen by the garbage
collector, so on the Linux and macOS flavours the finalizer occasionally ran
while another, still-live wrapper was reading from the same file. The checks
then failed with `Error in (function () : id is invalid`, `can't get ID ref
count`, or `Couldn't delete <id>`. The Windows flavours happened to collect
at different moments and reported `OK`, which made the problem look
platform-specific.

## What changed in 0.0.3

* `R/class.icesat2.h5_local.R`: no HDF5 wrapper registers a finalizer any
  more. Only the wrapper that actually opened a file can close it, and only
  through an explicit, idempotent `close()`. Identifiers that become
  unreachable are left to `hdf5r`, which owns them.
* `R/class.icesat2.h5ds_local.R`: dataset wrappers no longer close borrowed
  dataset identifiers.
* `R/predict_h5.R` and the regenerated `man/predict_h5*.Rd`: the
  `predict_h5()` examples now close their input ICESat-2 files instead of
  leaving HDF5 handles open for the remainder of the example session.
* `tests/testthat/test-integration-local.R`: added regression tests that
  abandon never-closed file wrappers, force `gc()`, then reopen the same
  paths and run the previously failing example with forced collection
  enabled (`gctorture2()`).

Package loading also no longer probes or initialises Python (Python
dependencies are resolved only when a Python-backed feature is used), which
removes the `pypy` warnings seen on the Fedora check system.

## Verification

Local (Windows, R 4.6.1, `hdf5r` 1.3.12):

* all examples executed in a single session, in `R CMD check` order: 71
  examples OK, no `id is invalid`, no `Couldn't delete`;
* complete `testthat` suite: passing (the credentialed Earthdata and Google
  Earth Engine tests are opt-in and skipped);
* the previously failing example run with forced garbage collection
  (`gctorture2()` at steps 300 and 3000) and with deliberately abandoned,
  never-closed wrappers: no failures.

R-hub (R-devel, commit `ce59bb8`), results in the submission thread:

* Fedora Linux 44, gcc 16 (r-devel r90586): `Status: OK`, including
  `checking examples ... OK` and `checking examples with --run-donttest ...
  OK` -- the latter runs every example, including the join example, in a
  single session;
* Windows (R-devel): `OK`;
* macOS Intel and macOS arm64 (R-devel): see the run linked in the
  submission thread.

No credentials or authentication files are included in the package.

Kind regards,
Carlos Alberto Silva
