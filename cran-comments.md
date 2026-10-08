ICESat2VegR 0.0.6 fixes the macOS PROJ database incompatibility reported for
0.0.4. The package no longer bundles inst/proj or copies a build-machine
PROJ data directory during compilation. At runtime it no longer overrides
PROJ's resource search paths, leaving the linked GDAL/PROJ libraries to use
their matching installation data. The GDALOpen example remains enabled, and
a new regression test creates and reopens an EPSG:4326 raster while checking
that the installed package contains no PROJ data directory.

The inappropriate startup greeting has also been removed.

Checks to date:
- Windows R 4.6.1: the corrected 0.0.5 candidate passed a full
  R CMD check --as-cran --no-manual with no ERRORs or WARNINGs. It had three
  NOTEs: one day since the preceding release, Pandoc absent from this local
  terminal's PATH, and one example taking 6.35 seconds on this machine.
  Regular examples, --run-donttest examples, and tests all passed.
- R-hub Ubuntu release and macOS arm64 R-devel passed the PROJ-fixed 0.0.5
  candidate. R-hub GCC 16 Linux and Intel macOS R-devel passed the distinct
  0.0.6 source.
- Verified the 0.0.6 source archive contains no proj.db or inst/proj files
  and all configure/cleanup scripts have LF line endings.
- Official win-builder R-devel checked 0.0.6 with 0 ERRORs, 0 WARNINGs,
  and 1 NOTE: "Days since last update: 2". Installation, shell-script line
  endings, compiled code, examples, tests, and PDF/HTML manuals all passed.
  Check time was 1,034 seconds, so CRAN's separate incoming pre-test may
  still flag the overall check time.
- R-hub's generated Windows Git checkout converted shell scripts to CRLF even
  though the 0.0.6 source archive and its stored Git blobs have LF endings.
  This R-hub-specific WARNING did not occur on official win-builder or in
  the local Windows source-package check.
- A prior MacBuilder arm64 check could not reach package installation because
  its runner lacked several dependencies, including sf and terra. Subsequent
  R-hub macOS arm64 and Intel checks completed successfully.
