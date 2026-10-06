ICESat2VegR 0.0.4 corrects the allocation/deallocation mismatch reported by
gcc-ASAN and Valgrind for version 0.0.3. ANNIndex arrays allocated with new[]
are now released with delete[]. The fixed-radius search now receives the square
of the requested distance, as required by ANN. Regression tests cover index
creation, search, finalization, and spaced-sampling distance behavior.

The PROJ configure test now skips execution when cross-compiling; native builds
continue to run the test, and the package still links to PROJ. There are no new
features or dependencies.

Checks for this source:
- R 4.6.1 on Windows: R CMD check --as-cran --no-manual, Status: OK, including
  regular examples, --run-donttest examples, and tests.
- R-hub R-devel: Ubuntu GCC 12, Ubuntu Clang, GCC 16, GCC-ASAN, macOS Intel,
  and macOS arm64 all passed.
- The preceding 0.0.4 archive, before the configure-only change, passed
  WinBuilder R-devel and MacBuilder arm64.
