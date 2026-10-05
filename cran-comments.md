## ICESat2VegR 0.0.4 - memory-safety correction

The gcc-ASAN check of version 0.0.3 reported an allocation/deallocation
mismatch in `ANNIndex::~ANNIndex()` while running the `spacedSampling` example.
Both index and distance arrays were allocated with `new[]` and released with
`delete`. They are now released with `delete[]`. The Valgrind example log
reported the same two mismatches. The ANN fixed-radius search also expected a
squared radius but was given the requested distance; it now receives the square
of that distance. Regression tests cover index creation, search, finalization,
and spaced-sampling distance behavior. No new features or dependencies are
included.

Version 0.0.3 was otherwise `Status: OK` in CRAN's regular checks on Linux,
Windows, and macOS. The Valgrind check itself was also `Status: OK`; its
detailed example output exposed the mismatched frees addressed here.
