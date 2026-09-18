Dear CRAN Team and Professor Ripley,

Thank you for reporting these issues. We identified and corrected both
problems in ICESat2VegR version 0.0.2.

The intermittent HDF5 `id is invalid` error was caused by delayed finalizers
on non-owning `hdf5r` group and dataset identifiers. These unsafe finalizers
have been removed, and only the wrapper that opens the HDF5 file is now
responsible for closing it.

Package loading also no longer probes or initializes Python. Python
dependencies are resolved only when a Python-backed feature is explicitly
used, avoiding the warnings caused by `/usr/sbin/pypy` on the Fedora check
system.

We verified the correction using the complete test suite, the previously
failing ATL03/ATL08 example, repeated forced-garbage-collection tests, and
live NASA Earthdata and Google Earth Engine tests. All tests passed. No
credentials or authentication files are included in the package.

We are submitting the corrected version 0.0.2 to CRAN.

Best regards,
Carlos Alberto Silva
