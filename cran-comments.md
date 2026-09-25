## ICESat2VegR 0.0.3 - resubmission

The example failures on the Linux and macOS flavours are fixed. Version 0.0.2
registered an R finalizer on the object wrapping an opened HDF5 file, and that
finalizer called `hdf5r`'s `H5File$close_all()`, which closes every identifier
of the file. When the garbage collector ran it while another wrapper was still
reading from the same file, the check stopped with `id is invalid`. No wrapper
finalizer is registered any more; only an explicit `close()` releases a file.
No new features and no new dependencies.

Verified with the complete example set and the test suite (Windows, R 4.6.1),
including runs with forced garbage collection, and with R-hub checks on Fedora
Linux (gcc 16), Ubuntu Linux and Windows (R-devel), all `Status: OK`.
