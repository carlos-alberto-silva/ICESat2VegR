root <- file.path(Sys.getenv("RUNNER_TEMP"), "icesat2-fresh-python")
dir.create(root, recursive = TRUE, showWarnings = FALSE)
Sys.unsetenv(c("RETICULATE_PYTHON", "RETICULATE_PYTHON_ENV",
  "RETICULATE_PYTHON_FALLBACK", "CONDA_DEFAULT_ENV", "VIRTUAL_ENV"))
Sys.setenv(R_USER_CACHE_DIR = file.path(root, "cache"),
  R_USER_DATA_DIR = file.path(root, "data"),
  RETICULATE_USE_MANAGED_VENV = "yes", RETICULATE_UV = "managed",
  RETICULATE_CONDA = file.path(root, "absent-conda"))
options(timeout = 600)
.libPaths(c(normalizePath("checks/ICESat2VegR.Rcheck"), .libPaths()))
library(ICESat2VegR)
stopifnot(as.character(packageVersion("ICESat2VegR")) == "0.0.7")
stopifnot(!reticulate::py_available(initialize = FALSE))
stopifnot(ICESat2VegR_configure(auto_restart = FALSE,
  auto_install_miniconda = FALSE, quick_probe = character()))
modules <- c("numpy", "h5py", "packaging", "ee", "earthaccess",
  "googleapiclient", "boto3", "botocore", "s3transfer")
stopifnot(all(vapply(modules, reticulate::py_module_available, logical(1))))
cfg <- reticulate::py_config()
print(cfg)
stopifnot(startsWith(normalizePath(cfg$python, winslash = "/"),
  normalizePath(root, winslash = "/")))
reticulate::py_run_string("import h5py, numpy as np, tempfile, os\np = os.path.join(tempfile.gettempdir(), 'icesat2-configure-smoke.h5')\nwith h5py.File(p, 'w') as f: f['values'] = np.arange(3)\nwith h5py.File(p, 'r') as f: assert list(f['values'][:]) == [0,1,2]\nos.remove(p)")
cat("FRESH PYTHON SMOKE TEST PASSED\n")
