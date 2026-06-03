# ICESat2VegR Test Case Inventory

This inventory is organized by test layer. The automated tests should cover the local, deterministic cases by default and keep Earthdata and Google Earth Engine checks opt-in because they require credentials, network access, or project-specific configuration.

## Unit Tests

- `fit_metrics()` returns the expected statistics for known observed/predicted vectors.
- `fit_metrics()` drops incomplete pairs before calculating metrics.
- `fit_metrics()` rejects vectors with different lengths.
- `fit_metrics()` rejects inputs with fewer than two complete pairs.
- `fit_metrics()` reports relative metrics as `NA` when the observed mean is zero.
- `to_vect()` infers `longitude`/`latitude` columns from a plain `data.frame`.
- `to_vect()` infers `lon_ph`/`lat_ph` columns from ICESat-2 photon tables.
- `to_vect()` honors explicit `lon` and `lat` column names.
- `to_vect()` removes non-finite coordinates before creating a vector.
- `to_vect()` errors clearly when coordinate columns cannot be inferred.
- `randomSampling()` returns the requested absolute sample size.
- `randomSampling()` supports fractional sample sizes.
- `gridSampling()` adds grid grouping columns and respects per-cell sample limits.
- `stratifiedSampling()` samples within histogram strata and adds the stratum column.
- Chained sampling methods execute in order and return a data.table-compatible result.
- `fit_model()` with `method = "none"` returns a fitted Random Forest model and training metrics.
- `fit_model()` with `method = "split"` returns train/test indexes, predictions, and test metrics.
- `varSel()` returns selected variables and importance tables for a small regression problem.

## Integration Tests

- Bundled ATL03 and ATL08 HDF5 fixtures can be opened and closed with `ATL03_read()` and `ATL08_read()`.
- ATL03 segment metadata extraction returns expected coordinate, beam, and requested attribute columns.
- ATL08 segment attribute extraction returns expected coordinate, beam, and requested attribute columns.
- ATL03 and ATL08 photon attributes can be joined from the bundled fixture pair.
- Joined photons can be segmented with `ATL03_ATL08_segment_create()` at a fixed segment length.
- Segment statistics can be computed from joined ATL03/ATL08 photons.
- Segment statistics can be converted to `terra::SpatVector`.
- Segment statistics can be clipped to the bundled AOI geometry.
- Extracted ATL08 attributes can be summarized to a raster grid with `ATL08_seg_attributes_dt_gridStat()`.
- Extracted vectors can be written to GeoPackage or GeoJSON with standard spatial writers.
- A simple local model can be fit with `fit_model()` and evaluated with `fit_metrics()`.

## End-to-End Tests

- README local workflow: read bundled ATL03/ATL08 granules, extract ATL03/ATL08 attributes, create vectors, summarize ATL08 canopy height, and write vector outputs.
- README complete script local phase: read AOI, derive bbox/date range, pair ATL03/ATL08 fixtures by timestamp, join photons, create 20 m segments, compute `rh98` metrics, clip to AOI, and export GeoJSON.
- README modelling phase: load bundled example segments, sample rows, select predictors from a small deterministic synthetic predictor table, fit a Random Forest, and compute train/test metrics.
- Earthdata workflow: search for ATL03/ATL08 granules and optionally download them when valid Earthdata credentials are present.
- Cloud workflow: initialize cloud reading through `earthaccess`/`h5py` when the reticulate environment and Earthdata credentials are present.
- Google Earth Engine workflow: initialize Earth Engine, create an AOI geometry, build predictor stacks, extract ancillary data, create a map, and export only when an EE project/token is configured.

## Manual QA / Frontend

- This repository is an R package and does not expose a browser frontend, so Playwright/Cypress coverage is not applicable here.
- Interactive map outputs from README examples should be manually inspected when generated locally: layer controls should render, map tiles should load, vector/raster overlays should appear in the expected AOI, and legends should match the selected attribute.
