# atl08 <- ATL08_read("../inst/extdata/ATL08_20191002161359_00880506_006_02.h5")
# library(ICESat2VegR)
devtools::load_all()
library(data.table)

atl03 <- ATL03_read("Z:/01_Projects/04_NASA_ICESat2/11_others/rICESat2Veg/inst/exdata/ATL03_20220401221822_01501506_006_02.h5")
atl08 <- ATL08_read("Z:/01_Projects/04_NASA_ICESat2/11_others/rICESat2Veg/inst/exdata/ATL08_20220401221822_01501506_006_02.h5")


atl03_clip <- ATL03_h5_clip(atl03, "atl03_clip.h5", "", function(x, ...) {
  ds <- x[["geolocation/segment_id"]]
  idx <- seq_len(ds$dims)
  idx[ds[] %in% 771236:771276]
}, beam = "gt1r", additional_groups = c("orbit_info"))

atl03_clip <- ATL03_read("ATL03_clip.h5")
atl03_seg_clip <- ATL03_seg_attributes_dt(atl03_clip)

atl08_clip <- ATL08_h5_clip(
  atl08,
  output = "atl08_clip.h5",
  clip_obj = "",
  landSegmentsMask_fn = function(x, ...) {
    mask <- x[["land_segments/segment_id_beg"]][] <= 771276 & x[["land_segments/segment_id_end"]][] >= 771236
    seq_len(x[["land_segments/segment_id_beg"]]$dims)[mask]
  },
  beam = "gt1r",
  additional_groups = c("orbit_info")
)
close(atl08_clip)
