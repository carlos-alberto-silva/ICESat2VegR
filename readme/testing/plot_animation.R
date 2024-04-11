library(gganimate)
library(ggplot2)

atl03_path <- system.file("extdata/atl03_clip.h5", package = "ICESat2VegR")
atl08_path <- system.file("extdata/atl08_clip.h5", package = "ICESat2VegR")

atl03 <- ATL03_read(atl03_path)
atl08 <- ATL08_read(atl08_path)

atl03_08_dt <- ATL03_ATL08_photons_attributes_dt_join(atl03, atl08)
atl03_08_dt <- na.omit(atl03_08_dt)

dt <- atl03_08_dt[ph_h < 100 & beam == "gt1r"]
dt <- dt[, classed_pc_flag := factor(classed_pc_flag)]


init <- 1
end <- 179
view_width <- 60

xmin <- seq(init, end - view_width, view_width)
xmax <- seq(init + view_width, end, view_width)

dt_clip <- dt[dist_ph_along >= init & dist_ph_along <= end]
res <- ggplot() +
  geom_point(data = dt_clip, aes(dist_ph_along, h_ph, color = classed_pc_flag)) +
  scale_color_manual(
    values = c("gray", "orange", "forestgreen", "green"),
    labels = c("Noise", "Ground", "Vegetation", "High Vegetation")
  ) +
  labs(
    color = "Photon classification",
    title = gettextf("Granule: %s", ATL03_name)
  ) +
  view_step_manual(
    1, 2,
    pause_first = FALSE, ease = "linear", wrap = FALSE,
    xmin = xmin, xmax = xmax, ymin = min(dt_clip$h_ph), ymax = max(dt_clip$h_ph)
  ) +
  theme_bw()

anim_save("output.gif", res)