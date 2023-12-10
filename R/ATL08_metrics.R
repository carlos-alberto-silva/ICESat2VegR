#ATL08_canopy.var.map
ATL08_canopy.var.map=list()
ATL08_canopy.var.map[["h_canopy"]]="h_canopy"
ATL08_canopy.var.map[["canopy_rh_conf"]]="canopy_rh_conf"
ATL08_canopy.var.map[["h_median_canopy_abs"]]="h_median_canopy_abs"
ATL08_canopy.var.map[["h_min_canopy"]]="h_min_canopy"
ATL08_canopy.var.map[["h_mean_canopy_abs"]]="h_mean_canopy_abs"
ATL08_canopy.var.map[["h_median_canopy"]]="h_median_canopy"
ATL08_canopy.var.map[["h_canopy_abs"]]="h_canopy_abs"
ATL08_canopy.var.map[["toc_roughness"]]="toc_roughness"
ATL08_canopy.var.map[["h_min_canopy_abs"]]="h_min_canopy_abs"
ATL08_canopy.var.map[["h_dif_canopy"]]="h_dif_canopy"
ATL08_canopy.var.map[["h_canopy_quad"]]="h_canopy_quad"
ATL08_canopy.var.map[["h_canopy_20m"]]="h_canopy_20m"
ATL08_canopy.var.map[["n_ca_photons"]]="n_ca_photons"
ATL08_canopy.var.map[["photon_rate_can"]]="photon_rate_can"
ATL08_canopy.var.map[["centroid_height"]]="centroid_height"
ATL08_canopy.var.map[["canopy_h_metrics_abs"]]="canopy_h_metrics_abs"
ATL08_canopy.var.map[["h_mean_canopy"]]="h_mean_canopy"
ATL08_canopy.var.map[["subset_can_flag"]]="subset_can_flag"
ATL08_canopy.var.map[["canopy_h_metrics"]]="canopy_h_metrics"
ATL08_canopy.var.map[["n_toc_photons"]]="n_toc_photons"
ATL08_canopy.var.map[["h_max_canopy_abs"]]="h_max_canopy_abs"
ATL08_canopy.var.map[["h_canopy_uncertainty"]]="h_canopy_uncertainty"
ATL08_canopy.var.map[["canopy_openness"]]="canopy_openness"
ATL08_canopy.var.map[["h_max_canopy"]]="h_max_canopy"
ATL08_canopy.var.map[["segment_cover"]]="segment_cover"

#ATL08_terrain.var.map
ATL08_terrain.var.map=list()
ATL08_terrain.var.map[["h_te_best_fit"]]="h_te_best_fit"
ATL08_terrain.var.map[["h_te_best_fit_20m"]]="h_te_best_fit_20m"
ATL08_terrain.var.map[["h_te_interp"]]="h_te_interp"
ATL08_terrain.var.map[["h_te_max"]]="h_te_max"
ATL08_terrain.var.map[["h_te_mean"]]="h_te_mean"
ATL08_terrain.var.map[["h_te_median"]]="h_te_median"
ATL08_terrain.var.map[["h_te_mode"]]="h_te_mode"
ATL08_terrain.var.map[["h_te_rh25"]]="h_te_rh25"
ATL08_terrain.var.map[["h_te_skew"]]="h_te_skew"
ATL08_terrain.var.map[["h_te_std"]]="h_te_std"
ATL08_terrain.var.map[["h_canopy"]]="h_canopy"
ATL08_terrain.var.map[["h_te_uncertainty"]]="h_te_uncertainty"
ATL08_terrain.var.map[["n_te_photons"]]="n_te_photons"
ATL08_terrain.var.map[["photon_rate_te"]]="photon_rate_te"
ATL08_terrain.var.map[["subset_te_flag"]]="subset_te_flag"
ATL08_terrain.var.map[["terrain_slope"]]="terrain_slope"



#'Extract ICESat-2 ATL08-derived land elevation and vegetation heights
#'
#'@description This function extracts user defined land elevation and vegetation heights metrics from ICESat-2 ATL08 data
#'and convert to data.table [data.table::data.table]
#'
#'@usage ATL08_metrics(atl08_h5, beam)
#'
#'@param atl08_h5 A ICESat-2 ATL08 object (output of [readATL08()] function).
#'An S4 object of class "icesat2.atl08".
#'@param beam Character vector indicating beams to process (e.g. "gt1l", "gt1r", "gt2l", "gt2r", "gt3l", "gt3r")
#'
#'@return Returns an S4 object of class [data.table::data.table]
#'containing the aTL08-derived terrain elevation and vegetation relative heights metrics.
#'
#'@seealso \url{https://icesat-2.gsfc.nasa.gov/sites/default/files/page_files/ICESat2_ATL08_ATBD_r006.pdf}
#'
#'
#'@examples
#'
#'# Specifying the path to ATL08 file (zip file)
#'outdir = tempdir()
#'atl08_zip <- system.file("extdata",
#'                   "ATL08_20220401221822_01501506_005_01.zip",
#'                   package="rICESat2Veg")
#'
#'# Unzipping ATL08 file
#'atl08_path <- unzip(atl08_zip,exdir = outdir)
#'
#'# Reading ATL08 data (h5 file)
#atl08_h5<-readATL08(ATL08path=atl08_path)
#'
#'# Extracting ATL08 Elevation and Canopy Height Metrics
#'ATL08_metrics<-ATL08_metrics(atl08_h5=atl08_h5)
#'head(ATL08_metrics$terrain)
#'
#'close(atl08_h5)
#'@export
ATL08_metrics <- function(atl08_h5,
                       beam = c("gt1l", "gt1r", "gt2l", "gt2r", "gt3l", "gt3r"),
                       canopy_attribute=c(
                         "h_canopy",
                         "canopy_rh_conf",
                         "h_median_canopy_abs",
                         "h_min_canopy",
                         "h_mean_canopy_abs",
                         "h_median_canopy",
                         "h_canopy_abs",
                         "toc_roughness",
                         "h_min_canopy_abs",
                         "h_dif_canopy",
                         "h_canopy_quad",
                         "h_canopy_20m",
                         "n_ca_photons",
                         "photon_rate_can",
                         "centroid_height",
                         "canopy_h_metrics_abs",
                         "h_mean_canopy",
                         "subset_can_flag",
                         "canopy_h_metrics",
                         "n_toc_photons",
                         "h_max_canopy_abs",
                         "h_canopy_uncertainty",
                         "canopy_openness",
                         "h_max_canopy",
                         "segment_cover"),
                      terrain_attribute=c(
                        "h_te_best_fit",
                        "h_te_best_fit_20m",
                        "h_te_interp",
                        "h_te_max",
                        "h_te_mean",
                        "h_te_median",
                        "h_te_min",
                        "h_te_mode",
                        "h_te_rh25",
                        "h_te_skew",
                        "h_te_std",
                        "h_te_uncertainty",
                        "n_te_photons",
                        "photon_rate_te",
                        "subset_te_flag",
                        "terrain_slope")
                       ) {

  # Check file input
  if (!class(atl08_h5)=="icesat2.atl08") {
    stop("atl08_h5 must be an object of class 'icesat2.atl08' - output of [readATL08()] function ")
  }

  #h5
  atl08_h5v2<-atl08_h5@h5

  # Check beams to select
  groups_id<-hdf5r::list.groups(atl08_h5v2, recursive = F)

  check_beams<-groups_id %in% beam
  beam<-groups_id[check_beams]

  terrain.dt <- data.table::data.table()
  canopy.dt <- data.table::data.table()

  if (length(terrain_attribute) > 1) max=length(beam)
  if (length(canopy_attribute) > 1) max=length(beam)
  if (length(canopy_attribute) > 1 & length(terrain_attribute)) max=length(beam)*2

  pb <- utils::txtProgressBar(min = 0, max = max, style = 3)

  i_s = 0

  if (length(terrain_attribute) > 1) {

    for (i in beam) {
      i_s = i_s + 1
      utils::setTxtProgressBar(pb, i_s)

      atl08_h5v2_i<-atl08_h5v2[[paste0(i,"/land_segments/terrain")]]

      lat_i<-atl08_h5v2[[paste0(i,"/land_segments/latitude")]][]
      lon_i<-atl08_h5v2[[paste0(i,"/land_segments/longitude")]][]

      m = data.table::data.table(latitude=lat_i,longitude=lon_i)

      for (col in terrain_attribute) {
        #print(col)
        metric_address = ATL08_terrain.var.map[[col]]

        if (is.null(metric_address)) {
          if (atl08_h5v2_i$exists(col)) {
            metric_address = col
          } else {
            if (!col %in% names(atl08_h5v2_i)) warning(
              sprintf(
                "The column '%s' is not available in the ATL08 product!",
                col
              )
            )
            m[, eval(col) := NA]
            next
          }
        }
        base_addr = gsub("^(.*)/.*", "\\1", metric_address)
        if (atl08_h5v2_i$exists(base_addr) && atl08_h5v2_i$exists(metric_address))
          if (metric_address %in% c("h_te_best_fit_20m","subset_te_flag")) {
            m<-cbind(m,t(atl08_h5v2_i[[metric_address]][,]))
          } else {
            m[, eval(col) := atl08_h5v2_i[[metric_address]][]]
          }
      }
      terrain_dt = data.table::rbindlist(list(terrain.dt, m), fill = TRUE)
    }
   }



      if (length(canopy_attribute) > 1) {

        for (i in beam) {
          i_s = i_s + 1
          utils::setTxtProgressBar(pb, i_s)
          atl08_h5v2_i<-atl08_h5v2[[paste0(i,"/land_segments/canopy")]]
          lat_i<-atl08_h5v2[[paste0(i,"/land_segments/latitude")]][]
          lon_i<-atl08_h5v2[[paste0(i,"/land_segments/longitude")]][]

          m = data.table::data.table(latitude=lat_i,longitude=lon_i)

          for (col in canopy_attribute) {
            #print(col)
            metric_address = ATL08_canopy.var.map[[col]]

            if (is.null(metric_address)) {
              if (atl08_h5v2_i$exists(col)) {
                metric_address = col
              } else {
                if (i.s == 1) warning(
                  sprintf(
                    "The column '%s' is not available in the ATL08 product!",
                    col
                  )
                )
                m[, eval(col) := NA]
                next
              }
            }
            base_addr = gsub("^(.*)/.*", "\\1", metric_address)
            if (atl08_h5v2_i$exists(base_addr) && atl08_h5v2_i$exists(metric_address))
              if (metric_address %in% c("h_canopy_20m","canopy_h_metrics_abs","subset_can_flag","canopy_h_metrics")) {
                m<-cbind(m,t(atl08_h5v2_i[[metric_address]][,]))
              } else {
                m[, eval(col) := atl08_h5v2_i[[metric_address]][]]
              }
          }
          canopy_dt = data.table::rbindlist(list(canopy.dt, m), fill = TRUE)
      }


      }

  close(pb)
  m_dt<-list(terrain=terrain_dt, canopy=canopy_dt)


  return(m_dt)
}




atl08_path<-"C:\\Users\\c.silva\\Documents\\rICESat2Veg\\inst\\exdata\\ATL08_20220401221822_01501506_005_01.h5"
atl08_h5<-readATL08(ATL08path="C:\\Users\\c.silva\\Documents\\rICESat2Veg\\inst\\exdata\\ATL08_20220826063041_09981605_005_01.h5")

readATL03 <-function(ATL03path) {
  ATL03_h5 <- hdf5r::H5File$new(ATL03path, mode = 'r')
  ATL03<- new("icesat2.ATL03", h5 = ATL03_h5)
  return(ATL03)
}

atl03_path<-"C:\\Users\\c.silva\\Documents\\rICESat2Veg\\inst\\exdata\\ATL03_20220401221822_01501506_005_01.h5"
atl03_h5<-readATL03(ATL03path="C:\\Users\\c.silva\\Documents\\rICESat2Veg\\inst\\exdata\\ATL03_20220826063041_09981605_005_01.h5")

names(atl03_h5@h5[["gt1r/heights"]])

length(atl03_h5@h5[["gt1r/heights/lat_ph"]][])

length(atl08_h5v2[["gt1r/signal_photons/ph_h"]][])

head(atl08_h5v2[["gt1r/signal_photons/ph_h"]][])

atl08_h5v2[[paste0(beam[1],"/signal_photons/")]]

hist(atl08_h5v2[["gt1r/signal_photons/ph_h"]][])


required <- c(paste0(beam[1],"/signal_photons/ph_segment_id"),
              paste0(beam[1],"/signal_photons/classed_pc_indx"),
              paste0(beam[1],"/signal_photons/classed_pc_flag"))

hhh<-atl08_h5v2[[required[3]]][]



ph_segment_id<-atl08_h5v2[["gt1r/land_segments/segment_id_beg"]][]
length(unique(ph_segment_id))

length(min(ph_segment_id):max(ph_segment_id))




atl08_h5v2[["gt1r/land_segments/segment_id_beg"]][]

atl08_h5v2[["gt1r/land_segments/segment_id_end"]][]

names(atl08_h5v2[["gt1r/land_segments"]])

names(atl08_h5v2[["gt1r/land_segments"]])

length(atl08_h5v2_i[["gt1r/signal_photons/ph_segment_id"]][])/5


length(hhh)


required <- c(paste0(beam,"/signal_photons/ph_segment_id"),
              paste0(beam,"/signal_photons/classed_pc_indx"),
              paste0(beam,"/signal_photons/classed_pc_flag"))

check <- required %in% atl08_h5v2list

canoly_list <- atl08_h5v2[[paste0(beam,"/land_segments/canopy")]]
terrain_list <- atl08_h5v2[[paste0(beam,"/land_segments/terrain")]]

m_dt <- list(terrain=data.table::data.table(),
             canopy=data.table::data.table())

canopy.dt <- data.table::data.table()

pb <- utils::txtProgressBar(min = 0, max = length(canopy_attribute), style = 3)
i_s = 0

if (length(canopy_attribute) > 1) { glist="canopy"}
if (length(terrain_attribute) > 1) { glist="terrain"}
if (length(canopy_attribute) > 1 & length(terrain_attribute) > 1) { glist=c("terrain","canopy")}



remotes::install_git("https://github.com/carlos-alberto-silva/rGEDI", dependencies = TRUE)

require(rGEDI)

# Set output dir for downloading the files
outdir=getwd()

# Downloading GEDI data
gediDownload(filepath=gLevel1B,outdir=outdir)
gediDownload(filepath=gLevel2A,outdir=outdir)
gediDownload(filepath=gLevel2B,outdir=outdir)

#######
# Herein, we are using only a GEDI sample dataset for this tutorial.
#######
# downloading zip file
download.file("https://github.com/carlos-alberto-silva/rGEDI/releases/download/datasets/examples.zip",destfile=file.path(outdir, "examples.zip"))

# unzip file
unzip(file.path(outdir,"examples.zip"))

# Reading GEDI data
gedilevel1b<-readLevel1B(level1Bpath = file.path(outdir,"GEDI01_B_2019108080338_O01964_T05337_02_003_01_sub.h5"))
gedilevel2a<-readLevel2A(level2Apath = file.path(outdir,"GEDI02_A_2019108080338_O01964_T05337_02_001_01_sub.h5"))
gedilevel2b<-readLevel2B(level2Bpath = file.path(outdir,"GEDI02_B_2019108080338_O01964_T05337_02_001_01_sub.h5"))


level2a<-gedilevel2a@h5
groups_id<-grep("BEAM\\d{4}$",gsub("/","",
                                   hdf5r::list.groups(level2a, recursive = F)), value = T)







#'Get ICESat-2 ATL08-derived land elevation and vegetation heights
#'
#'@description This function extracts land elevation and vegetation heights from ICESat-2 ATL08 data
#'and convert to data.table [data.table::data.table]
#'
#'@usage getATL08(atl08_class)
#'
#'@param atl08_class A ICESat-2 ATL08 object (output of [readATL08()] function).
#'An S4 object of class "icesat2.atl08".
#' @param beam Character vector indicating beams to process
#' @param beam_strength Character vector indicating the strength of beams to process
#'
#'@return Returns an S4 object of class [data.table::data.table]
#'containing the land elevation and vegetation relative heights metrics.
#'
#'@seealso \url{https://icesat-2.gsfc.nasa.gov/sites/default/files/page_files/ICESat2_ATL08_ATBD_r006.pdf}
#'
#'@details Characteristics. Flag indicating likely invalid waveform (1=valid, 0=invalid).
#'\itemize{
#'\item \emph{beam} Beam identifier
#'\item \emph{shot_number} Shot number
#'\item \emph{degrade_flag} Flag indicating degraded state of pointing and/or positioning information
#'\item \emph{quality_flag} Flag simplifying selection of most useful data
#'\item \emph{delta_time} Transmit time of the shot since Jan 1 00:00 2018
#'\item \emph{sensitivity} Maxmimum canopy cover that can be penetrated
#'\item \emph{solar_elevation} Solar elevation
#'\item \emph{lat_lowestmode} Latitude of center of lowest mode
#'\item \emph{lon_lowestmode} Longitude of center of lowest mode
#'\item \emph{elev_highestreturn} Elevation of highest detected return relative to reference ellipsoid Meters
#'\item \emph{elev_lowestmode} Elevation of center of lowest mode relative to reference ellipsoid
#'\item \emph{rh} Relative height metrics at 1% interval
#'}
#'
#'@examples
#'
#'# Specifying the path to GEDI level2A data (zip file)
#'outdir = tempdir()
#'level2A_fp_zip <- system.file("extdata",
#'                   "GEDI02_A_2019108080338_O01964_T05337_02_001_01_sub.zip",
#'                   package="rGEDI")
#'
#'# Unzipping GEDI level2A data
#'level2Apath <- unzip(level2A_fp_zip,exdir = outdir)
#'
#'# Reading GEDI level2A data (h5 file)
#'level2a<-readLevel2A(level2Apath=level2Apath)
#'
#'# Extracting GEDI Elevation and Height Metrics
#'level2AM<-getLevel2AM(level2a)
#'head(level2AM)
#'
#'close(level2a)
#'@export
get_atl08 <- function(atl08_class,
                      beam = c("gt1l", "gt1r", "gt2l", "gt2r", "gt3l", "gt3r"),
                      beam_strength = c("weak", "strong")) {

  # Check file input
  if (!is.character(atl08_h5) | !tools::file_ext(atl08_h5) == "h5") {
    stop("atl08_h5 must be a path to a h5 file")
  }

  # Filter beams to analyze
  if (any(!beam_strength %in% c("weak", "strong"))) {
    stop("beam_strength must contain weak, beam or both ")
  }

  # Check beams to select
  suppressWarnings(sc_orient <- as.numeric(rhdf5::h5read(file = atl08_h5, name = "orbit_info/sc_orient")))

  if(sc_orient == 0) { # Backward
    list_strength <- list(weak = c("gt1r", "gt2r", "gt3r"),
                          strong = c("gt1l", "gt2l", "gt3l"))
  }else{ #Forward
    list_strength <- list(weak = c("gt1l", "gt2l", "gt3l"),
                          strong = c("gt1r", "gt2r", "gt3r"))
  }

  if(length(beam_strength) == 1) {

    keep_beams <- beam %in% list_strength[[beam_strength]]

    if (sum(keep_beams) == 0) stop(sprintf("Selected beams are not %s",beam_strength))
    if (sum(keep_beams) != length(beam)) message(sprintf("Keeping only the following %s beams: %s", beam_strength, paste(beam[keep_beams], collapse = ",")))

    beam <- beam[keep_beams]

  }

  # Setting up min and max lat
  if(!missing(lat_range)) {
    if (length(lat_range) != 2) stop("lat_range must have a min and and max value")
    if(lat_range[1] >= lat_range[2]) stop("First element of lat_range must be lower than second element")
  }else{
    lat_range <- c(NA, NA)
  }

  # Check that the necessary data is available
  atl08_h5list <- suppressWarnings(rhdf5::h5ls(atl08_h5))

  required <- c(paste0(beam,"/signal_photons/ph_segment_id"),
                paste0(beam,"/signal_photons/classed_pc_indx"),
                paste0(beam,"/signal_photons/classed_pc_flag"))

  check <- paste0("/",required) %in% paste(atl08_h5list$group, atl08_h5list$name, sep = "/")


  if (any(check == FALSE)) {
    missing_beams <- beam[(beam %in% unique(stringr::str_extract(required[!check],"gt..")))]
    warning(sprintf("Missing datasets in file: %s. Ignoring following beams: %s ", paste(required[!check], collapse = ","), paste(missing_beams, collapse = ", ") ))
    beam <- beam[!beam %in% missing_beams]
  }

  # Retrieve RGT, region and cycle
  suppressWarnings(rgt <- as.integer(rhdf5::h5read(file = atl08_h5, name = "orbit_info/rgt")))
  suppressWarnings(cycle <- as.integer(rhdf5::h5read(file = atl08_h5, name = "orbit_info/cycle_number")))
  suppressWarnings(region <- as.integer(rhdf5::h5read(file = atl08_h5, name = "ancillary_data/start_region")))

  out <- list() #Initiate list to return
  n=
    for ( n in beam) {

      beam_strength_n <- names(list_strength)[unlist(lapply(list_strength, function(x) n %in% x))]

      beam_strength_n= "weak"
      print(sprintf("Working on beam %s (%s)", n, beam_strength_n))

      # Find all 100-m segment center points within lat_range
      lat_center <- as.numeric(rhdf5::h5read(file = atl08_h5, name = paste0(n,"/land_segments/latitude")))

      if (is.na(lat_range[1])) {
        lat_range[1] <- min(lat_center)
      }

      if (is.na(lat_range[2])) {
        lat_range[2] <- max(lat_center)
      }

      # Index of points to keep
      lat_idx <- which(lat_center >= lat_range[1] & lat_center <= lat_range[2])

      if (sum(lat_idx) != 0) {


        # Land segment flags and info
        suppressWarnings(trg_fields <-  atl08_h5list %>%
                           dplyr::filter(stringr::str_detect(group, stringr::str_c(n, "/land_segments$")),
                                         otype == "H5I_DATASET") %>%
                           dplyr::mutate(full_name = stringr::str_c(group, name, sep = "/")) %>%
                           tidyr::separate(dim, into = c("dim1", "dim2"), sep = " x "))

        trg_fields_list <- lapply(1:length(trg_fields$name), function(x) {

          if (any(is.na(c(trg_fields$dim1[x], trg_fields$dim2[x])))) {
            trg_index <- list(lat_idx)
          }else{
            trg_index <- list(1:trg_fields$dim1[x], lat_idx)
          }

          suppressWarnings(out <- rhdf5::h5read(file = atl08_h5,
                                                name = trg_fields$full_name[x],
                                                index = trg_index))

          if(trg_fields$dclass[x] == "FLOAT") {
            if(length(dim(out)) > 1) {
              out <- as.data.frame(apply(out, 1, as.numeric))
              colnames(out) <- stringr::str_c("X", 1:ncol(out))

            }else{
              out <- as.numeric(out)
            }
          }else if (trg_fields$dclass[x] == "INTEGER") {
            if(length(dim(out)) > 1) {
              out <- as.data.frame(apply(out, 1, as.integer))
              colnames(out) <- stringr::str_c("X", 1:ncol(out))
            }else{
              out <- as.integer(out)
            }
          }

          return(out)

        })

        names(trg_fields_list) <- trg_fields$name

        land_df <- as.data.frame(trg_fields_list)

        # Canopy products

        suppressWarnings(trg_fields <-  atl08_h5list %>%
                           dplyr::filter(stringr::str_detect(group, stringr::str_c(n, "/land_segments/canopy$")),
                                         otype == "H5I_DATASET") %>%
                           dplyr::mutate(full_name = stringr::str_c(group, name, sep = "/"))  %>%
                           tidyr::separate(dim, into = c("dim1", "dim2"), sep = " x "))

        trg_fields_list <- lapply(1:length(trg_fields$name), function(x) {

          if (any(is.na(c(trg_fields$dim1[x], trg_fields$dim2[x])))) {
            trg_index <- list(lat_idx)
          }else{
            trg_index <- list(1:trg_fields$dim1[x], lat_idx)
          }

          suppressWarnings(out <- rhdf5::h5read(file = atl08_h5,
                                                name = trg_fields$full_name[x],
                                                index = trg_index))

          if(trg_fields$dclass[x] == "FLOAT") {
            if(length(dim(out)) > 1) {
              out <- as.data.frame(apply(out, 1, as.numeric))
              colnames(out) <- stringr::str_c("X", 1:ncol(out))

            }else{
              out <- as.numeric(out)
            }
          }else if (trg_fields$dclass[x] == "INTEGER") {
            if(length(dim(out)) > 1) {
              out <- as.data.frame(apply(out, 1, as.integer))
              colnames(out) <- stringr::str_c("X", 1:ncol(out))
            }else{
              out <- as.integer(out)
            }
          }

          return(out)

        })

        names(trg_fields_list) <- trg_fields$name

        canopy_df <- as.data.frame(trg_fields_list)

        # Terrain products

        suppressWarnings(trg_fields <-  atl08_h5list %>%
                           dplyr::filter(stringr::str_detect(group, stringr::str_c(n, "/land_segments/terrain$")),
                                         otype == "H5I_DATASET") %>%
                           dplyr::mutate(full_name = stringr::str_c(group, name, sep = "/")) %>%
                           tidyr::separate(dim, into = c("dim1", "dim2"), sep = " x "))

        trg_fields_list <- lapply(1:length(trg_fields$name), function(x) {

          if (any(is.na(c(trg_fields$dim1[x], trg_fields$dim2[x])))) {
            trg_index <- list(lat_idx)
          }else{
            trg_index <- list(1:trg_fields$dim1[x], lat_idx)
          }

          suppressWarnings(out <- rhdf5::h5read(file = atl08_h5,
                                                name = trg_fields$full_name[x],
                                                index = trg_index))

          if(trg_fields$dclass[x] == "FLOAT") {
            if(length(dim(out)) > 1) {
              out <- as.data.frame(apply(out, 1, as.numeric))
              colnames(out) <- stringr::str_c("X", 1:ncol(out))

            }else{
              out <- as.numeric(out)
            }
          }else if (trg_fields$dclass[x] == "INTEGER") {
            if(length(dim(out)) > 1) {
              out <- as.data.frame(apply(out, 1, as.integer))
              colnames(out) <- stringr::str_c("X", 1:ncol(out))
            }else{
              out <- as.integer(out)
            }
          }

          return(out)

        })

        names(trg_fields_list) <- trg_fields$name

        terrain_df <- as.data.frame(trg_fields_list)

        # ATL08 output

        atl08_out <- cbind(data.frame(cycle = cycle,
                                      region = region,
                                      beam = n,
                                      beam_strength = beam_strength_n),
                           land_df,
                           canopy_df,
                           terrain_df)

        colnames(atl08_out) <- stringr::str_replace(colnames(atl08_out),
                                                    pattern = "\\.",
                                                    replacement = "_")

      }

    }
}




getLevel2AM<-function(level2a){
  level2a<-level2a@h5
  groups_id<-grep("BEAM\\d{4}$",gsub("/","",
                                     hdf5r::list.groups(level2a, recursive = F)), value = T)
  rh.dt<-data.table::data.table()
  pb <- utils::txtProgressBar(min = 0, max = length(groups_id), style = 3)
  i.s=0

  for ( i in groups_id){
    i.s<-i.s+1
    utils::setTxtProgressBar(pb, i.s)
    level2a_i<-level2a[[i]]

    if (any(hdf5r::list.datasets(level2a_i)=="shot_number")){

      if(length(level2a_i[["rh"]]$dims)==2) {
        rh=t(level2a_i[["rh"]][,])
      } else {
        rh=t(level2a_i[["rh"]][])
      }

      rhs<-data.table::data.table(
        beam<-rep(i,length(level2a_i[["shot_number"]][])),
        shot_number=level2a_i[["shot_number"]][],
        degrade_flag=level2a_i[["degrade_flag"]][],
        quality_flag=level2a_i[["quality_flag"]][],
        quality_flag=level2a_i[["delta_time"]][],
        sensitivity=level2a_i[["sensitivity"]][],
        solar_elevation=level2a_i[["solar_elevation"]][],
        lat_lowestmode=level2a_i[["lat_lowestmode"]][],
        lon_lowestmode=level2a_i[["lon_lowestmode"]][],
        elev_highestreturn=level2a_i[["elev_highestreturn"]][],
        elev_lowestmode=level2a_i[["elev_lowestmode"]][],
        rh)
      rh.dt<-rbind(rh.dt,rhs)
    }
  }

  colnames(rh.dt)<-c("beam","shot_number","degrade_flag","quality_flag","delta_time",
                     "sensitivity","solar_elevation","lat_lowestmode","lon_lowestmode",
                     "elev_highestreturn","elev_lowestmode",paste0("rh",seq(0,100)))
  close(pb)
  return(rh.dt)
}


