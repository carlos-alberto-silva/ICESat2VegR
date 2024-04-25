
polygon_filepath <- "C:/Users/c.silva/Documents/ICESat2VegR/polygon.shp"

require(sf)
pol = st_polygon(
  list(
    poly1=cbind(
      lon<-c(-107.7,-106.5,-106.5,-107.7,-107.7),
      lat<-c(32.75,32.75,42.75,42.75,32.75)),
    poly2=cbind(
      lon<-c(-107.7,-106.5,-106.5,-107.7,-107.7)+1.5,
      lat<-c(32.75,32.75,42.75,42.75,32.75)+12)
  )
)

pol1 = st_polygon(
  list(
    poly1=cbind(
      lon<-c(-107.7,-106.5,-106.5,-107.7,-107.7),
      lat<-c(32.75,32.75,42.75,42.75,32.75))
     )
)

pol2 = st_polygon(
  list(
       poly2=cbind(
      lon<-c(-107.7,-106.5,-106.5,-107.7,-107.7)+1.5,
      lat<-c(32.75,32.75,42.75,42.75,32.75)+12)
  )
)


ncSplit <- split(nc, f = nc[[splitField]])

class(pol)
polc = st_sfc(list(poly1=pol1,poly2=pol2))

st_geometry(polc)

polc

getwd()
st_write(polc, dsn = "polygon.shp", layer = "polygon.shp", driver = "ESRI Shapefile", overwrite=T)
## Writing layer `nc' to data source `nc1.shp' using driver `ESRI Shapefile'
## Writing 100 features with 14 fields and geometry type Multi Polygon.

plot(canopy_metrics$longitude,canopy_metrics$latitude)
plot(canopy_metrics_clip$longitude,canopy_metrics_clip$latitude, col="red")

plot(newFile$longitude,newFile$latitude, col="red")





if (hasLeaflet) {
  leaflet() %>%
    addCircleMarkers(canopy_metrics_clip$longitude,
                     canopy_metrics_clip$latitude,
                     radius = 1,
                     opacity = 1,
                     color = "red"
    ) %>%
    addScaleBar(options = list(imperial = FALSE)) %>%
    addPolygons(
      data = polygon, weight = 1, col = "white",
      opacity = 1, fillOpacity = 0
    ) %>%
    addProviderTiles(providers$Esri.WorldImagery)
}



points(canopy_metrics_clip$longitude,canopy_metrics_clip$latitude, col="red")

plot(pol, add=T)


polygon0 <- sf::st_read(polygon_filepath)

polygon0$geometry

require(hdf5r)

atl08_path<-"Z:\\01_Projects\\04_NASA_ICESat2\\10_others\\ICESat2VegR\\inst\\exdata\\ATL08_20220401221822_01501506_005_01.h5"
atl08_h5<-ATL08_read(ATL08path=atl08_path)

ATL03_read <-function(ATL03path) {
  ATL03_h5 <- hdf5r::H5File$new(ATL03path, mode = 'r')
  ATL03<- new("icesat2.ATL03", h5 = ATL03_h5)
  return(ATL03)
}

atl03_path<-"C:\\Users\\c.silva\\Documents\\ICESat2VegR\\inst\\exdata\\ATL03_20220401221822_01501506_005_01.h5"
atl03_h5<-ATL03_read(ATL03path="C:\\Users\\c.silva\\Documents\\ICESat2VegR\\inst\\exdata\\ATL03_20220826063041_09981605_005_01.h5")

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
if (length(canopy_attribute) > 1 && length(terrain_attribute) > 1) { glist=c("terrain","canopy")}



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
                                   hdf5r::list.groups(level2a, recursive = F)), value = TRUE)







#'Get ICESat-2 ATL08-derived land elevation and vegetation heights
#'
#'@description This function extracts land elevation and vegetation heights from ICESat-2 ATL08 data
#'and convert to data.table [data.table::data.table]
#'
#'@param atl08_class A ICESat-2 ATL08 object (output of [ATL08_read()] function).
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
  if (!is.character(atl08_h5) || !tools::file_ext(atl08_h5) == "h5") {
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
  suppressWarnings(cycle <- as.integer(rhdf5::h5read(file = atl08_h5, name = "orbit_info/cycle_number")))
  suppressWarnings(region <- as.integer(rhdf5::h5read(file = atl08_h5, name = "ancillary_data/start_region")))

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
                                     hdf5r::list.groups(level2a, recursive = F)), value = TRUE)
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


