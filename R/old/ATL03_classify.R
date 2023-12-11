#ATL08_photon.var.map
ATL08_photon.var.map=list()
ATL08_photon.var.map[["ph_segment_id"]]="ph_segment_id"
ATL08_photon.var.map[["classed_pc_indx"]]="classed_pc_indx"
ATL08_photon.var.map[["classed_pc_flag"]]="classed_pc_flag"
ATL08_photon.var.map[["ph_h"]]="ph_h"
ATL08_photon.var.map[["d_flag"]]="d_flag"
ATL08_photon.var.map[["delta_time"]]="delta_time"
#'
#'ATL08 computed photons attributes
#'
#'@description This function extracts computed photons attributes from ICESat-2 ATL08 data
#'
#'@usage ATL08_photons_attributes_dt(atl08_h5, beam)
#'
#'@param atl08_h5 A ICESat-2 ATL08 object (output of [ATL08read()] function).
#'An S4 object of class [rICESat2Veg::icesat2.atl08_dt].
#'@param beam Character vector indicating beams to process (e.g. "gt1l", "gt1r", "gt2l", "gt2r", "gt3l", "gt3r")
#'
#'@return Returns an S4 object of class [data.table::data.table]
#'containing the ATL08 computed photons attributes.
#'
#'@details These are the photons attributes extracted by default:
#'\itemize{
#'\item \emph{ph_segment_id} Georeferenced	bin	number (20-m) associated	with	each photon
#'\item \emph{classed_pc_indx} Indices of photons	tracking back	to ATL03	that	surface finding	software	identified and	used	within	the
#'creation of the	data products.
#'\item \emph{classed_pc_flag} The L2B algorithm is run if this flag is set to 1 indicating data have sufficient waveform fidelity for L2B to run
#'\item \emph{ph_h} Height of photon above interpolated ground surface
#'#'\item \emph{d_flag} Flag indicating	whether DRAGANN	labeled	the photon as noise or signal
#'\item \emph{delta_time} Mid-segment	GPS	time	in seconds past	an epoch. The epoch is provided	in the metadata	at the file	level
#'}
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
#atl08_h5<-ATL08read(atl08_path=atl08_path)
#'
#'# Extracting ATL08 classified photons and heights
#'ATL08_photons_attributes_dt<-ATL08_photons_attributes_dt(atl08_h5=atl08_h5)
#'head(ATL08_photons_attributes_dt)
#'
#'close(atl08_h5)
#'@export
ATL08_photons_attributes_dt <- function(atl03_h5,atl08_h5,
                       beam = c("gt1l", "gt1r", "gt2l", "gt2r", "gt3l", "gt3r")) {

  # Check file input
  if (!class(atl03_h5)=="icesat2.atl03") {
    stop("atl03_h5 must be an object of class 'icesat2.atl08' - output of [ATL03read()] function ")
  }

  # Check file input
  if (!class(atl08_h5)=="icesat2.atl08") {
    stop("atl08_h5 must be an object of class 'icesat2.atl08' - output of [ATL08read()] function ")
  }

  #h5
  atl03_h5v2<-atl03_h5@h5
  atl08_h5v2<-atl08_h5@h5

  # Check beams to select
  groups_id_atl03<-hdf5r::list.groups(atl03_h5v2, recursive = F)
  groups_id_atl08<-hdf5r::list.groups(atl08_h5v2, recursive = F)

  check_beams_atl03<-groups_id_atl03 %in% beam
  check_beams_atl08<-groups_id_atl08 %in% beam

  beam_atl03<-groups_id_atl03[check_beams_atl03]
  beam_atl08<-groups_id_atl08[check_beams_atl08]

  beam<-beam_atl03[beam_atl03 %in% beam_atl08]

  photon.dt <- data.table::data.table()

  pb <- utils::txtProgressBar(min = 0, max = length(beam), style = 3)

  i_s = 0

  for (i in beam) {
      i_s = i_s + 1
      utils::setTxtProgressBar(pb, i_s)

      dataTableATL03Photons <- data.table::data.table(data.frame(
        lon_ph = atl03_h5v2[[paste0(i,"/heights/lon_ph")]][],
        lat_ph = atl03_h5v2[[paste0(i,"/heights/lat_ph")]][],
        h_ph = atl03_h5v2[[paste0(i,"/heights/h_ph")]][]
      ))

      dataTableATL03Segs <- data.table::data.table(data.frame(
        segment_id = atl03_h5v2[[paste0(i,"/geolocation/segment_id")]][],
        ph_index_beg = atl03_h5v2[[paste0(i,"/geolocation/ph_index_beg")]][]
      ))


      dataTableATL08Photons <- data.table::data.table(data.frame(
        ph_segment_id = atl08_h5v2[[paste0(i,"/signal_photons/ph_segment_id")]][],
        classed_pc_indx = atl08_h5v2[[paste0(i,"/signal_photons/classed_pc_indx")]][],
        classed_pc_flag = atl08_h5v2[[paste0(i,"/signal_photons/classed_pc_flag")]][],
        ph_h = atl08_h5v2[[paste0(i,"/signal_photons/ph_h")]][]
      ))
      data.table::setindex(dataTableATL03Segs, "segment_id")
      data.table::setindex(dataTableATL08Photons, "ph_segment_id")

      idx = data.table::merge.data.table(dataTableATL03Segs, dataTableATL08Photons, by.x = "segment_id", by.y = "ph_segment_id")[, ph_index_beg + classed_pc_indx-1]



      dataTableATL03Photons[idx, c("ph_segment_id","classed_pc_indx", "classed_pc_flag","ph_h") := list(
                                      dataTableATL08Photons$ph_segment_id,
                                      dataTableATL08Photons$classed_pc_indx,
                                      dataTableATL08Photons$classed_pc_flag,
                                      dataTableATL08Photons$ph_h)]


      dataTableATL03Photons$beam<-i

      photon_dt = data.table::rbindlist(list(photon.dt, dataTableATL03Photons), fill = TRUE)
    }

  close(pb)

  return(photon_dt)
}

head(photon_dt)

plot(photon_dt$)
summary(photon_dt)


library(hdf5r)
atl03_path <- "Z:\\01_Projects\\04_NASA_ICESat2\\10_others\\rICESat2Veg\\inst\\exdata\\ATL03_20220401221822_01501506_005_01.h5"
atl08_path <- "Z:\\01_Projects\\04_NASA_ICESat2\\10_others\\rICESat2Veg\\inst\\exdata\\ATL08_20220401221822_01501506_005_01.h5"

atl03_h5<-ATL03read(atl03_path=atl03_path)
atl02_h5<-ATL08read(atl08_path=atl08_path)


atl03 <- hdf5r::H5File$new(atl03_path, mode = "r")
atl08 <- hdf5r::H5File$new(atl08_path, mode = "r")


sum(atl03[["gt1r/geolocation/segment_ph_cnt"]][])


library(data.table)
dataTableATL03Photons <- data.table::data.table(data.frame(
  x = atl03[["gt1r/heights/lon_ph"]][],
  y = atl03[["gt1r/heights/lat_ph"]][],
  z = atl03[["gt1r/heights/h_ph"]][]
))

dataTableATL03Segs <- data.table::data.table(data.frame(
  id = atl03[["gt1r/geolocation/segment_id"]][],
  phIndexBeg = atl03[["gt1r/geolocation/ph_index_beg"]][]
))


dataTableATL08Photons <- data.table::data.table(data.frame(
  segmentId = atl08[["gt1r/signal_photons/ph_segment_id"]][],
  phIdx = atl08[["gt1r/signal_photons/classed_pc_indx"]][],
  class = atl08[["gt1r/signal_photons/classed_pc_flag"]][]
))
data.table::setindex(dataTableATL03Segs, "id")
data.table::setindex(dataTableATL08Photons, "segmentId")

idx = data.table::merge.data.table(dataTableATL03Segs, dataTableATL08Photons, by.x = "id", by.y = "segmentId")[, phIndexBeg + phIdx]


dataTableATL03Photons[idx, class := dataTableATL08Photons$class]

