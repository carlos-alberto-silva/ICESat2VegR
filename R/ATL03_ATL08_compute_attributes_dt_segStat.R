#' Statistics of ATL03 and ATL08 labeled photons at the segment level
#'
#' @description Computes a series of statistics from ATL03 and ATL08 labeled photons
#' within a given segment length
#'
#' @usage ATL03_ATL08_compute_seg_attributes_dt_segStat(atl03_atl08_dt, func, seg_length,ph_class,beam,quality_ph,night_flag)
#'
#' @param atl03_atl08_dt  An S4 object of class "icesat2.atl03atl08_dt" containing ATL03 and ATL08 data
#' (output of [ATL03_ATL08_photons_attributes_dt_join()] function).
#' @param func The function to be applied for computing the defined statistics
#' @param seg_length Segment length. Default is 30 m
#' @param ph_class Character vector indicating photons to process based on the classification (1=ground, 2=canopy, 3=top canopy),
#' Default is c(2,3)
#' @param beam Character vector indicating beams to process. Default is c("gt1l", "gt1r", "gt2l", "gt2r", "gt3l", "gt3r")
#' @param quality_ph Indicates the quality of the associated photon. 0=nominal, 1=possible_afterpulse, 2=possible_impulse_response_
#'effect, 3=possible_tep. Default is 0
#' @param night_flag Flag indicating the data were acquired in night conditions: 0=day, 1=night. Default is 1
#'
#' @return Returns an S4 object of class [rICESat2Veg::icesat2.atl08_dt]
#' Containing Statistics of ATL03 and ATL08 labeled photons
#'
#' @examples
#'
#'# Specifying the path to ATL03 and ATL08 file (zip file)
#'outdir = tempdir()
#'atl03_zip <- system.file("extdata",
#'                   "ATL03_20220401221822_01501506_005_01.zip",
#'                   package="rICESat2Veg")
#'
#'atl08_zip <- system.file("extdata",
#'                   "ATL08_20220401221822_01501506_005_01.zip",
#'                   package="rICESat2Veg")
#'
#'# Unzipping ATL03 file
#'atl03_path <- unzip(atl03_zip,exdir = outdir)
#'
#'# Unzipping ATL08 file
#'atl08_path <- unzip(atl08_zip,exdir = outdir)
#'
#'# Reading ATL03 data (h5 file)
#atl03_h5<-ATL08_read(atl03_path=atl03_path)
#'
#'# Reading ATL08 data (h5 file)
#atl08_h5<-ATL08_read(atl08_path=atl08_path)
#'
#'# Extracting ATL03 and ATL08 labeled photons
#'atl03_atl08_dt<-ATL03_ATL08_photons_attributes_dt_join(atl03_h5,atl08_h5)
#'head(atl03_atl08_dt)
#'
#'# Computing the max canopy height at 30 m segments
#'max_canopy <-ATL03_ATL08_compute_seg_attributes_dt_segStat(atl03_atl08_dt, func=max(ph_h),
#'                                seg_length = 30,
#'                                ph_class=c(2,3),
#'                                beam=c("gt1l", "gt1r", "gt2l", "gt2r", "gt3l", "gt3r"),
#'                                quality_ph=0,
#'                                night_flag=1)
#'
#'head(max_canopy)
#'
#'# Define your own function
#' mySetOfMetrics <- function(x) {
#'   metrics <- list(
#'     min = min(x), # Min of x
#'     max = max(x), # Max of x
#'     mean = mean(x), # Mean of x
#'     sd = sd(x) # Sd of x
#'     h_canopy = quantile(x,0.98)
#'   )
#'   return(metrics)
#' }
#'
#' # Computing a series of canopy height statistics from customized function
#'canopy_metrics <-ATL03_ATL08_compute_seg_attributes_dt_segStat(atl03_atl08_dt, func=mySetOfMetrics(ph_h),
#'                                seg_length = 30, # 30 m
#'                                ph_class=c(2,3),
#'                                beam=c("gt1l", "gt1r", "gt2l", "gt2r", "gt3l", "gt3r"),
#'                                quality_ph=0,
#'                                night_flag=1)
#'
#'head(canopy_metrics)
#'
#'close(atl03_h5)
#'close(atl08_h5)
#' @import data.table lazyeval
#' @export
ATL03_ATL08_compute_seg_attributes_dt_segStat <- function(atl03_atl08_dt,
                                                          func,
                                                          seg_length = 30,
                                                          ph_class=c(2,3),
                                 beam=c("gt1l", "gt1r", "gt2l", "gt2r", "gt3l", "gt3r"),
                                 quality_ph=0,
                                 night_flag=1) {

  if (!class(atl03_atl08_dt)[1]=="icesat2.atl03atl08_dt"){
    stop("atl03_atl08_dt needs to be an object of class 'icesat2.atl03atl08_dt' ")
  }

  atl03_atl08_dt2<-atl03_atl08_dt@dt[atl03_atl08_dt@dt$classed_pc_flag %in% ph_class &
                                       atl03_atl08_dt@dt$quality_ph %in% quality_ph &
                                       atl03_atl08_dt@dt$beam %in% beam &
                                       atl03_atl08_dt@dt$night_flag %in% night_flag,]

  segment<-cut(atl03_atl08_dt2$dist_ph_along, breaks=seq(0,max(atl03_atl08_dt2$dist_ph_along),seg_length),labels = FALSE)

  atl03_atl08_dt2$segment<-segment


  # Add data.table operator
  `:=` <- data.table::`:=`

  call <- lazy_call(func)

  call2 <- lazy_call(seg_position(lon_ph,lat_ph))

  if (is.null(segment)) {
    metrics <- lazy_apply_dt_call(dt = atl03_atl08_dt2, call = call)
    metrics <- as.data.table(metrics)
    if (ncol(metrics) < 2) {
      colnames(metrics) <- paste0(call)[1]
    }

  } else {
    metrics <- lazy_apply_dt_call(dt = atl03_atl08_dt2, call = call, group.by = "by = segment")
    latlon <- lazy_apply_dt_call(dt = atl03_atl08_dt2, call = call2, group.by = "by = segment")

    metrics$latitude<-latlon$V1
    metrics$longitude<-latlon$V2

    if (ncol(metrics) == 4) {
      colnames(metrics)[2] <- paste0(call)[1]
    }

    metrics$dist_along<-metrics$segment*seg_length - seg_length/2
    metrics<-metrics[-nrow(metrics),]
    colnames(metrics)[2:c(ncol(metrics)-3)]<-paste0(names(metrics)[2:c(ncol(metrics)-3)],"_",paste0(call)[2])
  }

  metrics<- new("icesat2.atl08_dt", dt = metrics)
  return(metrics)
}
