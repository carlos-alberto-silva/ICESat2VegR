atl08_path<-"C:\\Users\\c.silva\\Documents\\rICESat2Veg\\inst\\exdata\\ATL08_20220401221822_01501506_005_01.h5"
atl08_h5<-ATL08read(atl08_path="C:\\Users\\c.silva\\Documents\\rICESat2Veg\\inst\\exdata\\ATL08_20220826063041_09981605_005_01.h5")

ATL03read <-function(atl03_path) {
  ATL03_h5 <- hdf5r::H5File$new(atl03_path, mode = 'r')
  ATL03<- new("icesat2.ATL03", h5 = ATL03_h5)
  return(ATL03)
}

atl03_path<-"C:\\Users\\c.silva\\Documents\\rICESat2Veg\\inst\\exdata\\ATL03_20220401221822_01501506_005_01.h5"
atl03_h5<-ATL03read(atl03_path="C:\\Users\\c.silva\\Documents\\rICESat2Veg\\inst\\exdata\\ATL03_20220826063041_09981605_005_01.h5")

names(atl03_h5@h5[["gt1r/heights"]])

length(atl03_h5@h5[["gt1r/heights/lat_ph"]][])

length(atl08_h5v2[["gt1r/signal_photons/ph_segment_id"]][])

head(atl08_h5v2[["gt1r/signal_photons/classed_pc_indx"]][])

atl08_h5v2[[paste0(beam[1],"/signal_photons/")]]

hist(atl08_h5v2[["gt1r/signal_photons/ph_h"]][])



col<-atl08_h5v2[["gt1r/signal_photons/classed_pc_flag"]][][1:50000]

col[col==0]<-"lightgray"
col[col==1]<-"chocolate4"
col[col==2]<-"dark green"
col[col==3]<-"green"

windows()
plot(1:50000,atl08_h5v2[["gt1r/signal_photons/ph_h"]][][1:50000], col=col, ylim=c(0,30))


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
