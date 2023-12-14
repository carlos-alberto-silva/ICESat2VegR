![](https://github.com/carlos-alberto-silva/rICESat2Veg/blob/master/readme/cover.png)<br/>
[![R-CMD-check](https://github.com/carlos-alberto-silva/rICESat2Veg/actions/workflows/r.yml/badge.svg?branch=master)](https://github.com/carlos-alberto-silva/rICESat2Veg/actions/workflows/r.yml)
[![CRAN](https://www.r-pkg.org/badges/version/rICESat2Veg)](https://cran.r-project.org/package=rICESat2Veg)
![Github](https://img.shields.io/badge/Github-0.1.12-green.svg)
![licence](https://img.shields.io/badge/Licence-GPL--3-blue.svg) 
![Downloads](https://cranlogs.r-pkg.org/badges/grand-total/rICESat2Veg)
[![Build Status](https://travis-ci.com/carlos-alberto-silva/rICESat2Veg.svg?token=Jqizwyc6gBxNafNccTdU&branch=master)](https://travis-ci.com/carlos-alberto-silva/rICESat2Veg)

**rICESat2Veg: An R Package for NASA's Ice, Cloud, and Elevation Satellite (ICESat-2) Data Processing and Visualization for Land and Vegetation Applications.**

Authors: Carlos Alberto Silva and Caio Hamamura  

The rICESat2Veg package provides functions for downloading, processing, visualizing and exporting ICESat-2 ATL03 and ATL08 data in R environment.

# Getting Started

## Installation
```r
#The CRAN version:
install.packages("rICESat2Veg")

# The development version:
#install.packages("remotes")
remotes::install_github("https://github.com/carlos-alberto-silva/rICESat2Veg", dependencies = TRUE)

# loading rGEDI package
library(rICESat2Veg)

```    
## Downloading ICESat-2 data
```r
#######
# Herein, we are using only a ICESat-2 sample dataset for this tutorial.
#######
# downloading zip file
download.file("https://github.com/carlos-alberto-silva/rICESat2Veg/releases/download/datasets/examples.zip",destfile=file.path(outdir, "examples.zip"))

# unzip file 
unzip(file.path(outdir,"examples.zip"))

```

## Read ATL03 and ATL08 files
```r
atl03_h5<-ATL03_read(atl03_path=atl03_path)
atl08_h5<-ATL08_read(atl08_path=atl08_path)
```

## Extracting ATL03 and ATL08 photon attributes
```r
atl03_photons_dt<-ATL03_photons_attributes_dt(atl03_h5=atl03_h5, beam="gt1l")
head(atl03_photons_dt)

atl08_photons_dt<-ATL08_photons_attributes_dt(atl08_h5=atl08_h5, beam="gt1l")
head(atl08_photons_dt)
```

## Clipping ATL03 photons attributes

```r
# Clipping by  bounding box

# Define the bounding box
xmin <- -107.7
xmax <- -106.5
ymin <- 32.75
ymax <- 42.75

# Clip
atl03_photons_dt_clip <- ATL03_photons_attributes_dt_clipBox(atl03_photons_dt, xmin, xmax, ymin, ymax)
head(atl03_photons_dt_clip) # print the first six observations
```

## Clipping by geometry
```r
# Specify the path to shapefile
polygon_filepath <- system.file("extdata", "polygon.shp", package = "rICESat2Veg")

# Read shapefile
polygon <- terra::vect(polygon_filepath)

# Clip
atl03_photons_dt_clipg <- ATL03_photons_attributes_dt_clipGeometry(atl03_photons_dt, polygon, split_by = "FID")
head(atl03_photons_dt_clipg) # print the first six observations
```

## Merging ATL03 and ATL08 photon attributes
```r
atl03_atl08_dt<-ATL03_ATL08_photons_attributes_dt_join(atl03_h5,atl08_h5, beam = "gt1l")
head(atl03_atl08_dt)
```

## Plotting ATL08 and ATL08 photon cloud
```r
# plot by "ph_h"
plot(atl03_atl08_dt, y="ph_h",colors, beam = "gt1l")

# plot by "h_ph"
plot(atl03_atl08_dt, y="h_ph",colors, beam = "gt1l")
```


# Acknowledgements
We gratefully acknowledge funding from NASA’s ICESat-2 (ICESat-2, grant 22-ICESat2_22-0006), Carbon Monitoring System (CMS, grant 22-CMS22-0015) and Commercial Smallsat Data Scientific Analysis(CSDSA, grant 22-CSDSA22_2-0080). 

# Reporting Issues 
Please report any issue regarding the rICESat2Veg package to Dr. Silva (c.silva@ufl.edu)

# Citing rICESat2Veg
Silva,C.A; Hamamura,C.rICESat2Veg: An R Package for NASA's Ice, Cloud, and Elevation Satellite (ICESat-2) Data Processing and Visualization for Terrestrial Applications.version 0.0.1, accessed on November. 22 2023, available at: <https://CRAN.R-project.org/package=rICESat2Veg>

# Disclaimer
**rICESat2Veg package has not been developted by the GEDI team. It comes with no guarantee, expressed or implied, and the authors hold no responsibility for its use or reliability of its outputs.**

