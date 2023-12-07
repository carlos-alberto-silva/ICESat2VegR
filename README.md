![](https://github.com/carlos-alberto-silva/rICESat2Veg/blob/master/readme/cover.png)<br/>

[![CRAN](https://www.r-pkg.org/badges/version/rICESat2Veg)](https://cran.r-project.org/package=rICESat2Veg)
![Github](https://img.shields.io/badge/Github-0.0.1-green.svg)
![licence](https://img.shields.io/badge/Licence-GPL--3-blue.svg) 
![Downloads](https://cranlogs.r-pkg.org/badges/grand-total/rICESat2Veg)
[![Build Status](https://travis-ci.com/carlos-alberto-silva/rICESat2Veg.svg?token=Jqizwyc6gBxNafNccTdU&branch=master)](https://travis-ci.com/carlos-alberto-silva/rICESat2Veg)

**rICESat2Veg: An R Package For NASA's Ice, Cloud, and Elevation Satellite (ICESat-2) Data Processing and Visualization For Terrestrial Applications.**

Authors: Carlos Alberto Silva and Caio Hamamura  

The rICESat2Veg package provides functions for downloading, processing, exporting and visualizing ICESat-2 ATL03 and ATL08 data.

# Getting Started

## Installation
```r
#The CRAN version:
install.packages("rICESat2Veg")

# The development version:
#install.packages("remotes")
library(remotes)
install_github("https://github.com/carlos-alberto-silva/rICESat2Veg", dependencies = TRUE)

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

## Reading ICESat-2 data
```r
atl08_h5<-readATL08(ATL08path=atl08_path)

```
## Extracting ATL08-derived Terrain Metrics
```r
terrain_metrics<-ATL08_terrain(atl08_h5=atl08_h5)
```

## Extracting ATL08-derived Canopy Metrics
```r
canopy_metrics<-ATL08_canopy(atl08_h5=atl08_h5)
```

# Acknowledgements
We gratefully acknowledge funding from NASA’s ICESat-2 (ICESat-2, grant 22-ICESat2_22-0006), Carbon Monitoring System (CMS, grant 22-CMS22-0015) and Commercial Smallsat Data Scientific Analysis(CSDSA, grant 22-CSDSA22_2-0080). 

# Reporting Issues 
Please report any issue regarding the rICESat2Veg package to Dr. Silva (c.silva@ufl.edu)

# Citing rICESat2Veg
Silva,C.A; Hamamura,C.rICESat2Veg: An R Package for NASA's Ice, Cloud, and Elevation Satellite (ICESat-2) Data Processing and Visualization for Terrestrial Applications.version 0.0.1, accessed on November. 22 2023, available at: <https://CRAN.R-project.org/package=rICESat2Veg>

# Disclaimer
**rICESat2Veg package comes with no guarantee, expressed or implied, and the authors hold no responsibility for its use or reliability of its outputs.**

