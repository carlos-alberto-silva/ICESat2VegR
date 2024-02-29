library(Rcpp)

setwd("C:/Users/caioh/src/cpp/ANN/")
system("make realclean")
system("make linux-g++")

Sys.setenv(PATH=paste0(Sys.getenv("PATH"), ";C:/Users/caioh/src/cpp/ann/lib"))
Sys.setenv(PKG_LIBS="-LC:/Users/caioh/src/cpp/ann/lib -lANN")
Sys.setenv(PKG_CPPFLAGS="-IC:/Users/caioh/src/cpp/ann/include")

setwd("C:/Users/caioh/src/r/rICESat2Veg/")

Rcpp::cppFunction(
  paste(readLines("src/functions.cpp"), collapse='\n'),
  rebuild = TRUE
)


FindNearestNeighbors(1:10, 1:10, 2, 6)
