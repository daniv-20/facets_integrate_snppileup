library(Rhtslib)
library(tidyverse)

Sys.setenv(
  PKG_CPPFLAGS = "-I'C:/Users/vaithid1/AppData/Local/R/cache/R/renv/library/facets_integrate_snppileup-a88b4cc8/windows/R-4.4/x86_64-w64-mingw32/Rhtslib/include'",
  PKG_LIBS = "C:/Users/vaithid1/AppData/Local/R/cache/R/renv/library/facets_integrate_snppileup-a88b4cc8/windows/R-4.4/x86_64-w64-mingw32/Rhtslib/libs/x64/Rhtslib.dll"
)


Rcpp::sourceCpp("inst/extcode/rcpp_vers_snp-pileup-rev.cpp", verbose = TRUE)


