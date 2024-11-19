library(tidyverse)
#library(link2Rhtslib)
Sys.setenv(PATH = paste("G:/Projects/FACETS/htslib-1.21", Sys.getenv("PATH"), sep = ";"))
Sys.setenv(
  PKG_CPPFLAGS = "-IG:/Projects/FACETS/htslib-1.21",
  PKG_LIBS = "-LG:/Projects/FACETS/htslib-1.21 -lhts -LG:/RBuildTools/4.4/mingw_64/lib -lz -lbz2 -llzma -lcurl"
)


  # Compile and load the C++ code
Rcpp::sourceCpp("test.cpp", verbose = TRUE, rebuild = TRUE)

  # Call the function from your C++ file
test_rhtslib()

  ### builds test successfully with actual htslib!!

Rcpp::sourceCpp("C:/Users/vaithid1/OneDrive - Memorial Sloan Kettering Cancer Center/Repos/facets_integrate_snppileup/inst/extcode/rcpp_vers_snp-pileup-rev.cpp", verbose = TRUE, rebuild = TRUE)

snp_plp_test_rhtslib()
#snp-pileup <vcf file> <output file> <sequence files...>

datapath ="C:/Users/vaithid1/OneDrive - Memorial Sloan Kettering Cancer Center/Repos/facets_integrate_snppileup/inst/extdata"

humanvcf = "C:/Users/vaithid1/OnLaptop/00-common_all.vcf"

input_args = c(humanvcf, file.path(datapath, "snp-pileup-test-r-output.csv"), file.path(datapath,"sorted_c11.bam"))

start.time = proc.time()

rcpp_snp_pileup(input_args)

end.time = proc.time()

runtime = start.time - end.time

print(paste("Snp_pileup runtime = ", runtime))




## ----------------
# 
# 
# # Sys.setenv(
# #   PKG_CPPFLAGS = "-IG:/Projects/FACETS/htslib-1.21",
# #   PKG_LIBS = "-LG:/Projects/FACETS/htslib-1.21 -lhts -lz -lbz2 -llzma"
# # )
# 
# Sys.setenv(
#   PKG_CPPFLAGS = "-IG:/Projects/FACETS/htslib-1.21",
#   PKG_LIBS = "-LG:/Projects/FACETS/htslib-1.21 -lhts -LG:/RBuildTools/4.4/mingw_64/lib -lz -lbz2 -llzma -lcurl"
# )
# 
# 
# # Compile and load the C++ code
# Rcpp::sourceCpp("test.cpp", verbose = TRUE, rebuild = TRUE)
# 
# # Call the function from your C++ file
# test_rhtslib()
# 
# 
# Sys.setenv(PATH = paste("G:/Projects/FACETS/htslib-1.21", Sys.getenv("PATH"), sep = ";"))
# 
# 
# ### builds test successfully with actual htslib!!
# 
# Rcpp::sourceCpp("C:/Users/vaithid1/OneDrive - Memorial Sloan Kettering Cancer Center/Repos/facets_integrate_snppileup/inst/extcode/rcpp_vers_snp-pileup-rev.cpp", verbose = TRUE, rebuild = TRUE)
# 
# 
# 
# 
# 
# 
# ## ignore rhtslib for now -------------------
# 
# 
# # Set environment variables for compilation
# # Set environment variables for linking and compilation
# Sys.setenv(
#   PKG_CPPFLAGS = "-I'C:/Users/vaithid1/AppData/Local/R/cache/R/renv/library/facets_integrate_snppileup-a88b4cc8/windows/R-4.4/x86_64-w64-mingw32/Rhtslib/include'",
#   PKG_LIBS = "-L'C:/Users/vaithid1/AppData/Local/R/cache/R/renv/library/facets_integrate_snppileup-a88b4cc8/windows/R-4.4/x86_64-w64-mingw32/Rhtslib/usrlib/x64' -lhts"
# )
# 
# 
# Sys.getenv("PKG_CPPFLAGS")
# Sys.getenv("PKG_LIBS")
# 
# Sys.getenv("PATH")
# 
# 
# ##---------
# Sys.setenv(
#   PKG_CPPFLAGS = paste(
#     "-I'C:/Users/vaithid1/AppData/Local/R/cache/R/renv/library/facets_integrate_snppileup-a88b4cc8/windows/R-4.4/x86_64-w64-mingw32/Rhtslib/include'",
#     "-I'C:/RBuildTools/4.4/mingw_64/include'",
#     sep = " "
#   ),
#   PKG_LIBS = paste(
#     "-L'C:/Users/vaithid1/AppData/Local/R/cache/R/renv/library/facets_integrate_snppileup-a88b4cc8/windows/R-4.4/x86_64-w64-mingw32/Rhtslib/usrlib/x64' -lhts",
#     "-L'C:/RBuildTools/4.4/mingw_64/lib' -lz -lbz2 -llzma -lcurl",
#     sep = " "
#   )
# )
# 
# 
# 
# 
# # Compile and load the C++ code
# Rcpp::sourceCpp("test.cpp", verbose = TRUE)
# 
# # Call the function from your C++ file
# test_rhtslib()
# 
# 
# Sys.setenv(
#   PKG_CPPFLAGS = "-I'C:/Users/vaithid1/AppData/Local/R/cache/R/renv/library/facets_integrate_snppileup-a88b4cc8/windows/R-4.4/x86_64-w64-mingw32/Rhtslib/include'",
#   PKG_LIBS = paste(
#     "-L'C:/Users/vaithid1/AppData/Local/R/cache/R/renv/library/facets_integrate_snppileup-a88b4cc8/windows/R-4.4/x86_64-w64-mingw32/Rhtslib/usrlib/x64'",
#     "-lhts -lz -lbz2 -llzma -lcurl -lws2_32",
#     sep = " "
#   )
# )
# 
# Rcpp::sourceCpp("test.cpp", verbose = TRUE)
# test_rhtslib()
# 
# 
# 
# 
# 
# 
# ## installed link2Rhtslib
# ##devtools::install("C:\\Users\\vaithid1\\OneDrive - Memorial Sloan Kettering Cancer Center\\Projects\\FACETS\\link2Rhtslib")
# 
# # Sys.setenv(
# #   PKG_CPPFLAGS = "-I'C:/Users/vaithid1/AppData/Local/R/cache/R/renv/library/facets_integrate_snppileup-a88b4cc8/windows/R-4.4/x86_64-w64-mingw32/Rhtslib/include'",
# #   PKG_LIBS = "C:/Users/vaithid1/AppData/Local/R/cache/R/renv/library/facets_integrate_snppileup-a88b4cc8/windows/R-4.4/x86_64-w64-mingw32/Rhtslib/libs/x64/Rhtslib.dll"
# # )
# # Rcpp::sourceCpp("test.cpp", verbose = TRUE)
# # test_rhtslib()
# 
# #^^^ works but is missing bcf_sr*
