#include <Rcpp.h>
// #include "C:/Users/vaithid1/OneDrive - Memorial Sloan Kettering Cancer Center/Repos/facets_integrate_snppileup/inst/extcode/snp-pileup-rev.h"
#include <htslib/synced_bcf_reader.h>
#include <htslib/hts.h>
#include <stdio.h>
#include <iostream>
#include <cstring>
#include <vector>
#include <cstdlib> // For exit()
#include <cstdio>  // For printf, FILE
#include <sstream>
#include <ctime>

using namespace Rcpp;

// [[Rcpp::export]]
void test_rhtslib()
{
  Rcpp::Rcout << "Rhtslib linked successfully!" << std::endl;
  printf("HTSlib version: %s\n", hts_version());
}
