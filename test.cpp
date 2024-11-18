#include <Rcpp.h>
//#include "C:/Users/vaithid1/OneDrive - Memorial Sloan Kettering Cancer Center/Repos/facets_integrate_snppileup/inst/extcode/snp-pileup-rev.h"
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
void test_rhtslib() {
  Rcpp::Rcout << "Rhtslib linked successfully!" << std::endl;
  printf("HTSlib version: %s\n", hts_version());
}

uint64_t get_snp_count(char *file)
{
  uint64_t count = 0;
  // load vcf file
  bcf_srs_t *vcfReader = bcf_sr_init();
  if (!bcf_sr_add_reader(vcfReader, file))
  {
    printf("Failed to read VCF file: %s\n", bcf_sr_strerror(vcfReader->errnum));
    return 1;
  }
  while (bcf_sr_next_line(vcfReader))
  {
    count++;
  }

  bcf_sr_destroy(vcfReader);
  return count;
}

