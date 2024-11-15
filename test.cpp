#include <Rcpp.h>


#include "htslib/bgzf.h"
#include "htslib/hfile.h"
#include "htslib/hts.h"
#include "htslib/sam.h"
#include "htslib/faidx.h" 
#include "htslib/kstring.h"
#include "htslib/khash_str2int.h"
#include "htslib/synced_bcf_reader.h"
#include "htslib/vcf.h"

using namespace Rcpp;

// [[Rcpp::export]]
void test_rhtslib() {
  Rcpp::Rcout << "Rhtslib linked successfully!" << std::endl;
}

