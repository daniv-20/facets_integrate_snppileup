#ifndef _MAIN_H_
#define _MAIN_H_

#include <ctime>
#include <iostream>
#include <sstream>
#include <vector>

#include "htslib/bgzf.h"
#include "htslib/hfile.h"
#include "htslib/hts.h"
#include "htslib/sam.h"
#include "htslib/faidx.h" 
#include "htslib/kstring.h"
#include "htslib/khash_str2int.h"
#include "htslib/synced_bcf_reader.h"
#include "htslib/vcf.h"

using namespace std;

// struct arguments {
// 	vector<char*> args;
// 	bool count_orphans;
// 	bool gzipped;
// 	BGZF* gzippedPointer;
// 	bool ignore_overlaps;
// 	int min_base_quality;
// 	int min_map_quality;
// 	vector<int> min_read_counts;
// 	int max_depth;
// 	bool progress;
// 	int pseudo_snps;
// 	bool verbose;
// 	int rflag_filter;
// 	void (*outFunc)(arguments, string, FILE *);
// };

struct file_info {
	int refs;
	int alts;
	int errors;
	int deletions;
};

typedef struct { 
    int min_mq, flag, min_baseQ, capQ_thres, max_depth, max_indel_depth, fmt_flag, all; 
    int rflag_require, rflag_filter; 
    int openQ, extQ, tandemQ, min_support; // for indels 
    double min_frac; // for indels 
    char *reg, *pl_list, *fai_fname, *output_fname; 
    faidx_t *fai; 
    void *bed, *rghash; 
    int argc; 
    char **argv;
} mplp_conf_t; 


#define MPLP_REF_INIT {{NULL,NULL},{-1,-1},{0,0}}
typedef struct { 
    char *ref[2]; 
    int ref_id[2]; 
    int ref_len[2]; 
} mplp_ref_t; 


typedef struct { 
    samFile *fp; 
	struct arguments *conf; 
    bam_hdr_t *h;
    mplp_ref_t *ref;
} mplp_aux_t; 


#endif
