#include "snp-pileup.h"
#include <htslib/hts.h>
#include <stdio.h>
#include <argp.h>

int main()
{
    printf("HTSlib version: %s\n", hts_version());
    return 0;
}