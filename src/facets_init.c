#include <R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>

// Declare external functions
void htslib_version();
SEXP _facets_htslib_version(SEXP);

// Create a table of R-callable routines
static const R_CallMethodDef CallEntries[] = {
    {"_facets_htslib_version", (DL_FUNC) &_facets_htslib_version, 1},
    {NULL, NULL, 0}
};

// Initialize the package
void R_init_facets(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, TRUE); // Enable dynamic lookup
}
