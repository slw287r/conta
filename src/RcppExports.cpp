// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

// intersect
void intersect(const char* tsvFileName, const char* outTsvFileName, const char* vcfFileName, bool nonDbSnp, bool DEBUG);
RcppExport SEXP conta_intersect(SEXP tsvFileNameSEXP, SEXP outTsvFileNameSEXP, SEXP vcfFileNameSEXP, SEXP nonDbSnpSEXP, SEXP DEBUGSEXP) {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const char* >::type tsvFileName(tsvFileNameSEXP);
    Rcpp::traits::input_parameter< const char* >::type outTsvFileName(outTsvFileNameSEXP);
    Rcpp::traits::input_parameter< const char* >::type vcfFileName(vcfFileNameSEXP);
    Rcpp::traits::input_parameter< bool >::type nonDbSnp(nonDbSnpSEXP);
    Rcpp::traits::input_parameter< bool >::type DEBUG(DEBUGSEXP);
    intersect(tsvFileName, outTsvFileName, vcfFileName, nonDbSnp, DEBUG);
    return R_NilValue;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"conta_intersect", (DL_FUNC) &conta_intersect, 5},
    {NULL, NULL, 0}
};

RcppExport void R_init_conta(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
