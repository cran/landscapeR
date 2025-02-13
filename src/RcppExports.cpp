// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// contigCells_cpp
IntegerVector contigCells_cpp(int pt, int bgr, NumericMatrix mtx);
RcppExport SEXP _landscapeR_contigCells_cpp(SEXP ptSEXP, SEXP bgrSEXP, SEXP mtxSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type pt(ptSEXP);
    Rcpp::traits::input_parameter< int >::type bgr(bgrSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type mtx(mtxSEXP);
    rcpp_result_gen = Rcpp::wrap(contigCells_cpp(pt, bgr, mtx));
    return rcpp_result_gen;
END_RCPP
}
// assignValues_cpp
NumericMatrix assignValues_cpp(int val, IntegerVector ad, NumericMatrix mtx);
RcppExport SEXP _landscapeR_assignValues_cpp(SEXP valSEXP, SEXP adSEXP, SEXP mtxSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type val(valSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type ad(adSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type mtx(mtxSEXP);
    rcpp_result_gen = Rcpp::wrap(assignValues_cpp(val, ad, mtx));
    return rcpp_result_gen;
END_RCPP
}
// indexTranspose_cpp
IntegerVector indexTranspose_cpp(IntegerVector id, int dim1, int dim2);
RcppExport SEXP _landscapeR_indexTranspose_cpp(SEXP idSEXP, SEXP dim1SEXP, SEXP dim2SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< IntegerVector >::type id(idSEXP);
    Rcpp::traits::input_parameter< int >::type dim1(dim1SEXP);
    Rcpp::traits::input_parameter< int >::type dim2(dim2SEXP);
    rcpp_result_gen = Rcpp::wrap(indexTranspose_cpp(id, dim1, dim2));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_landscapeR_contigCells_cpp", (DL_FUNC) &_landscapeR_contigCells_cpp, 3},
    {"_landscapeR_assignValues_cpp", (DL_FUNC) &_landscapeR_assignValues_cpp, 3},
    {"_landscapeR_indexTranspose_cpp", (DL_FUNC) &_landscapeR_indexTranspose_cpp, 3},
    {NULL, NULL, 0}
};

RcppExport void R_init_landscapeR(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
