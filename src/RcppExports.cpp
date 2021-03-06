// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

// massFisherTest
NumericMatrix massFisherTest(const LogicalMatrix& experiments, const IntegerVector& sums, const List& elements, std::string altString);
RcppExport SEXP metaRE_massFisherTest(SEXP experimentsSEXP, SEXP sumsSEXP, SEXP elementsSEXP, SEXP altStringSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const LogicalMatrix& >::type experiments(experimentsSEXP);
    Rcpp::traits::input_parameter< const IntegerVector& >::type sums(sumsSEXP);
    Rcpp::traits::input_parameter< const List& >::type elements(elementsSEXP);
    Rcpp::traits::input_parameter< std::string >::type altString(altStringSEXP);
    rcpp_result_gen = Rcpp::wrap(massFisherTest(experiments, sums, elements, altString));
    return rcpp_result_gen;
END_RCPP
}
// quickFisherTest
NumericVector quickFisherTest(NumericVector eff1, NumericVector n1, NumericVector eff2, NumericVector n2, std::string alternative);
RcppExport SEXP metaRE_quickFisherTest(SEXP eff1SEXP, SEXP n1SEXP, SEXP eff2SEXP, SEXP n2SEXP, SEXP alternativeSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type eff1(eff1SEXP);
    Rcpp::traits::input_parameter< NumericVector >::type n1(n1SEXP);
    Rcpp::traits::input_parameter< NumericVector >::type eff2(eff2SEXP);
    Rcpp::traits::input_parameter< NumericVector >::type n2(n2SEXP);
    Rcpp::traits::input_parameter< std::string >::type alternative(alternativeSEXP);
    rcpp_result_gen = Rcpp::wrap(quickFisherTest(eff1, n1, eff2, n2, alternative));
    return rcpp_result_gen;
END_RCPP
}
// enumerateMotifsCpp
SEXP enumerateMotifsCpp(List parameters, Function createGCS, Function logDebug);
RcppExport SEXP metaRE_enumerateMotifsCpp(SEXP parametersSEXP, SEXP createGCSSEXP, SEXP logDebugSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< List >::type parameters(parametersSEXP);
    Rcpp::traits::input_parameter< Function >::type createGCS(createGCSSEXP);
    Rcpp::traits::input_parameter< Function >::type logDebug(logDebugSEXP);
    rcpp_result_gen = Rcpp::wrap(enumerateMotifsCpp(parameters, createGCS, logDebug));
    return rcpp_result_gen;
END_RCPP
}
