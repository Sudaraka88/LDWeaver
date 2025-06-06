// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// ACGTN2num
void ACGTN2num(NumericMatrix nv, StringVector cv, int ncores);
RcppExport SEXP _LDWeaver_ACGTN2num(SEXP nvSEXP, SEXP cvSEXP, SEXP ncoresSEXP) {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type nv(nvSEXP);
    Rcpp::traits::input_parameter< StringVector >::type cv(cvSEXP);
    Rcpp::traits::input_parameter< int >::type ncores(ncoresSEXP);
    ACGTN2num(nv, cv, ncores);
    return R_NilValue;
END_RCPP
}
// fastHadamard
void fastHadamard(NumericMatrix MIt, NumericMatrix den, NumericMatrix uq_t, NumericMatrix pxy_t, NumericMatrix pxpy_t, NumericMatrix RXY, NumericMatrix pXrX, NumericMatrix pYrY, int ncores);
RcppExport SEXP _LDWeaver_fastHadamard(SEXP MItSEXP, SEXP denSEXP, SEXP uq_tSEXP, SEXP pxy_tSEXP, SEXP pxpy_tSEXP, SEXP RXYSEXP, SEXP pXrXSEXP, SEXP pYrYSEXP, SEXP ncoresSEXP) {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type MIt(MItSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type den(denSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type uq_t(uq_tSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type pxy_t(pxy_tSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type pxpy_t(pxpy_tSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type RXY(RXYSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type pXrX(pXrXSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type pYrY(pYrYSEXP);
    Rcpp::traits::input_parameter< int >::type ncores(ncoresSEXP);
    fastHadamard(MIt, den, uq_t, pxy_t, pxpy_t, RXY, pXrX, pYrY, ncores);
    return R_NilValue;
END_RCPP
}
// compareToRow
LogicalVector compareToRow(NumericMatrix x, NumericVector y);
RcppExport SEXP _LDWeaver_compareToRow(SEXP xSEXP, SEXP ySEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type x(xSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type y(ySEXP);
    rcpp_result_gen = Rcpp::wrap(compareToRow(x, y));
    return rcpp_result_gen;
END_RCPP
}
// vecPosMatch
NumericVector vecPosMatch(NumericVector x, NumericVector y);
RcppExport SEXP _LDWeaver_vecPosMatch(SEXP xSEXP, SEXP ySEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type x(xSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type y(ySEXP);
    rcpp_result_gen = Rcpp::wrap(vecPosMatch(x, y));
    return rcpp_result_gen;
END_RCPP
}
// compareTriplet
bool compareTriplet(NumericVector MI0X, NumericVector MI0Z, double MI0);
RcppExport SEXP _LDWeaver_compareTriplet(SEXP MI0XSEXP, SEXP MI0ZSEXP, SEXP MI0SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type MI0X(MI0XSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type MI0Z(MI0ZSEXP);
    Rcpp::traits::input_parameter< double >::type MI0(MI0SEXP);
    rcpp_result_gen = Rcpp::wrap(compareTriplet(MI0X, MI0Z, MI0));
    return rcpp_result_gen;
END_RCPP
}
// fast_intersect
std::vector<int> fast_intersect(std::vector<int> A, std::vector<int> B);
RcppExport SEXP _LDWeaver_fast_intersect(SEXP ASEXP, SEXP BSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< std::vector<int> >::type A(ASEXP);
    Rcpp::traits::input_parameter< std::vector<int> >::type B(BSEXP);
    rcpp_result_gen = Rcpp::wrap(fast_intersect(A, B));
    return rcpp_result_gen;
END_RCPP
}
// extractAlnParam
List extractAlnParam(std::string file, int filter, double gap_thresh, double maf_thresh);
RcppExport SEXP _LDWeaver_extractAlnParam(SEXP fileSEXP, SEXP filterSEXP, SEXP gap_threshSEXP, SEXP maf_threshSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< std::string >::type file(fileSEXP);
    Rcpp::traits::input_parameter< int >::type filter(filterSEXP);
    Rcpp::traits::input_parameter< double >::type gap_thresh(gap_threshSEXP);
    Rcpp::traits::input_parameter< double >::type maf_thresh(maf_threshSEXP);
    rcpp_result_gen = Rcpp::wrap(extractAlnParam(file, filter, gap_thresh, maf_thresh));
    return rcpp_result_gen;
END_RCPP
}
// extractSNPs
List extractSNPs(std::string file, int n_seq, int n_snp, std::vector<int> POS);
RcppExport SEXP _LDWeaver_extractSNPs(SEXP fileSEXP, SEXP n_seqSEXP, SEXP n_snpSEXP, SEXP POSSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< std::string >::type file(fileSEXP);
    Rcpp::traits::input_parameter< int >::type n_seq(n_seqSEXP);
    Rcpp::traits::input_parameter< int >::type n_snp(n_snpSEXP);
    Rcpp::traits::input_parameter< std::vector<int> >::type POS(POSSEXP);
    rcpp_result_gen = Rcpp::wrap(extractSNPs(file, n_seq, n_snp, POS));
    return rcpp_result_gen;
END_RCPP
}
// extractRef
List extractRef(std::string file);
RcppExport SEXP _LDWeaver_extractRef(SEXP fileSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< std::string >::type file(fileSEXP);
    rcpp_result_gen = Rcpp::wrap(extractRef(file));
    return rcpp_result_gen;
END_RCPP
}
// test_openmp
void test_openmp();
RcppExport SEXP _LDWeaver_test_openmp() {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    test_openmp();
    return R_NilValue;
END_RCPP
}
// readFasta
List readFasta(std::string file, int pos_len);
RcppExport SEXP _LDWeaver_readFasta(SEXP fileSEXP, SEXP pos_lenSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< std::string >::type file(fileSEXP);
    Rcpp::traits::input_parameter< int >::type pos_len(pos_lenSEXP);
    rcpp_result_gen = Rcpp::wrap(readFasta(file, pos_len));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_LDWeaver_ACGTN2num", (DL_FUNC) &_LDWeaver_ACGTN2num, 3},
    {"_LDWeaver_fastHadamard", (DL_FUNC) &_LDWeaver_fastHadamard, 9},
    {"_LDWeaver_compareToRow", (DL_FUNC) &_LDWeaver_compareToRow, 2},
    {"_LDWeaver_vecPosMatch", (DL_FUNC) &_LDWeaver_vecPosMatch, 2},
    {"_LDWeaver_compareTriplet", (DL_FUNC) &_LDWeaver_compareTriplet, 3},
    {"_LDWeaver_fast_intersect", (DL_FUNC) &_LDWeaver_fast_intersect, 2},
    {"_LDWeaver_extractAlnParam", (DL_FUNC) &_LDWeaver_extractAlnParam, 4},
    {"_LDWeaver_extractSNPs", (DL_FUNC) &_LDWeaver_extractSNPs, 4},
    {"_LDWeaver_extractRef", (DL_FUNC) &_LDWeaver_extractRef, 1},
    {"_LDWeaver_test_openmp", (DL_FUNC) &_LDWeaver_test_openmp, 0},
    {"_LDWeaver_readFasta", (DL_FUNC) &_LDWeaver_readFasta, 2},
    {NULL, NULL, 0}
};

RcppExport void R_init_LDWeaver(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
