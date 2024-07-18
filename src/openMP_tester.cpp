#include <Rcpp.h>
#ifdef _OPENMP
#include <omp.h>
// [[Rcpp::plugins(openmp)]]
#endif


// [[Rcpp::export]]
void test_openmp() {
#ifdef _OPENMP
  int nthreads = omp_get_max_threads();
  Rcpp::Rcout << "with OpenMP support - system supports up to: " << nthreads << " cores" << std::endl;
#else
  Rcpp::Rcout << "without OpenMP support - WARNING! Some functions will run slower..." << std::endl;
#endif
}
