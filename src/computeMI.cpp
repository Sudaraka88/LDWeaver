#include <Rcpp.h>

#ifdef _OPENMP
#include <omp.h>
// [[Rcpp::plugins(openmp)]]
#endif

using namespace Rcpp;

// [[Rcpp::export(name = '.fastHadamard')]]
void fastHadamard(NumericMatrix MIt, NumericMatrix den, NumericMatrix uq_t, NumericMatrix pxy_t, NumericMatrix pxpy_t, NumericMatrix RXY, NumericMatrix pXrX, NumericMatrix pYrY, int ncores){
  // function to perform the final Hadamard operation (faster than the base-R version, but can it be done better?)
  const int nrow = MIt.nrow();
  const int ncol = MIt.ncol();
  int counter_n = nrow*ncol;
#ifdef _OPENMP
#pragma omp parallel for num_threads(ncores)
#endif
  for(int counter = 0; counter < counter_n; counter++) MIt[counter] += uq_t[counter]*pxy_t[counter]/den[counter]*log(pxy_t[counter]/(pxpy_t[counter]+ RXY[counter] + pXrX[counter] + pYrY[counter])*den[counter]);

}


// [[Rcpp::export(name = '.compareToRow')]]
LogicalVector compareToRow(NumericMatrix x, NumericVector y) {
  const int nr = x.nrow();
  const int nc = x.ncol();
  const int nt = y.length();
  LogicalVector ret(nr, false);
  for (int j=0; j < nr; ++j) {
    for (int k=0; k < nc; ++k) {
      for(int m=0; m<nt; ++m){
        if (x(j, k) == y[m]) {
          ret[j] = true;
          break;
        }
      }
    }
  }
  return ret;
}

// [[Rcpp::export(name = '.vecPosMatch')]]
NumericVector vecPosMatch(NumericVector x, NumericVector y) {
  // find pos of x in y
  int nx = x.length();
  int ny = y.length();

  NumericVector ret(nx);
  for(int ii=0; ii < nx; ++ii){
    for(int jj=0; jj<ny; ++jj){
      if(y[jj] == x[ii]){
        ret[ii] = jj+1;
        break;
      }
    }
  }
  return ret;
}


// [[Rcpp::export(name = '.compareTriplet')]]
bool compareTriplet(NumericVector MI0X, NumericVector MI0Z, double MI0){
  // This function still stands for ARACNE, we probably need to redo this using sparseMx
  bool ARACNE = true;
  const int nl = MI0X.length();
  for(int i=0; i<nl; ++i){
    //  Rprintf("%d", i);
    if(MI0 < MI0X[i]){
      if(MI0 < MI0Z[i]){
        ARACNE = false;
        break;
      }
    }
  }
  return ARACNE;
}
