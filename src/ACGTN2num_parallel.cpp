#include <Rcpp.h>
#ifdef _OPENMP
#include <omp.h>
// [[Rcpp::plugins(openmp)]]
#endif

using namespace Rcpp;

// [[Rcpp::export(name = '.ACGTN2num')]]
void ACGTN2num(NumericMatrix nv, StringVector cv, int ncores){
  // This function is used to mask the reference allele of each snp in a []5xnsnp] matrix
  const int len = cv.length();
  //NumericMatrix nv(5, len);
  //Rcout << len  << '\n';
#ifdef _OPENMP
#pragma omp parallel for num_threads(ncores)
#endif
  for(int c=0; c<len; c++){
    char cc = as<char>(cv[c]);
    //Rcout << cc << '\n';
    if(cc == 'A'){
      nv[c*5] = 0;
      continue;
    } else if(cc == 'C'){
      nv[c*5 + 1] = 0;
      continue;
    } else if(cc == 'G'){
      nv[c*5 + 2] = 0;
      continue;
    } else if(cc == 'T'){
      nv[c*5 + 3] = 0;
      continue;
    } else if(cc == 'N'){
      nv[c*5 + 4] = 0;
      //continue;
    } else if(cc == '-'){
      nv[c*5 + 4] = 0;
      //continue;
    }

  }
  //return(nv);
}
