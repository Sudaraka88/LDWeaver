#include <Rcpp.h>
using namespace Rcpp;
// Function written to speed up ARACNE and to drop dependency on Rfast2

// [[Rcpp::export(name = '.fast_intersect')]]
std::vector<int> fast_intersect(std::vector<int> A, std::vector<int> B) {

  std::vector<int> Av(A);
  std::vector<int> Bv(B);
  std::sort(std::begin(Av), std::end(Av));
  std::sort(std::begin(Bv), std::end(Bv));

  std::vector<int> intersect;
  intersect.reserve( std::min(Av.size(), Bv.size()) );

  int i = 0; // Index for vector1
  int j = 0; // Index for vector2
  while (i < Av.size() && j < Bv.size()) {
    if (Av[i] < Bv[j]) {
      i++;
    } else if (Av[i] > Bv[j]) {
      j++;
    } else { // Found an intersection element
      intersect.push_back(Av[i]);
      i++;
      j++;
    }
  }

  return(intersect);

}
