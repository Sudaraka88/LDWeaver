#include <Rcpp.h>
using namespace Rcpp;
// Function written to speed up ARACNE and to drop dependency on Rfast2

// [[Rcpp::plugins("cpp11")]]

// [[Rcpp::export(name = '.fast_intersect')]]
std::vector<int> fast_intersect(std::vector<int> A, std::vector<int> B) {

  std::vector<int> Av(A);
  std::vector<int> Bv(B);
  std::sort(std::begin(Av), std::end(Av));
  std::sort(std::begin(Bv), std::end(Bv));
  // std::vector<int> match;
  // match.reserve( std::min(Av.size(), Bv.size()) );
  std::vector<int> intersect;
  intersect.reserve( std::min(Av.size(), Bv.size()) );

  for (auto it1 = Av.begin(), it2 = Bv.begin();
       it1 != Av.end() && it2 != Bv.end();
       ++it2) {
    while (it1 != Av.end() && *it1 < *it2) ++it1;
    if (it1 != Av.end() && *it1 == *it2) {
      // match.push_back(it2 - Bv.begin());
      intersect.push_back(Bv[it2 - Bv.begin()]);
    }
  }

  return(intersect);

}

