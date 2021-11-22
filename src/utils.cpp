#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericMatrix colwiseBinaryDistance(IntegerMatrix mat) {
  NumericMatrix res(mat.cols(), mat.cols());
  for (int c1 = 0; c1 < mat.cols(); ++c1) {
    for (int c2 = c1 + 1; c2 < mat.cols(); ++c2) {
      double n_comp = 0, n_eq = 0;
      for (int r = 0; r < mat.rows(); ++r) {
        int v1 = mat(r, c1), v2 = mat(r, c2);
        if ((v1 == NA_INTEGER) || (v2 == NA_INTEGER))
          continue;

        n_eq += (v1 == v2);
        n_comp++;
      }

      if (n_comp < 0.1) {
          n_comp = 1;
      }

      res(c1, c2) = res(c2, c1) = 1 - n_eq / n_comp;
    }
  }

  rownames(res) = colnames(res) = colnames(mat);
  return res;
}