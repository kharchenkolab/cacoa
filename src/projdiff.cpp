// [[Rcpp::depends(RcppArmadillo)]]

#include <RcppArmadillo.h>



// determine mean difference between groups of columns, distances on that projection
// g1, g2 - vectors of indices for the two groups (0-based)
// [[Rcpp::export]]
arma::rowvec projdiff(const arma::mat & mat, const arma::ivec & g1, const arma::ivec & g2) {
  // mat: columns - samples; rows - genes;
  // determine consensus difference
  arma::mat dm(mat.n_rows,g1.n_elem * g2.n_elem);

  for(int i=0;i<g1.n_elem;i++) {
    for(int j=0;j<g2.n_elem;j++) {
      dm.col(i*g2.n_elem + j) = mat.col(g1[i])-mat.col(g2[j]);
    }
  }
  arma::vec dmm=mean(dm,1); // would be nice to do trimming here
  dmm/=sqrt(sum(dmm % dmm));
  return dmm.t() * dm;

}

