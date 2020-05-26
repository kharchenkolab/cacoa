#include <RcppArmadillo.h>

using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
arma::sp_mat calcZScores(const arma::sp_mat& nonref_local_mean, const arma::sp_mat& ref_local_mean, arma::sp_mat sds, const double minsd) {
  arma::mat c = arma::mat(nonref_local_mean) - arma::mat(ref_local_mean);
  arma::mat d = arma::mat(sds).transform( [& minsd](double val) {
    if(val > minsd) { return val; } else { return minsd; };
  } );
  arma::mat res = c / d;
  return arma::sp_mat(res);
}

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
arma::sp_mat matsub(const arma::sp_mat& a, const arma::sp_mat& b) {
  return a - b;
}
