// [[Rcpp::depends(RcppArmadillo)]]

#include <RcppArmadillo.h>



// determine mean difference between groups of columns, distances on that projection
// g1, g2 - vectors of indices for the two groups (0-based)
// [[Rcpp::export]]
arma::rowvec projdiff(const arma::mat & mat, const arma::ivec & g1, const arma::ivec & g2) {
  // TODO: remove with common expression shifts

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

// fit linear model for density differences
struct fit_density_lm_result {
  arma::mat  KP;
  arma::mat  KPge;
  arma::cube KP_perm;
  int        n_randomizations;
};

inline fit_density_lm_result fit_density_lm_impl(const arma::mat& M,
                                                 const arma::mat& P,
                                                 int n_randomizations) {
  arma::mat K  = arma::solve(M.t() * M, M.t());
  arma::mat KP = K * P;

  arma::mat  KPge(KP.n_rows, KP.n_cols, arma::fill::zeros);
  arma::cube KP_perm(KP.n_rows, KP.n_cols, n_randomizations);

  for (int r = 0; r < n_randomizations; ++r) {
    arma::mat P_prime = P.rows(arma::randperm(P.n_rows));
    arma::mat KP_prime = K * P_prime;
    KP_perm.slice(r) = KP_prime;
    KPge += arma::conv_to<arma::mat>::from(KP_prime >= KP);
  }
  return {KP, KPge, KP_perm, n_randomizations};
}

// [[Rcpp::export]]
Rcpp::List fit_density_lm(const arma::mat& M,
                          const arma::mat& P,
                          int n_randomizations) {
  auto res = fit_density_lm_impl(M, P, n_randomizations);
  return Rcpp::List::create(
    Rcpp::Named("KP") = res.KP,
    Rcpp::Named("KPge") = res.KPge,
    Rcpp::Named("KP_perm") = res.KP_perm,
    Rcpp::Named("n_randomizations") = res.n_randomizations
  );
}