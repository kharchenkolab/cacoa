#include <RcppArmadillo.h>  
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;

// =======================================================
// 1. Pairwise distances
// =======================================================
// [[Rcpp::export]]
DataFrame computePairwiseDistances(NumericMatrix dist_mat,
                                             NumericMatrix model_mat,
                                             CharacterVector sample_names) {
  int n = dist_mat.nrow();
  int num_covs = model_mat.ncol();
  
  std::vector<std::string> Sample1, Sample2;
  std::vector<double> Distance;
  std::vector<std::vector<double>> cov_diff_cols(num_covs);
  
  CharacterVector cov_names = colnames(model_mat);
  
  for (int i = 0; i < n; ++i) {
    for (int j = i + 1; j < n; ++j) {
      Sample1.push_back(Rcpp::as<std::string>(sample_names[i]));
      Sample2.push_back(Rcpp::as<std::string>(sample_names[j]));
      Distance.push_back(dist_mat(i, j));
      
      for (int c = 0; c < num_covs; ++c) {
        double diff = std::abs(model_mat(i, c) - model_mat(j, c));
        cov_diff_cols[c].push_back(diff);
      }
    }
  }
  
  List out = List::create(
    _["Sample1"] = Sample1,
    _["Sample2"] = Sample2,
    _["Distance"] = Distance
  );
  
  for (int c = 0; c < num_covs; ++c) {
    std::string diff_name = Rcpp::as<std::string>(cov_names[c]) + "_diff";
    out[diff_name] = cov_diff_cols[c];
  }
  
  return DataFrame(out);
}

// =======================================================
// 2. Linear model fit
// =======================================================
// [[Rcpp::export]]
Rcpp::List lm_fit(const arma::mat& X, const arma::colvec& y) {
  int n = X.n_rows;
  int p = X.n_cols;
  
  arma::mat XtX = X.t() * X;
  if (arma::rank(XtX) < p) {
    stop("Design matrix is singular or rank deficient");
  }
  
  arma::mat XtX_inv = inv(XtX);
  arma::colvec beta = XtX_inv * X.t() * y;
  
  arma::colvec fitted = X * beta;
  arma::colvec residuals = y - fitted;
  
  int df = n - p;
  if (df <= 0) {
    stop("Degrees of freedom <= 0: not enough observations for predictors");
  }
  
  double sigma2 = arma::dot(residuals, residuals) / df;
  arma::colvec std_err = arma::sqrt(sigma2 * XtX_inv.diag());
  arma::colvec t_stats = beta / std_err;
  
  NumericVector p_values(beta.n_elem);
  for (int i = 0; i < beta.n_elem; ++i) {
    double t_val = t_stats(i);
    double p = 2 * R::pt(-std::abs(t_val), df, 1, 0);
    p_values[i] = p;
  }

  double rss = arma::dot(residuals, residuals);
  
  return Rcpp::List::create(
    Rcpp::Named("coefficients") = beta,
    Rcpp::Named("std_errors") = std_err,
    Rcpp::Named("t_values") = t_stats,
    Rcpp::Named("p_values") = p_values,
    Rcpp::Named("residuals") = residuals,
    Rcpp::Named("sigma2") = sigma2,
    Rcpp::Named("rss") = rss,
    Rcpp::Named("df") = df
  );
}

// =======================================================
// 3. PCA projection
// =======================================================
// [[Rcpp::export]]
arma::mat pca_project(const arma::mat& cm_norm, int n_pcs) {
  int n_rows = cm_norm.n_rows;
  int n_cols = cm_norm.n_cols;
  int min_dim = std::min(n_rows, n_cols) - 1;
  
  if (n_pcs > min_dim) {
    Rcpp::Rcout << "Warning: n_pcs is too large. Setting it to maximal allowed value " << min_dim << std::endl;
    n_pcs = min_dim;
  }
  
  arma::mat U;
  arma::vec s;
  arma::mat V;
  
  bool status = arma::svd_econ(U, s, V, cm_norm, "right");
  if (!status) {
    Rcpp::stop("SVD failed.");
  }
  
  arma::mat pcs = V.cols(0, n_pcs - 1);
  arma::mat projected = cm_norm * pcs;
  
  return projected;
}