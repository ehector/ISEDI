// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;

// [[Rcpp::export]]
arma::mat matrix_inv(const arma::mat& X){
  arma::mat X_inv = inv(X);
  return(X_inv);
}

// [[Rcpp::export]]
double binom_ll_2(const arma::vec& par, const arma::vec& y_2, const arma::mat& X_2, const arma::vec& Z_2){
  double sum = 0;
  for(int i=0; i < y_2.size(); i++){
    sum += -log(1 + as_scalar(exp(X_2.row(i)*par))) ;
  }
  sum += as_scalar(Z_2.t()*par);
  return(sum);
}

// [[Rcpp::export]]
arma::vec binom_ll_2_deriv(const arma::vec& par, const arma::vec& y_2, const arma::mat& X_2, const arma::vec& Z_2){
  arma::vec score = zeros<vec>(X_2.n_cols);
  const arma::vec eta = exp(X_2*par)/(1+exp(X_2*par));
  
  for(int q=0; q < X_2.n_cols; q++){
    score(q) = as_scalar(X_2.col(q).t() * eta);
  }
  
  return(Z_2 - score);
}

// [[Rcpp::export]]
double binom_KLdiv(const arma::vec& par, const arma::vec& hat_beta_1, const arma::vec& y_1, const arma::mat& X_1, 
                   const arma::vec& hat_eta_1, const arma::vec& hat_mu_1){
  const arma::vec eta_1 = X_1*par;
  const arma::vec mu_1 = exp(eta_1)/(1+exp(eta_1));
  
  double sum = 0;
  for(int i=0; i<y_1.size(); i++){
    sum += hat_mu_1(i) * (hat_eta_1(i) - eta_1(i)) + log(1+exp(eta_1(i))) - log(1+exp(hat_eta_1(i)));
  }
  
  return(sum);
}

// [[Rcpp::export]]
arma::vec binom_KLdiv_deriv(const arma::vec& par, const arma::vec& hat_beta_1, const arma::vec& y_1, const arma::mat& X_1, 
                            const arma::vec& hat_eta_1, const arma::vec& hat_mu_1){
  arma::vec deriv = zeros<vec>(X_1.n_cols);
  const arma::vec eta_1 = X_1*par;
  const arma::vec mu_1 = exp(eta_1)/(1+exp(eta_1));
  
  for(int q=0; q<X_1.n_cols; q++){
    deriv(q) = as_scalar(X_1.col(q).t() * (mu_1 - hat_mu_1 ));
  }
  
  return(deriv);
}

// [[Rcpp::export]]
double binom_min_func(const arma::vec& par, const arma::vec& hat_beta_1, const arma::mat& X_1, const arma::mat& X_2,
                      const arma::vec& y_1, const arma::vec& y_2, const arma::vec& Z_2, const arma::vec& hat_eta_1,
                      const arma::vec& hat_mu_1, const int& n_1, const int& n_2, const double& lambda){
  double sum = -binom_ll_2(par, y_2, X_2, Z_2) / n_2 + lambda*binom_KLdiv(par, hat_beta_1, y_1, X_1, hat_eta_1, hat_mu_1)/n_1;
  return(sum);
}

// [[Rcpp::export]]
arma::vec binom_min_func_deriv(const arma::vec& par, const arma::vec& hat_beta_1, const arma::mat& X_1, const arma::mat& X_2,
                      const arma::vec& y_1, const arma::vec& y_2, const arma::vec& Z_2, const arma::vec& hat_eta_1,
                      const arma::vec& hat_mu_1, const int& n_1, const int& n_2, const double& lambda){
  arma::vec deriv = -binom_ll_2_deriv(par, y_2, X_2, Z_2) / n_2 + lambda*binom_KLdiv_deriv(par, hat_beta_1, y_1, X_1, hat_eta_1, hat_mu_1)/n_1;
  return(deriv);
}
