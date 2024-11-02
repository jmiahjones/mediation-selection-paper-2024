#include <RcppArmadillo.h>
using namespace arma;

// [[Rcpp::depends(RcppArmadillo)]]

// Function soft_thresh("soft_thresh");
// [[Rcpp::export]]
double soft_thresh(double rho, double lam) {
  if (fabs(rho) < lam) {
    return 0;
  } else {
    return (rho - (rho > 0 ? lam : -lam));
  }
}

// [[Rcpp::export]]
Rcpp::NumericVector coord_desc_cpp(const mat &X, const vec &Y, double lam, const vec &w) {
  double w_tol_zero = 1e8;
  int max_iter = 2e9;
  double coef_tol = 1e-5;
  int n = X.n_rows;
  int p = X.n_cols;
  vec theta_old = arma::zeros(p);
  vec theta = arma::zeros(p);
  int iter = 0;
  bool theta_conv = false;
  double conv_crit = INFINITY;
  while (iter < max_iter && !theta_conv) {
    iter++;
    theta_old = theta;
    for (int j = 0; j < p; j++) {
      if(w(j) <= w_tol_zero){
	vec residj = Y - X * theta + theta(j) * X.col(j);
	double rhoj = arma::mean(X.col(j) % residj);
        theta(j) = soft_thresh(rhoj, lam * w(j));
      }
    }
    conv_crit = arma::max(arma::abs(theta - theta_old));
    theta_conv = (conv_crit < coef_tol);
  }
  return Rcpp::wrap(theta);
}


// [[Rcpp::export]]
Rcpp::NumericVector coord_desc_cpp_obsweight(const mat &X, const vec &Y, double lam, const vec &w, const vec &obsW, const vec &thetaInit) {
  double w_tol_zero = 1e8;
  int max_iter = 2e9;
  double coef_tol = 1e-5;
  int n = X.n_rows;
  int p = X.n_cols;
  vec theta_old = thetaInit;
  vec theta = thetaInit;
  int iter = 0;
  bool theta_conv = false;
  double conv_crit = INFINITY;
  while (iter < max_iter && !theta_conv) {
    iter++;
    theta_old = theta;
    for (int j = 0; j < p; j++) {
      if(w(j) <= w_tol_zero){
        vec residj = obsW % (Y - X * theta + theta(j) * X.col(j));
        double rhoj = arma::mean(X.col(j) % residj);
        theta(j) = soft_thresh(rhoj, lam * w(j));
      }
    }
    conv_crit = arma::max(arma::abs(theta - theta_old));
    theta_conv = (conv_crit < coef_tol);
  }
  return Rcpp::wrap(theta);
}

/*** R
coord_desc <- function(X,Y,lam,w,obsWeight=NULL,thetaInit=NULL){
  X <- scale(X,T,T)
  y <- scale(Y,T,F)
  z <- attributes(X)$`scaled:scale`
  internal_lam <-  lam / n

  p <- ncol(X)
  if(is.null(thetaInit)) {
    thetaInit <- rep(0, p)
  }
  # else {
  #   # stopifnot(is.numeric(thetaInit),
  #   #           length(thetaInit)==p)
  # }

  if(!is.null(obsWeight)){
    stopifnot(is.numeric(obsWeight),
              length(obsWeight) == nrow(X),
              !any(is.na(obsWeight)))
    res <- coord_desc_cpp_obsweight(X,y,internal_lam,w,obsWeight,thetaInit)
  } else {
    res <- coord_desc_cpp(X,y,internal_lam,w)
  }
  theta <- res
  theta <- theta / z
  return(theta)
}
*/

/*
 *
# old R version that generated this on ChatGPT

 coord_desc_R <- function(X,Y,lam,w,
 w_tol_zero=1e8,
 max_iter=2e9L,
 coef_tol=1e-5){
 p <- ncol(X)
 theta_old <- rep(Inf, p)
 theta <- rnorm(p)
 X <- scale(X,T,T)
 y <- scale(Y,T,F)
 z <- attributes(X)$`scaled:scale`
 iter=0L
 theta_conv <- F
 conv_crit <- Inf
 force_zero <- (w > w_tol_zero)
 j_idxs <- (1:p)[!force_zero]
 theta[force_zero] <- 0
# browser()
 while(iter < max_iter && !theta_conv){
 iter = iter + 1L
 theta_old <- theta
 for(j in j_idxs){ # only iterate over those not forced to zero
# if(!force_zero[j]) {
 residj <- (y - as.numeric(X[,-j] %*% theta[-j]))
 rhoj <- mean(X[,j] * residj)
 if(is.na(rhoj) || is.na(lam*w[j])){
 browser()
 }
 theta[j] <- soft_thresh(rhoj, lam * w[j])
# }
 }

 conv_crit <- max(abs(theta - theta_old))
 theta_conv <- (conv_crit < coef_tol)
 }
# browser()
 theta <- theta / z
 return(list(coef=theta, convergence=theta_conv, iters=iter))
 }
 *
 */
