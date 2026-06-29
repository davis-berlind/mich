#ifndef SCP_H
#define SCP_H

#include <Rcpp.h>
#include <Rmath.h>

Rcpp::List mean_scp(
    Rcpp::NumericVector y,
    Rcpp::NumericVector lambda,
    double omega,
    double log_omega,
    Rcpp::NumericVector log_pi
);

Rcpp::List multi_mean_scp(
    Rcpp::NumericMatrix QTr,
    Rcpp::NumericVector lambda,
    Rcpp::NumericMatrix mean_weights,
    Rcpp::NumericMatrix sandwich_weights,
    Rcpp::NumericVector log_prob_weights
);

Rcpp::List var_scp(
    Rcpp::NumericVector y, 
    Rcpp::NumericVector omega,
    Rcpp::NumericVector u_bar,
    Rcpp::NumericVector lgamma_u_bar,
    Rcpp::NumericVector v, 
    double log_v,
    Rcpp::NumericVector log_pi
);

Rcpp::List meanvar_scp(
    Rcpp::NumericVector y,
    Rcpp::NumericVector lambda,
    double omega,
    double log_omega,
    Rcpp::NumericVector u_bar,
    Rcpp::NumericVector lgamma_u_bar, 
    Rcpp::NumericVector v,
    double log_v,
    Rcpp::NumericVector log_pi
);

#endif
