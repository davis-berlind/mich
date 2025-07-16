#ifndef SCP_H
#define SCP_H

#include <Rcpp.h>
#include <Rmath.h>

Rcpp::List mean_scp(Rcpp::NumericVector y, Rcpp::NumericVector lambda,
                    double omega, Rcpp::NumericVector log_pi);

Rcpp::List multi_mean_scp(Rcpp::NumericMatrix y, Rcpp::NumericVector omega_bar,
                          Rcpp::NumericVector log_omega_bar,
                          Rcpp::NumericVector log_pi);

Rcpp::List var_scp(Rcpp::NumericVector y, Rcpp::NumericVector omega,
                   Rcpp::NumericVector u_bar, Rcpp::NumericVector lgamma_u_bar,
                   Rcpp::NumericVector v, Rcpp::NumericVector log_pi);

Rcpp::List meanvar_scp(Rcpp::NumericVector y, Rcpp::NumericVector lambda,
                       double omega, Rcpp::NumericVector u_bar,
                       Rcpp::NumericVector lgamma_u_bar, Rcpp::NumericVector v,
                       Rcpp::NumericVector log_pi);

#endif
