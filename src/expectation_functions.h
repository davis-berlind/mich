#ifndef EXPECTATION_FUNCTIONS_H
#define EXPECTATION_FUNCTIONS_H

#include <Rcpp.h>
#include <Rmath.h>

Rcpp::NumericMatrix multi_mu_bar_fn(Rcpp::NumericMatrix b,
                                    Rcpp::NumericVector prob);

Rcpp::NumericMatrix multi_mu2_bar_fn(Rcpp::NumericMatrix b,
                                     Rcpp::NumericVector omega,
                                     Rcpp::NumericVector prob);

Rcpp::NumericVector mu_bar_fn(Rcpp::NumericVector b,
                              Rcpp::NumericVector prob);

Rcpp::NumericVector mu2_bar_fn(Rcpp::NumericVector b,
                               Rcpp::NumericVector omega,
                               Rcpp::NumericVector prob);

Rcpp::NumericVector lambda_bar_fn(Rcpp::NumericVector u,
                                  Rcpp::NumericVector v,
                                  Rcpp::NumericVector prob);

Rcpp::NumericVector mu_lambda_fn(Rcpp::NumericVector b,
                                 Rcpp::NumericVector u,
                                 Rcpp::NumericVector v,
                                 Rcpp::NumericVector prob);

Rcpp::NumericVector mu2_lambda_fn(Rcpp::NumericVector b,
                                  Rcpp::NumericVector omega,
                                  Rcpp::NumericVector u,
                                  Rcpp::NumericVector v,
                                  Rcpp::NumericVector prob);

#endif
