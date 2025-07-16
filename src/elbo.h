#ifndef ELBO_H
#define ELBO_H

#include <Rcpp.h>
#include <Rmath.h>

double elbo_fn(int T, double mu_0, double lambda_0, Rcpp::NumericVector r_tilde,
               Rcpp::NumericVector lambda_bar, Rcpp::NumericVector delta,
               Rcpp::NumericMatrix b_bar_j, Rcpp::NumericMatrix omega_bar_j,
               Rcpp::NumericVector u_bar_j, Rcpp::NumericMatrix v_bar_j,
               Rcpp::NumericMatrix pi_bar_j, Rcpp::NumericMatrix log_pi_bar_j,
               Rcpp::NumericVector lgamma_u_bar_j,
               Rcpp::NumericVector digamma_u_bar_j, double omega_j, double u_j,
               double v_j, double log_omega_j, double log_u_j,
               double lgamma_u_j, double log_v_j, Rcpp::NumericMatrix log_pi_j,
               Rcpp::NumericMatrix b_bar_l, Rcpp::NumericMatrix omega_bar_l,
               Rcpp::NumericMatrix pi_bar_l, Rcpp::NumericMatrix log_pi_bar_l,
               double omega_l, double log_omega_l, Rcpp::NumericMatrix log_pi_l,
               Rcpp::NumericVector u_bar_k, Rcpp::NumericMatrix v_bar_k,
               Rcpp::NumericMatrix pi_bar_k, Rcpp::NumericVector lgamma_u_bar_k,
               Rcpp::NumericVector digamma_u_bar_k,
               Rcpp::NumericMatrix log_pi_bar_k, double u_k, double v_k,
               double log_u_k, double lgamma_u_k, double log_v_k,
               Rcpp::NumericMatrix log_pi_k);

double multi_elbo_fn(Rcpp::NumericVector mu_0, Rcpp::NumericMatrix r_bar,
                     double omega_l, double log_omega_l,
                     Rcpp::NumericMatrix log_pi_l, Rcpp::NumericVector mu_var,
                     Rcpp::List post_params, Rcpp::NumericVector omega_bar_l,
                     Rcpp::NumericVector log_omega_bar_l);

#endif
