#include <Rcpp.h>
#include <Rmath.h>
#include "scp.h"
#include "elbo.h"
#include "expectation_functions.h"
using namespace Rcpp;

//' Multivariate MICH Algorithm
//'
//' @description
//' Implementation of Algorithm 3 from Berlind, Cappello, and Madrid Padilla
//' (2025). This algorithm takes a sequence of \eqn{d}-dimensional observations
//' \eqn{\mathbf{y}_{1:T}}, and iteratively fits \eqn{L} mean-scp models
//' resulting in a variational approximation to the posterior distribution of
//' the \eqn{L} mean change-points. The algorithm terminates once the
//' percentage increase in the ELBO falls below `tol`.
//'
//' @param y A numeric matrix. \eqn{T \times N} matrix of observations.
//' @param mu_0 A numeric vector. Vector of intercept parameters.
//' @param fit_intercept A logical. If `fit_intercept == TRUE`, then an
//'   intercept is estimated and `mu_0` gets updated.
//' @param refit A logical. If `refit == TRUE`, then the MICH algorithm is
//'   initialized by the fit provided in `post_params`, otherwise the null
//'   model \eqn{\boldsymbol{\mu}_t = \mathbf{0}} and
//'   \eqn{\bar{\pi}_{\ell t} = 1/T} is used as the initialization.
//' @param max_iter An integer. Maximum number of iterations. If ELBO does not
//'   converge before `max_iter` is reached, then `converged == FALSE` in the
//'   returned fit object.
//' @param tol A scalar. Convergence tolerance for relative increase in ELBO.
//' @param verbose A logical. If `verbose == TRUE`, then value of the ELBO is
//'   printed every 5000th iteration.
//' @param omega_bar_l,log_omega_bar_l Numeric vectors. Vector of posterior
//'   precision parameters \eqn{\{\bar{\omega}_t\}_{t=1}^T} and log evaluations
//'   such that \eqn{V(\mathbf{b}_\ell|\tau=t) = \bar{\omega}_t\mathbf{I}_d}.
//' @param post_params A list. A length \eqn{L} list of the posterior parameters
//'   of each mean-scp component. Each element of the list is a list containing
//'   a \eqn{T\times L} matrix of posterior mean parameters and length \eqn{T}
//'   vectors of posterior change-point location probabilities and their log
//'   evaluations.
//'
//' @return A List. Parameters of the variational approximation the MICH
//' posterior distribution, including:
//'   * `residual`: A numeric vector. Residual \eqn{\mathbf{r}_{1:T}} after
//'     subtracting out each \eqn{E[\boldsymbol{\mu}_{\ell t}]} from
//'     \eqn{\mathbf{y}_{1:T}}.
//'   * `mu`: A numeric matrix. Posterior estimate of
//'     \eqn{\Sigma_{\ell=1}^L E[\boldsymbol{\mu}_{\ell,1:T}|\mathbf{y}_{1:T}]}.
//'   * `mu_0`: A numeric vector. Estimate of the intercept.
//'   * `post_params`: A list. List of posterior parameters for each mean-scp
//'     component.
//'   * `elbo`: A numeric vector. Value of the ELBO after each iteration.
//'   * `converged`: A logical. Indicates whether relative increase in the ELBO
//'     is less than `tol`.
//'
// [[Rcpp::export]]
List multi_mich_cpp(NumericMatrix y, NumericVector mu_0,
                    bool fit_intercept, bool refit,
                    double max_iter, double tol, bool verbose,
                    double omega_l, NumericMatrix log_pi_l,
                    NumericVector omega_bar_l, NumericVector log_omega_bar_l,
                    List post_params) {

  int T = y.nrow();
  int d = y.ncol();
  int L = post_params.length();

  // log parameters
  double log_omega_l = std::log(omega_l);

  // initialize residual
  NumericMatrix r_bar = clone(y);
  for (int t = 0; t < T; t++) {
    for (int i = 0; i < d; i++) {
      r_bar(t,i) -= mu_0[i];
    }
  }

  // initialize ELBO
  std::vector<double> elbo;
  elbo.reserve(std::min(max_iter, 1e6));
  elbo.push_back(R_NegInf);

  // initialize L expected mean parameters
  List post_params_l, mu_bar, mu2_bar, mean_scp_fit;
  NumericVector mu_var (T), mu_0_new (d);
  NumericMatrix mu_bar_l (T, d), mu2_bar_l (T, d), mu_var_l (T, L);

  for (int l = 0; l < L; l++) {
    post_params_l = post_params[l];
    if (refit) {
      mu_bar_l = multi_mu_bar_fn(post_params_l["b_bar"], post_params_l["pi_bar"]);
      mu2_bar_l = multi_mu2_bar_fn(post_params_l["b_bar"], omega_bar_l, post_params_l["pi_bar"]);
    }
    mu_bar.push_back(mu_bar_l);
    mu2_bar.push_back(mu2_bar_l);

    // initialize residual and variance terms
    for (int t = 0; t < T ; t++) {
      for (int i = 0; i < d; i++) {
        if (refit) {
          r_bar(t,i) -= mu_bar_l(t,i);
          mu_var_l(t,l) += mu2_bar_l(t,i) - mu_bar_l(t,i) * mu_bar_l(t,i);
        }
      }
    }
  }

  // vb algorithm
  int iter = 0;

  while (iter < max_iter) {

    // updating q(b_l, gamma_l)
    for (int l = 0; l < L; l++) {
      post_params_l = post_params[l];
      mu_bar_l = as<NumericMatrix>(mu_bar[l]);
      mu2_bar_l = as<NumericMatrix>(mu2_bar[l]);

      for (int t = 0; t < T; t++) {
        for (int i = 0; i < d; i++) {
          // add back l^th component of residual
          r_bar(t,i) += mu_bar_l(t,i);
          // add back l^th component of variance of mu
          mu_var_l(t,l) -= mu2_bar_l(t,i) - mu_bar_l(t,i) * mu_bar_l(t,i);
        }
      }

      // fit single multi_mean_scp model on residuals
      mean_scp_fit = multi_mean_scp(r_bar, omega_bar_l, log_omega_bar_l, log_pi_l(_,l));

      // store updated prior parameters
      post_params_l["b_bar"] = as<NumericMatrix>(mean_scp_fit["b_bar"]);
      post_params_l["pi_bar"] = as<NumericVector>(mean_scp_fit["pi_bar"]);
      post_params_l["log_pi_bar"] = as<NumericVector>(mean_scp_fit["log_pi_bar"]);
      post_params[l] = post_params_l;

      // update l^th component of mean
      mu_bar_l = multi_mu_bar_fn(post_params_l["b_bar"], post_params_l["pi_bar"]);
      mu_bar[l] = mu_bar_l;
      mu2_bar_l = multi_mu2_bar_fn(post_params_l["b_bar"], omega_bar_l, post_params_l["pi_bar"]);
      mu2_bar[l] = mu2_bar_l;

      for (int t = T - 1; t >= 0; t--) {
        for (int i = 0; i < d; i++){
          // subtract out l^th component of residual
          r_bar(t, i) -= mu_bar_l(t,i);
          // subtract out l^th component of variance of mu
          mu_var_l(t, l) += mu2_bar_l(t,i) - mu_bar_l(t,i) * mu_bar_l(t,i);
        }
      }
    }

    // update mu_0
    if (fit_intercept) {
      mu_0_new.fill(0.0);
      for (int i = 0; i < d; i++) {
        for (int t = 0; t < T; t++) {
          r_bar(t, i) += mu_0[i];
          mu_0_new[i] += r_bar(t, i);
        }

        mu_0[i] = mu_0_new[i] / T;
        for (int t = 0; t < T; t++) {
          r_bar(t, i) -= mu_0[i];
        }
      }
    }

    // calculate var(mu_t)
    mu_var = Rcpp::rowSums(mu_var_l);

    iter++;
    elbo.push_back(0.0);
    elbo[iter] = multi_elbo_fn(mu_0, r_bar, omega_l, log_omega_l, log_pi_l,
                               mu_var, post_params, omega_bar_l, log_omega_bar_l);
    if (verbose & (iter % 5000 == 0)) Rcout << "Iteration: " << iter << " elbo: " << elbo[iter] << "\n";
    if (std::abs((elbo[iter] - elbo[iter - 1]) / elbo[iter - 1]) < tol) break;
  }

  // construct mean signal
  NumericMatrix mu (T, d);
  for (int t = 0; t < T; t++) {
    mu(t,_) = y(t,_) - r_bar(t,_);
  }

  // creating lists of posterior parameters
  List result = List::create(_["residual"] = r_bar,
                             _["mu"] = mu,
                             _["mu_0"] = mu_0,
                             _["L"] = L,
                             _["post_params"] = post_params,
                             _["elbo"] = elbo,
                             _["converged"] = (max_iter > iter));
  return result;
}
