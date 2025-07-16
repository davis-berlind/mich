#include <Rcpp.h>
#include <Rmath.h>
#include "scp.h"
#include "elbo.h"
#include "expectation_functions.h"
using namespace Rcpp;

//' MICH Algorithm
//'
//' @description
//' Implementation of Algorithm 1 & 2 from Berlind, Cappello, and Madrid Padilla
//' (2025). This algorithm takes a sequence of \eqn{T} observations
//' \eqn{\mathbf{y}_{1:T}}, and iteratively fits \eqn{L} mean-scp models,
//' \eqn{K} var-scp models, and \eqn{J} meanvar-scp models resulting in a
//' variational approximation to the posterior distribution of
//' the \eqn{L} mean, \eqn{K} variance, and \eqn{J} mean and variance
//' change-points. The algorithm terminates once the percentage increase in the
//' ELBO falls below `tol`.
//'
//' @param y A numeric vector. Length \eqn{T} vector of observations.
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
List mich_cpp(NumericVector y,
              int J, int L, int K, double mu_0, double lambda_0,
              bool fit_intercept, bool fit_scale,
              bool refit, double max_iter, bool verbose, double tol,
              double omega_j, double u_j, double v_j, NumericMatrix log_pi_j,
              NumericMatrix pi_bar_j, NumericMatrix log_pi_bar_j,
              NumericMatrix b_bar_j, NumericMatrix omega_bar_j,
              NumericVector u_bar_j, NumericMatrix v_bar_j,
              NumericVector lgamma_u_bar_j, NumericVector digamma_u_bar_j,
              double omega_l, NumericMatrix log_pi_l,
              NumericMatrix pi_bar_l, NumericMatrix log_pi_bar_l,
              NumericMatrix b_bar_l, NumericMatrix omega_bar_l,
              double u_k, double v_k, NumericMatrix log_pi_k,
              NumericMatrix pi_bar_k, NumericMatrix log_pi_bar_k,
              NumericVector u_bar_k, NumericMatrix v_bar_k,
              NumericVector lgamma_u_bar_k, NumericVector digamma_u_bar_k) {

  int T = y.length();

  // log parameters
  double log_omega_j = std::log(omega_j);
  double log_v_j = std::log(v_j); double log_u_j = std::log(u_j);
  double lgamma_u_j = std::lgamma(u_j);
  double log_omega_l = std::log(omega_l);
  double log_u_k = std::log(u_k); double log_v_k = std::log(v_k);
  double lgamma_u_k = std::lgamma(u_k);

  // initialize residual
  NumericVector r_tilde = clone(y);

  // initialize precision
  NumericVector lambda_bar (T, lambda_0);

  // initialize correction term
  NumericVector delta (T, 0.0);

  // initialize J expected mean-variance parameters
  NumericMatrix mu_lambda_j (T, J), mu2_lambda_j (T, J), lambda_bar_j (T, J);
  double mu_bar_jt;
  if (refit) {
    for (int j = 0; j < J; j++) {
      mu_lambda_j(_,j) = mu_lambda_fn(b_bar_j(_,j), u_bar_j, v_bar_j(_,j), pi_bar_j(_,j));
      mu2_lambda_j(_,j) = mu2_lambda_fn(b_bar_j(_,j), omega_bar_j(_,j), u_bar_j, v_bar_j(_,j), pi_bar_j(_,j));
      lambda_bar_j(_,j) = lambda_bar_fn(u_bar_j, v_bar_j(_,j), pi_bar_j(_,j));
    }
  } else {
    lambda_bar_j.fill(1.0);
  }

  // initialize L expected mean parameters
  NumericMatrix mu_bar_l (T, L), mu2_bar_l (T, L);
  if (refit) {
    for (int l = 0; l < L; l++) {
      mu_bar_l(_,l) = mu_bar_fn(b_bar_l(_,l), pi_bar_l(_,l));
      mu2_bar_l(_,l) = mu2_bar_fn(b_bar_l(_,l), omega_bar_l(_,l), pi_bar_l(_,l));
    }
  }

  // initialize K expected variance parameters
  NumericMatrix lambda_bar_k (T, K);
  if (refit) {
    for (int k = 0; k < K; k++) {
      lambda_bar_k(_,k) = lambda_bar_fn(u_bar_k, v_bar_k(_,k), pi_bar_k(_,k));
    }
  } else {
    lambda_bar_k.fill(1.0);
  }

  // initialize combined residual, variance, and correction term
  for (int t = 0; t < T; t++) {
    r_tilde[t] -= mu_0;

    if (refit) {
      for (int j = 0; j < J; j++) {
        mu_bar_jt = mu_lambda_j(t,j) / lambda_bar_j(t,j);
        r_tilde[t] -= mu_bar_jt;
        lambda_bar[t] *= lambda_bar_j(t,j);
        delta[t] -= mu_bar_jt * mu_bar_jt - mu2_lambda_j(t,j) / lambda_bar_j(t,j);
      }

      for (int l = 0; l < L; l++) {
        r_tilde[t] -= mu_bar_l(t,l);
        delta[t] -= mu_bar_l(t,l) * mu_bar_l(t,l) - mu2_bar_l(t,l);
      }

      for (int k = 0; k < K; k++) {
        lambda_bar[t] *= lambda_bar_k(t,k);
      }
    }
  }

  // initialize ELBO
  std::vector<double> elbo;
  elbo.reserve(std::min(max_iter, 1e6));
  elbo.push_back(R_NegInf);

  // initialize model fits
  List mean_scp_fit, var_scp_fit, meanvar_scp_fit;

  // initialize correction priors
  NumericVector v_tilde (T, 0.0), log_pi_tilde (T, 0.0);
  double v_tilde_sum, log_pi_tilde_sum, log_pi_tilde_max;

  // vb algorithm
  int iter = 0;

  while (iter < max_iter) {

    // updating q(b_l, gamma_l)
    for (int l = 0; l < L; l++) {
      for (int t = T - 1; t >= 0; t--) {
        // add back l^th component of residual
        r_tilde[t] += mu_bar_l(t,l);
        // subtract l^th component of variance correction
        delta[t] += mu_bar_l(t,l) * mu_bar_l(t,l) - mu2_bar_l(t,l);
        delta[t] = std::max(0.0, delta[t]);
      }

      // fit single mean_scp model on residuals
      mean_scp_fit = mean_scp(r_tilde, lambda_bar, omega_l, log_pi_l(_,l));

      // store updated posterior parameters
      pi_bar_l(_,l) =  as<NumericVector>(mean_scp_fit["pi_bar"]);
      log_pi_bar_l(_,l) =  as<NumericVector>(mean_scp_fit["log_pi_bar"]);
      b_bar_l(_,l) =  as<NumericVector>(mean_scp_fit["b_bar"]);
      omega_bar_l(_,l) =  as<NumericVector>(mean_scp_fit["omega_bar"]);

      // update l^th component of mean
      mu_bar_l(_,l) = mu_bar_fn(b_bar_l(_,l), pi_bar_l(_,l));
      mu2_bar_l(_,l) = mu2_bar_fn(b_bar_l(_,l), omega_bar_l(_,l), pi_bar_l(_,l));

      for (int t = 0; t < T; t++) {
        // subtract out l^th component of residual
        r_tilde[t] -= mu_bar_l(t,l);
        // add back l^th component of variance correction
        delta[t] -= mu_bar_l(t,l) * mu_bar_l(t,l) - mu2_bar_l(t,l);
        delta[t] = std::max(0.0, delta[t]);
      }
    }

    // updating q(s_k, gamma_k)
    for (int k = 0; k < K; k++) {
      v_tilde_sum = 0.0;
      log_pi_tilde_sum = log_pi_k(0,k);
      log_pi_tilde_max = R_NegInf;

      for (int t = T - 1; t >= 0; t--) {
        // divide out k^th component of precision
        lambda_bar[t] /= lambda_bar_k(t,k);

        // calculate corrected priors
        v_tilde_sum += 0.5 * lambda_bar[t] * delta[t];
        v_tilde[t] = v_k + v_tilde_sum;
      }

      for (int t = 1; t < T; t++) {
        log_pi_tilde_sum += -0.5 * lambda_bar[t-1] * delta[t-1];
        log_pi_tilde[t] = log_pi_k(t,k) + log_pi_tilde_sum;
        log_pi_tilde_max = std::max(log_pi_tilde_max, log_pi_tilde[t]);
      }

      for (int t = 1; t < T; t++) {
        log_pi_tilde[t] -= log_pi_tilde_max;
      }

      // fit single var_scp model on residuals
      var_scp_fit = var_scp(r_tilde, lambda_bar, u_bar_k, lgamma_u_bar_k,
                            v_tilde, log_pi_tilde);

      // store updated posterior parameters
      pi_bar_k(_,k) = as<NumericVector>(var_scp_fit["pi_bar"]);
      log_pi_bar_k(_,k) = as<NumericVector>(var_scp_fit["log_pi_bar"]);
      v_bar_k(_,k) = as<NumericVector>(var_scp_fit["v_bar"]);

      // updated k^th component precision
      lambda_bar_k(_,k) = lambda_bar_fn(u_bar_k, v_bar_k(_,k), pi_bar_k(_,k));

      for (int t = 0; t < T; t++) {
        // multiply back k^th component or precision
        lambda_bar[t] *= lambda_bar_k(t,k);
      }
    }

    // updating q(b_j, s_j, gamma_j)
    for (int j = 0; j < J; j++) {
      v_tilde_sum = 0.0;
      log_pi_tilde_sum = log_pi_j(0,j);
      log_pi_tilde_max = R_NegInf;

      for (int t = T-1; t >= 0; t--) {
        // add back j^th component of residual
        mu_bar_jt = mu_lambda_j(t,j) / lambda_bar_j(t,j);
        r_tilde[t] += mu_bar_jt;
        // divide out j^th component of precision
        lambda_bar[t] /= lambda_bar_j(t,j);
        // subtract j^th component of variance correction
        delta[t] += mu_bar_jt * mu_bar_jt - mu2_lambda_j(t,j) / lambda_bar_j(t,j);
        delta[t] = std::max(0.0, delta[t]);

        // calculate corrected priors
        v_tilde_sum += 0.5 * lambda_bar[t] * delta[t];
        v_tilde[t] = v_j + v_tilde_sum;
      }

      for (int t = 1; t < T; t++) {
        log_pi_tilde_sum += -0.5 * lambda_bar[t-1] * delta[t-1];
        log_pi_tilde[t] = log_pi_j(t,j) + log_pi_tilde_sum;
        log_pi_tilde_max = std::max(log_pi_tilde_max, log_pi_tilde[t]);
      }

      for (int t = 1; t < T; t++) {
        log_pi_tilde[t] -= log_pi_tilde_max;
      }

      // fit single meanvar_scp model on residuals
      meanvar_scp_fit = meanvar_scp(r_tilde, lambda_bar, omega_j, u_bar_j,
                                    lgamma_u_bar_j, v_tilde, log_pi_tilde);

      // store updated posterior parameters
      pi_bar_j(_,j) =  as<NumericVector>(meanvar_scp_fit["pi_bar"]);
      log_pi_bar_j(_,j) =  as<NumericVector>(meanvar_scp_fit["log_pi_bar"]);
      b_bar_j(_,j) =  as<NumericVector>(meanvar_scp_fit["b_bar"]);
      omega_bar_j(_,j) =  as<NumericVector>(meanvar_scp_fit["omega_bar"]);
      v_bar_j(_,j) =  as<NumericVector>(meanvar_scp_fit["v_bar"]);

      // update j^th component of mean and precision
      mu_lambda_j(_,j) = mu_lambda_fn(b_bar_j(_,j), u_bar_j, v_bar_j(_,j), pi_bar_j(_,j));
      mu2_lambda_j(_,j) = mu2_lambda_fn(b_bar_j(_,j), omega_bar_j(_,j), u_bar_j, v_bar_j(_,j), pi_bar_j(_,j));
      lambda_bar_j(_,j) = lambda_bar_fn(u_bar_j, v_bar_j(_,j), pi_bar_j(_,j));

      for (int t = 0; t < T; t++) {
        mu_bar_jt = mu_lambda_j(t,j) / lambda_bar_j(t,j);
        // subtract out j^th component of residual
        r_tilde[t] -= mu_bar_jt;
        // multiply back j^th component or precision
        lambda_bar[t] *= lambda_bar_j(t,j);
        // add back j^th component of variance correction
        delta[t] -= mu_bar_jt * mu_bar_jt - mu2_lambda_j(t,j) / lambda_bar_j(t,j);
        delta[t] = std::max(0.0, delta[t]);
      }
    }

    // updating mu_0 and lambda_0
    if (fit_intercept) {
      r_tilde += Rcpp::rep(mu_0, T);
      mu_0 = Rcpp::sum(lambda_bar * r_tilde) / Rcpp::sum(lambda_bar);
      r_tilde += Rcpp::rep(-mu_0, T);
    }

    if (fit_scale) {
      for (int t = 0; t < T; t++) {
        lambda_bar[t] /= lambda_0;
      }

      lambda_0 = T / Rcpp::sum(lambda_bar * (Rcpp::pow(r_tilde, 2) + delta));

      for (int t = 0; t < T; t++) {
        lambda_bar[t] *= lambda_0;
      }
    }

    iter++;
    elbo.push_back(0.0);

    // calculate elbo
    elbo[iter] = elbo_fn(T, mu_0, lambda_0,
                         r_tilde, lambda_bar, delta,
                         b_bar_j, omega_bar_j, u_bar_j, v_bar_j, pi_bar_j, log_pi_bar_j,
                         lgamma_u_bar_j, digamma_u_bar_j,
                         omega_j, u_j, v_j, log_omega_j, log_u_j, lgamma_u_j, log_v_j, log_pi_j,
                         b_bar_l, omega_bar_l, pi_bar_l, log_pi_bar_l,
                         omega_l, log_omega_l, log_pi_l,
                         u_bar_k, v_bar_k, pi_bar_k,
                         lgamma_u_bar_k, digamma_u_bar_k, log_pi_bar_k,
                         u_k, v_k, log_u_k, lgamma_u_k, log_v_k, log_pi_k);

    if (std::isnan(elbo[iter])) throw std::runtime_error("NaN in elbo");
    if (verbose & (iter % 5000 == 0)) Rcout << "Iteration: " << iter << "; elbo: " << elbo[iter] << ";\n";
    // terminate if relative increase in elbo falls below tol
    if (std::abs((elbo[iter] - elbo[iter - 1]) / elbo[iter - 1]) < tol) break;
  }

  // construct mean signal
  NumericVector mu (T, mu_0);
  for (int j = 0; j < J; j++) {
    mu += mu_bar_fn(b_bar_j(_,j), pi_bar_j(_,j));
  }
  for (int l = 0; l < L; l++) {
    mu += mu_bar_fn(b_bar_l(_,l), pi_bar_l(_,l));
  }

  // creating lists of posterior parameters
  List J_model, L_model, K_model;

  List result = List::create(_["y"] = y,
                             _["residual"] = r_tilde,
                             _["mu"] = mu,
                             _["lambda"] = lambda_bar,
                             _["delta"] = delta,
                             _["converged"] = (max_iter > iter),
                             _["elbo"] = elbo,
                             _["mu_0"] = mu_0,
                             _["lambda_0"] = lambda_0,
                             _["J"] = J,
                             _["L"] = L,
                             _["K"] = K);

  if (J > 0) {
    J_model =  List::create(_["pi_bar"] = pi_bar_j,
                            _["b_bar"] = b_bar_j,
                            _["omega_bar"] = omega_bar_j,
                            _["v_bar"] = v_bar_j,
                            _["u_bar"] = u_bar_j,
                            _["mu_lambda_bar"] = mu_lambda_j,
                            _["mu2_lambda_bar"] = mu2_lambda_j,
                            _["lambda_bar"] = lambda_bar_j);
    result["meanvar_model"] = J_model;
  }

  if (L > 0) {
    L_model =  List::create(_["pi_bar"] = pi_bar_l,
                            _["b_bar"] = b_bar_l,
                            _["omega_bar"] = omega_bar_l,
                            _["mu_bar"] = mu_bar_l,
                            _["mu2_bar"] = mu2_bar_l);
    result["mean_model"] = L_model;
  }

  if (K > 0) {
    K_model =  List::create(_["pi_bar"] = pi_bar_k,
                            _["v_bar"] = v_bar_k,
                            _["u_bar"] = u_bar_k,
                            _["lambda_bar"] = lambda_bar_k);
    result["var_model"] = K_model;
  }

  return result;
}
