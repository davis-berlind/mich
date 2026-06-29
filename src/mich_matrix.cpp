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
//' @keywords internal
//'
//' @param y A numeric matrix. \eqn{T \times d} matrix of observations.
//' @param mu_0 A numeric vector. Vector of intercept parameters.
//' @param lambda A numeric vector. Eigenvalues of \eqn{\Lambda}.
//' @param log_lambda A numeric vector. Log eigenvalues of \eqn{\Lambda}.
//' @param Q A numeric matrix. Eigenvectors of \eqn{\Lambda}.
//' @param Lambda_bar_log_det Numeric vector. Log determinant of
//'   \eqn{\bar{\Lambda}}
//' @param inv_Lambda_bar_trace Numeric vector. Reciprocal of trace of
//'   \eqn{\bar{\Lambda}}
//' @param mean_weights A numeric matrix. Weights for calculating b_bar.
//' @param sandwich_weights A numeric matrix. Weights for calculating norm of
//'   b_bar.
//' @param log_prob_weights A numeric matrix. Constant part of log_pi_bar.
//' @param fit_intercept A logical. If `fit_intercept == TRUE`, then an
//'   intercept is estimated and `mu_0` gets updated.
//' @param fit_scale A logical. If `fit_scale == TRUE`, then \eqn{\Lambda} is
//'   updated on each vb iteration.
//' @param max_iter An integer. Maximum number of iterations. If ELBO does not
//'   converge before `max_iter` is reached, then `converged == FALSE` in the
//'   returned fit object.
//' @param tol A scalar. Convergence tolerance for relative increase in ELBO.
//' @param verbose A logical. If `verbose == TRUE`, then the value of the ELBO
//'   is printed every 5000th iteration.
//' @param omega_l, d_log_omeal_l Scalars. Prior precision parameters for
//'   mean-scp components
//'   of model.
//' @param log_pi_l A numeric matrix. A \eqn{T \times L} matrix of prior log
//'   change-point location probabilities for each of the \eqn{L} mean
//'   change-points.
//' @param post_params A list. A length \eqn{L} list of the posterior parameters
//'   of each mean-scp component. Each element of the list is a list containing
//'   a \eqn{T\times L} matrix of posterior mean parameters and length \eqn{T}
//'   vectors of posterior change-point location probabilities and their log
//'   evaluations.
//'
//' @return A List. Parameters of the variational approximation the MICH
//' posterior distribution, including:
//'   * `L`: An integer. Number of components included in model.
//'   * `lambda`: A numeric vector. Eigenvalues of \eqn{\Lambda}.
//'   * `Q`: A numeric matrix. Eigenvectors of \eqn{\Lambda}.
//'   * `residual`: A numeric matrix. Residual \eqn{\mathbf{r}_{1:T}} after
//'     subtracting out each \eqn{E[\boldsymbol{\mu}_{\ell t}]} from
//'     \eqn{\mathbf{y}_{1:T}}.
//'   * `mu_0`: A numeric vector. Estimate of the intercept.
//'   * `post_params`: A list. List of posterior parameters for each mean-scp
//'     component.
//'   * `elbo`: A numeric vector. Value of the ELBO after each iteration.
//'   * `converged`: A logical. Indicates whether relative increase in the ELBO
//'     is less than `tol`.
//'
// [[Rcpp::export]]
List mich_matrix_cpp(
    NumericMatrix y,
    NumericVector mu_0,
    NumericVector lambda,
    NumericVector log_lambda,
    NumericMatrix Q,
    NumericVector Lambda_bar_log_det,
    NumericVector inv_Lambda_bar_trace,
    NumericMatrix mean_weights,
    NumericMatrix sandwich_weights,
    NumericMatrix log_prob_weights,
    bool fit_intercept,
    bool fit_scale,
    double max_iter,
    double tol,
    bool verbose,
    double omega_l,
    double log_omega_l,
    double d_log_omega_l,
    NumericMatrix log_pi_l,
    List post_params
  ) {

  int T = y.nrow();
  int d = y.ncol();
  int L = post_params.length();

  // initialize residual
  NumericVector QTmu_0 (d);
  NumericMatrix r_bar = clone(y);
  NumericMatrix QTr_bar (T, d); // initialize decorrelated residual

  // subtract out mu_0
  for (int t = 0; t < T; t++) {
    for (int i = 0; i < d; i++) {
      for (int j = 0; j < d; j++) {
        if (t == 0) QTmu_0[i] += mu_0[j] * Q(j, i);
        QTr_bar(t, i) += r_bar(t, j) * Q(j, i);
      }
      QTr_bar(t, i) -= QTmu_0[i];
    }
  }

  // subtract out mu_l
  for (int l = 0; l < L; l++) {
    List post_params_l = post_params[l];
    NumericMatrix QTmu_bar_l = post_params_l["QTmu_bar"];
    for (int t = 0; t < T; t++) {
      for (int i = 0; i < d; i++) {
        QTr_bar(t, i) -= QTmu_bar_l(t, i);
      }
    }
  }

  // initialize ELBO
  std::vector<double> elbo;
  elbo.reserve(std::min(max_iter, 1e6));
  elbo.push_back(R_NegInf);

  // vb algorithm
  int iter = 1;

  while (iter <= max_iter) {

    // update covariance matrix
    if (fit_scale){
      // todo: refit
      // NumericVector lambda,
      // NumericVector log_lambda,
      // NumericMatrix Q,
      // NumericVector Lambda_bar_log_det,
      // NumericVector inv_Lambda_bar_trace,
      // NumericMatrix mean_weights,
      // NumericMatrix sandwich_weights,
      // NumericMatrix log_prob_weights,
    }

    // updating q(b_l, gamma_l)
    for (int l = 0; l < L; l++) {
      // add back lth partial mean
      List post_params_l = post_params[l];
      NumericMatrix QTmu_bar_l = post_params_l["QTmu_bar"];
      for (int t = 0; t < T; t++) {
        for (int i = 0; i < d; i++) {
          QTr_bar(t, i) += QTmu_bar_l(t, i);
        }
      }

      // update posterior parameters
      post_params_l = multi_mean_scp(
        QTr_bar, lambda, mean_weights, sandwich_weights, log_prob_weights(_,l)
      );

      QTmu_bar_l = as<NumericMatrix>(post_params_l["QTmu_bar"]);
      post_params[l] = post_params_l;

      for (int t = 0; t < T; t++) {
        for (int i = 0; i < d; i++) {
          QTr_bar(t, i) -= QTmu_bar_l(t, i);
        }
      }
    }

    // update mu_0
    if (fit_intercept) {
      for (int i = 0; i < d; i++) {
        double mn = 0.0;
        for (int t = 0; t < T; t++) {
          QTr_bar(t, i) += QTmu_0[i];
          mn += QTr_bar(t, i);
        }
        QTmu_0[i] = mn / T;
        for (int t = 0; t < T; t++) {
          QTr_bar(t, i) -= QTmu_0[i];
        }
      }
    }

    elbo.push_back(0.0);
    elbo[iter] = multi_elbo_fn(
      QTr_bar, post_params, lambda, log_lambda,
      Lambda_bar_log_det, inv_Lambda_bar_trace,
      log_pi_l, omega_l, d_log_omega_l
    );

    if (verbose & (iter % 5000 == 0)) Rcout << "Iteration: " << iter << " elbo: " << elbo[iter] << "\n";
    if (std::abs((elbo[iter] - elbo[iter - 1]) / elbo[iter - 1]) < tol) break;
    iter++;
  }

  // correlated intercept
  for (int i = 0; i < d; i++) {
    mu_0[i] = 0.0;
    for (int j = 0; j < d; j++) {
      mu_0[i] += Q(i, j) * QTmu_0[j];
    }
  }

  // creating lists of posterior parameters
  List result = List::create(
    _["L"] = L,
    _["mu_0"] = mu_0,
    _["lambda"] = lambda,
    _["Q"] = Q,
    _["post_params"] = post_params,
    _["QTr"] = QTr_bar,
    _["elbo"] = elbo,
    _["converged"] =  (max_iter > iter)
  );

  return result;
}
