#include <Rcpp.h>
#include <Rmath.h>
using namespace Rcpp;

// const double log_2_pi = std::log(2 * M_PI);

//' MICH Multivariate ELBO
//'
//' @description
//' Computes the Evidence Lower Bound (ELBO) of the MICH model for a
//' multivariate sequence \eqn{\mathbf{y}_{1:T}}. See Appendix of Berlind,
//' Cappello, and Madrid Padilla (2025) for detailed derivations.
//'
//' @param QTr A numeric matrix. Decorrelated residual terms.
//' @param post_params A list. A length \eqn{L} list of the posterior parameters
//'   of each mean-scp component. Each element of the list is a list containing
//'   a \eqn{T\times L} matrix of posterior mean parameters and length \eqn{T}
//'   vectors of posterior change-point location probabilities and their
//'   log evaluations.
//' @param lambda,log_lambda Numeric vectors. Eigenvalues of
//'   \eqn{\Lambda}.
//' @param trace_omega_omega_bar_inv Numeric vector. Trace of 
//'   \eqn{\Omega_{\ell t} \bar{\Omega}^{-1}_{\ell t}}.
//' @param log_det_omega_omega_bar Numeric vector. log determinant of 
//'   \eqn{\Omega_{\ell t} \bar{\Omega}^{-1}_{\ell t}}.
//' @param omega A numeric Matrix. Eigen values of \eqn{\Omega_{\ell t}}.
//' @param log_pi_l A numeric matrix. A \eqn{T\times L} matrix of log
//'   prior change-point location probabilities for mean-scp components of
//'   model.
//'
//' @return A scalar. Value of ELBO for MICH model at current parameter values.
//'
// [[Rcpp::export]]
double multi_elbo_fn(
  NumericMatrix QTr,
  List post_params,
  NumericVector lambda,
  NumericVector log_lambda,
  NumericVector trace_omega_omega_bar_inv,
  NumericVector log_det_omega_omega_bar,
  NumericMatrix omega,
  NumericMatrix log_pi_l
) {

  int T = QTr.nrow();
  int d = QTr.ncol();
  int L = post_params.length();

  // E[log p]
  double elbo = 0.0;
  for (int i = 0; i < d; i++) {
    elbo += T * log_lambda[i];
    for (int t = 0; t < T; t++) {
      elbo -= lambda[i] * QTr(t, i) * QTr(t, i);
    }
  }

  for (int l = 0; l < L; l++) {
    List post_params_l = post_params[l];
    NumericMatrix QTmu_bar_l = post_params_l["QTmu_bar"];
    NumericVector muTLmu_bar_l = post_params_l["muTLmu_bar"];
    for (int t = 0; t < T; t++) {
      elbo -= muTLmu_bar_l[t];
      for (int i = 0; i < d; i++) {
        elbo += lambda[i] * QTmu_bar_l(t, i) * QTmu_bar_l(t, i);
      }
    }

    // KL(q || p)
    NumericVector pi_bar_l = post_params_l["pi_bar"];
    NumericVector log_pi_bar_l = post_params_l["log_pi_bar"];
    NumericMatrix QTb_bar_l = post_params_l["QTb_bar"];
    if (pi_bar_l[T] > 1e-20) elbo -= 2.0 * pi_bar_l[T] * (log_pi_bar_l[T] - log_pi_l(T,l));
    for (int t = 0; t < T; t++) {
      if (pi_bar_l[t] > 1e-20) elbo -= 2.0 * pi_bar_l[t] * (log_pi_bar_l[t] - log_pi_l(t,l));
      elbo += pi_bar_l[t] * (log_det_omega_omega_bar[t] + d - trace_omega_omega_bar_inv[t]);
      for (int i = 0; i < d; i++) {
        elbo -= pi_bar_l[t] * omega(t,i) * QTb_bar_l(t, i) * QTb_bar_l(t, i);
      }
    }
  }

  return elbo;
}
