#include <Rcpp.h>
#include <Rmath.h>
using namespace Rcpp;

// const double log_2_pi = std::log(2 * M_PI);

//' MICH Multivariate ELBO
//'
//' @description
//' Computes the Evidence Lower Bound (ELBO) of the MICH model for a
//' multivariate sequence \eqn{\mathbf{y}_{1:T}}. See Appendix B of Berlind,
//' Cappello, and Madrid Padilla (2025) for detailed derivations.
//'
//' @param QTr A numeric matrix. Expected residual terms. See Proposition 3 of
//'   Berlind, Cappello, and Madrid Padilla (2025).
//' @param post_params A list. A length \eqn{L} list of the posterior parameters
//'   of each mean-scp component. Each element of the list is a list containing
//'   a \eqn{T\times L} matrix of posterior mean parameters and length \eqn{T}
//'   vectors of posterior change-point location probabilities and their
//'   log evaluations.
//' @param lambda log_lambda Numeric vectors. Eigenvalues of
//'   \eqn{\Lambda}.
//' @param Lambda_bar_log_det Numeric vector. Log determinant of
//'   \eqn{\bar{\Lambda}}
//' @param inv_Lambda_bar_trace Numeric vector. Reciprocal of trace of
//'   \eqn{\bar{\Lambda}}
//' @param log_pi_l A numeric matrix. A \eqn{T\times L} matrix of log
//'   prior change-point location probabilities for mean-scp components of
//'   model.
//' @param omega,d_log_omega Scalars. Prior precision parameters.
//'
//' @return A scalar. Value of ELBO for MICH model at current parameter values.
//'
double multi_elbo_fn(
  NumericMatrix QTr,
  List post_params,
  NumericVector lambda,
  NumericVector log_lambda,
  NumericVector Lambda_bar_log_det,
  NumericVector inv_Lambda_bar_trace,
  NumericMatrix log_pi_l,
  double omega,
  double d_log_omega
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
      elbo -= pi_bar_l[t] * (inv_Lambda_bar_trace[t] - d + Lambda_bar_log_det[t] - d_log_omega);
      for (int i = 0; i < d; i++) {
        elbo -= pi_bar_l[t] * omega * QTb_bar_l(t, i) * QTb_bar_l(t, i);
      }
    }
  }

  return elbo;
}
