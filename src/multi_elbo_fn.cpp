#include <Rcpp.h>
#include <Rmath.h>
using namespace Rcpp;

const double log_2_pi = std::log(2 * M_PI);

//' MICH Multivariate ELBO
//'
//' @description
//' Computes the Evidence Lower Bound (ELBO) of the MICH model for a
//' multivariate sequence \eqn{\mathbf{y}_{1:T}}. See Appendix B of Berlind,
//' Cappello, and Madrid Padilla (2025) for detailed derivations.
//'
//' @param mu_0 A numeric vector. Intercept parameters.
//' @param r_bar A numeric matrix. Expected residual terms. See Proposition 3 of
//'   Berlind, Cappello, and Madrid Padilla (2025).
//' @param omega_l,log_omega_l Scalars. Prior precision parameter and log
//'   evaluation for mean-scp components of model.
//' @param log_pi_l A numeric matrix. A \eqn{T\times L} matrix of log
//'   prior change-point location probabilities for mean-scp components of
//'   model.
//' @param mu_var A numeric vector. Variance of mean signal, i.e.
//'   \eqn{E\lVert\boldsymbol{\mu}_t\rVert_2^2}.
//' @param post_params A list. A length \eqn{L} list of the posterior parameters
//'   of each mean-scp component. Each element of the list is a list containing
//'   a \eqn{T\times L} matrix of posterior mean parameters and length \eqn{T}
//'   vectors of posterior change-point location probabilities and their
//'   log evaluations.
//' @param omega_bar_l,log_omega_bar_l Numeric vectors. Vector of posterior
//'   precision parameters \eqn{\{\bar{\omega}_t\}_{t=1}^T} and log evaluations
//'   such that \eqn{V(\mathbf{b}_\ell|\tau=t) = \bar{\omega}_t\mathbf{I}_d}.
//'
//' @return A scalar. Value of ELBO for MICH model at current parameter values.
//'
double multi_elbo_fn(NumericVector mu_0, NumericMatrix r_bar,
                     double omega_l, double log_omega_l, NumericMatrix log_pi_l,
                     NumericVector mu_var, List post_params,
                     NumericVector omega_bar_l, NumericVector log_omega_bar_l) {

  int T = r_bar.nrow();
  int d = r_bar.ncol();
  int L = post_params.length();

  double elbo = -0.5 * d * T * log_2_pi;
  List post_params_l;
  NumericMatrix b_bar_l (T, d);
  NumericVector pi_bar_l (T), log_pi_bar_l (T);

  for (int t = 0; t < T; t++) {
    elbo += mu_var[t];
    for (int i = 0; i < d; i++) {
      elbo -= 0.5 * r_bar(t,i) * r_bar(t,i);
    }

    for (int l = 0; l < L; l++) {
      post_params_l = post_params[l];
      pi_bar_l = as<NumericVector>(post_params_l["pi_bar"]);
      log_pi_bar_l = as<NumericVector>(post_params_l["log_pi_bar"]);
      b_bar_l = as<NumericMatrix>(post_params_l["b_bar"]);

      // E[log p_l] - E[log q_l]
      // adding log odds with numerical stability for small pi_bar_l values
      if (pi_bar_l[t] > 1e-20) elbo += pi_bar_l[t] * (log_pi_l(t,l) - log_pi_bar_l[t]);
      elbo += 0.5 * d * pi_bar_l[t]  * (log_omega_l - log_omega_bar_l[t]);

      for (int i = 0; i < d; i++) {
        elbo += 0.5 * pi_bar_l[t] * (1 - omega_l * (b_bar_l(t,i) * b_bar_l(t,i) + 1 / omega_bar_l[t]));
      }
    }
  }
  return elbo;
}
