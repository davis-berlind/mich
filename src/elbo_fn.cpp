#include <Rcpp.h>
#include <Rmath.h>
using namespace Rcpp;

const double log_2_pi = std::log(2 * M_PI);

//' MICH Univariate ELBO
//'
//' @description
//' Computes the Evidence Lower Bound (ELBO) of the MICH model for a univariate
//' sequence \eqn{\mathbf{y}_{1:T}}. See Appendix B of Berlind, Cappello, and
//' Madrid Padilla (2025) for detailed derivations.
//'
//' @param T An integer. The length of \eqn{\mathbf{y}_{1:T}}.
//' @param mu_0 A scalar. Intercept parameter.
//' @param lambda_0 A scalar. Initial precision parameter.
//' @param r_tilde A numeric vector. Expected residual term. See (B.4) of
//'   Berlind, Cappello, and  Madrid Padilla (2025).
//' @param lambda_bar A numeric vector. Expected precision term. See (B.5) of
//'   Berlind, Cappello, and  Madrid Padilla (2025).
//' @param delta A numeric vector. Variance correction term. See (B.6) of
//'   Berlind, Cappello, and  Madrid Padilla (2025).
//' @param b_bar_j A numeric matrix. A \eqn{T \times J} matrix of posterior mean
//'   parameters for meanvar-scp components of model.
//' @param omega_bar_j A numeric matrix. A \eqn{T\times J} matrix of posterior
//'   variance parameters for meanvar-scp components of model.
//' @param u_bar_j,lgamma_u_bar_j,digamma_u_bar_j Numeric vectors. A length
//'   \eqn{T} vector of posterior shape parameters and log gamma and digamma
//'   evaluations for meanvar-scp components of model.
//' @param v_bar_j A numeric matrix. A \eqn{T\times J} matrix of posterior
//'   rate parameters for meanvar-scp components of model.
//' @param pi_bar_j,log_pi_bar_j Numeric matrices. A \eqn{T\times J} matrix of
//'   posterior change-point location probabilities and the log evaluations for
//'   the meanvar-scp components of model.
//' @param omega_j,log_omega_j Scalars. Prior precision parameter and log
//'   evaluation for meanvar-scp components of model.
//' @param u_j,log_u_j,lgamma_u_j Scalars. Prior shape parameter, log and log
//'   gamma evaluations for meanvar-scp components of model.
//' @param v_j,log_v_j Scalars. Prior rate parameter and log evaluation for
//'   meanvar-scp components of model.
//' @param log_pi_j A numeric matrix. A \eqn{T\times J} matrix of log
//'   prior change-point location probabilities for meanvar-scp components
//'   of model.
//' @param b_bar_l A numeric matrix. A \eqn{T \times L} matrix of posterior mean
//'   parameters for mean-scp components of model.
//' @param omega_bar_l A numeric matrix. A \eqn{T\times L} matrix of posterior
//'   variance parameters for mean-scp components of model.
//' @param pi_bar_l,log_pi_bar_l Numeric matrices. A \eqn{T\times L} matrix of
//'   posterior change-point location probabilities and the log evaluations for
//'   the mean-scp components of model.
//' @param omega_l,log_omega_l Scalars. Prior precision parameter and log
//'   evaluation for mean-scp components of model.
//' @param log_pi_l A numeric matrix. A \eqn{T\times L} matrix of log
//'   prior change-point location probabilities for mean-scp components of
//'   model.
//' @param u_bar_k,lgamma_u_bar_k,digamma_u_bar_k Numeric vectors. A length
//'   \eqn{T} vector of posterior shape parameters and log gamma and digamma
//'   evaluations for var-scp components of model.
//' @param v_bar_k A numeric matrix. A \eqn{T\times K} matrix of posterior
//'   rate parameters for var-scp components of model.
//' @param pi_bar_k,log_pi_bar_k Numeric matrices. A \eqn{T\times K} matrix of
//'   posterior change-point location probabilities and the log evaluations for
//'   the var-scp components of model.
//' @param u_k,log_u_k,lgamma_u_k Scalars. Prior shape parameter, log and log
//'   gamma evaluations for var-scp components of model.
//' @param v_k,log_v_k Scalars. Prior rate parameter and log evaluation for
//'   var-scp components of model.
//' @param log_pi_k A numeric matrix. A \eqn{T\times K} matrix of log
//'   prior change-point location probabilities for nvar-scp components of
//'   model.
//'
//' @return A scalar. Value of ELBO for MICH model at current parameter values.
//'
double elbo_fn(int T, double mu_0, double lambda_0,
               NumericVector r_tilde,
               NumericVector lambda_bar,
               NumericVector delta,
               NumericMatrix b_bar_j,
               NumericMatrix omega_bar_j,
               NumericVector u_bar_j,
               NumericMatrix v_bar_j,
               NumericMatrix pi_bar_j,
               NumericMatrix log_pi_bar_j,
               NumericVector lgamma_u_bar_j,
               NumericVector digamma_u_bar_j,
               double omega_j, double u_j, double v_j, double log_omega_j,
               double log_u_j, double lgamma_u_j, double log_v_j,
               NumericMatrix log_pi_j,
               NumericMatrix b_bar_l,
               NumericMatrix omega_bar_l,
               NumericMatrix pi_bar_l,
               NumericMatrix log_pi_bar_l,
               double omega_l, double log_omega_l,
               NumericMatrix log_pi_l,
               NumericVector u_bar_k,
               NumericMatrix v_bar_k,
               NumericMatrix pi_bar_k,
               NumericVector lgamma_u_bar_k,
               NumericVector digamma_u_bar_k,
               NumericMatrix log_pi_bar_k,
               double u_k, double v_k, double log_u_k,
               double lgamma_u_k, double log_v_k,
               NumericMatrix log_pi_k) {

  int J = pi_bar_j.ncol(), L = pi_bar_l.ncol(), K = pi_bar_k.ncol();

  // calculate E[log p(y)]
  double elbo = 0.5 * T * (std::log(lambda_0) - log_2_pi);

  for (int t = 0; t < T; t++) {
    elbo += -0.5 * lambda_bar[t] * (r_tilde[t] * r_tilde[t] + delta[t]);

    for (int j = 0; j < J; j++) {
      // E[log lambda_j] component of E[log p(y)]
      elbo += 0.5 * (T - t) * (digamma_u_bar_j[t] - std::log(v_bar_j(t,j))) * pi_bar_j(t,j);
      // E[log p_j] - E[log q_j]
      // adding log odds with numerical stability for small pi_bar_j values
      if (pi_bar_j(t,j) > 1e-20) elbo += pi_bar_j(t,j) * (log_pi_j(t,j) - log_pi_bar_j(t,j));
      // mean component
      elbo += 0.5 * pi_bar_j(t,j) * (log_omega_j - std::log(omega_bar_j(t,j)));
      elbo += 0.5 * pi_bar_j(t,j) * (1 - omega_j * (b_bar_j(t,j) * b_bar_j(t,j) * (u_bar_j[t] / v_bar_j(t,j)) + 1 / omega_bar_j(t,j)));
      // // precision component
      elbo += pi_bar_j(t,j) * (u_j * (log_v_j - std::log(v_bar_j(t,j))) + lgamma_u_bar_j[t] - lgamma_u_j);
      elbo += pi_bar_j(t,j) * ((u_j - u_bar_j[t]) * digamma_u_bar_j[t] + u_bar_j[t] * (1 - v_j / v_bar_j(t,j)));
    }
    for (int l = 0; l < L; l++) {
      // E[log p_l] - E[log q_l]
      // adding log odds with numerical stability for small pi_bar_l values
      if (pi_bar_l(t,l) > 1e-20) elbo += pi_bar_l(t,l) * (log_pi_l(t,l) - log_pi_bar_l(t,l));
      elbo += 0.5 * pi_bar_l(t,l) * (log_omega_l - std::log(omega_bar_l(t,l)));
      elbo += 0.5 * pi_bar_l(t,l) * (1 - omega_l * (b_bar_l(t,l) * b_bar_l(t,l) + 1 / omega_bar_l(t,l)));
    }
    for (int k = 0; k < K; k++) {
      // E[log lambda_k] component of E[log p(y)]
      elbo += 0.5 * (T - t) * (digamma_u_bar_k[t] - std::log(v_bar_k(t,k))) * pi_bar_k(t,k);
      // E[log p_j] - E[log q_j]
      // adding log odds with numerical stability for small pi_bar_l values
      if (pi_bar_k(t,k) > 1e-20) elbo += pi_bar_k(t,k) * (log_pi_k(t,k) - log_pi_bar_k(t,k));
      elbo += pi_bar_k(t,k) * (u_k * (log_v_k - std::log(v_bar_k(t,k))) + lgamma_u_bar_k[t] - lgamma_u_k);
      elbo += pi_bar_k(t,k) * ((u_k - u_bar_k[t]) * digamma_u_bar_k[t] + u_bar_k[t] * (1 - v_k / v_bar_k(t,k)));
    }
  }
  return elbo;
}
