#include <Rcpp.h>
#include <Rmath.h>
using namespace Rcpp;

//' Multivariate Expected Mean Signal
//'
//' @description
//' Given that \eqn{E[\mathbf{y}_t|\mathbf{b},\tau] = \mathbf{b}I(t\geq \tau)},
//' \eqn{P(\tau = t)=\pi_t}, and \eqn{\bar{\mathbf{b}}_t=E[\mathbf{b}|\tau=t]},
//' `multi_mu_bar_fn()` calculates \eqn{E[\mathbf{b}I(t\geq \tau)]} as
//' \eqn{\Sigma_{t'=1}^t \bar{\mathbf{b}}_{t'}\pi_{t'}}.
//'
//' @keywords internal
//'
//' @param b A numeric matrix. \eqn{T\times d} matrix of conditional mean
//'   parameters.
//' @param prob A numeric vector. Vector of change-point location probabilities.
//'
//' @return A numeric matrix. A \eqn{T\times d} matrix of
//' \eqn{E[\mathbf{b}I(t\geq \tau)]}.
//'
// [[Rcpp::export]]
NumericMatrix multi_mu_bar_fn(NumericMatrix b, NumericVector prob) {
  int T = prob.length();
  int d = b.ncol();
  NumericVector fwd_sum (d);
  NumericMatrix mu_bar (T, d);
  for (int t = 0; t < T; t++) {
    for (int i = 0; i < d; i++) {
      fwd_sum[i] += b(t, i) * prob[t];
      mu_bar(t, i) = fwd_sum[i];
    }
  }
  return mu_bar;
}

//' Multivariate Expected Squared-Mean Signal
//'
//' @description
//' Given that \eqn{E[\mathbf{y}_t|\mathbf{b},\tau] = \mathbf{b}I(t\geq \tau)},
//' \eqn{P(\tau = t)=\pi_t}, and for some \eqn{1\leq i \leq d},
//' \eqn{\bar{b}_{it}=E[b_i|\tau=t]} and \eqn{\bar{\omega}_{it}=V(b_i|\tau=t)},
//' `multi_mu2_bar_fn()` calculates \eqn{E[b_i^2I(t\geq \tau)]} as
//' \eqn{\Sigma_{t'=1}^t (\bar{b}^2_{it'} + 1/\bar{\omega}_{it'}) \pi_{t'}}.
//'
//' @keywords internal
//'
//' @param b A numeric matrix. \eqn{T\times d} matrix of conditional mean
//'   parameters.
//' @param omega A numeric vector. Length \eqn{T} vector of conditional variance
//'   parameters.
//' @param prob A numeric vector. Vector of change-point location probabilities.
//'
//' @return A numeric matrix. A \eqn{T\times d} matrix of
//' \eqn{E[b^2_iI(t\geq \tau)]}.
//'
// [[Rcpp::export]]
NumericMatrix multi_mu2_bar_fn(NumericMatrix b, NumericVector omega, NumericVector prob) {
  int T = prob.length();
  int d = b.ncol();
  NumericVector fwd_sum (d);
  NumericMatrix mu2_bar (T, d);
  for (int t = 0; t < T; t++) {
    for (int i = 0; i < d; i++) {
      fwd_sum[i] += (b(t,i) * b(t,i) + 1 / omega[t]) * prob[t];
      mu2_bar(t, i) = fwd_sum[i];
    }
  }
  return mu2_bar;
}

//' Expected Mean Signal
//'
//' @description
//' Given that \eqn{E[y_t|b,\tau]=\mu_t=bI(t\geq \tau)}, \eqn{P(\tau=t)=\pi_t},
//' and \eqn{\bar{b}_t=E[b|\tau=t]}, `mu_bar_fn()` calculates \eqn{E[\mu_t]} as
//' \eqn{\Sigma_{t'=1}^t \bar{b}_{t'}\pi_{t'}}.
//'
//' @keywords internal
//'
//' @param b A numeric vector. Length \eqn{T} vector of conditional mean
//'   parameters.
//' @param prob A numeric vector. Vector of change-point location probabilities.
//'
//' @return A numeric vector. A length \eqn{T} vector of \eqn{E[\mu_t]}.
//'
// [[Rcpp::export]]
NumericVector mu_bar_fn(NumericVector b, NumericVector prob) {
  int T = prob.length();
  double fwd_sum = 0.0;
  NumericVector mu_bar (T, 0.0);
  for (int t = 0; t < T; t++) {
    fwd_sum += b[t] * prob[t];
    mu_bar[t] = fwd_sum;
  }
  return mu_bar;
}

//' Expected Squared-Mean Signal
//'
//' @description
//' Given that \eqn{E[y_t|b,\tau]=\mu_t=bI(t\geq \tau)}, \eqn{P(\tau=t)=\pi_t},
//' \eqn{\bar{b}_t=E[b|\tau=t]}, and \eqn{\bar{\omega}_t=V(b|\tau=t)},
//' `mu2_bar_fn()` calculates \eqn{E[\mu^2_t]} as
//' \eqn{\Sigma_{t'=1}^t (\bar{b}^2_{t'} + 1 / \bar{\omega}_{t'})\pi_{t'}}.
//'
//' @keywords internal
//'
//' @param b A numeric vector. Length \eqn{T} vector of conditional mean
//'   parameters.
//' @param omega A numeric vector. Length \eqn{T} vector of conditional variance
//'   parameters.
//' @param prob A numeric vector. Vector of change-point location probabilities.
//'
//' @return A numeric vector. A length \eqn{T} vector of \eqn{E[\mu_t]}.
//'
// [[Rcpp::export]]
NumericVector mu2_bar_fn(NumericVector b, NumericVector omega, NumericVector prob) {
  int T = prob.length();
  double fwd_sum = 0.0;
  NumericVector mu2_bar (T, 0.0);
  for (int t = 0; t < T; t++) {
    fwd_sum += (b[t] * b[t] + 1 / omega[t]) * prob[t];
    mu2_bar[t] = fwd_sum;
  }
  return mu2_bar;
}

//' Expected Precision Signal
//'
//' @description
//' Given that \eqn{V(y_t|s,\tau)=1/\lambda_t=1/s^{I(t\geq \tau)}},
//' \eqn{P(\tau=t)=\pi_t}, and \eqn{\bar{u}_t/\bar{v}_t=E[s|\tau=t]},
//' `lambda_bar_fn()` calculates \eqn{E[\lambda_t]} as
//' \eqn{1 - \Sigma_{t'=1}^t \pi_{t'}(1 - \bar{v}_{t'} / \bar{u}_{t'})}.
//'
//' @keywords internal
//'
//' @param u A numeric vector. Length \eqn{T} vector of conditional shape
//'   parameters
//' @param v A numeric vector. Length \eqn{T} vector of conditional rate
//'   parameters
//' @param prob A numeric vector. Vector of change-point location probabilities.
//'
//' @return A numeric vector. A length \eqn{T} vector of \eqn{E[\lambda_t]}.
//'
// [[Rcpp::export]]
NumericVector lambda_bar_fn(NumericVector u, NumericVector v, NumericVector prob) {
  int T = prob.length();
  double fwd_sum = 0.0, rev_sum = 1.0;
  NumericVector lambda_bar (T, 0.0);
  for (int t = 0; t < T; t++) {
    rev_sum -= prob[t];
    fwd_sum += (u[t] / v[t]) * prob[t];
    lambda_bar[t] = fwd_sum + rev_sum;
  }
  return lambda_bar;
}

//' Expected Mean-Precision Signal
//'
//' @description
//' Given that \eqn{E[y_t|b,\tau]=\mu_t=bI(t\geq \tau)},
//' \eqn{\bar{b}_t=E[b|\tau=t]}, \eqn{\bar{\omega}_t=V(b|\tau=t)},
//' \eqn{V(y_t|s,\tau)=1/\lambda_t=1/s^{I(t\geq \tau)}},
//' \eqn{\bar{u}_t/\bar{v}_t=E[s|\tau=t]}, and \eqn{P(\tau=t)=\pi_t},
//' `mu_lambda_fn()` calculates \eqn{E[\mu_t\lambda_t]} as
//' \eqn{\Sigma_{t'=1}^t \bar{b}_{t'}\bar{v}_{t'} \pi_{t'}/ \bar{u}_{t'}}.
//'
//' @keywords internal
//'
//' @param b A numeric vector. Length \eqn{T} vector of conditional mean
//'   parameters.
//' @param u A numeric vector. Length \eqn{T} vector of conditional shape
//'   parameters
//' @param v A numeric vector. Length \eqn{T} vector of conditional rate
//'   parameters
//' @param prob A numeric vector. Vector of change-point location probabilities.
//'
//' @return A numeric vector. A length \eqn{T} vector of
//' \eqn{E[\mu_t\lambda_t]}.
//'
// [[Rcpp::export]]
NumericVector mu_lambda_fn(NumericVector b, NumericVector u, NumericVector v,
                           NumericVector prob) {
  int T = prob.length();
  double fwd_sum = 0.0;
  NumericVector mu_lambda (T, 0.0);
  for (int t = 0; t < T; t++) {
    fwd_sum += b[t] * (u[t] / v[t]) * prob[t];
    mu_lambda[t] = fwd_sum;
  }
  return mu_lambda;
}

//' Expected Squared-Mean-Precision Signal
//'
//' @description
//' Given that \eqn{E[y_t|b,\tau]=\mu_t=bI(t\geq \tau)},
//' \eqn{\bar{b}_t=E[b|\tau=t]},
//' \eqn{V(y_t|s,\tau)=1/\lambda_t=1/s^{I(t\geq \tau)}},
//' \eqn{\bar{u}_t/\bar{v}_t=E[s|\tau=t]}, and \eqn{P(\tau=t)=\pi_t},
//' `mu_lambda_fn()` calculates \eqn{E[\mu^2_t\lambda_t]} as
//' \eqn{\Sigma_{t'=1}^t \pi_{t'}(\bar{b}^2_{t'}\bar{v}_{t'} / \bar{u}_{t'} + 1/\bar{\omega}_{t'})}.
//'
//' @keywords internal
//'
//' @param b A numeric vector. Length \eqn{T} vector of conditional mean
//'   parameters.
//' @param u A numeric vector. Length \eqn{T} vector of conditional shape
//'   parameters
//' @param omega A numeric vector. Length \eqn{T} vector of conditional variance
//'   parameters.
//' @param v A numeric vector. Length \eqn{T} vector of conditional rate
//'   parameters
//' @param prob A numeric vector. Vector of change-point location probabilities.
//'
//' @return A numeric vector. A length \eqn{T} vector of
//' \eqn{E[\mu^2_t\lambda_t]}.
//'
// [[Rcpp::export]]
NumericVector mu2_lambda_fn(NumericVector b, NumericVector omega, NumericVector u,
                            NumericVector v, NumericVector prob) {
  int T = prob.length();
  double fwd_sum = 0.0;
  NumericVector mu2_lambda (T, 0.0);
  for (int t = 0; t < T; t++) {
    fwd_sum += (b[t] * b[t] * (u[t] / v[t]) + 1 / omega[t]) * prob[t];
    mu2_lambda[t] = fwd_sum;
  }
  return mu2_lambda;
}
