#include <Rcpp.h>
#include <Rmath.h>
using namespace Rcpp;

//' Mean Single Change-Point Model
//'
//' @description
//' Implementation of the univariate Mean-SCP model from Berlind, Cappello, and
//' Madrid Padilla (2025). The function `mean_scp()` takes a length \eqn{T}
//' vector \eqn{y_{1:T}} with a single mean change and returns the posterior
//' distribution of the change-point.
//'
//' @param y A numeric vector. \eqn{T} observations with a single variance
//'   change.
//' @param lambda A numeric vector. Precision vector of `y`.
//' @param omega A positive scalar. Prior precision parameter.
//' @param log_pi A numeric vector. Vector of log prior probabilities for the
//'   location of the change-point.
//'
//' @return A list. A list of posterior parameters including the mean `b_bar`,
//' precision `omega_bar`, and posterior probabilities of the change-point
//' location `pi_bar`.
//'
//' @export
//'
// [[Rcpp::export]]
List mean_scp(
  NumericVector y,
  NumericVector lambda,
  double omega,
  double log_omega,
  NumericVector log_pi
){
  int T = y.length();
  NumericVector b_bar (T+1, 0.0), pi_bar (T+1, 0.0), log_pi_bar (T+1, 0.0), omega_bar (T+1, 0.0);
  double tot = 0, ls = 0, ys = 0, log_pi_max = log_pi[T] - 0.5 * log_omega;

  // no change params
  log_pi_bar[T] = log_pi_max;
  omega_bar[T] = omega;

  for (int t = T - 1; t >= 0; t--) {
    ys += y[t] * lambda[t];
    ls += lambda[t];
    omega_bar[t] = omega + ls;
    b_bar[t] = ys / omega_bar[t];
    log_pi_bar[t] = log_pi[t] - 0.5 * (std::log(omega_bar[t]) - omega_bar[t] * b_bar[t] * b_bar[t]);
    log_pi_max = std::max(log_pi_max, log_pi_bar[t]);
  }

  // normalize posterior probability
  for (int t = 0; t < T+1; t++) {
    log_pi_bar[t] -= log_pi_max;
    if (log_pi[t] > R_NegInf) pi_bar[t] = std::exp(log_pi_bar[t]);
    else pi_bar[t] = 0.0;
    tot += pi_bar[t];
  }

  return List::create(
    _["b_bar"] = b_bar,
    _["omega_bar"] = omega_bar,
    _["pi_bar"] = pi_bar / tot,
    _["log_pi_bar"] = log_pi_bar - std::log(tot)
  );
}

//' Mulitvariate Mean Single Change-Point Model
//'
//' @description
//' Implementation of the multivariate Mean-SCP model from Berlind, Cappello,
//' and Madrid Padilla (2025). The function `mean_scp()` takes a \eqn{T\times d}
//' matrix \eqn{\mathbf{y}_{1:T}} with a single mean change and returns the
//' posterior distribution of the change-point.
//'
//' @param y A numeric matrix. \eqn{T\times d} matrix of observations with a
//'   single mean change.
//' @param omega_bar A numeric vector. Posterior precision parameters
//'   \eqn{\bar{\omega}_t} such that
//'   \eqn{\bar{\Omega}_t = \bar{\omega}_t\mathbf{I}_d}.
//' @param log_omega_bar A numeric vector. Log of `omega_bar`.
//' @param log_pi A numeric vector. Vector of log prior probabilities for the
//'   location of the change-point.
//'
//' @return A list. A list of posterior parameters including the mean `b_bar`,
//' and posterior probabilities of the change-point location `pi_bar`.
//'
//' @export
//'
// [[Rcpp::export]]
List multi_mean_scp(
  NumericMatrix QTr,
  NumericVector lambda,
  NumericMatrix mean_weights,
  NumericMatrix sandwich_weights,
  NumericVector log_prob_weights
){
  int T = QTr.nrow();
  int d = QTr.ncol();
  NumericVector QTr_sum (d, 0.0);
  NumericMatrix QTmu_bar (T, d);
  NumericVector muTLmu_bar (T);
  NumericVector log_pi_bar (T+1, 0.0);
  NumericMatrix QTb_bar (T+1, d);

  log_pi_bar = clone(log_prob_weights);
  double log_pi_max = log_prob_weights[T];

  for (int t = T - 1; t >= 0; t--) {
    for (int i = 0; i < d; i++) {
      QTr_sum[i] += QTr(t, i);
      log_pi_bar[t] += 0.5 * sandwich_weights(t, i) * QTr_sum[i] * QTr_sum[i];
      QTb_bar(t,i) = mean_weights(t, i) * QTr_sum[i];
    }
    log_pi_max = std::max(log_pi_max, log_pi_bar[t]);
  }

  double tot = 0.0;
  NumericVector pi_bar (T+1, 0.0);
  for (int t = 0; t < T+1; t++) {
    log_pi_bar[t] -= log_pi_max;
    pi_bar[t] = std::exp(log_pi_bar[t]);
    tot += pi_bar[t];
  }

  double log_tot = std::log(tot);

  NumericVector QTmu_sum (d, 0.0);
  double muTLmu_sum = 0.0;

  pi_bar[T] /= tot;
  log_pi_bar[T] -= log_tot;
  for (int t = 0; t < T; t++) {
    pi_bar[t] /= tot;
    log_pi_bar[t] -= log_tot;
    for (int i = 0; i < d; i++) {
      QTmu_sum[i] += QTb_bar(t, i) * pi_bar[t];
      QTmu_bar(t, i) = QTmu_sum[i];

      muTLmu_sum += (mean_weights(t,i) + lambda[i] * QTb_bar(t,i) * QTb_bar(t,i)) * pi_bar[t];
      muTLmu_bar[t] = muTLmu_sum;
    }
  }

  return List::create(
    _["QTb_bar"] = QTb_bar,
    _["QTmu_bar"] = QTmu_bar,
    _["muTLmu_bar"] = muTLmu_bar,
    _["pi_bar"] = pi_bar,
    _["log_pi_bar"] = log_pi_bar
  );
}

//' Variance Single Change-Point Model
//'
//' @description
//' Implementation of the Var-SCP model from Berlind, Cappello, and Madrid
//' Padilla (2025). The function `var_scp()` takes a length \eqn{T} vector
//' \eqn{y_{1:T}} with a single variance change and returns the posterior
//' distribution of the change-point.
//'
//' @param y A numeric vector. \eqn{T} observations with a single variance
//'   change.
//' @param omega A numeric vector. Known trend component of precision of `y`.
//' @param u_bar A numeric vector. Posterior shape parameters equal to
//'   \eqn{u_0 + T - t + 1} for each \eqn{t}.
//' @param lgamma_u_bar A numeric vector. Log gamma function evaluated at u_bar.
//' @param v A numeric vector. Vector of prior rate parameters for each \eqn{t}.
//' @param log_pi A numeric vector. Vector of log prior probabilities for the
//'   location of the change-point.
//'
//' @return A list. A list of posterior parameters including the rate `v_bar`,
//' and posterior probabilities of the change-point location `pi_bar`.
//'
//' @export
//'
// [[Rcpp::export]]
List var_scp(
  NumericVector y,
  NumericVector omega,
  NumericVector u_bar,
  NumericVector lgamma_u_bar,
  NumericVector v,
  double log_v,
  NumericVector log_pi
){
  int T = y.length();
  NumericVector pi_bar (T+1, 0.0), log_pi_bar (T+1, 0.0), v_bar (T+1, 0.0);
  double tot = 0, fs = 0, rs = 0, log_pi_max = R_NegInf;

  for (int t = T - 1; t >= 0; t--) {
    rs += y[t] * y[t] * omega[t];
    v_bar[t] = v[t] + 0.5 * rs;
    log_pi_bar[t] += log_pi[t] + lgamma_u_bar[t] - u_bar[t] * std::log(v_bar[t]);
    log_pi_bar[T-t-1] -= 0.5 * fs;
    if (t <= T-t-1) {
      log_pi_max = std::max(log_pi_max, log_pi_bar[T-t-1]);
      log_pi_max = std::max(log_pi_max, log_pi_bar[t]);
    }
    fs += y[T-t-1] * y[T-t-1] * omega[T-t-1];
  }

  // no change params
  log_pi_bar[T] = log_pi[T] + lgamma_u_bar[T] - u_bar[T] * log_v - 0.5 * fs;
  log_pi_max = std::max(log_pi_max, log_pi_bar[T]);
  v_bar[T] = v[T];

  // normalize posterior probability
  for (int t = 0; t < T+1; t++) {
    log_pi_bar[t] -= log_pi_max;
    if (log_pi[t] > R_NegInf) pi_bar[t] = std::exp(log_pi_bar[t]);
    else pi_bar[t] = 0.0;
    tot += pi_bar[t];
  }

  return List::create(
    _["v_bar"] = v_bar,
    _["pi_bar"] = pi_bar / tot,
    _["log_pi_bar"] = log_pi_bar - std::log(tot)
  );
}

//' Mean-Variance Single Change-Point Model
//'
//' @description
//' Implementation of the MeanVar-SCP model from Berlind, Cappello, and Madrid
//' Padilla (2025). The function `meanvar_scp()` takes a length \eqn{T} vector
//' \eqn{y_{1:T}} with a single joint mean and variance change and returns the
//' posterior distribution of the change-point.
//'
//' @param y A numeric vector. \eqn{T} observations with a single joint mean and
//' variance change.
//' @param lambda A numeric vector. Known trend component of precision of `y`.
//' @param omega A positive scalar. Prior precision parameter.
//' @param u_bar A numeric vector. Posterior shape parameters equal to
//'   \eqn{u_0 + T - t + 1} for each \eqn{t}.
//' @param lgamma_u_bar A numeric vector. Log gamma function evaluated at u_bar.
//' @param v A numeric vector. Vector of prior rate parameters for each \eqn{t}.
//' @param log_pi A numeric vector. Vector of log prior probabilities for the
//'   location of the change-point.
//'
//' @return A list. A list of posterior parameters including the  mean `b_bar`,
//' precision `omega_bar`, rate `v_bar`, and posterior probabilities of the
//' change-point location `pi_bar`.
//'
//' @export
//'
// [[Rcpp::export]]
List meanvar_scp(
  NumericVector y,
  NumericVector lambda,
  double omega,
  double log_omega,
  NumericVector u_bar,
  NumericVector lgamma_u_bar,
  NumericVector v,
  double log_v,
  NumericVector log_pi
){

  int T = y.length();
  NumericVector b_bar (T+1, 0.0), pi_bar (T+1, 0.0), log_pi_bar (T+1, 0.0), v_bar (T+1, 0.0), omega_bar (T+1, 0.0);
  double tot = 0, fs = 0, rs = 0, ls = 0, ys = 0, log_pi_max = R_NegInf;

  for (int t = T - 1; t >= 0; t--) {
    ys += y[t] * lambda[t];
    ls += lambda[t];
    rs += y[t] * y[t] * lambda[t];
    omega_bar[t] = omega + ls;
    b_bar[t] = ys / omega_bar[t];
    v_bar[t] = v[t] + 0.5 * (rs - b_bar[t] * b_bar[t] * omega_bar[t]);
    log_pi_bar[t] += log_pi[t] + lgamma_u_bar[t] - u_bar[t] * std::log(v_bar[t]) - 0.5 * std::log(omega_bar[t]);
    log_pi_bar[T-t-1] -= 0.5 * fs;
    if (t <= T-t-1) {
      log_pi_max = std::max(log_pi_max, log_pi_bar[T-t-1]);
      log_pi_max = std::max(log_pi_max, log_pi_bar[t]);
    }
    fs += y[T-t-1] * y[T-t-1] * lambda[T-t-1];
  }

  // no change params
  log_pi_bar[T] = log_pi[T] - 0.5 * log_omega + lgamma_u_bar[T] - u_bar[T] * log_v - 0.5 * fs;
  log_pi_max = std::max(log_pi_max, log_pi_bar[T]);
  v_bar[T] = v[T];

  for (int t = 0; t < T+1; t++) {
    log_pi_bar[t] -= log_pi_max;
    if (log_pi[t] > R_NegInf) pi_bar[t] = std::exp(log_pi_bar[t]);
    else pi_bar[t] = 0.0;
    tot += pi_bar[t];
  }

  return List::create(
    _["b_bar"] = b_bar,
    _["omega_bar"] = omega_bar,
    _["v_bar"] = v_bar,
    _["pi_bar"] = pi_bar / tot,
    _["log_pi_bar"] = log_pi_bar - std::log(tot)
  );
}
