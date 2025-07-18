% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/priors.R
\name{log_var_prior}
\alias{log_var_prior}
\title{Log Weighted Var-SCP Prior}
\usage{
log_var_prior(T)
}
\arguments{
\item{T}{An integer. Number of observations in \eqn{y_{1:T}}.}
}
\value{
A numeric vector. Log prior probabilities.
}
\description{
Log weighted prior for the Var-SCP model as described in Appendix C.2 of
Berlind, Cappello, and Madrid Padilla (2025). This prior ensures that in the
absence of a change-point, the posterior probabilities are approximately
uniform.
}
