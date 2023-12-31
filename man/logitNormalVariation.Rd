% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/proportionality.R
\name{logitNormalVariation}
\alias{logitNormalVariation}
\title{Logit Normal Variation}
\usage{
logitNormalVariation(
  mu,
  Sigma,
  type = c("standard", "phi", "phis", "rho"),
  order = c("second", "first")
)
}
\arguments{
\item{mu}{The mle estimate of the mu matrix}

\item{Sigma}{The mle estimate of the Sigma matrix}

\item{type}{Type of variation metric to be calculated: \code{standard}, \code{phi},
\code{phis} (a symmetrical version of \code{phi}), \code{rho}, or \code{logp} (the variance-covariance matrix of log-transformed proportions)}

\item{order}{The order of the Taylor-series approximation to be used in the
estimation}
}
\value{
An estimate of the requested metric of proportionality.
}
\description{
Estimates the variation matrix of count-compositional data
based on a multinomial logit-Normal distribution. Estimation is performed using
only the parameters of the distribution.
}
\examples{
data(singlecell)
mle <- mleLR(singlecell)
mu.hat <- mle$mu
Sigma.hat <- mle$Sigma

logitNormalVariation(mu.hat, Sigma.hat)
logitNormalVariation(mu.hat, Sigma.hat, type="phi")
logitNormalVariation(mu.hat, Sigma.hat, type="rho")

}
