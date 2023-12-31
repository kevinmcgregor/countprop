% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mle.R
\name{ebic}
\alias{ebic}
\title{Extended Bayesian Information Criterion}
\usage{
ebic(l, n, d, df, gamma)
}
\arguments{
\item{l}{Log-likelihood estimates of the model}

\item{n}{Number of rows of the data set for which the log-likelihood has been
calculated}

\item{d}{The size of the (k-1) by (k-1) covariance matrix of a
k by k count-compositional data matrix}

\item{df}{Degrees of freedom}

\item{gamma}{A tuning parameter. Larger values means more penalization}
}
\value{
The value of the EBIC.
}
\description{
Calculates the Extended Bayesian Information Criterion (EBIC) of a model.
Used for model selection to asses the fit of the multinomial logit-Normal
model which includes a graphical lasso penalty.
}
\note{
The graphical lasso penalty
is the sum of the absolute value of the elements of the covariance matrix \code{Sigma}.
The penalization parameter lambda controls the sparsity of Sigma.
}
\examples{
data(singlecell)
mle <- mleLR(singlecell, lambda.gl=0.5)
log.lik_1 <- mle$est[[1]]$log.lik
n <- NROW(singlecell)
k <- NCOL(singlecell)
df_1 <- mle$est[[1]]$df

ebic(log.lik_1, n, k, df_1, 0.1)

}
