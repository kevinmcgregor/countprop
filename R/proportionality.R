#' Logit Normal Variation
#'
#' Estimates the variation matrix of count-compositional data
#' based on a multinomial logit-Normal distribution. Estimate is performed using
#' only the parameters of the distribution.
#'
#' @param mu The mle estimate of the mu matrix
#' @param Sigma The mle estimate of the Sigma matrix
#' @param type Type of variation metric to be calculated: \code{standard}, \code{phi},
#' \code{phis} (a symmetrical version of \code{phi}), or \code{rho}
#' @param order The order of the Taylor-series approximation to be used in the
#' estimation
#'
#' @return An estimation of the variation matrix, \code{V}.
#'
#' @examples
#' data(singlecell)
#' mle <- mleLR(singlecell)
#' mu.hat <- mle$mu
#' Sigma.hat <- mle$Sigma
#'
#' logitNormalVariation(mu.hat, Sigma.hat)
#' logitNormalVariation(mu.hat, Sigma.hat, type="phi")
#' logitNormalVariation(mu.hat, Sigma.hat, type="rho")
#'
#' @export
logitNormalVariation <- function(mu, Sigma, type=c("standard","phi", "phis","rho"),
                                 order=c("second", "first")) {

  if (!is.vector(mu) | !is.numeric(mu)) stop("mu must be a numeric vector")
  if (!is.matrix(Sigma) | !is.numeric(Sigma)) stop("Sigma must be a numeric matrix")
  if (NROW(Sigma)!=NCOL(Sigma) | !isSymmetric(Sigma)) stop("Sigma must be a valid covariance matrix")
  if (length(mu)!=NROW(Sigma)) stop("Dimension mismatch between mu and Sigma")

  type <- match.arg(type)
  J <- length(mu)

  ones <- rep(1, J)
  ones.long <- rep(1, J+1)
  d.S <- diag(Sigma)
  V <- tcrossprod(d.S, ones) + tcrossprod(ones, d.S) - 2*Sigma
  V <- rbind(cbind(V, d.S), c(d.S, 0))
  rownames(V) <- colnames(V) <- NULL

  if (type=="phi") {
    lv <- logVarTaylorFull(mu, Sigma, order=order)
    lv.row <- tcrossprod(diag(lv), ones.long)
    V <- V/lv.row
  } else if (type=="phis") {
    lv <- logVarTaylorFull(mu, Sigma, order=order)
    lv.d <- diag(lv)
    den <- outer(lv.d, lv.d, "+") + 2*lv
    V <- V/den
  } else if (type=="rho") {
    lv <- logVarTaylorFull(mu, Sigma, order=order)
    lv.d <- diag(lv)
    den <- outer(lv.d, lv.d, "+")
    V <- 2*lv/den
  }

  return(V)
}

#' Plugin Variation
#'
#' Estimates the variation matrix of count-compositional data
#' based on a the same approximation used in logitNormalVariation()
#' only for this function it uses empirical estimates of mu and Sigma.
#' Also performs zero-imputation using \code{cmultRepl()} from the
#'
#' \code{zCompositions} package.
#'
#' @param counts Matrix of counts; samples are rows and features are columns.
#' @param type Type of variation metric to be calculated: \code{standard}, \code{phi},
#' \code{phis} (a symmetrical version of \code{phi}), or \code{rho}
#' @param order The order of the Taylor-series approximation to be used in the
#' estimation
#' @param impute.zeros If TRUE, then \code{cmultRepl()} from the \code{zCompositions} package is used to impute zero values in the counts matrix.
#' @param ... Optional arguments passed to zero-imputation function \code{cmultRepl()}
#'
#' @return An estimation of the variation matrix, \code{V}.
#'
#' @importFrom zCompositions cmultRepl
#' @examples
#' data(singlecell)
#'
#' pluginVariation(singlecell)
#' pluginVariation(singlecell, type="phi")
#' pluginVariation(singlecell, type="rho")
#'
#' @export
pluginVariation <- function(counts, type=c("standard","phi", "phis","rho"),
                                 order=c("second", "first"), impute.zeros=TRUE, ...) {
  if (!is.matrix(counts) | !is.numeric(counts)) stop("counts must be a numeric matrix")
  if (!is.logical(impute.zeros)) stop("impute.zeros must be TRUE or FALSE")

  type <- match.arg(type)
  J <- NCOL(counts) - 1

  # Imputing zeros
  if (any(counts==0) & impute.zeros) {
    y.no0 <- as.matrix(cmultRepl(counts, output = "p-count", ...))
  } else {
    y.no0 <- counts
  }
  y.alr <- log(y.no0[,-(J+1)]) - log(y.no0[,(J+1)])

  # Empirical estimates of mu and Sigma
  mu <- colMeans(y.alr)
  Sigma <- cov(y.alr)

  V <- logitNormalVariation(mu, Sigma, type, order)
  colnames(V) <- rownames(V) <- colnames(counts)
  return(V)
}

#' Naive (Empirical) Variation
#'
#' Naive (empirical) estimates of proportionality metrics using only the
#' observed counts.
#'
#' @param counts Matrix of counts; samples are rows and features are columns
#' @param pseudo.count Positive count to be added to all elements of count matrix.
#' @param type Type of variation metric to be calculated: \code{standard}, \code{phi},
#'  \code{phis} (a symmetric version of \code{phi}), \code{rho}, or \code{logp} (the variance-covariance matrix of log-transformed proportions)
#' @param impute.zeros If TRUE, then \code{cmultRepl()} from the \code{zCompositions} package is used to impute zero values in the counts matrix.
#' @param ... Optional arguments passed to zero-imputation function \code{cmultRepl()}
#'
#' @return A matrix containing the proportionality metric of interest calculated naively (empirically).
#'
#' @examples
#' #' data(singlecell)
#'
#' naiveVariation(singlecell)
#' naiveVariation(singlecell, type="phi")
#' naiveVariation(singlecell, type="rho")
#'
#' @export
#'
naiveVariation <- function(counts, pseudo.count=0, type=c("standard","phi", "phis","rho", "logp"),
                           impute.zeros=TRUE, ...) {

  if (!is.matrix(counts) | !is.numeric(counts)) stop("counts must be a numeric matrix")
  if (!is.logical(impute.zeros)) stop("impute.zeros must be TRUE or FALSE")
  if (!is.numeric(pseudo.count)) stop("pseudo.count must be numeric")
  if (pseudo.count<0) stop("pseudo.count must be non-negative")

  type <- match.arg(type)
  J <- NCOL(counts)
  l <- counts

  l <- l + pseudo.count
  l <- l/rowSums(l)
  l <- log(l)
  get.inf <- is.infinite(l)
  if (any(get.inf)) {
    stop("There are infinities after taking log.  Consider setting impute.zeros=TRUE")
  }


  v <- matrix(0,J,J)
  for (i in 1:J) {
    for (j in 1:J){
      if (type=="standard") {
        v[i,j] <- compositions::var(l[,i]-l[,j])
      } else if (type=="phi") {
        v[i,j] <- compositions::var(l[,i]-l[,j])/compositions::var(l[,i])
      } else if (type=="phis") {
        v[i,j] <- compositions::var(l[,i]-l[,j])/(compositions::var(l[,i]+l[,j]))
      } else if (type=="rho") {
        v[i,j] <- 2*compositions::cov(l[,i],l[,j])/(compositions::var(l[,i])+compositions::var(l[,j]))
      }
    }
  }

  if (type=="logp") v <- compositions::cov(l)

  colnames(v) <- rownames(v) <- colnames(counts)
  return(v)
}


#' Full logp Variance-Covariance
#'
#' Estimates the variance-covariance of the log of the proportions using a
#' Taylor-series approximation.
#'
#' @param mu The mean vector of the log-ratio-transformed data (ALR or CLR)
#' @param Sigma The variance-covariance matrix of the log-ratio-transformed data (ALR or CLR)
#' @param transf The desired transformation. If \code{transf="alr"} the inverse
#' additive log-ratio transformation is applied. If \code{transf="clr"} the
#' inverse centered log-ratio transformation is applied.
#' @param order The desired order of the Taylor Series approximation
#'
#' @return The estimated variance-covariance matrix for \code{log p}.
#'
#' @examples
#' data(singlecell)
#' mle <- mleLR(singlecell)
#' mu <- mle$mu
#' Sigma <- mle$Sigma
#'
#' logVarTaylorFull(mu, Sigma)
#'
#' @export
#'
logVarTaylorFull <- function(mu, Sigma, transf=c("alr", "clr"), order=c("second", "first")) {
  transf <- match.arg(transf)
  order <- match.arg(order)

  D <- length(mu)
  ones <- rep(1, D+1)
  emu <- exp(mu)
  if (transf=="alr") {
    ainv <- emu/(1+sum(emu))
  } else {
    ainv <- emu/sum(emu)
  }
  M <- rbind(diag(D),0)-tcrossprod(ones, ainv)
  t2 <- 0
  if (order=="second") {
    mat <- Sigma%*%(tcrossprod(ainv)-diag(ainv))
    t2 <- sum(diag(mat%*%mat))
  }
  M%*%tcrossprod(Sigma, M) + 0.5*t2
}
