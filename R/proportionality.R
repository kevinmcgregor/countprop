#' Logit Normal Variation
#'
#' Estimates the variation matrix of count-compositional data
#' based on a multinomial logit-Normal distribution. Estimation is performed using
#' only the parameters of the distribution.
#'
#' @param mu The mle estimate of the mu matrix
#' @param Sigma The mle estimate of the Sigma matrix (input Sigma on the ALR scale, even if requesting metrics on the CLR scale)
#' @param type Type of variation metric to be calculated: \code{standard}, \code{phi},
#' \code{phis} (a symmetrical version of \code{phi}), or \code{rho}
#' @param lr Which scale to calculate the proportionality metric on, either alr or clr.
#' @param order Deprecated: The order of the Taylor-series approximation to be used in the
#' estimation
#'
#' @return An estimate of the requested metric of proportionality.
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
                                 lr=c("alr", "clr"), order=c("second", "first")) {

  if (!is.vector(mu) | !is.numeric(mu)) stop("mu must be a numeric vector")
  if (!is.matrix(Sigma) | !is.numeric(Sigma)) stop("Sigma must be a numeric matrix")
  if (NROW(Sigma)!=NCOL(Sigma) | !isSymmetric(Sigma)) stop("Sigma must be a valid covariance matrix")
  if (length(mu)!=NROW(Sigma)) stop("Dimension mismatch between mu and Sigma")

  if (!missing(order)) stop("order argument is deprecated")

  type <- match.arg(type)
  lr <- match.arg(lr)

  J <- length(mu)
  ones <- rep(1, J)
  d.S <- diag(Sigma)
  V <- tcrossprod(d.S, ones) + tcrossprod(ones, d.S) - 2*Sigma
  V <- cbind(rbind(V, d.S), c(d.S, 0))
  dimnames(V) <- NULL

  if (lr=="alr") {
    if (type=="phi") {
      den <- tcrossprod(d.S, ones)
      V[1:J,1:J] <- V[1:J,1:J]/den
      V[J+1,1:J] <- rep(Inf, J)
      V[,J+1] <- c(rep(1, J), NaN)
    } else if (type=="phis" | type=="rho") {
      # Calculate rho.  If type is "phis" then rho is still calculated
      # but gets transformed to phis later.
      V <- cbind(rbind(2*Sigma/outer(d.S, d.S, "+"), 0), 0)
      V[J+1,J+1] <- 1
    }
  } else {
    H.inv <- qr.solve(diag(J) + matrix(1, J, J))
    Fm <- cbind(diag(J), -1)
    HiF <- H.inv%*%Fm
    Sigma.clr <- crossprod(HiF, Sigma)%*%HiF
    d.Sc <- diag(Sigma.clr)
    if (type=="phi") {
      V <- V/tcrossprod(d.Sc, c(ones, 1))
    } else if (type=="phis" | type=="rho") {
      # Calculate rho.  If type is "phis" then rho is still calculated
      # but gets transformed to phis later.
      V <- 2*Sigma.clr/outer(d.Sc, d.Sc, "+")
    }
  }

  if (type=="phis") {
    V <- (1-V)/(1+V)
    if (lr=="alr") {
      V[J+1,1:J] <- rep(1, J)
      V[,J+1] <- c(rep(1, J), 0)
    }
  }

  return(V)
}

#' Plugin Variation (Deprecated)
#'
#' Estimates the variation matrix of count-compositional data
#' based on a the same approximation used in logitNormalVariation()
#' only for this function it uses empirical estimates of mu and Sigma.
#' Also performs zero-imputation using \code{cmultRepl()} from the \code{zCompositions} package.
#'
#' @param counts Matrix of counts; samples are rows and features are columns.
#' @param type Type of variation metric to be calculated: \code{standard}, \code{phi},
#' \code{phis} (a symmetrical version of \code{phi}), \code{rho}, or \code{logp} (the variance-covariance matrix of log-transformed proportions).
#' @param order The order of the Taylor-series approximation to be used in the
#' estimation
#' @param impute.zeros If TRUE, then \code{cmultRepl()} from the \code{zCompositions} package is used to impute zero values in the counts matrix.
#' @param ... Optional arguments passed to zero-imputation function \code{cmultRepl()}
#'
#' @return An estimate of the requested metric of proportionality.
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
  .Deprecated("naiveVariation()")
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
#' @param lr Which scale to calculate the proportionality metric on, either alr or clr.
#' @param impute.zeros If TRUE, then \code{cmultRepl()} from the \code{zCompositions} package is used to impute zero values in the counts matrix.
#' @param ... Optional arguments passed to zero-imputation function \code{cmultRepl()}
#'
#' @return An estimate of the requested metric of proportionality.
#'
#' @examples
#' #' data(singlecell)
#'
#' naiveVariation(singlecell)
#' naiveVariation(singlecell, type="phi")
#' naiveVariation(singlecell, type="rho")
#' naiveVariation(singlecell, type="rho", lr="clr")
#'
#' @importFrom zCompositions cmultRepl
#' @importFrom stats var
#' @export
#'
naiveVariation <- function(counts, pseudo.count=0, type=c("standard","phi", "phis","rho"),
                           lr=c("alr", "clr"), impute.zeros=TRUE, ...) {

  if (!is.matrix(counts) | !is.numeric(counts)) stop("counts must be a numeric matrix")
  if (!is.logical(impute.zeros)) stop("impute.zeros must be TRUE or FALSE")
  if (!is.numeric(pseudo.count)) stop("pseudo.count must be numeric")
  if (pseudo.count<0) stop("pseudo.count must be non-negative")

  type <- match.arg(type)
  lr <- match.arg(lr)

  if (impute.zeros & any(counts==0)) {
    l <- as.matrix(zCompositions::cmultRepl(counts, output = "p-counts", z.warning = 0.9999,
                                            suppress.print=TRUE))
  } else {
    l <- counts + pseudo.count
  }

  if (lr=="alr") {
    l <- compositions::alr(l)
  } else {
    l <- compositions::clr(l)
  }

  J <- NCOL(l)

  v <- matrix(0,J,J)
  for (i in 1:J) {
    for (j in 1:J){
      if (type=="standard") {
        v[i,j] <- var(l[,i]-l[,j])
      } else if (type=="phi") {
        v[i,j] <- var(l[,i]-l[,j])/var(l[,i])
      } else if (type=="phis") {
        v[i,j] <- var(l[,i]-l[,j])/(var(l[,i]+l[,j]))
      } else if (type=="rho") {
        v[i,j] <- 2*cov(l[,i],l[,j])/(var(l[,i])+var(l[,j]))
      }
    }
  }

  return(v)
}


#' Convert between CLR and ALR covariance matrices
#'
#' @param S Covariance matrix to be converted
#' @param direction Which direction to convert between alr and clr.
#'
#' @returns A covariance matrix on the requested scale.
#' @export
#'
#' @examples convertSigma(diag(3), "alr2clr")
convertSigma <- function(S, direction=c("alr2clr", "clr2alr")) {
  direction <- match.arg(direction)
  k <- ifelse(direction=="alr2clr", NCOL(S), NCOL(S)-1)

  Fm <- cbind(diag(k), -1)

  if (direction=="alr2clr") {
    Hinv <- qr.solve(diag(k)+matrix(1,k,k))
    HiF <- Hinv%*%Fm
    crossprod(HiF, S)%*%HiF
  } else{
    Fm%*%tcrossprod(S, Fm)
  }
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
