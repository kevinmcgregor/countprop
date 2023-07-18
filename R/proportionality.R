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
#' mu.hat <- mle$est.min$mu
#' Sigma.hat <- mle$est.min$Sigma
#'
#' logitNormalVariation(mu.hat, Sigma.hat) #Standard logit-Normal estimates of variance
#' logitNormalVariation(mu.hat, Sigma.hat, type="phi", order="second") #Logit-Normal based estimate of phi
#' logitNormalVariation(mu.hat, Sigma.hat, type="phis", order="second") #Logit-Normal based estimate of phis
#' logitNormalVariation(mu.hat, Sigma.hat, type="rho", order="second") #Logit-Normal based estimate of rho
#'
#' @export
logitNormalVariation <- function(mu, Sigma, type=c("standard","phi", "phis","rho"),
                                 order=c("second", "first")) {
  type <- match.arg(type)
  J <- length(mu)

  ones <- rep(1, J)
  ones.long <- rep(1, J+1)
  d.S <- diag(Sigma)
  V <- tcrossprod(d.S, ones) + tcrossprod(ones, d.S) - 2*Sigma
  V <- rbind(cbind(V, d.S), c(d.S, 0))
  colnames(V) <- NULL

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

#' Naive Variation
#'
#' Estimates the variation matrix of count compositional data using the proportions
#' of a count compositional dataset.
#'
#' @param counts Count-compositional dataset
#' @param pseudo.count Scaler value added to the data matrix to prevent infinite
#' values caused by taking the log of the counts
#' @param type Type of variation metric to be calculated: \code{standard}, \code{phi},
#'  \code{phis} (a symmetrical version of \code{phi}), \code{rho}, or \code{logx}
#' @param use If equal to \code{"everything"} and there are no infinite values after
#' taking the log the calculation will use all data. If equal to \code{"everything"}
#' and there are infinite values after taking the log it is recommended to  run
#' the function again, instead setting the \code{use} parameter to \code{"pairwise.complete.obs"}
#' @param set.inf.na If \code{TRUE}, sets any infinite values in \code{counts} to \code{NA}
#' @param already.log If \code{FALSE}, the counts have not been transformed by
#' by the log. This transformation is of the form \eqn{log(frac{X_{ij}}{s_{i}})}, where
#' \eqn{s_{i} = \sum{n=1}^{j} X_{in}}, where \eqn{X_{ij}} is element \eqn{counts[i,j]}
#'
#' @return The naive variation matrix, \code{v}.
#'
#' @examples
#' n.g <- ncol(dat.ss)
#'
#' naiveVariation(dat.ss)[-n.g,-n.g] #Standard naive estimate of the variance
#' naiveVariation(dat.ss, type="phi")[-n.g,-n.g] #Logit-Normal based naive estimate of phi
#' naiveVariation(dat.ss, type="phis")[-n.g,-n.g] #Logit-Normal based naive estimate of phis
#' naiveVariation(dat.ss, type="rho")[-n.g,-n.g] #Logit-Normal based naive estimate of rho
#'
#' @export
#'
naiveVariation <- function(counts, pseudo.count=0, type=c("standard","phi", "phis","rho", "logx"),
                           use="everything",
                           set.inf.na=TRUE, already.log=FALSE) {
  type <- match.arg(type)
  J <- NCOL(counts)
  l <- counts
  if (!already.log) {
    l <- l + pseudo.count
    l <- l/rowSums(l)
    l <- log(l)
    get.inf <- is.infinite(l)
    l[get.inf] <- NA
    if (any(get.inf) & use=="everything") {
      warning("There are infinities after taking log.  Consider setting paramter use='pairwise.complete.obs'")
    }
  }

  v <- matrix(0,J,J)
  for (i in 1:J) {
    for (j in 1:J){
      if (type=="standard") {
        v[i,j] <- compositions::var(l[,i]-l[,j], use=use)
      } else if (type=="phi") {
        v[i,j] <- compositions::var(l[,i]-l[,j], use=use)/compositions::var(l[,i], use=use)
      } else if (type=="phis") {
        v[i,j] <- compositions::var(l[,i]-l[,j], use=use)/(compositions::var(l[,i]+l[,j], use=use))
      } else if (type=="rho") {
        v[i,j] <- 2*compositions::cov(l[,i],l[,j], use=use)/(compositions::var(l[,i], use=use)+compositions::var(l[,j], use=use))
      }
    }
  }

  if (type=="logx") v <- compositions::cov(l, use=use)

  return(v)
}


#' Full Logx Variance-Covariance
#'
#' Estimates the variance-covariance of the log of the data, using a
#' Taylor-series approximation. This function differs from the function
#' \code{Logx Variance-Covariance} in that the resultant matrix includes a reference category.
#'
#' @param mu The mean vector of a dataset following the multinomial logit-Normal model
#' @param Sigma The sigma matrix of a dataset following the multinomial logit-Normal model
#' @param transf The desired transformation. If \code{transf="alr"} the inverse
#' additive logratio transformation is applied. If \code{transf="clr"} the
#' inverse centered logratio transformation is applied.
#' @param order The desired order of the Taylor Series approximation
#'
#' @return The estimated variance-covariance matrix, \code{logx}.
#'
#' @examples
#' mu <- mle$est.min$mu
#' Sigma <- mle$est.min$Sigma
#'
#' #Second order approximation of the variance-covariance matrix of the log of the data,
#' #with an alr transformation.
#' logVarTaylorFull(mu, Sigma, transf="alr", order="second")
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
