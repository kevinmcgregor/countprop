globalVariables(c("n", "mu", "Sigma", "J"))

#' Maximum Likelihood Estimate
#' 
#' Returns the maximum likelihood estimates of multinomial logit-Normal model 
#' parameters given a count-compositional dataset. The MLE procedure is based on the 
#' multinomial logit-Normal distribution, using the EM algorithm from Hoff (2003).
#'
#' @param y Count-compositional dataset
#' @param max.iter Maximum number of iterations 
#' @param max.iter.nr Maximum number of Newton-Raphson iterations
#' @param tol Stopping rule
#' @param tol.nr Stopping rule for the Newton-Raphson algorithm
#' @param lambda.gl Penalization parameter lambda, for the graphical lasso penalty. Controls
#' the sparsity of Sigma
#' @param gamma Gamma value for EBIC calculation of the log-likelihood
#'
#' @return The additive log-ratio of y (\code{v}); maximum likelihood estimates of 
#' \code{mu}, \code{Sigma}, and \code{Sigma.inv};
#' the log-likelihood (\code{log.lik}); the EBIC (extended Bayesian information criterion) 
#' of the log-likelihood of the multinomial logit-Normal model with the 
#' graphical lasso penalty (\code{ebic}); degrees of freedom of the \code{Sigma.inv}
#' matrix (\code{df}).
#' 
#' @note The graphical lasso penalty 
#' is the sum of the absolute value of the elements of the covariance matrix \code{Sigma}.
#' The penalization parameter lambda controls the sparsity of Sigma.
#' 
#' @note This function is also used within the \code{mlePath()} function.
#' 
#' @examples
#' n <- 100
#' mle_dat <- mleLR(dat.ss, lambda.gl=0.5)
#' 
#' mle_dat$mu #mle mu of dat.ss
#' mle_dat$Sigma #mle Sigma of dat.ss
#' mle_dat$ebic #ebic of the fitted model
#' 
#' 
#' @export
#' 
mleLR <- function(y, max.iter=10000, max.iter.nr=100, tol=1e-6, tol.nr=1e-6, lambda.gl=0, gamma=0.1) {
  n <- NROW(y)
  k <- NCOL(y)
  ni <- rowSums(y)
  # observed proportions
  #p.hat <- y/rowSums(y)
  # initialize 
  pseudo.count <- 0.1
  v <- unclass(compositions::alr(y+pseudo.count))
  attr(v, "orig") <- NULL
  attr(v, "V") <- NULL
  colnames(v) <- NULL
  mu <- rep(0, k-1)
  mu.old <- mu+1
  
  # Sigma <- cov(v)
  # Sigma.old <- Sigma+0.1
  # i.inv <- array(0, dim=c(k-1,k-1,n))
  
  Sigma.inv <- diag(k-1) #MASS::ginv(cov(v))
  Sigma.inv.old <- Sigma.inv+0.1
  i.inv <- array(0, dim=c(k-1,k-1,n))
  
  count <- 1
  #while(sum(abs((mu-mu.old)/mu.old))>tol & sum(abs((diag(Sigma.inv)-diag(Sigma.inv.old))/diag(Sigma.inv.old)))>tol) {
  while(max(abs((mu-mu.old)/mu.old))>tol & max(abs((diag(Sigma.inv)-diag(Sigma.inv.old))/diag(Sigma.inv.old)))>tol) { 
    if (count%%10==0) {
      cat("Iter=", count, "\n")
      cat("Error=", max(abs((mu-mu.old)/mu.old)), "\n")
      #print(Sigma[1,1])
      print(mu)
      #print(eigen(h)$values)
    }
    if (count>max.iter) stop("Maximum number of iterations reached")
    
    #Sigma.inv <- qr.solve(Sigma)
    
    for (i in 1:n) {
      # cat("i=", i, "\n")
      g <- rep(1, k-1)
      count.nr <- 1
      #while(abs(sum(g))>tol.nr) {
      while(max(abs(g))>tol.nr) {
        if (count.nr>max.iter.nr) {
          print(g)
          stop("Maximum number of Newton-Raphson iterations reached")
        }
        g <- grad(v[i,], y[i,-k], ni[i], mu, Sigma.inv)
        h <- hess(v[i,], ni[i], Sigma.inv)
        #if (any(is.complex(eigen(h)$values))) browser()
        i.inv[,,i] <- qr.solve(-h)
        # Newton-Raphson step
        v[i,] <- v[i,] + i.inv[,,i]%*%g
        count.nr <- count.nr + 1
      }
    }
    
    mu.old <- mu
    mu <- colMeans(v)
    #Sigma.old <- Sigma
    Sigma.inv.old <- Sigma.inv
    #Sigma.inv <- huge::huge(v, lambda.gl, method="glasso", verbose=FALSE)$icov[[1]]
    gl <- suppressWarnings(glasso::glasso(compositions::cov(v), lambda.gl))
    Sigma.inv <- gl$wi
    Sigma <- gl$w
    # Sigma <- matrix(0, k-1, k-1)
    # for (i in 1:n){
    #   # Need to update Hessian to final value of v[i,]
    #   h <- hess(v[i,], ni[i], Sigma.inv)
    #   neg.hess.inv <- qr.solve(-h)
    #   Sigma <- Sigma + (tcrossprod(v[i,]-mu) + neg.hess.inv)/n
    # }
    count <- count+1
    
  }
  
  S <- compositions::cov(v)
  
  log.lik <- logLik(v, y, ni, S, Sigma.inv)
  # Currently only using Gaussian part for loglik
  #log.lik <- logLikG(v, S, Sigma.inv)
  
  df <- sum(Sigma.inv[lower.tri(Sigma.inv)]!=0)
  eb <- ebic(log.lik, n, k-1, df, gamma)
  
  # return(list(v=v, mu=mu, Sigma=Sigma, var.v=neg.hess.inv))
  return(list(v=v, mu=mu, Sigma.inv=Sigma.inv, Sigma=Sigma, log.lik=log.lik, 
              ebic=eb, df=df))
}



#' Wrapper for MleLR()
#' 
#' Executes a function call to \code{mlePath()}. Helps with the \code{mclapply()} 
#' within \code{mlePath()}
#'
#' @param x Input argument for the \code{MLE Function}
#' 
#' 
wrapMLE <- function(x) {
  do.call(mleLR, x)
}


#' Maximum Likelihood Estimator Paths
#' 
#' Calculates the maximum likelihood estimates of the parameters for the 
#' mutlinomial logit-Normal distribution under various values
#' of the penalization parameter \code{lambda}. Parameter \code{lambda} controls
#' the sparsity of the covariance matrix \code{Sigma}, and penalizes the false 
#' large correlations that may arise in microbiome data when a large
#' number of OTU-associations are being measured. 
#'
#' @param y Count-compositional dataset
#' @param max.iter Maximum number of iterations
#' @param max.iter.nr Maximum number of Newton-Raphson iterations
#' @param tol Stopping rule
#' @param tol.nr Stopping rule for the Newton Raphson algorithm
#' @param lambda.gl Vector of penalization parameters lambda, for the graphical lasso penalty
#' @param lambda.min.ratio Minimum lambda ratio of the maximum lambda, 
#' used for the sequence of lambdas
#' @param n.lambda Number of lambda to evaluate different paths for
#' @param n.cores Number of cores to use (for parallel computation)
#' @param gamma Gamma value for EBIC calculation of the log-likelihood
#'
#' @return The MLE estimates of \code{y} for each element lambda of lambda.gl, (\code{est}); 
#' the value of the estimates which produce the minimum EBIC, (\code{est.min}); 
#' the vector of lambdas used for graphical lasso, (\code{lambda.gl}); the index of 
#' the minimum EBIC (extended Bayesian information criterion), (\code{min.idx}); 
#' vector containing the EBIC for each lambda, (\code{ebic}).
#' 
#' @note If using parallel computing, consider setting \code{n.cores} to be equal
#' to the number of lambdas being evaluated for, \code{n.lambda}.
#' 
#' @note The graphical lasso penalty 
#' is the sum of the absolute value of the elements of the covariance matrix \code{Sigma}.
#' The penalization parameter lambda controls the sparsity of Sigma. 
#' 
#' @examples
#' n <- 100
#' mle.sim <- mlePath(dat.ss, tol=1e-4, tol.nr=1e-4, n.lambda = 8, lambda.min.ratio = 0.01, gamma = 0.1, n.cores = 1)
#' 
#' mu.hat <- mle.sim$est.min$mu #optimal mle of mu based on ebic 
#' Sigma.hat <- mle.sim$est.min$Sigma #optimal mle of Sigma based on ebic 
#' Sigma.hat <- mle.sim$est.min$ebic #minimum ebic 
#' 
#' 
#' @export
#' 
mlePath <- function(y, max.iter=10000, max.iter.nr=100, tol=1e-6, tol.nr=1e-6, lambda.gl=NULL,
                    lambda.min.ratio=0.1, n.lambda=1,
                    n.cores=NULL, gamma=0.1) {
  
  k <- NCOL(y)
  
  # Set lambda.gl if it's not specified by the user
  if (is.null(lambda.gl)) {
    pc <- 0.1
    alr.y <- log(y[,-k]+pc) - log(y[,k]+pc)
    d <- NCOL(alr.y)
    S <- compositions::cov(alr.y)
    lmax <- max(max(S - diag(d)), -min(S - diag(d)))
    lmin <- lambda.min.ratio*lmax
    lambda.gl <- exp(seq(log(lmin), log(lmax), length.out=n.lambda))
  }
  
  n.lam <- length(lambda.gl)
  
  m.pars <- vector("list", length=n.lam)
  for (i in 1:n.lam) {
    m.pars[[i]] <- list(y, max.iter, max.iter.nr, tol, tol.nr, lambda.gl[i], gamma)
  }
  
  if (is.null(n.cores)) n.cores <- parallel::detectCores()
  
  est <- parallel::mclapply(m.pars, wrapMLE, mc.cores = n.cores)
  
  ebic.vec <- unlist(lapply(est, function(x){x$ebic}))
  wm <- which.min(ebic.vec)
  
  est.min <- est[[wm]]
  
  return(list(est=est, est.min=est.min,
              lambda.gl=lambda.gl, min.idx=wm, ebic=ebic.vec))
}

#' Gradient of Normal Random Variables
#' 
#' Calculates the gradient of the normal random variables, on the logit scale.
#'
#' @param v Additive logratio (alr) transform of y
#' @param y Count compositional dataset
#' @param ni Row sums of y
#' @param mu Mu vector of y
#' @param Sigma.inv Sigma inverse matrix of y
#'
#' @return The gradient vector of Normal random variables.
#'
#' 
#' 
grad <- function(v, y, ni, mu, Sigma.inv) {
  ev <- exp(v)
  c(y - ni*ev/(1+sum(ev)) - Sigma.inv%*%(v-mu))
}

#' Hessian 
#' 
#' Calculates the hessian matrix.
#'
#' @param v Count-compositional dataset which has been transformed by the additive logratio
#' @param ni Row sums of the raw data 
#' @param Sigma.inv Inverse of the Sigma matrix
#'
#' @return The hessian matrix.
#'
#' 
#' 
hess <- function(v, ni, Sigma.inv) {
  ev <- exp(v)
  p.hat <- ev/(1+sum(ev))
  ni*(tcrossprod(p.hat)-diag(p.hat)) - Sigma.inv
}


#' Gaussian Log-Likelihood
#' 
#' Calculates the Gaussian log-likelihood, under the multinomial logit-Normal model.
#'
#' @param v Count-compositional dataset which has been transformed by the additive
#' logratio
#' @param S Covariance of \code{v}
#' @param invSigma Inverse of the Sigma matrix
#'
#' @return The estimated Gaussian log-likelihood under the Multinomial logit-Normal distribution.
#' 
#' 
#'
logLikG <- function(v, S, invSigma) {
  n <- NROW(v)
  ldet <- determinant(invSigma, logarithm = TRUE)$modulus
  n*0.5*(ldet - sum(diag(S%*%invSigma)))
}

#' Log-Likelihood
#' 
#' Calculates the log-likelihood, under the multinomial logit-Normal model.
#'
#' @param v The additive logratio transform of y
#' @param y Compositional dataset
#' @param ni The row sums of y
#' @param S Covariance of \code{v}
#' @param invSigma The inverse of the Sigma matrix
#'
#' @return The estimated log-likelihood under the Multinomial logit-Normal distribution.
#'
#'
#'
logLik <- function(v, y, ni, S, invSigma) {
  n.sp <- NCOL(y)
  rs <- log(rowSums(exp(v)+1))
  ldet <- determinant(invSigma, logarithm = TRUE)$modulus
  sum(y[,-n.sp]*v) - sum(ni*rs) + n*0.5*(ldet - sum(diag(S%*%invSigma)))
}


#' Extended Bayesian Information Criterion
#' 
#' Calculates the Extended Bayesian Information Criterion (EBIC) of a model.
#' Used for model selection to asses the fit of the multinomial logit-Normal 
#' model which includes a graphical lasso penalty.
#'
#' @param l Log-likelihood estimates of the model
#' @param n Number of rows of the data set for which the log-likelihood has been 
#' calculated
#' @param d The size of the (k-1) by (k-1) covariance matrix of a 
#' k by k count-compositional data matrix
#' @param df Degrees of freedom 
#' @param gamma A tuning parameter. Larger values means more penalization
#'
#' @return The value of the EBIC.
#' 
#' @note The graphical lasso penalty 
#' is the sum of the absolute value of the elements of the covariance matrix \code{Sigma}.
#' The penalization parameter lambda controls the sparsity of Sigma. 
#' 
#' @examples 
#' log.lik_1 <- mle$est[[1]]$log.lik
#' n <- NROW(dat.ss)
#' k <- NCOL(dat.ss)
#' df_1 <- mle$est[[1]]$df
#' 
#' ebic(log.lik_1, n, k, df_1, gamma=0.1)
#' 
#' @export
#'
ebic <- function(l, n, d, df, gamma) {
  -2*l + log(n)*df + 4*gamma*log(d)*df
}

#' Extended Bayesian Information Criterion Plot
#'
#' Plots the extended Bayesian information criterion (EBIC) of the model fit for
#' various penalization parameters \code{lambda}. 
#'
#' @param fit The model fit 
#' @param xlog TRUE or FALSE. Renders plot with the x-axis in the log-scale if \code{TRUE}
#'
#' @return Plot of the EBIC (y-axis) against each lambda (x-axis).
#'
#' @examples 
#' mle$ebic #y-axis of ebic plot
#' log(mle$lambda.gl) #x-axis of ebic plot
#' 
#' ebicPlot(mle, xlog = TRUE) 
#'
#' @export
#' 
#'
ebicPlot <- function(fit, xlog=FALSE) {
  if (xlog) {
    x.ax <- log(fit$lambda.gl)
  } else {
    x.ax <- fit$lambda.gl
  }
  plot(x.ax, fit$ebic, type="o", pch=19, cex=0.5,
       col="darkred", xlab="log(lambda)", ylab="EBIC")
}

#' Root Mean Square Error
#' 
#' Calculates the root mean square error (RMSE) between two values (i.e. scalers, 
#' vectors, or matrices).
#'
#' @param x Value one for comparison (a scaler, vector, or matrix)
#' @param y Value two for comparison (a scaler, vector, or matrix). Must be the 
#' same form as \code{x}
#'
#' @return A single scaler value, the RMSE.
#'
#' @examples
#' x <- sample(1:20, size = 4)
#' y <- sample(1:20, size = 4)
#' 
#' rmse(x,y)
#' 
#' @export
#' 
rmse <- function(x,y) {
  sqrt(mean((x-y)^2))
}

#' Row Root Mean Square Error
#' 
#' Calculates the root mean square error (RMSE) between each row of two matrices.
#'
#' @param x First matrix for comparison
#' @param y Second matrix for comparison. Must be the same dimensions as \code{x}
#'
#' @return A vector. Element \code{i} of the resultant vector is the RMSE which 
#' compares row \code{i} of matrix \code{x}, and row \code{i} of matrix \code{y}. 
#'
#' @examples
#' x <- matrix(sample(1:500, size = 21), nrow=3)
#' y <- matrix(sample(1:500, size = 21), nrow=3)
#' 
#' rmse_by_row(x,y)
#' 
#' @export
#' 
rmse_by_row <- function(x,y) {
  n <- NROW(x)
  r <- rep(0, n)
  for (i in 1:n) {
    r[i] <- rmse(x[i,], y[i,])
  }
  return(r)
}