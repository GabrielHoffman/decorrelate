# Gabriel Hoffman May 24, 2021 Estimate shrinkage parameter by emprical Bayes
# Simple R adaptation of the Rcpp code of CRAN's beam package Current code only
# works when target matrix is identity

# Updated Apr 23, 2024 to handle truncated SVD, which is different than n < p case




# double alphaToDelta(double alpha, int n, int p){ return
# (alpha*n+(1-alpha)*p+(1-alpha))/(1-alpha); }
alphaToDelta <- function(alpha, n, p) {
  (alpha * n + (1 - alpha) * p + (1 - alpha)) / (1 - alpha)
}

# double deltaToAlpha(double delta, int n, int p){ return
# (delta-p-1)/(n+delta-p-1); }
deltaToAlpha <- function(delta, n, p) {
  ifelse(is.finite(delta), (delta - p - 1) / (n + delta - p - 1), 1)
}

# Get nu so that (1-lambda) * Sigma_mle + nu*lambda*I
# has a diagonal closest to identity.
getNu <- function(ecl, k, a.value) {
  if (ecl$lambda == 0 | ecl$lambda == 1) {
    nu <- 1
  } else {
    # estimate nu to satisfy
    # p = (1-ecl$lambda)*a.value + ecl$lambda*nu*p
    nu <- with(ecl, (p - (1 - lambda) * a.value) / (p * lambda))
  }

  nu
}


# Integrated log-likelihood for empirical Bayes
#' @importFrom CholWishart lmvgamma
logML <- function(delta, ecl, k, a.value) {
  p <- ecl$p
  n <- ecl$n

  ecl$lambda <- deltaToAlpha(delta, n, p)
  nu <- getNu(ecl, k, a.value)

  eigs <- ecl$dSq[seq(k)] * ecl$n

  out <- -0.5 * n * p * log(pi)
  out <- out + lmvgamma((delta + n) * 0.5, p)
  out <- out - lmvgamma(delta * 0.5, p)
  out <- out + 0.5 * delta * p * log(delta - p - 1)

  # eigs.app <- append(eigs, rep(0, p - length(eigs)))
  # out <- out - 0.5 * (delta + n) * sum(log(nu*(delta - p - 1) + eigs.app))

  out <- out - 0.5 * (delta + n) * (sum(log(nu * (delta - p - 1) + eigs)) + (p - k) * log(nu * (delta - p - 1)))

  # out - 0.5 * n * logdetD
  # out - 0.5 * n * 2*p*log(sqrt(nu))
  out - 0.5 * n * p * log(nu)
}


#' Estimate shrinkage parameter by empirical Bayes
#'
#' Estimate shrinkage parameter by empirical Bayes
#'
#' @param ecl \code{eclairs()} decomposition
#' @param k number of singular vectors to use
#' @param minimum minimum value of lambda
#' @param lambda (default: NULL) If NULL, estimate lambda from data. Else evaluate logML using specified lambda value.
#'
#' @details Estimate shrinkage parameter for covariance matrix estimation using empirical Bayes method (Hannart and Naveau, 2014; Leday and Richardson 2019).  The shrinage estimate of the covariance matrix is \eqn{(1-\lambda)\hat\Sigma + \lambda \nu I}, where \eqn{\hat\Sigma} is the sample covariance matrix, given a value of \eqn{\lambda}.  A large value of \eqn{\lambda} indicates more weight on the prior.
#'
#' @return value \eqn{\lambda} and \eqn{\nu} indicating the shrinkage between sample and prior covariance matrices.
#'
#' @examples
#' library(Rfast)
#'
#' n <- 800 # number of samples
#' p <- 200 # number of features
#'
#' # create correlation matrix
#' Sigma <- autocorr.mat(p, .9)
#'
#' # draw data from correlation matrix Sigma
#' Y <- rmvnorm(n, rep(0, p), sigma = Sigma * 5.1, seed = 1)
#' rownames(Y) <- paste0("sample_", seq(n))
#' colnames(Y) <- paste0("gene_", seq(p))
#'
#' # eclairs decomposition: covariance
#' ecl <- eclairs(Y, compute = "correlation")
#'
#' # For full SVD
#' getShrinkageParams(ecl)
#'
#' # For truncated SVD at k = 20
#' getShrinkageParams(ecl, k = 20)
#'
#' @references
#' Hannart, A., & Naveau, P. (2014). Estimating high dimensional covariance matrices: A new look at the Gaussian conjugate framework. Journal of Multivariate Analysis, 131, 149-162.
#'
#' Leday, G. G., & Richardson, S. (2019). Fast Bayesian inference in large Gaussian graphical models. Biometrics, 75(4), 1288-1298.
#
#' @export
getShrinkageParams <- function(ecl, k = ecl$k, minimum = 1e-04, lambda = NULL) {
  n <- ecl$n
  p <- ecl$p

  # factors of Sigma
  # sum of diagonals of Sigma
  ecl$k <- k
  A <- with(ecl, dmult(U[, seq(k), drop = FALSE], sqrt(dSq[seq(k)]), "right"))
  a.value <- sum(colSums(A^2))

  # estimate lambda
  if (is.null(lambda)) {
    # get upper and lower values of range
    lowerVal <- alphaToDelta(minimum, n, p)
    upperVal <- alphaToDelta(1 - minimum, n, p)

    # get optimal estimate of delta using log-likelihood
    res <- optimize(function(x) logML(x, ecl, k, a.value), c(lowerVal, upperVal), maximum = TRUE)

    ll.value <- res$objective
    ecl$lambda <- deltaToAlpha(res$maximum, n, p)
  } else {
    # compute logLik
    ecl$lambda <- lambda
    delta <- alphaToDelta(lambda, n, p)
    ll.value <- logML(delta, ecl, k, a.value)
  }

  nu <- getNu(ecl, k, a.value)

  data.frame(lambda = ecl$lambda, logLik = ll.value, nu = nu)
}
