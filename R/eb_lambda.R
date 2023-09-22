# Gabriel Hoffman May 24, 2021 Estimate shrinkage parameter by emprical Bayes
# Simple R adaptation of the Rcpp code of CRAN's beam package Current code only
# works when target matrix is identity




# double alphaToDelta(double alpha, int n, int p){ return
# (alpha*n+(1-alpha)*p+(1-alpha))/(1-alpha); }
alphaToDelta <- function(alpha, n, p) {
  (alpha * n + (1 - alpha) * p + (1 - alpha)) / (1 - alpha)
}

# double deltaToAlpha(double delta, int n, int p){ return
# (delta-p-1)/(n+delta-p-1); }
deltaToAlpha <- function(delta, n, p) {
  (delta - p - 1) / (n + delta - p - 1)
}

# Integrated log-likelihood for empirical Bayes eigs stores non-zero
# eigen-values of p total eigen values
#' @importFrom CholWishart lmvgamma
logML <- function(delta, p, n, eigs, logdetD) {
  out <- -0.5 * n * p * log(pi)
  out <- out + lmvgamma((delta + n) * 0.5, p)
  out <- out - lmvgamma(delta * 0.5, p)
  out <- out + 0.5 * delta * p * log(delta - p - 1)

  # in bream Rcpp code assums eigs are zero after rank n eigs[(n+1):p] =0
  if (n > p) {
    out <- out - 0.5 * (delta + n) * sum(log((delta - p - 1) + eigs))
  } else {
    # if full svd, length(eigs) is 0 but for low rank svd, use lenght(eigs)
    # out = out - 0.5*(delta+n)*(sum(log((delta-p-1)+eigs)))+
    # (p-n)*sum(log(delta-p-1)))
    out <- out - 0.5 * (delta + n) * (sum(log((delta - p - 1) + eigs)) + (p -
      length(eigs)) * sum(log(delta - p - 1)))
  }

  out <- out - 0.5 * n * logdetD

  out
}


getShrinkageValue <- function(n, p, eigs, logdetD, minimum = 1e-04, lambda = NULL) {
  # if lambda is NULL, estimate it
  if (is.null(lambda)) {
    # get upper and lower values of range
    lowerVal <- alphaToDelta(minimum, n, p)
    upperVal <- alphaToDelta(1 - minimum, n, p)

    # get optimal estiamte of delta using log-likelihood
    res <- optimize(logML,
      lower = lowerVal, upper = upperVal, p = p, n = n, eigs = eigs,
      logdetD = logdetD, maximum = TRUE
    )

    # convert delta to alpha between 0 and 1
    lambda <- deltaToAlpha(res$maximum, n, p)
    value <- res$objective
  } else {
    # compute logML
    delta <- alphaToDelta(lambda, n, p)
    value <- logML(delta, p = p, n = n, eigs = eigs, logdetD = logdetD)
  }

  list(lambda = lambda, logML = value)
}




#' Estimate shrinkage parameter by empirical Bayes
#'
#' Estimate shrinkage parameter by empirical Bayes
#'
#' @param ev array of eigen values
#' @param n number of samples
#' @param p number of features
#' @param nu scale of prior covariance matrix
#' @param lambda (default: NULL) If NULL, estimate lambda from data. Else evaluate logML using specified lambda value.
#'
# @seealso \link{getShrinkageValue}
#'
#' @details Estimate shrinkage parameter for covariance matrix estimation using empirical Bayes method \insertCite{leday2019fast,hannart2014estimating}{decorrelate}.  The shrinage estimate of the covariance matrix is \eqn{(1-\lambda)\hat\Sigma + \lambda \nu I}, where \eqn{\hat\Sigma} is the sample covariance matrix, given a value of \eqn{\lambda}.  A large value of \eqn{\lambda} indicates more weight on the prior.
#'
#' @return value \eqn{\lambda} indicating the shrinkage between sample and prior covariance matrices.
#'
#' @examples
#' ev <- c(10, 2, 1) # eigen values
#' n <- 12 # samples
#' p <- 3 # features
#' nu <- 2 # scale of target covariance
#'
#' estimate_lambda_eb(ev, n, p, nu)
#'
#' @references \insertAllCited{}
#'
#' @import Rdpack
#' @export
estimate_lambda_eb <- function(ev, n, p, nu, lambda = NULL) {
  # when D = diag(nu,p), logdet(D) is 2*p*log(sqrt(nu)) and the eigen values
  # of X %*% Dinv %*% t(X) equal eigen(tcrossprod(X)) / nu \tso divide the
  # eigen-values by nu

  # D = diag(nu, p) logdetD = 2*sum(log(diag(chol(D))))\t

  # if is low rank
  if (length(ev) < min(n, p)) {
    # if eigen-values are truncated, include additial eigen-values so that
    # sum(ev) equals the total variance, regardless of rank Note that
    # eclairs calls \testimate_lambda_eb( n*dcmp$d^2, n, p, nu) so
    # eigen-values are already scaled by n
    totalVar <- p * n * nu

    # if the emprical sum of variance is larger then the theoretical then
    # the specified n is likely mis-specified
    if (sum(ev) > totalVar) {
      warning("The sample size n is likely mis-specified")
    }

    # min(n,p) so works regardless of n > p or m < p can give over-estimate
    # of lambda for small k \tbut the deviation from full k can be huge idx
    # = seq(length(ev)+1, min(n,p)) ev[idx] = (totalVar-sum(ev)) /
    # length(idx)

    # can give under-estimate of lambda for small k but much smaller
    # difference than above
    len <- min(n, p) - length(ev)
    remainingVar <- totalVar - sum(ev)
    ev <- append(ev, series_start_total(min(ev), remainingVar, len))
  }

  # estimate optimal lambda (i.e. alpha) value
  getShrinkageValue(n, p, ev / nu,
    logdetD = 2 * p * log(sqrt(nu)), minimum = 1e-04,
    lambda = lambda
  )
}




#' Create decreasing series of values
#'
#' Create series of \eqn{n} values decreasing from \eqn{start} with a constant ratio between terms where the values sum to \eqn{totalSum}
#'
#' @param start maximum value
#' @param totalSum specify target sum of values produces
#' @param n number of terms
#'
# @examples start = 10 totalSum = 40 n = 4 values = series_start_total(start,
# totalSum, n)
series_start_total <- function(start, totalSum, n) {
  # get alpha so that: sum(start*alpha^seq(1,n)) == totalSum

  f <- function(alpha, start, n, totalSum) {
    (sum(start * alpha^seq(1, n)) - totalSum)^2
  }

  fit <- optimize(f, c(1e-12, 1 - 1e-12),
    start = start, n = n, totalSum = totalSum,
    tol = 1e-10
  )
  alpha <- fit$minimum

  start * alpha^seq(1, n)
}
