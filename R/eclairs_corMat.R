# Gabriel Hoffman June 1, 2021 extend eclairs() to run on correlation matrix

#' Estimate covariance/correlation with low rank and shrinkage
#'
#' Estimate covariance/correlation with low rank and shrinkage from the correlation matrix
#'
#' @param C sample correlation matrix between features
#' @param n number of samples used to estimate the sample correlation matrix
#' @param k the rank of the low rank component.  Defaults to min of sample size and feature number, \code{min(n, p)}.
#' @param lambda shrinkage parameter. If not specified, it is estimated from the data.
#'
#' @return \link{eclairs} object storing:
#' \describe{
#'  \item{U: }{orthonormal matrix with k columns representing the low rank component}
#'  \item{dSq: }{eigen-values so that \eqn{U diag(d^2) U^T} is the low rank component}
#'  \item{lambda: }{shrinkage parameter \eqn{\lambda} for the scaled diagonal component}
#'  \item{sigma: }{standard deviations of input columns}
#'  \item{nu: }{diagonal value, \eqn{\nu}, of target matrix in shrinkage}
#'  \item{n: }{number of samples (i.e. rows) in the original data}
#'  \item{p: }{number of features (i.e. columns) in the original data}
#'  \item{k: }{rank of low rank component}
#'  \item{rownames: }{sample names from the original matrix}
#'  \item{colnames: }{features names from the original matrix}
#'  \item{method: }{method used for decomposition}
#'  \item{call: }{the function call}
#' }
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
#' rownames(Y) <- paste0("sample_", 1:n)
#' colnames(Y) <- paste0("gene_", 1:p)
#'
#' # eclairs decomposition
#' eclairs(Y, compute = "correlation")
#'
#' # eclairs decomposition from correlation matrix
#' eclairs_corMat(cor(Y), n = n)
#'
#' @importFrom Rfast eachrow
#' @importFrom irlba partial_eigen
#' @importFrom methods new
#' @importFrom Matrix diag isSymmetric
#'
#' @export
eclairs_corMat <- function(C, n, k = min(n, nrow(C)), lambda = NULL) {
  # check that diagonals are all equal, up to a tolerance
  if (!all.equal(as.numeric(diag(C)), rep(1, ncol(C)))) {
    stop("All diagonal entries of C must be 1")
  }

  # check that C is symmetric
  if (!isSymmetric(C)) {
    stop("Matrix C is not symmetric")
  }

  if (missing(n)) {
    stop("n must be specified")
  }

  p <- nrow(C) # number of features
  nu <- mean(diag(C)) # scale of identity matrix

  # k is the min of the number of samples and features
  k <- min(n, p)

  if (k > p / 3) {
    dcmp <- eigen(C)
  } else {
    # partial eigen decomposition
    suppressWarnings({
      dcmp <- partial_eigen(C, k)
    })
  }

  # keep first k eigen values and vectors
  dcmp$values <- dcmp$values[seq_len(k)]
  dcmp$vectors <- dcmp$vectors[, seq_len(k), drop = FALSE]

  # keep only positive eigen-values
  if (any(dcmp$values <= .Machine$double.eps)) {
    keep <- which(dcmp$values > .Machine$double.eps)

    dcmp$values <- dcmp$values[keep]
    dcmp$vectors <- dcmp$vectors[, keep, drop = FALSE]
  }

  # Modify sign of dcmp$v and dcmp$u so principal components are consistant
  # This is motivated by whitening:::makePositivDiagonal() but here adjust
  # both U and V so reconstructed data is correct
  values <- sign0(diag(dcmp$vectors))

  # faster version
  dcmp$vectors <- eachrow(dcmp$vectors, values, "*")

  ecl <- list(
    U = dcmp$vectors,
    dSq = dcmp$values,
    V = NULL,
    lambda = NA,
    sigma = rep(1, ncol(C)),
    nu = NA,
    n = n,
    p = p,
    k = k,
    logLik = NA,
    rownames = NULL,
    colnames = colnames(C),
    method = "eigen",
    call = match.call()
  )

  ecl <- new("eclairs", ecl)

  # estimate lambda and nu values
  res <- getShrinkageParams(ecl, lambda = lambda)
  ecl$lambda <- res$lambda
  ecl$nu <- res$nu
  ecl$logLik <- res$logLik

  ecl
}
