# Gabriel Hoffman


#' Summarize correlation matrix
#'
#' Summarize correlation matrix as a scalar scalar value, given its SVD and shrinkage parameter. \code{averageCorr()} computes the average correlation, and \code{averageCorrSq()} computes the average squared correlation, where both exclude the diagonal terms. \code{sumInverseCorr()} computes the sum of entries in the inverse correlation matrix to give the 'effective number of independent features'. \code{effVariance()} evaluates effective variance of the correlation (or covariance) matrix.  These values can be computed using the correlation matrix using standard MLE, or EB shrinkage.
#'
#' @param ecl estimate of correlation matrix from \code{eclairs()} storing \eqn{U}, \eqn{d_1^2}, \eqn{\lambda} and \eqn{\nu}
#' @param method compute average correlation for either the empirical Bayes (EB) shinken correlation matrix or the MLE correlation matrix
#'
#' @return value of summary statistic
#' @details
#' \code{tr()}: trace of the matrix. Sum of diagonals is the same as the sum of the eigen-values.
#'
#' \code{averageCorr()}: The average correlation is computed by summing the off-diagonal values in the correlation matrix. The sum of all elements in a matrix is \eqn{g = \sum_{i,j} C_{i,j} = 1^T C 1 }, where \eqn{1} is a vector of \eqn{p} elements with all entries 1. This last term is a quadratic form of the correlation matrix that can be computed efficiently using the SVD and shrinkage parameter from \code{eclairs()}.  Given the value of \eqn{g}, the average is computed by subtracting the diagonal values and dividing by the number of off-diagonal values: \eqn{(g - p) / (p(p-1))}.
#'
#' \code{averageCorrSq()}: The average squared correlation is computed using only the eigen-values. Surprisingly, this is a function of the variance of the eigen-values.  The is reviewed by \insertCite{watanabe2022statistics;textual}{decorrelate} and \insertCite{durand2017linkage;textual}{decorrelate}.  Letting \eqn{\lambda_i} be the \eqn{i^{th}} sample or shrunk eigen-value, and \eqn{\tilde{\lambda}} be the mean eigen-value, then \eqn{\sum_i (\lambda_i - \tilde{\lambda})^2 / p(p-1)\tilde{\lambda}^2}.
#'
#' \code{sumInverseCorr()}: The 'effective number of independent features' is computed by summing the entires of the inverse covariance matrix.  This has the form \eqn{\sum_{i,j} C^{-1}_{i,j} = 1^T C^{-1} 1}. This last term is a quadratic form of the correlation matrix that can be computed efficiently using the SVD and shrinkage parameter from \code{eclairs()} as described above.
#'
#' \code{effVariance()}: Compute a metric of the amount of variation represented by a correlation (or covariance) matrix that is comparable across matrices of difference sizes. Proposed by \insertCite{pena2003descriptive;textual}{decorrelate}, the 'effective variance' is \eqn{|C|^\frac{1}{p}} where \eqn{C} is a correlation (or covariance matrix) between \eqn{p} variables. The effective variance is the mean of the log eigen-values.
#'
#'
#' @references
#' \insertAllCited{}
#'
#' @examples
#' library(Rfast)
#' n <- 200 # number of samples
#' p <- 800 # number of features
#'
#' # create correlation matrix
#' Sigma <- matrix(.2, p, p)
#' diag(Sigma) <- 1
#'
#' # draw data from correlation matrix Sigma
#' Y <- rmvnorm(n, rep(0, p), sigma = Sigma, seed = 1)
#' rownames(Y) <- paste0("sample_", seq(n))
#' colnames(Y) <- paste0("gene_", seq(p))
#'
#' # eclairs decomposition
#' ecl <- eclairs(Y, compute = "cor")
#'
#' # Average correlation value
#' averageCorr(ecl)
#'
#' # Average squared correlation value
#' averageCorrSq(ecl)
#'
#' # Sum elements in inverse correlation matrix
#' # Gives the effective number of independent features
#' sumInverseCorr(ecl)
#'
#' # Effective variance
#' effVariance(ecl)
#'
#' @rdname summarizeCorr
#' @export
averageCorr <- function(ecl, method = c("EB", "MLE")) {
  method <- match.arg(method)

  # check if scale of all features is 1
  if (!all(ecl$sigma == 1)) {
    stop("Only computable from a corrrelation matrix")
  }

  if (method == "MLE") {
    ecl$lambda <- 0
  }

  p <- ecl$p

  # Compute Grand sum: sum of all elements
  one <- rep(1, p)
  gs <- quadForm(ecl, one, alpha = 1 / 2)

  # subtract diagonal and divide by number of elements
  (gs - p) / (p * (p - 1))
}


# Also see https://doi.org/10.1111/evo.14382
#' @rdname summarizeCorr
#' @export
averageCorrSq <- function(ecl, method = c("EB", "MLE")) {
  method <- match.arg(method)

  # check if scale of all features is 1
  if (!all(ecl$sigma == 1)) {
    stop("Only computable from a corrrelation matrix")
  }

  if (method == "MLE") {
    ecl$lambda <- 0
  }

  p <- ecl$p

  # get eigen-values
  dSq <- c(ecl$dSq, rep(0, max(0, ecl$p - length(ecl$dSq))))

  # get shrunken eigen-values
  dSq.shrink <- dSq * (1 - ecl$lambda) + ecl$nu * ecl$lambda

  # mean eigen-value
  m <- mean(dSq.shrink)

  # mean r^2 value
  sum((dSq.shrink - m)^2) / (p * (p - 1) * m^2)
}

#' @rdname summarizeCorr
#' @export
sumInverseCorr <- function(ecl, method = c("EB", "MLE")) {
  method <- match.arg(method)

  # check if scale of all features is 1
  if (!all(ecl$sigma == 1)) {
    stop("Only computable from a corrrelation matrix")
  }

  if (method == "MLE") {
    ecl$lambda <- 0
  }

  p <- ecl$p
  one <- rep(1, p)

  quadForm(ecl, one)
}


#' @rdname summarizeCorr
#' @export
effVariance <- function(ecl, method = c("EB", "MLE")) {
  method <- match.arg(method)

  if (method == "MLE") {
    ecl$lambda <- 0
  }

  # Can be either a covariance or correlation matrix det(C)^(1/p) Compute in
  # log space then exponentiate
  exp(logDet(ecl) / ecl$p)
}


#' @rdname summarizeCorr
#' @export
tr <- function(ecl, method = c("EB", "MLE")) {
  method <- match.arg(method)

  # check if scale of all features is 1
  if (!all(ecl$sigma == 1)) {
    stop("Only computable from a corrrelation matrix")
  }

  if (method == "MLE") {
    ecl$lambda <- 0
  }

  # sum of eigen-values
  p <- ecl$p

  # get eigen-values
  dSq <- c(ecl$dSq, rep(0, max(0, ecl$p - length(ecl$dSq))))

  # get shrunken eigen-values
  dSq.shrink <- dSq * (1 - ecl$lambda) + ecl$nu * ecl$lambda

  sum(dSq.shrink)
}
