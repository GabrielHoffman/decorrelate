# Gabriel Hoffman April 8, 2021 Estimate covariance matrix after applying
# transformation



#' Estimate covariance matrix after applying transformation
#'
#' Given covariance between features in the original data, estimate the covariance matrix after applying a transformation to each feature.  Here we use the eclairs decomposition of the original covariance matrix, perform a parametric bootstrap and return the eclairs decomposition of the covariance matrix of the transformed data.
#'
#' @param ecl covariance/correlation matrix as an \link{eclairs} object
#' @param f function specifying the transformation.
# The default function is \eqn{\log(x^2 + 1e-3)} which offsets \eqn{x^2} away
# from zero since \eqn{\log(x^2)} is not defined for \eqn{x=0}.
#' @param n.boot number of parametric bootstrap samples.  Increasing n gives more precise estimates.
#' @param lambda shrinkage parameter.  If not specified, it is estimated from the data.
#' @param compute evaluate either the \code{"covariance"} or \code{"correlation"} of \code{X}
#' @param seed If you want the same to be generated again use a seed for the generator, an integer number
#'
#' @details
#' When the transformation is linear, these covariance matrices are the same.
#'
#' @return \link{eclairs} decomposition representing correlation/covariance on the transformed data
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
#' # sample matrix from MVN with covariance Sigma
#' Y <- rmvnorm(n, rep(0, p), sigma = Sigma, seed = 1)
#'
#' # perform eclairs decomposition
#' ecl <- eclairs(Y)
#'
#' # Parametric boostrap to estimate covariance
#' # after transformation
#'
#' # transformation function
#' f <- function(x) log(x^2 + 1e-3)
#'
#' # number of bootstrap samples
#' n_boot <- 5000
#'
#' # Evaluate eclairs decomposition on boostrap samples
#' ecl2 <- cov_transform(ecl, f = f, n_boot, lambda = 1e-4)
#'
#' # Get full covariance matrix from eclairs decomposition
#' C1 <- getCov(ecl2)
#'
#' # Parametric boostrap samples directly from full covariance matrix
#' X <- rmvnorm(n_boot, rep(0, p), getCov(ecl))
#'
#' # get covariance of transformed data
#' C2 <- cov(f(X))
#'
#' # Evaluate differences
#' # small differences are due to Monte Carlo error from boostrap sampling
#' range(C1 - C2)
#'
#' # Plot entries from two covariance estimates
#' par(pty = "s")
#' plot(C1, C2, main = "Concordance between covariances")
#' abline(0, 1, col = "red")
#'
#' # Same above but compute eclairs for correlation matrix
#' #-------------------------------------------------------
#'
#' # Evaluate eclairs decomposition on boostrap samples
#' ecl2 <- cov_transform(ecl, f = f, n_boot, compute = "correlation", lambda = 1e-4)
#'
#' # Get full covariance matrix from eclairs decomposition
#' C1 <- getCor(ecl2)
#'
#' # Parametric boostrap samples directly from full covariance matrix
#' X <- rmvnorm(n_boot, rep(0, p), getCov(ecl))
#'
#' # get correlation of transformed data
#' C2 <- cor(f(X))
#'
#' # Evaluate differences
#' # small differences are due to Monte Carlo error from boostrap sampling
#' range(C1 - C2)
#'
#' # Plot entries from two correlation estimates
#' oldpar <- par(pty = "s")
#' plot(C1, C2, main = "Correlation between covariances")
#' abline(0, 1, col = "red")
#'
#' par(oldpar)
#' @export
cov_transform <- function(ecl, f, n.boot, lambda = NULL, compute = c(
                            "covariance",
                            "correlation"
                          ), seed = NULL) {
  stopifnot(is(ecl, "eclairs"))
  stopifnot(is(f, "function"))
  compute <- match.arg(compute)

  mu <- rep(0, ecl$p)

  # draw samples from the multivariate normal with covariance ecl
  X <- rmvnorm_eclairs(n.boot, mu, ecl, seed = seed)

  # Perform eclairs decomposition of sampled data after applying transform f
  # set scale=FALSE to get covariance
  eclairs(f(X), lambda = lambda, compute = compute)
}
