# Gabriel Hoffman
# April 8, 2021
#
# Estimate covariance matrix after applying transformation



#' Estimate covariance matrix after applying transformation
#'
#' Given covariance between features in the original data, estimate the covariance matrix after applying a transformation to each feature.  Here we use the eclairs decomposition of the original covariance matrix, perform a parametric bootstrap and return the eclairs decomposition of the covariance matrix of the transformed data.  
#'
#' @param Sigma.eclairs covariance/correlation matrix as an \link{eclairs} object
#' @param n number of parametric bootstrap samples.  Increasing n gives more precise estimates.
#' @param f function specifying the transformation.  The default function is \eqn{2\log(|x| + 1e-4)} which is a stable approximation of \eqn{\log(x^2)}
#' @param lambda shrinkage parameter.  If not specified, it is estimated from the data.
#' @param compute compute the 'covariance' (default) or 'correlation'
#'
#' @details
#' When the transformation is linear, these covariance matrices are the same.   
#'
#' @examples
#' library(Matrix)
#' library(Rfast)
#' set.seed(1)
#' n = 2000 # number of samples
#' p = 8*6 # number of features
#' 
#' # Create correlation matrix with autocorrelation
#' autocorr.mat <- function(p = 100, rho = 0.9) {
#'     mat <- diag(p)
#'     return(rho^abs(row(mat)-col(mat)))
#' }
#' 
#' # create correlation matrix
#' Sigma = autocorr.mat(p/8, .9)
#' Sigma = bdiag(Sigma, Sigma)
#' Sigma = bdiag(Sigma, Sigma)
#' Sigma = bdiag(Sigma, Sigma)
#' 
#' # sample matrix from MVN with covariance Sigma
#' Y = rmvnorm(n, rep(0, p), sigma=Sigma)
#' 
#' # perform eclairs decomposition
#' ecl = eclairs(Y)
#' 
#' # Parametric boostrap to estimate covariance
#' # after transformation
#' 
#' # transformation function is stable approximation of log(x^2)
#' f = function(x) 2*log(abs(x)+1e-4)
#' 
#' # number of bootstrap samples
#' n_boot = 50000
#' 
#' # Evaluate eclairs decomposition on boostrap samples
#' ecl2 = cov_transform( ecl, n_boot, f=f, lambda = 1e-4)
#' 
#' # Get full covariance matrix from eclairs decomposition
#' C1 = getCov(ecl2)
#' 
#' # Parametric boostrap samples directly from full covariance matrix 
#' X = rmvnorm(n_boot, rep(0, p), getCov(ecl))
#' 
#' # get covariance of transformed data
#' C2 = cov(f(X))
#' 
#' # Evaluate differences
#' # small differences are due to Monte Carlo error from boostrap sampling
#' range(C1-C2)
#' 
#' # Plot entries from two covariance estimates
#' par(pty="s")
#' plot(C1, C2, main="Concordance between covariances")
#' abline(0, 1, col='red')
#' 
#'  
#' # Save above but compute eclairs for correlation matrix
#' 
#' # Evaluate eclairs decomposition on boostrap samples
#' ecl2 = cov_transform( ecl, n_boot, f=f, compute="correlation", lambda = 1e-4)
#' 
#' # Get full covariance matrix from eclairs decomposition
#' C1 = getCor(ecl2)
#' 
#' # Parametric boostrap samples directly from full covariance matrix 
#' X = rmvnorm(n_boot, rep(0, p), getCov(ecl))
#' 
#' # get correlation of transformed data
#' C2 = cor(f(X))
#' 
#' # Evaluate differences
#' # small differences are due to Monte Carlo error from boostrap sampling
#' range(C1-C2)
#' 
#' # Plot entries from two correlation estimates
#' par(pty="s")
#' plot(C1, C2, main="Correlation between covariances")
#' abline(0, 1, col='red')
#' 
#' @export
cov_transform = function(Sigma.eclairs, n, f = function(x) 2*log(abs(x)+1e-4), lambda=NULL, compute=c("covariance", "correlation")){

	stopifnot(is(Sigma.eclairs, "eclairs"))
	stopifnot(is(f, "function"))
	compute = match.arg(compute)

	# f1 = function(x) 2*log(abs(x)+1e-4) 
	# f2 = function(x) log(x^2)

	mu = rep(0, Sigma.eclairs$p)

	# draw samples from the multivariate normal with covariance Sigma.eclairs
	X = rmvnorm_eclairs(n, mu, Sigma.eclairs)

	# Perofrm eclairs decomposition of sampled data after applying transform f
	# set scale=FALSE to get covariance
	eclairs( f(X), lambda=lambda, compute=compute, warmStart=Sigma.eclairs )
}













