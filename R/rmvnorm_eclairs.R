# Gabriel Hoffman
# April 7, 2021
#
# Linear time drawing from multivariate normal and t distributions


#' Draw from multivariate normal and t distributions
#'
#' Draw from multivariate normal and t distributions using eclairs decomposition
#'
#' @param n sample size
#' @param mu mean vector
#' @param Sigma.eclairs covariance matrix as an \link{eclairs} object
#' @param v degrees of freedom, defaults to Inf.  If finite, uses a multivariate t distribution
#'
#' @details 
#' Draw from multivariate normal and t distributions using eclairs decomposition.  If the (implied) covariance matrix is \eqn{p \times p}, the standard approach is \eqn{O(p^3)}. Taking advantage of the previously computed eclairs decomposition of rank \eqn{k}, this can be done in \eqn{O(pk^2)}.
#'
#' @return matrix where rows are samples from multivariate normal or t distribution where columns have covariance specified by \code{Sigma.eclairs}
#'
#' @examples
#' library(Rfast)
#' set.seed(1)
#' n = 800 # number of samples
#' p = 200 # number of features
#' 
#' # create correlation matrix
#' Sigma = autocorr.mat(p, .9)
#' 
#' # draw data from correlation matrix Sigma
#' Y = rmvnorm(n, rep(0, p), sigma=Sigma)
#' 
#' # perform eclairs decomposition
#' ecl = eclairs(Y)
#' 
#' # draw from multivariate normal
#' n = 100000
#' mu = rep(0, ncol(Y))
#' 
#' # using eclairs decomposition
#' X.draw1 = rmvnorm_eclairs(n, mu, ecl)
#' 
#' # using full covariance matrix implied by eclairs model
#' X.draw2 = rmvnorm(n, mu, getCov(ecl))
#' 
#' # assess difference betweeen covariances from two methods
#' range( cov(X.draw1) - cov(X.draw2))
#' 
#' # compare covariance to the covariance matrix used to simulated the data
#' range( cov(X.draw1) - getCov(ecl))
#'
#' @importFrom stats rchisq
#' @importFrom methods is
#' @importFrom Rfast matrnorm
#' @export
rmvnorm_eclairs = function(n, mu, Sigma.eclairs, v=Inf){

	stopifnot(is(Sigma.eclairs, "eclairs"))

	p = length(mu)

	# simulate from standard normal
	X = matrnorm(n,p)

	# get full covariance matrix from eclairs
	# sigma = getCov(Sigma.eclairs)

	# transform by covariance and add mean	
	# res1 =  X %*% chol(sigma) + rep(mu, rep(n, p))

	# alpha = 1/2 induces correlation
	# res2 = mult_eclairs(X, Sigma.eclairs$U, Sigma.eclairs$dSq, Sigma.eclairs$lambda, alpha = 1/2) + rep(mu, rep(n, p))


	# cov(res1)[1:3, 1:3]
	# cov(res2)[1:3, 1:3]

	# hist(cov(res1) - cov(res2))
	# hist(cov(res2) - sigma)

	# Create X %*% sqrt of Sigma
	# alpha = 1/2 induces correlation, -1/2 removes correlation
	X_transform = mult_eclairs(X, Sigma.eclairs$U, Sigma.eclairs$dSq, Sigma.eclairs$lambda, Sigma.eclairs$nu, alpha = 1/2)

	if( is.infinite(v) ){

		# Multivariate Normal
		X_values = X_transform + rep(mu, rep(n, p))

	}else{
		# Multivariate t 

		# create multivariate t from inverse chi-square weights, X_transform and mean
		w <- sqrt(v/rchisq(n, v))
		X_values = w * X_transform + rep(mu, rep(n, p))
	}

	# return matrix
	X_values
}
