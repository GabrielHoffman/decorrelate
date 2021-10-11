

#' Evaluate effective variance
#'
#' Evaluate effective variance of the covariance/correlation matrix
#'
#' @param Sigma.eclairs estimate of covariance/correlation matrix from \code{eclairs} storing \eqn{U}, \eqn{d_1^2}, \eqn{\lambda} and \eqn{\nu}
#'
#' @references{
#'   \insertRef{pena2003descriptive}{decorrelate}
#' }
#'
#' @details Compute a metric of the amount of variation represented by a covariance/correlation matrix that is comparable across matrices of difference sizes. Proposed by Pena and Rodriguez (2003), the 'effective variance' is \eqn{|C|^\frac{1}{p}} where \eqn{C} is a correlation or covariance matrix between \eqn{p} variables. The effective variance is the mean of the log eigen-values.
#'
#' @examples
#' library(Rfast)
#' set.seed(1)
#' n = 2000 # number of samples
#' p = 1000 # number of features
#' rho = .9
#' 
#' # create correlation matrix
#' Sigma = matrix(rho, p, p)
#' diag(Sigma) = 1
#' 
#' # draw data from correlation matrix Sigma
#' Y = rmvnorm(n, rep(0, p), sigma=Sigma*5.1)
#' rownames(Y) = paste0("sample_", 1:n)
#' colnames(Y) = paste0("gene_", 1:p)
#' 
#' # eclairs decomposition
#' Sigma.eclairs = eclairs(Y, compute="cor")
#' 
#' # effective variance
#' effVariance( Sigma.eclairs )
#' 
#' @export
effVariance = function( Sigma.eclairs ){

	# Can be either a covariance or correlation matrix
	# det(C)^(1/p)
	# Compute in log space then exponentiate
	exp(logDet(Sigma.eclairs) / Sigma.eclairs$p)
}

# #' Evaluate effective dependence
# #'
# #' Evaluate effective dependence of the correlation matrix
# #'
# #' @param Sigma.eclairs estimate of correlation matrix from \code{eclairs} storing \eqn{U}, \eqn{d_1^2}, \eqn{\lambda} and \eqn{\nu}
# #'
# #' @details Compute a metric of the average degree of dependence between all pairs of variables. Proposed by Pena and Rodriguez (2003), the 'effective dependence' is \eqn{1-|C|^\frac{1}{p}} where \eqn{C} is a correlation matrix between \eqn{p} variables. In the case where all pairwise correlations are \eqn{\rho}, the effective dependence converges to \eqn{\rho} for sufficiently large \eqn{p}.  Like the squared correlation, the value is bounded between 0 and 1.
# #'
# #' @references{
# #'   \insertRef{pena2003descriptive}{decorrelate}
# #' }
# #
# #' @examples
# #' library(Rfast)
# #' set.seed(1)
# #' n = 2000 # number of samples
# #' p = 1000 # number of features
# #' rho = .9
# #' 
# #' # create correlation matrix
# #' Sigma = matrix(rho, p, p)
# #' diag(Sigma) = 1
# #' 
# #' # draw data from correlation matrix Sigma
# #' Y = rmvnorm(n, rep(0, p), sigma=Sigma*5.1)
# #' rownames(Y) = paste0("sample_", 1:n)
# #' colnames(Y) = paste0("gene_", 1:p)
# #' 
# #' # eclairs decomposition
# #' Sigma.eclairs = eclairs(Y, compute="cor")
# #' 
# #' # effective variance
# #' effVariance( Sigma.eclairs )
# #' 
# #' # effective dependence
# #' effDependence( Sigma.eclairs )
# #' 
# #' # expected value with constant correlation
# #' rho
# #' 
# #' @export
# effDependence = function( Sigma.eclairs ){

# 	if( Sigma.eclairs$nu != 1 ){
# 		stop("Can only compute effective dependence on a correlation matrix")
# 	}

# 	# where C is a correlation matrix
# 	# 1 - det(C)^(1/p)
# 	# Compute in log space then exponentiate
# 	1 - exp(logDet(Sigma.eclairs) / Sigma.eclairs$p)
# }


