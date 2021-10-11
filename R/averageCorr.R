
#' Compute average correlation of implied correlation matrix
#' 
#' Compute average correlation of implied correlation matrix
#' 
#' @param Sigma.eclairs estimate of correlation matrix from \code{eclairs} storing \eqn{U}, \eqn{d_1^2}, \eqn{lambda} and \eqn{nu}
#'
#' @details If correlation matrix \eqn{C} has all correlations set to \eqn{rho}, and \eqn{dSq} is the vector of eigen-values of C, then all but the first eigen-value equals \eqn{1-rho}.  When \eqn{C} is an sample correlation matrix, \eqn{rho} can be estimated as the mean of all but the first eigen-value.  Since each variable is scaled to have variance 1 to create the correlation matrix, the total variance is \eqn{p}.  Therefore the mean of the remaining eigen-values is \eqn{(p - dSq[1])/p}. 
#'
#' @examples
#' library(Rfast)
#' set.seed(1)
#' n = 200 # number of samples
#' p = 800 # number of features
#' 
#' # create correlation matrix
#' Sigma = matrix(.2, p, p)
#' diag(Sigma) = 1
#' 
#' # draw data from correlation matrix Sigma
#' Y = rmvnorm(n, rep(0, p), sigma=Sigma)
#' rownames(Y) = paste0("sample_", 1:n)
#' colnames(Y) = paste0("gene_", 1:p)
#' 
#' # eclairs decomposition
#' Sigma.eclairs = eclairs(Y, compute="cor")
#' 
#' averageCorr( Sigma.eclairs )
#' 
#' @export
averageCorr = function(Sigma.eclairs){

    if (Sigma.eclairs$nu != 1) {
        stop("Can only compute effective dependence on a correlation matrix")
    }

    with(Sigma.eclairs, 1 - (p - dSq[1])/p)
}