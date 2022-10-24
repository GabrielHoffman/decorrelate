
#' Compute average correlation of implied correlation matrix
#' 
#' Compute average correlation of implied correlation matrix
#' 
#' @param Sigma.eclairs estimate of correlation matrix from \code{eclairs()} storing \eqn{U}, \eqn{d_1^2}, \eqn{\lambda} and \eqn{\nu}
#' @param method compute average correlation for either the MLE correlation matrix, or the empirical Bayes (EB) correlation matrix 
#' @param algorithm \code{"closedApprox"} uses a closed form approximation that works well for large \eqn{n}. \code{"lossFxn"} uses optimizes a quadratic loss function between observed eigen-values and evalues of constant correlation matrix.
#'  
#' @details Consider a correlation matrix \eqn{C} with all correlations set to \eqn{\rho}, where \eqn{d^2} is the vector of eigen-values of \eqn{C}.  Then the first eigen-value is \eqn{1+(p-1)\rho}, and the next \eqn{p-1} eigen-values are \eqn{1-\rho}.  (See \url{https://statisticaloddsandends.wordpress.com/2018/07/20/} for clear explanation.)  Now consider \eqn{C} as a sample correlation matrix, and estimate the average correlation based on the eigen-values.  In this case, all additional eigen-values are zero when \eqn{p > n}.  
#' 
#' Consider estimating \eqn{\rho} using each eigen-value separately. The first eigen-value gives the estimate \eqn{(d^2_1 - 1) / (p-1)}, the next \eqn{p-1} eigen-values give the estimate \eqn{1 - d^2_i}, and all remaining estimates are 0.  Finally, take the mean of these single estimates for the final estimate of \eqn{\rho}.  This gives the average correlation for the MLE correlation matrix.
#' 
#' Since the empirical Bayes estimate is a convex combination of the MLE and the identity matrix, the average correlation of the EB correlation matrix is \eqn{(1-\lambda)\rho}
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
averageCorr = function(Sigma.eclairs, method = c("MLE", "EB"), algorithm = c("closedApprox", "lossFxn")){

    method = match.arg(method)
    algorithm = match.arg(algorithm)

    if (Sigma.eclairs$nu != 1) {
        stop("Can only compute effective dependence on a correlation matrix")
    }

    if( algorithm == "closedApprox"){
        message(algorithm)
        # this is approximate for small n and large p, but is very accurate
        # i.e. 2 decimal places

        # MLE
        rho1 = with(Sigma.eclairs, (dSq[1] - 1) / (p-1))
        rho_else_sum = with(Sigma.eclairs, sum(1 - dSq[-1]) + max(p-n, 0))
        rho_hat = (rho1 + rho_else_sum) / Sigma.eclairs$p

        # EB
        if( method == "EB"){
            rho_hat = (1-Sigma.eclairs$lambda) * rho_hat
        }
    }else{
        message(algorithm)
        # estimate rho using loss function
        p = Sigma.eclairs$p
        dSq = Sigma.eclairs$dSq

        # if high dimensional data, pad with zeros
        if( with(Sigma.eclairs, p > n) ){
            dSq = c(dSq, rep(0, with(Sigma.eclairs, p -n)))
        }

        # if EB, scale by lambda
        if( method == "EB" ){
            dSq = dSq*(1-Sigma.eclairs$lambda)
        }

        # define quadratic loss function between
        # observed eigen-values and evalues of constant
        # correlation matrix
        loss = function(rho){

            # eigen-value given rho
            value = 1 + (p-1)*rho
            value = c(value, rep(1-rho, p-1))

            sum((value - dSq)^2)
        }

        # minimize loss function
        fit = optimize(loss, interval=c(0 + 1e-4, 1-1e-4))

        rho_hat = fit$minimum
    }

    rho_hat
}