

# lambda values are not affected by nu
# decorrelate:::shrinkcovmat.equal_lambda( t(Y) )

# decorrelate:::shrinkcovmat.equal_lambda( t(Y)*10 )

# Estimate lambda and nu from the dataset
# 
# Estimate lambda and nu from the dataset using code adapted from ShrinkCovMat (https://github.com/AnestisTouloumis/ShrinkCovMat/blob/master/R/shrinkcovmat.identity.R)
#
#
# @reference
# Touloumis, A. 2015. Nonparametric Stein-type shrinkage covariance matrix estimators in high-dimensional settings. Computational Statistics & Data Analysis. doi:10.1016/j.csda.2014.10.018
#
#' @importFrom stats var cov
#' @import ShrinkCovMat
shrinkcovmat.equal_lambda = function (data, centered = FALSE){
    if (!is.matrix(data))
        data <- as.matrix(data)
    p <- nrow(data)
    n <- ncol(data)
    centered <- as.logical(centered)
    if (centered != TRUE && centered != FALSE) {
        stop("'centered' must be either 'TRUE' or 'FALSE'")
    }
    if (!centered) {
        if (n < 4)
            stop("The number of columns should be greater than 3")
        # sample_covariance_matrix <- cov(t(data))
        trace_statistics <- ShrinkCovMat:::trace_stats_uncentered(data)
        trace_sigma_hat <- trace_statistics[1]
        nu_hat <- trace_sigma_hat/p
        trace_sigma_squared_hat <- trace_statistics[2]
        lambda_hat <- (trace_sigma_hat^2 + trace_sigma_squared_hat)/(n *
            trace_sigma_squared_hat + (p - n + 1)/p * trace_sigma_hat^2)
        lambda_hat <- min(lambda_hat, 1)
    }
    else {
        if (n < 2)
            stop("The number of columns should be greater than 1")
        # sample_covariance_matrix <- tcrossprod(data)/n
        trace_statistics <- ShrinkCovMat:::trace_stats_centered(data)
        trace_sigma_hat <- trace_statistics[1]
        nu_hat <- trace_sigma_hat/p
        trace_sigma_squared_hat <- trace_statistics[2]
        lambda_hat <- (trace_sigma_hat^2 + trace_sigma_squared_hat)/((n +
            1) * trace_sigma_squared_hat + (p - n)/p * trace_sigma_hat^2)
        lambda_hat <- min(lambda_hat, 1)
    }
    data.frame(lambda_hat = lambda_hat, nu_hat = nu_hat)
    # if (lambda_hat < 1) {
    #     sigmahat <- (1 - lambda_hat) * sample_covariance_matrix +
    #         diag(nu_hat * lambda_hat, p)
    # }
    # else {
    #     sigmahat <- diag(lambda_hat * nu_hat, p)
    # }
    # target <- diag(nu_hat, p)
    # ans <- list(Sigmahat = sigmahat, lambdahat = lambda_hat,
    #     Sigmasample = sample_covariance_matrix, Target = target,
    #     centered = centered)
    # class(ans) <- "shrinkcovmathat"
    # ans
}


# Estimate lambda from eigen-values
#
# Estimate \eqn{lambda] from the squared singular values of the data matrix.  This is much faster than \code{shrinkcovmat.equal_lambda}, but less precise.  This uses the equation for \eqn{\lambda^star} in Section 3 of Touloumis (2015).
#
# @param ev squared singular values of the data matrix
# @param n number of samples
# @param p number of features
#
# @references Touloumis, A. 2015. Nonparametric Stein-type shrinkage covariance matrix estimators in hgh-dimensional settings. Computational Statistics & Data Analysis. doi:10.1016/j.csda.2014.10.018
#
est_lambda_ev = function(ev, n, p){
    (sum(ev^2) + sum(ev)^2)/ (n*sum(ev^2) + (p-n+1)/p * sum(ev)^2)
}








