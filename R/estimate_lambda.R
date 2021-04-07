
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

