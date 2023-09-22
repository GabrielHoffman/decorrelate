
#' Fit linear model after decorrelating 
#' 
#' Fit linear model after applying decorrelation projection to response and predictors. 
#'
#' @param formula an object of class 'formula' (or one that can be coerced to that class): a symbolic description of the model to be fitted.

#' @param data a matrix or data.frame containing the variables in the model
#' @param Sigma.eclairs estimate of covariance/correlation matrix from \link{eclairs} storing \eqn{U}, \eqn{d_1^2}, \eqn{\lambda} and \eqn{\nu}
#' @param subset same as for \link{lm}
#' @param weights same as for \link{lm}
#' @param na.action same as for \link{lm}
#' @param method same as for \link{lm}
#' @param model same as for \link{lm}
#' @param x same as for \link{lm}
#' @param y same as for \link{lm}
#' @param qr same as for \link{lm}
#' @param singular.ok same as for \link{lm} 
#' @param contrasts same as for \link{lm}
#' @param offset same as for \link{lm}
#' @param ... same as for \link{lm}
#'
#' @details This function fit a linear regression to the transformed response, and transformed design matrix.  Note that the design matrix, not just the data.frame of variables is transformed so that 1) factors are transformed and 2) the intercept term is transformed.  
#'
#' @return Object of class \code{lm} returned by function \link{lm}
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
#' Y = rmvnorm(n, rep(0, p), sigma=Sigma*5.1)
#' 
#' # eclairs decomposition
#' ecl = eclairs(Y)
#' 
#' # simulate covariates
#' data = data.frame(matrnorm(p,2))
#' colnames(data) = paste0('v', 1:2)
#' 
#' # simulate response
#' y = rnorm(p)
#' 
#' # fit linear model on transformed data
#' lm_eclairs(y ~ v1 + v2, data, ecl )
#'
#' @import stats
#' @export
lm_eclairs = function (formula, data, Sigma.eclairs, subset, weights, na.action, method = "qr", model = TRUE, x = FALSE, y = FALSE, qr = TRUE, singular.ok = TRUE, 
    contrasts = NULL, offset, ...){
    ret.x <- x
    ret.y <- y
    cl <- match.call()
    mf <- match.call(expand.dots = FALSE)
    m <- match(c("formula", "data", "subset", "weights", "na.action", 
        "offset"), names(mf), 0L)
    mf <- mf[c(1L, m)]
    mf$drop.unused.levels <- TRUE
    mf[[1L]] <- quote(stats::model.frame)
    mf <- eval(mf, parent.frame())
    if (method == "model.frame") 
        return(mf)
    else if (method != "qr") 
        warning(gettextf("method = '%s' is not supported. Using 'qr'", 
            method), domain = NA)
    mt <- attr(mf, "terms")

    # get response from above scope and then apply transformation
    y <- decorrelate(model.response(mf, "numeric"), Sigma.eclairs, transpose=TRUE)

    w <- as.vector(model.weights(mf))
    if (!is.null(w) && !is.numeric(w)) 
        stop("'weights' must be a numeric vector")
    offset <- model.offset(mf)
    mlm <- is.matrix(y)
    ny <- if (mlm) 
        nrow(y)
    else length(y)
    if (!is.null(offset)) {
        if (!mlm) 
            offset <- as.vector(offset)
        if (NROW(offset) != ny) 
            stop(gettextf("number of offsets is %d, should equal %d (number of observations)", 
                NROW(offset), ny), domain = NA)
    }
    if (is.empty.model(mt)) {
        x <- NULL
        z <- list(coefficients = if (mlm) matrix(NA_real_, 0, 
            ncol(y)) else numeric(), residuals = y, fitted.values = 0 * 
            y, weights = w, rank = 0L, df.residual = if (!is.null(w)) sum(w != 
            0) else ny)
        if (!is.null(offset)) {
            z$fitted.values <- offset
            z$residuals <- y - offset
        }
    }
    else {
        x <- model.matrix(mt, mf, contrasts)

        # get predictors and apply transformation
        x = decorrelate( x, Sigma.eclairs, transpose=TRUE)

        z <- if (is.null(w)) 
            lm.fit(x, y, offset = offset, singular.ok = singular.ok, 
                ...)
        else lm.wfit(x, y, w, offset = offset, singular.ok = singular.ok, 
            ...)
    }
    class(z) <- c(if (mlm) "mlm", "lm")
    z$na.action <- attr(mf, "na.action")
    z$offset <- offset
    z$contrasts <- attr(x, "contrasts")
    z$xlevels <- .getXlevels(mt, mf)
    z$call <- cl
    z$terms <- mt
    if (model) 
        z$model <- mf
    if (ret.x) 
        z$x <- x
    if (ret.y) 
        z$y <- y
    if (!qr) 
        z$qr <- NULL
    z
}





#' Fit linear model on each feature after decorrelating 
#' 
#' Fit linear model on each feature after applying decorrelation projection to response and predictors. 
#'
#' @param formula an object of class 'formula' (or one that can be coerced to that class): a symbolic description of the model to be fitted.

#' @param data a matrix or data.frame containing the variables in the model
#' @param X matrix or data.frame where each column stores a predictor to be evaluated by the regression model one at a time.  The \eqn{i^{th}} model includes \code{X[,i]} as a predictor.
#' @param Sigma.eclairs estimate of covariance/correlation matrix from \link{eclairs} storing \eqn{U}, \eqn{d_1^2}, \eqn{\lambda} and \eqn{\nu}
#' @param subset same as for \link{lm}
#' @param weights same as for \link{lm}
#' @param na.action same as for \link{lm}
#' @param method same as for \link{lm}
#' @param model same as for \link{lm}
#' @param x same as for \link{lm}
#' @param y same as for \link{lm}
#' @param qr same as for \link{lm}
#' @param singular.ok same as for \link{lm} 
#' @param contrasts same as for \link{lm}
#' @param offset same as for \link{lm}
#' @param ... same as for \link{lm}

#' @param ... other arguments passed to \code{lm()}
#'
#' @return data.frame with columns \code{beta}, \code{se}, \code{tsat}, \code{pvalue} storing results for regression model fit for each feature
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
#' Y = rmvnorm(n, rep(0, p), sigma=Sigma*5.1)
#' 
#' # eclairs decomposition
#' ecl = eclairs(Y)
#' 
#' # simulate covariates
#' data = data.frame(matrnorm(p,2))
#' colnames(data) = paste0('v', 1:2)
#' 
#' # simulate response
#' y = rnorm(p)
#' 
#' # Simulate 1000 features to test
#' X = matrnorm(p, 1000)
#' colnames(X) = paste0('set_', seq(ncol(X)))
#' 
#' # Use linear model to test each feature stored as columns in X
#' res = lm_each_eclairs(y ~ v1 + v2, data, X, ecl )
#' 
#' head(res)
#' 
#' # Analysis after non-linear transform
#' #------------------------------------
#' 
#' # Apply function to transforme data
#' f = function(x) log(x^2 + 0.001)
#' 
#' # evaluate covariance of transformed data
#' ecl_transform = cov_transform(ecl, f, 100)
#' 
#' # Use linear model to test each feature stored as columns in X
#' # in data transformed by f()
#' res2 = lm_each_eclairs( f(y) ~ v1 + v2, data, X, ecl_transform )
#' 
#' head(res)
#'
#' @import stats
#' @export
lm_each_eclairs = function(formula, data, X, Sigma.eclairs, subset, weights, na.action, method = "qr", model = TRUE, x = FALSE, y = FALSE, qr = TRUE, singular.ok = TRUE, 
    contrasts = NULL, offset, ...){

    # method = "qr"
    # singular.ok = TRUE
    # contrasts = NULL

    cl <- match.call()
    mf <- match.call(expand.dots = FALSE)
    m <- match(c("formula", "data", "subset", "weights", "na.action", 
        "offset"), names(mf), 0L)
    mf <- mf[c(1L, m)]
    mf$drop.unused.levels <- TRUE
    mf[[1L]] <- quote(stats::model.frame)
    mf <- eval(mf, parent.frame())
    if (method == "model.frame") 
        return(mf)
    else if (method != "qr") 
        warning(gettextf("method = '%s' is not supported. Using 'qr'", 
            method), domain = NA)
    mt <- attr(mf, "terms")

    # get response from above scope and then apply transformation
    y <- model.response(mf, "numeric")
    y.transform <- decorrelate(y, Sigma.eclairs, transpose=TRUE)

    # w <- as.vector(model.weights(mf))
    # if (!is.null(w) && !is.numeric(w)) 
    #     stop("'weights' must be a numeric vector")
    offset <- model.offset(mf)
    mlm <- is.matrix(y)
    ny <- if (mlm) 
        nrow(y)
    else length(y)
    if (!is.null(offset)) {
        if (!mlm) 
            offset <- as.vector(offset)
        if (NROW(offset) != ny) 
            stop(gettextf("number of offsets is %d, should equal %d (number of observations)", 
                NROW(offset), ny), domain = NA)
    }

    # get covariates
    x <- model.matrix(mt, mf, contrasts)

    # get predictors and apply transformation
    X.transform = decorrelate( x, Sigma.eclairs, transpose=TRUE)

    # transform matrix of features
    X.features.transform = decorrelate(X, Sigma.eclairs, transpose=TRUE)

    # Perform pre-projection
    obj = lm.projection(y.transform, X.transform)

    # for each feature
    res = t(apply(X.features.transform, 2, function(x_feat){

        # Perform hypothesis test on 1 additional feature using pre-projection
        lm.test( obj, x_feat)
    }))
    rownames(res) = colnames(X)
    colnames(res) = c("beta", "se", "tstat", "pvalue")

    as.data.frame(res)
}
   











