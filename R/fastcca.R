# Sept 13, 2021



#' Class fastcca
#'
#' Class \code{fastcca} 
#'
#' @details Object storing:
#' \itemize{
#'  \item{n.comp: }{number of canonical components}
#'  \item{cors: }{canonical correlations}
#'  \item{x.coefs: }{canonical coefficients for X}
#'  \item{x.vars: }{canonical variates for X}
#'  \item{y.coefs: }{canonical coefficients for Y}
#'  \item{y.vars: }{canonical variates for Y}
#'  \item{lambdas: }{shrinkage parameters from \code{eclairs}}
#' }
#' @name fastcca-class
#' @rdname fastcca-class
#' @exportClass fastcca
setClass("fastcca", contains= "list")


setMethod("show", 'fastcca',
      function(object){
	      print(object)
})


setMethod("print", 'fastcca',
      function(x){

	cat("       Fast regularized canonical correlation analysis\n\n")

	k = min(3, x$n.comp)

	cat("  Original data rows:", x$dims['n1'],"\n")
	cat("  Original data cols: ", x$dims['p1'], ', ',x$dims['p2'],"\n", sep='')
	cat("  Num components:    ", x$n.comp, "\n")
	cat("  Cor:               ", round(x$cor[seq_len(k)], digits=3), "...\n")
	cat("  rho.mod:           ", round(x$rho.mod[seq_len(k)], digits=3), "...\n")
	cat("  Cramer's V:        ", round(x$cramer.V, digits=3), "\n")
	cat("  lambda.x:          ", format(x$lambdas['x'], digits=3), "\n")
	cat("  lambda.y:          ", format(x$lambdas['y'], digits=3), "\n")
})




#' Fast canonical correlation analysis 
#' 
#' Fast Canonical correlation analysis that is scalable to high dimensional data.  Uses covariance shrinkage and algorithmic speed ups to be linear time in p when p > n.
#' 
#' @param X first matrix (n x p1)
#' @param Y first matrix (n x p2)
#' @param k number of canonical components to return
#' @param lambda.x optional shrinkage parameter for estimating covariance of X. If NULL, estimate from data.
#' @param lambda.y optional shrinkage parameter for estimating covariance of Y. If NULL, estimate from data.
#'
#' @details
#' Results from standard CCA are based on the SVD of \eqn{\Sigma_{xx}^{-\frac{1}{2}} \Sigma_{xy} \Sigma_{yy}^{-\frac{1}{2}}}.
#'
#' Uses \code{eclairs()} and empirical Bayes covariance regularization, and applies speed up of RCCA \insertCite{tuzhilina2020canonical;textual}{decorrelate} to perform CCA on n PCs and instead of p features.  Memory usage is \eqn{\mathcal{O}(np)} instead of \eqn{\mathcal{O}(p^2)}.  Computation is \eqn{\mathcal{O}(n^2p)} instead of \eqn{\mathcal{O}(p^3)} or \eqn{\mathcal{O}(np^2)}
#Computation is n^2p instead of p^3 of np^2
#'
#' @references \insertAllCited{}
#' 
#' @examples
#' pop <- LifeCycleSavings[, 2:3]
#' oec <- LifeCycleSavings[, -(2:3)]
#' fastcca(pop, oec)
#' 
#' @importFrom Rfast standardise
#' @export
fastcca = function(X, Y, k=min(dim(X), dim(Y)), lambda.x = NULL, lambda.y = NULL){

	if( ! is.matrix(X) ) X = as.matrix(X)
    if( ! is.matrix(Y) ) Y = as.matrix(Y)

    if( anyNA(X) | anyNA(Y) ){
    	stop("No NA values are allowed in data")
    }

    # mean-center columns
    # X = standardise(X, scale = FALSE)
    # Y = standardise(Y, scale = FALSE)
    X = scale(X, scale = FALSE)
    Y = scale(Y, scale = FALSE)

	n1 = nrow(X)
	n2 = nrow(Y)
	p1 = ncol(X)
	p2 = ncol(Y)

	if( n1 != n2 ){
		stop("X and Y must have same number of rows")
	}

	k = min(dim(X), dim(Y), k)

	if( k <= 0 ){
		stop("k must be >= 1")
	}

	# if( (p1 < n1) & (p2 < n2) ){
	# 	# Low dimensional
	# 	res = cca(X, Y, lambda.x, lambda.y)
	# 	return( res )
	# }

	# transform to lower dimensions
	X.tr = transform(X, lambda.x)
	Y.tr = transform(Y, lambda.y)

	n.comp = min(ncol(X), ncol(Y), nrow(X), k)

	# eigen decomp of Cxx^-.5 Cxy Cyy^-.5
	# geigen can be made faster if using low rank eclairs
	Cxx = X.tr$cor
	Cyy = Y.tr$cor
	# Cxy = cov(X.tr$mat, Y.tr$mat, use = "pairwise")
	# get covariance as product of centered matricies
	Cxy = crossprod(X.tr$mat, Y.tr$mat) / (n1-1) 

	# scale by shrinkage parameters
	Cxy = Cxy * sqrt(1-X.tr$lambda) * sqrt(1-Y.tr$lambda)

	sol = geigen2(Cxy, Cxx, Cyy, k) # compute SVD of matrix products
	names(sol) = c("rho", "alpha", "beta")
	
	# modified canonical correlation
	rho.mod = sol$rho[seq_len(n.comp)]
	names(rho.mod) = paste('can.comp', seq_len(n.comp), sep = '')

	# inverse transform
	X.inv.tr = inverseTransform(X.tr$mat, sol$alpha, X.tr$tr)

	x.coefs = X.inv.tr$coefs[,seq_len(n.comp),drop=FALSE]
	rownames(x.coefs) = colnames(X)

	x.vars = X.inv.tr$vars[,seq_len(n.comp),drop=FALSE]
	rownames(x.vars) = rownames(X)

	Y.inv.tr = inverseTransform(Y.tr$mat, sol$beta, Y.tr$tr)

	y.coefs = Y.inv.tr$coefs[,seq_len(n.comp),drop=FALSE]
	rownames(y.coefs) = colnames(Y)

	y.vars = Y.inv.tr$vars[,seq_len(n.comp),drop=FALSE]
	rownames(y.vars) = rownames(Y)

    rho = diag(cor(X.inv.tr$vars, Y.inv.tr$vars))[seq(n.comp)]
    names(rho) = paste("can.comp", seq(n.comp), sep = "")

    if(k < min( dim(X), dim(Y)) ){
		idx = seq(1, k)
	}else{
		idx = seq(1, k-1)
	}

	cramer.V = sqrt(mean(rho.mod[idx]^2)) # Cramer's V-statistic for CCA

	# based on CRAN::yacca
	# compute loadings
	load.x = cor(X, x.vars)
	load.y = cor(Y, y.vars)

	# compute redundancy index
	ri.x = rho.mod^2 * apply(load.x, 2, function(x) mean(x^2))
	ri.y = rho.mod^2 * apply(load.y, 2, function(x) mean(x^2))

	names(ri.x) = paste('can.comp', seq_len(n.comp), sep = '') 
	names(ri.y) = paste('can.comp', seq_len(n.comp), sep = '')

	res = list(	n.comp 	= n.comp, 
				rho.mod = rho.mod, 
				cor 	= rho,
				cramer.V= cramer.V,
				x.coefs = x.coefs,
				x.vars 	= x.vars, 
				x.ri 	= ri.x,
				y.coefs = y.coefs, 
				y.vars 	= y.vars, 
				y.ri 	= ri.y,
				lambdas = c(x = X.tr$lambda, y = Y.tr$lambda),
				dims 	= c(n1=n1, n2=n2, p1=p1, p2=p2))

	new('fastcca', res)
}


# Adapted from RCCA:::RCCA.tr
transform = function(X, lambda=NULL){

	# get shrinkage estimate to 1) compute SVD and 2) estimate lambda
	ecl = eclairs(X, lambda=lambda, compute="corr")

	# if more features that samples
	if( ncol(X) > nrow(X) ){

	    # scale to be consistent with PCA
		mat = with(ecl, V %*% diag(sqrt(dSq))) * sqrt(nrow(X)-1)
		V = ecl$U
	}else{
		V = NULL
		mat = X
	}

	# get shrinkage estimate of covariance matrix
	# C = with(ecl, (1-lambda)*var(mat, use = "pairwise") + diag(lambda*nu, ncol(mat)))
	C = with(ecl, (1-lambda)*crossprod(mat)/(nrow(mat)-1) + diag(lambda*nu, ncol(mat)))

	list(mat = mat, tr = V, cor=C, lambda = ecl$lambda, Sigma.eclairs=ecl)
}





# Adapted from RCCA:::RCCA.inv.tr
inverseTransform = function(X, alpha, V){

	n.comp = ncol(alpha)

	# find canonical variates
	u = X %*% alpha
	colnames(u) = paste('can.comp', seq_len(n.comp), sep = '')

	# inverse transfrom canonical coefficients 
	if(!is.null(V)){
		alpha = V %*% alpha
	}

	colnames(alpha) = paste('can.comp', seq_len(n.comp), sep = '')

	list(coefs = alpha, vars = u)
}


geigen = function(Amat, Bmat, Cmat, k){
    Bdim <- dim(Bmat)
    Cdim <- dim(Cmat)
    if (Bdim[1] != Bdim[2]) 
        stop("BMAT is not square")
    if (Cdim[1] != Cdim[2]) 
        stop("CMAT is not square")
    p <- Bdim[1]
    q <- Cdim[1]
    s <- min(c(p, q))
    if (max(abs(Bmat - t(Bmat)))/max(abs(Bmat)) > 1e-10) 
        stop("BMAT not symmetric.")
    if (max(abs(Cmat - t(Cmat)))/max(abs(Cmat)) > 1e-10) 
        stop("CMAT not symmetric.")
    Bmat <- (Bmat + t(Bmat))/2
    Cmat <- (Cmat + t(Cmat))/2
    Bfac <- chol(Bmat)
    Cfac <- chol(Cmat)
    Bfacinv <- solve(Bfac)
    Cfacinv <- solve(Cfac)
    Dmat <- crossprod(Bfacinv, Amat) %*% Cfacinv
    if (p >= q) {

    	if( k < min(p,q) / 3){
			result = irlba(Dmat, k) # should be faster thatn PRIMME::svds
        }else{
	        result <- svd(Dmat)
	    }

        values <- result$d
        Lmat <- Bfacinv %*% result$u
        Mmat <- Cfacinv %*% result$v
    }
    else {
        result <- svd(t(Dmat))
        values <- result$d
        Lmat <- Bfacinv %*% result$v
        Mmat <- Cfacinv %*% result$u
    }
    geigenlist <- list(values, Lmat, Mmat)
    names(geigenlist) <- c("values", "Lmat", "Mmat")
    return(geigenlist)
}


geigen2 = function(Amat, Bmat, Cmat, k){

	minvsqrt = function(S){
	  dcmp = eigen(S)
	  with(dcmp, vectors %*% diag(1/sqrt(values)) %*% t(vectors))
	}

	D = minvsqrt(Bmat) %*% Amat %*% minvsqrt(Cmat)
	
	p = nrow(Bmat)
	q = nrow(Cmat)

	if( q > p){
		D = t(D)
	}

	if( k > min(dim(D)) / 3 ){
		dcmp = svd(D)

		if( k < min(dim(D)) ){
			dcmp$d = dcmp$d[seq_len(k)]
			dcmp$u = dcmp$u[,seq_len(k), drop=FALSE]
			dcmp$v = dcmp$v[,seq_len(k), drop=FALSE]
		}
	}else{
		dcmp = irlba(D, nv=k)
	}

	# set diagonals to be positive
	dcmp$v = eachrow(dcmp$v, sign(diag(dcmp$v)), "*")
    dcmp$u = eachrow(dcmp$u, sign(diag(dcmp$u)), "*")
	
    values <- dcmp$d

    if( p >=q ){
	    Lmat <- minvsqrt(Bmat) %*% dcmp$u
	    Mmat <- minvsqrt(Cmat) %*% dcmp$v
	}else{		
	    Lmat <- minvsqrt(Bmat) %*% dcmp$v
	    Mmat <- minvsqrt(Cmat) %*% dcmp$u
	}

    geigenlist <- list(values, Lmat, Mmat)
    names(geigenlist) <- c("values", "Lmat", "Mmat")

    return(geigenlist)
}


