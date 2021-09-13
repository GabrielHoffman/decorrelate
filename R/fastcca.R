# Sept 13, 2021




setMethod("t", 'eclairs',
function(x){
	ecl.t = x

	ecl.t$U = x$V
	ecl.t$V = x$U

	ecl.t$n = x$p
	ecl.t$p = x$n

	ecl.t$rownames = x$colnames
	ecl.t$colnames = x$rownames

	ecl.t
})




#' Fast canonical correlation analysis 
#' 
#' Fast Canonical correlation analysis that is scalable to high dimensional data.  Uses covariance shrinkage and algorithmic speed ups to be linear time in p when p > n.
#' 
#' @param X first matrix (n x p1)
#' @param Y first matrix (n x p2)
#' @param lambda.x optional shrinkage parameter for estimating covariance of X. If NULL, estimate from data.
#' @param lambda.y optional shrinkage parameter for estimating covariance of Y. If NULL, estimate from data.
#' @param k number of canonical components to return
#'
#' @details
#' Results from standard CCA are based on the SVD of \eqn{Sig_{xx}^{-.5} Sig_{xy} Sig_{yy}^{-.5}}.
#'
#' Uses eclairs and EB cov regularization, uses speed up of RCCA (Tuzhilina, et al, 2020) to perform CCA on n PCs and instead of p features.  Memory usage is n*p instead of p*p.  Computation is n^2p instead of p^3 of np^2
#'
#' @references{
#'   \insertRef{tuzhilina2020canonical}{decorrelate}
#' }
#' @importFrom fda geigen
#' @export
fastcca = function(X, Y, lambda.x = NULL, lambda.y = NULL){


	n1 = nrow(X)
	n2 = nrow(Y)
	p1 = ncol(X)
	p2 = ncol(Y)

	if( n1 != n2 ){
		stop("X and Y must have same number of rows")
	}

	if( (p1 < n1) & (p2 < n2) ){
		# Low dimensional
		res = cca(X, Y, lambda.x, lambda.y)
		return( res )
	}

	# transform to lower dimensions
	X.tr = transform(X, lambda.x)
	Y.tr = transform(Y, lambda.y)

	n.comp = min(ncol(X), ncol(Y), nrow(X))

	if( 1){
		# eigen decomp of Cxx^-.5 Cxy Cyy^-.5
		# geigen can be made faster if using low rank eclairs
		Cxx = X.tr$cor
		Cyy = Y.tr$cor
		Cxy = cov(X.tr$mat, Y.tr$mat, use = "pairwise")
		sol = fda::geigen(Cxy, Cxx, Cyy)
		names(sol) = c("rho", "alpha", "beta")
	}else{

		# full decomposition
		S = t(solve(chol(X.tr$cor))) %*% cov(X.tr$mat,Y.tr$mat) %*% solve(chol(Y.tr$cor))
		S[1:3, 1:3]


		# n = nrow(X.tr$mat)
		# C = t(solve(chol(X.tr$cor))) %*%(crossprod(X.tr$mat,Y.tr$mat)/(n-1)) %*% solve(chol(Y.tr$cor))
		# C[1:4, 1:4]

		# a = X.tr$mat %*% solve(chol(X.tr$cor))
		# b = Y.tr$mat %*% solve(chol(Y.tr$cor))

		# S = crossprod(a,b) / (n-1)

		# dcmp = svd(S)
		# dcmp$d[1:4]

		# S[1:3, 1:3]

		# This will allow low rank approximation using irlba instead of SVD.

		n = nrow(X.tr$mat)
		X1 = decorrelate(X.tr$mat, t(X.tr$Sigma.eclairs), transpose=TRUE)
		Y1 = decorrelate(Y.tr$mat, t(Y.tr$Sigma.eclairs), transpose=TRUE)

		S = crossprod(X1,Y1)  / (n-1)
		S[1:3, 1:3]

		dcmp = svd(S)
		dcmp$d[1:4]

		k_svd = 20

		values <- dcmp$d
	    Lmat <- decorrelate(dcmp$u[,seq_len(k_svd)], t(X.tr$Sigma.eclairs), transpose=TRUE)
	    Mmat <- decorrelate(dcmp$v[,seq_len(k_svd)], t(Y.tr$Sigma.eclairs), transpose=TRUE)

	    sol = list(values = dcmp$d, Lmat = Lmat, Mmat = Mmat)
	    names(sol) = c("rho", "alpha", "beta")
	 }


    # sol$alpha[1:3, 1:3]
    # Lmat[1:3, 1:3]

    # sol$beta[1:3, 1:3]
    # Mmat[1:3, 1:3]

    # diag(can(sol$Lmat, Lmat)$cor)
    # diag(cor(sol$Mmat, Mmat))


    # A = solve(chol(X.tr$cor)) %*% dcmp$u
    # A[1:3, 1:3]

	# modified canonical correlation
	rho.mod = sol$rho[seq_len(n.comp)]
	names(rho.mod) = paste('can.comp', seq_len(n.comp), sep = '')

	# inverse transform
	X.inv.tr = RCCA:::RCCA.inv.tr(X.tr$mat, sol$alpha, X.tr$tr)

	x.coefs = as.matrix(X.inv.tr$coefs[,seq_len(n.comp),drop=FALSE])
	rownames(x.coefs) = colnames(X)

	x.vars = X.inv.tr$vars[,1:n.comp]
	rownames(x.vars) = rownames(X)

	Y.inv.tr = RCCA:::RCCA.inv.tr(Y.tr$mat, sol$beta, Y.tr$tr)

	y.coefs = as.matrix(Y.inv.tr$coefs[,seq_len(n.comp),drop=FALSE])
	rownames(y.coefs) = colnames(Y)

	y.vars = Y.inv.tr$vars[,seq_len(n.comp),drop=FALSE]
	rownames(y.vars) = rownames(Y)

	# canonical correlation
	# rho = diag(cor(X.inv.tr$vars,  Y.inv.tr$vars))[1:n.comp]
	rho = sapply(seq_len(n.comp), function(i) cor(X.inv.tr$vars[,i],  Y.inv.tr$vars[,i]))
	names(rho) = paste('can.comp', seq_len(n.comp), sep = '')

	list('n.comp' = n.comp, 'cors' = rho, 'mod.cors' = rho.mod, 'x.coefs' = x.coefs, 'x.vars' = x.vars, 'y.coefs' = y.coefs, 'y.vars' = y.vars, lambdas = c(x = X.tr$lambda, y = Y.tr$lambda))
}

transform = function(X, k=ncol(X), lambda=NULL){

	# get shrinkage estimate to 1) compute SVD and 2) estimate lambda
	ecl = eclairs(X, k = k, lambda=lambda, compute="corr")

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
	C = with(ecl, (1-lambda)*var(mat, use = "pairwise") + diag(lambda*nu, ncol(mat)))

	list(mat = mat, tr = V, cor=C, lambda = ecl$lambda, Sigma.eclairs=ecl)
}




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
