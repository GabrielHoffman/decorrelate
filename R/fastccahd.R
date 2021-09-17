# # Sept 17, 2021

# # Try to use RCCA algorithm with decorrelate
# # But accuracy does down if rank of SVD is decreased
# # something is wrong with the code here

# #' Fast canonical correlation analysis 
# #' 
# #' Fast Canonical correlation analysis that is scalable to high dimensional data.  Uses covariance shrinkage and algorithmic speed ups to be linear time in p when p > n.
# #' 
# #' @param X first matrix (n x p1)
# #' @param Y first matrix (n x p2)
# #' @param k number of canonical components to return
# #' @param lambda.x optional shrinkage parameter for estimating covariance of X. If NULL, estimate from data.
# #' @param lambda.y optional shrinkage parameter for estimating covariance of Y. If NULL, estimate from data.
# #'
# #' @details
# #' Results from standard CCA are based on the SVD of \eqn{Sig_{xx}^{-.5} Sig_{xy} Sig_{yy}^{-.5}}.
# #'
# #' Uses eclairs and EB cov regularization, uses speed up of RCCA (Tuzhilina, et al, 2020) to perform CCA on n PCs and instead of p features.  Memory usage is n*p instead of p*p.  Computation is n^2p instead of p^3 of np^2
# #'
# #' @references{
# #'   \insertRef{tuzhilina2020canonical}{decorrelate}
# #' }
# #' 
# #' @examples
# #' pop <- LifeCycleSavings[, 2:3]
# #' oec <- LifeCycleSavings[, -(2:3)]
# #' fastcca(pop, oec)
# #' 
# #' @importFrom Rfast standardise
# #' @export
# fastccahd = function(X, Y, k=min(dim(X), dim(Y)), k.svd=min(dim(X), dim(Y)), lambda.x = NULL, lambda.y = NULL, std=FALSE){

# 	if( ! is.matrix(X) ) X = as.matrix(X)
#     if( ! is.matrix(Y) ) Y = as.matrix(Y)

#     # mean-center columns
#     # X = standardise(X, scale = FALSE)
#     # Y = standardise(Y, scale = FALSE)
#     X = scale(X, scale = FALSE)
#     Y = scale(Y, scale = FALSE)

# 	n1 = nrow(X)
# 	n2 = nrow(Y)
# 	p1 = ncol(X)
# 	p2 = ncol(Y)

# 	if( n1 < p1 | n2 < p2 ){
# 		stop("Only works with n >> p")
# 	}

# 	if( n1 != n2 ){
# 		stop("X and Y must have same number of rows")
# 	}

# 	k = min(dim(X), dim(Y), k)

# 	if( k < 0 ){
# 		stop("k must be > 1")
# 	}
	
# 	n.comp = min(ncol(X), ncol(Y), nrow(X), k)

# 	if( std ){
# 		X.tr = transform(X, lambda.x)
# 		Y.tr = transform(Y, lambda.y)
# 		# eigen decomp of Cxx^-.5 Cxy Cyy^-.5
# 		# geigen can be made faster if using low rank eclairs
# 		Cxx = X.tr$cor
# 		Cyy = Y.tr$cor
# 		# Cxy = cov(X.tr$mat, Y.tr$mat, use = "pairwise")
# 		# get covariance as product of centered matricies
# 		Cxy = crossprod(X.tr$mat, Y.tr$mat) / (n1-1) 

# 		sol = geigen2(Cxy, Cxx, Cyy, k) # compute SVD of matrix products
# 		names(sol) = c("rho", "alpha", "beta")

# 	}else{

# 		X.tr = transformHD(X, lambda.x, k=k.svd)
# 		Y.tr = transformHD(Y, lambda.y, k=k.svd)

# 		sol = geigen3( X.tr, Y.tr, k=k)
# 		names(sol) = c("rho", "alpha", "beta")
# 	}


# 	# print( sol$beta[1:3])

# 	# modified canonical correlation
# 	rho.mod = sol$rho[seq_len(n.comp)]
# 	names(rho.mod) = paste('can.comp', seq_len(n.comp), sep = '')

# 	# inverse transform
# 	X.inv.tr = inverseTransform(X.tr$mat, sol$alpha, X.tr$tr)

# 	x.coefs = X.inv.tr$coefs[,seq_len(n.comp),drop=FALSE]
# 	rownames(x.coefs) = colnames(X)

# 	x.vars = X.inv.tr$vars[,1:n.comp]
# 	rownames(x.vars) = rownames(X)

# 	Y.inv.tr = inverseTransform(Y.tr$mat, sol$beta, Y.tr$tr)

# 	y.coefs = Y.inv.tr$coefs[,seq_len(n.comp),drop=FALSE]
# 	rownames(y.coefs) = colnames(Y)

# 	y.vars = Y.inv.tr$vars[,seq_len(n.comp),drop=FALSE]
# 	rownames(y.vars) = rownames(Y)

#     rho = diag(cor(X.inv.tr$vars, Y.inv.tr$vars))[1:n.comp]
#     names(rho) = paste("can.comp", 1:n.comp, sep = "")

# 	res = list(	n.comp 	= n.comp, 
# 				rho.mod = rho.mod, 
# 				cor 	= rho,
# 				x.coefs = x.coefs,
# 				x.vars 	= x.vars, 
# 				y.coefs = y.coefs, 
# 				y.vars 	= y.vars, 
# 				lambdas = c(x = X.tr$lambda, y = Y.tr$lambda),
# 				dims 	= c(n1=n1, n2=n2, p1=p1, p2=p2))

# 	new('fastcca', res)
# }


# # Adapted from RCCA:::RCCA.tr
# transformHD = function(X, lambda=NULL, k = min(dim(X))){

# 	# get shrinkage estimate to 1) compute SVD and 2) estimate lambda
# 	ecl = eclairs(X, k = k, lambda=lambda, compute="corr")

# 	# if more features that samples
# 	if( ncol(X) > nrow(X) ){

# 	    # scale to be consistent with PCA
# 		mat = with(ecl, V %*% diag(sqrt(dSq))) * sqrt(nrow(X)-1)
# 		V = ecl$U
# 	}else{
# 		V = NULL
# 		mat = X
# 	}

# 	list(mat = mat, tr = V, lambda = ecl$lambda, Sigma.eclairs=ecl)
# }





# geigen3 = function(X.tr, Y.tr, k){

# 	# D = minvsqrt(Bmat) %*% Amat %*% minvsqrt(Cmat)
# 	tmp = decorrelate(X.tr$mat, X.tr$Sigma.eclairs)
# 	tmp2 = decorrelate(Y.tr$mat, Y.tr$Sigma.eclairs)

# 	n = X.tr$Sigma.eclairs$n
# 	D = crossprod(tmp, tmp2) / (n-1) 

# 	# print(X.tr$mat[1:3, 1:3])
# 	# print(X.tr$Sigma.eclairs$U[1:3,1:3])
# 	# print(tmp2[1:3, 1:3])
# 	print(sum(X.tr$Sigma.eclairs$dSq))

# 	p = nrow(D)
# 	q = ncol(D)

# 	if( q > p){
# 		D = t(D)
# 	}

# 	if( k > min(dim(D)) / 3 ){
# 		dcmp = svd(D)

# 		if( k < min(dim(D)) ){
# 			dcmp$d = dcmp$d[seq_len(k)]
# 			dcmp$u = dcmp$u[,seq_len(k), drop=FALSE]
# 			dcmp$v = dcmp$v[,seq_len(k), drop=FALSE]
# 		}
# 	}else{
# 		dcmp = irlba(D, nv=k)
# 	}

#     values <- dcmp$d

#     if( p > q ){
# 	    # Lmat <- minvsqrt(Bmat) %*% dcmp$u
# 	    # Mmat <- minvsqrt(Cmat) %*% dcmp$v
# 	    Lmat = decorrelate(dcmp$u, X.tr$Sigma.eclairs, transpose=TRUE)
# 	    Mmat = decorrelate(dcmp$v, Y.tr$Sigma.eclairs, transpose=TRUE)
# 	}else{		
# 	    # Lmat <- minvsqrt(Bmat) %*% dcmp$v
# 	    # Mmat <- minvsqrt(Cmat) %*% dcmp$u	    
# 	    Lmat = decorrelate(dcmp$v, X.tr$Sigma.eclairs, transpose=TRUE)
# 	    Mmat = decorrelate(dcmp$u, Y.tr$Sigma.eclairs, transpose=TRUE)
# 	}

#     geigenlist <- list(values, Lmat, Mmat)
#     names(geigenlist) <- c("values", "Lmat", "Mmat")

#     return(geigenlist)
# }


