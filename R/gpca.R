

#' Generalized PCA
#'
#' Generalized PCA using eclairs estimate of correlation betweeen features
#'
#' @export
gpca_eclairs = function(Y, cor.est, k){

	n = nrow(Y)
	p = ncol(Y)

	if( ! is.null(cor.est) ){
		# whiten
		# Y_tilde = t(lrmult_right(t(Y), cor.est$U, cor.est$dSq *(1-cor.est$lambda), cor.est$lambda, -.5))

		# Y_tilde = lrmult_left(Y, cor.est$U, cor.est$dSq *(1-cor.est$lambda), cor.est$lambda, -.5)
		Y_tilde = decorrelate(Y, cor.est)

		# range(Y_tilde-Y_tilde2)

		# browser()

		# SVD in whitened space
		if( k < min(p, n)/2){
			dcmp.tilde = irlba(Y_tilde, k)
		}else{
			dcmp.tilde = svd(Y_tilde)
		}

		# back transform
		V.star = t(mult_eclairs(t(dcmp.tilde$v), cor.est$U, cor.est$dSq, cor.est$lambda, 0.5))
		U.star = dcmp.tilde$u
		d.star = dcmp.tilde$d
	}else{
		if( k < min(p, n)/2){
			dcmp.tilde = irlba(Y, k)	
		}else{
			dcmp.tilde = svd(Y)
		}
		V.star = dcmp.tilde$v
		U.star = dcmp.tilde$u
		d.star = dcmp.tilde$d
	}

	list(U = U.star[,1:k],
		V = V.star[,1:k],
		D = d.star[1:k])
}



#' PCA learning correlation structure
#'
#' PCA learning correlation structure
#'
#' @export
pca_cor2 = function(Y, k, lambda, k.features, tol=1e-5, maxit=10000){

	Y.scale = scale(Y)
	n = nrow(Y)
	p = ncol(Y)

	# initialize to identify correlation matrix
	cor.est = NULL

	for(i in 1:maxit){

	 	# Perform Generalized PCA based on ecliars covariance structure
		res.gpca = gpca_eclairs(Y.scale, cor.est, ifelse(i>1, k, min(dim(Y.scale))) )

		if( i == 1){
			k.features = sv_threshold( n, p, dcmp$D)
			message("k.features: ", k.features)
		}

		# Estimate residuals
		# resid = Y - with(res.gpca, U %*% diag(D) %*% t(V))
		resid = Y.scale - with(res.gpca, U %*% (D * t(V)))

		# Estimate correlation structure based on residuals
	 	cor.est = eclairs( resid, lambda = lambda, k=k.features)

	 	if(i > 2){
			delta = mean(abs(d.prev - res.gpca$D))
			if(i %% 20 == 0) message(i, ": ", delta)

			if( delta < tol) break
		}
		d.prev = res.gpca$D
	}

	res.gpca = gpca_eclairs(Y.scale, cor.est, k)

	res.gpca$cor.est = cor.est

	res.gpca
}

