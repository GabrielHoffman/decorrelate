



# #' Generalized PCA
# #'
# #' Generalized PCA using eclairs estimate of correlation betweeen features
# #'
# #' @importFrom irlba irlba
# #' @export
# gpca = function(Y, cor.est, k){

# 	n = nrow(Y)
# 	p = ncol(Y)

# 	# k cannot exceed n or p
# 	k = min(c(k, p, n))

# 	if( ! is.null(cor.est) & ! missing(cor.est)){
# 		# whiten
# 		# Y_tilde = t(lrmult_right(t(Y), cor.est$U, cor.est$dSq *(1-cor.est$lambda), cor.est$lambda, -.5))

# 		# Y_tilde = lrmult_left(Y, cor.est$U, cor.est$dSq *(1-cor.est$lambda), cor.est$lambda, -.5)
# 		Y_tilde = decorrelate(Y, cor.est)

# 		# range(Y_tilde-Y_tilde2)

# 		# browser()

# 		# SVD in whitened space
# 		if( k < min(p, n)/2){
# 			dcmp.tilde = irlba(Y_tilde, k)
# 		}else{
# 			dcmp.tilde = svd(Y_tilde)
# 			dcmp.tilde$v = dcmp.tilde$v[,seq_len(k), drop=FALSE]
# 			dcmp.tilde$u = dcmp.tilde$u[,seq_len(k), drop=FALSE]
# 			dcmp.tilde$d = dcmp.tilde$d[seq_len(k)]
# 		}

# 		# back transform
# 		# V.star = t(mult_eclairs(t(dcmp.tilde$v), cor.est$U, cor.est$dSq, cor.est$lambda, cor.est$nu, 0.5))
# 		# April 14
# 		V.star = mult_eclairs(dcmp.tilde$v, cor.est$U, cor.est$dSq, cor.est$lambda, cor.est$nu, 0.5, transpose=TRUE)
# 		U.star = dcmp.tilde$u
# 		d.star = dcmp.tilde$d
# 	}else{
# 		if( k < min(p, n)/3){
# 			dcmp.tilde = irlba(Y, k)	
# 		}else{
# 			dcmp.tilde = svd(Y)
# 		}
# 		V.star = dcmp.tilde$v
# 		U.star = dcmp.tilde$u
# 		d.star = dcmp.tilde$d
# 	}

# 	# name dimensions
# 	colnames(U.star) = paste0("PC", seq_len(k))
# 	rownames(U.star) = rownames(Y)
# 	colnames(V.star) = paste0("PC", seq_len(k))
# 	rownames(V.star) = colnames(Y)

# 	# Modify sign of dcmp$v and dcmp$u so principal components are consistant
# 	# This is motivated by whitening:::makePositivDiagonal()
# 	# but here adjust both U and V so reconstructed data is correct
# 	values = sign(diag(U.star))
# 	U.star = sweep(U.star, 2, values, "*")
# 	V.star = sweep(V.star, 2, values, "*")

# 	list(	U 		= U.star[,seq_len(k), drop=FALSE],
# 	 		V 		= V.star[,seq_len(k), drop=FALSE],
# 			D 		= d.star[seq_len(k)])#,
# 			# decomp 	= dcmp.tilde)
# }



# #' Empirically Reweighted Generalized PCA
# #'
# #' Empirically Reweighted Generalized PCA learns correlation structure between features in order to improve PCA projection
# #'
# #' @importFrom Rfast standardise
# #' @export
# ergpca = function(Y, k, lambda, k.features, cor.est = NULL, tol=1e-5, maxit=10000){

# 	Y.scale = standardise(Y, scale=FALSE)
# 	n = nrow(Y)
# 	p = ncol(Y)

# 	# warmStart = NULL
# 	lambda = NULL

# 	if( maxit > 0){
# 		for(i in 1:maxit){

# 			# if( i ==30) browser()
# 		 	# Perform Generalized PCA based on ecliars covariance structure
# 		 	# if cor.est == NULL, it is equivalent to identity correlation matrix
# 			res.gpca = gpca(Y.scale, cor.est, ifelse(i>1, k, min(dim(Y.scale))) )

# 			res.gpca$iter = i

# 			if( (i == 1) & missing(k.features) ){
# 				k.features = sv_threshold( n, p, res.gpca$D)
# 				message("k.features: ", k.features)
# 			}

# 			# Estimate residuals
# 			# resid = Y - with(res.gpca, U %*% diag(D) %*% t(V))
# 			resid = Y.scale - with(res.gpca, U %*% (D * t(V)))

# 			# if( i > 2){
# 			# 	warmStart = cor.est$decomp
# 			# }

# 			# Estimate correlation structure based on residuals
# 		 	cor.est = eclairs( resid, lambda = lambda, k=k.features, compute="corr")#, warmStart=warmStart)
# 		 	lambda = cor.est$lambda

# 		 	if(i > 2){
# 				delta = mean(abs(d.prev - res.gpca$D))
# 				if(i %% 1 == 0) message(i, ": ", delta)
# 				if( delta < tol) break			
# 			}
# 			d.prev = res.gpca$D
# 		}
# 	}

# 	res.gpca = gpca(Y.scale, cor.est, k)

# 	res.gpca$cor.est = cor.est
# 	res.gpca$lambda = lambda

# 	list(res.gpca, cor.est)
# }











