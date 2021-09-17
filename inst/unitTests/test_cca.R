# Sept 9, 2021

test_cca = function(){

	library(Rfast)
	set.seed(1)
	n = 100 # number of samples
	p1 = 3 # number of features
	p2 = 4

	# draw data from correlation matrix Sigma
	X = rmvnorm(n, rep(0, p1), autocorr.mat(p1, .9))
	Y = rmvnorm(n, rep(0, p2), autocorr.mat(p2, .9))

	res1 = decorrelate::cca(X,Y, lambda.x=0, lambda.y=0)
	res2 = CCA::rcc(X,Y, 0,0)
	res3 = RCCA::RCCA(X,Y, 0,0)

	checkEqualsNumeric(res1$rho.mod, res2$cor)
	checkEqualsNumeric(abs(res1$x.coefs), abs(res2$xcoef))
	checkEqualsNumeric(abs(res1$y.coefs), abs(res2$ycoef))

	checkEqualsNumeric(res1$rho.mod, res3$cor)
	checkEqualsNumeric(abs(res1$x.coefs), abs(res3$x.coefs))
	checkEqualsNumeric(abs(res1$y.coefs), abs(res3$y.coefs))

	# set.seed(1)
	# n = 100 # number of samples
	# p1 = 300 # number of features
	# p2 = 400

	# # draw data from correlation matrix Sigma
	# X = rmvnorm(n, rep(0, p1), autocorr.mat(p1, .9))
	# Y = rmvnorm(n, rep(0, p2), autocorr.mat(p2, .9))

	# res1 = cca(X,Y)
	# res2 = CCA::rcc(X,Y, 0.1, 0.1)
	# res3 = RCCA::RCCA(X,Y, 0.1, 0.1)

	# checkEqualsNumeric(res1$cor, res2$cor)
	# checkEqualsNumeric(abs(res1$xcoef), abs(res2$xcoef))
	# checkEqualsNumeric(abs(res1$ycoef), abs(res2$ycoef))

	# checkEqualsNumeric(res2$cor, res3$cor)
	# checkEqualsNumeric(abs(res2$xcoef), abs(res3$x.coefs))
	# checkEqualsNumeric(abs(res2$ycoef), abs(res3$y.coefs))


	# checkEqualsNumeric(res1$cor, res3$cor)
	# checkEqualsNumeric(abs(res1$xcoef), abs(res3$x.coefs))
	# checkEqualsNumeric(abs(res1$ycoef), abs(res3$y.coefs))



}	


# tr1 = function(X, lambda=NULL){

# 	ecl = eclairs(X, lambda=lambda, compute="corr")

# 	if( ncol(X) > nrow(X) ){

# 		# SVD = svd(X)
# 	 #    A = SVD$u %*% diag(SVD$d)

# 	    # scale to be consistent with PCA
# 		mat = with(ecl, V %*% diag(sqrt(dSq))) * sqrt(nrow(X)-1)

# 		# A[1:4, 1:4]
# 		# mat[1:4, 1:4]

# 		V = ecl$U
# 	}else{

# 		V = NULL
# 		mat = X
# 	}

# 	C = with(ecl, (1-lambda)*var(mat, use = "pairwise") + diag(lambda*nu, ncol(mat)))

# 	list(mat = mat, tr = V, cor=C, lambda = ecl$lambda)
# }

# n = 300
# p = 301
# X = rmvnorm(n, rep(0, p), autocorr.mat(p, .9))
# X = scale(X)

# res = RCCA:::RCCA.tr(X, 1e-20)
# str(res)

# res2 = tr1(X, 0)
# str(res2)




# n = 100
# p = 3
# X = rmvnorm(n, rep(0, p), autocorr.mat(p, .9))
	
# res = RCCA:::RCCA.tr(X, 0)
# str(res)

# res2 = tr1(X, 0)
# str(res2)

test_fastcca_big = function(){
	
	library(Rfast)
	n = 100 # number of samples
	p1 = 800 # number of features
	p2 = 800

	# draw data from correlation matrix Sigma
	# X = rmvnorm(n, rep(0, p1), autocorr.mat(p1, .9))
	# Y = rmvnorm(n, rep(0, p2), autocorr.mat(p2, .9))

	X = matrnorm(n,p1)
	Y = matrnorm(n,p2)

	X = scale(X)
	Y = scale(Y)

	resa = RCCA::RCCA(X, Y, 1e-5, 1e-5)

	system.time(res <- fastcca(X, Y))

	system.time(res2 <- decorrelate::cca(X, Y))

	checkEqualsNumeric(res$rho.mod[1:(n-1)], res2$rho.mod[1:(n-1)])

	checkEqualsNumeric(abs(res$x.vars)[,1:10], abs(X %*% res2$x.coefs)[,1:10])

	# cor(res$x.vars[,1:5], (X %*% res2$xcoef))
	checkEqualsNumeric( sum(cancor(res$x.vars[,1:5], (X %*% res2$x.coefs))$cor), 5)

	checkEqualsNumeric(abs(res$x.coefs[1:3, 1:3]), abs(res2$x.coefs[1:3, 1:3]))

	checkEqualsNumeric(res$lambdas, res2$lambdas)
}


# source("/Users/gabrielhoffman/workspace/repos/decorrelate/R/fastcca.R")


test_fastcca = function(){
	
	library(Rfast)
	# check in simple, low dimensional analysis
	n = 100 # number of samples
	p1 = 3 # number of features
	p2 = 4

	# draw data from correlation matrix Sigma
	X = rmvnorm(n, rep(0, p1), autocorr.mat(p1, .9))
	Y = rmvnorm(n, rep(0, p2), autocorr.mat(p2, .9))

	X = scale(X)
	Y = scale(Y)

	res = fastcca(X, Y, lambda.x=0, lambda.y=0)

	res2 = cca(X, Y, lambda.x=0, lambda.y=0)

	checkEqualsNumeric(res$rho.mod, res2$rho.mod)

	checkEqualsNumeric(abs(res$x.vars),abs(X %*% res2$x.coefs))

	checkEqualsNumeric(abs(res$x.coefs), abs(res2$x.coefs))

	checkEqualsNumeric(res$lambdas, res2$lambdas)
}










