

test_logDet = function(){

	library(Rfast)
	set.seed(1)
	n = 800 
	p = 12

	# create correlation matrix
	Sigma = autocorr.mat(p, .9)

	# draw data from correlation matrix Sigma
	Y = rmvnorm(n, rep(0, p), sigma=Sigma)

	# eclairs decomposition
	Sigma.eclairs = eclairs( Y )

	a = logDet(Sigma.eclairs)
	b = determinant(getCov(Sigma.eclairs))$modulus
	checkEqualsNumeric(a,b)

	a = logDet(Sigma.eclairs, 0.5)
	b = determinant(chol(getCov(Sigma.eclairs)))$modulus
	checkEqualsNumeric(a,b)

	a = logDet(Sigma.eclairs, -1)
	b = determinant(solve(getCov(Sigma.eclairs)))$modulus
	checkEqualsNumeric(a,b)
}


test_tr = function(){

	library(Rfast)
	set.seed(1)
	n = 800 
	p = 12

	# create correlation matrix
	Sigma = autocorr.mat(p, .9)

	# draw data from correlation matrix Sigma
	Y = rmvnorm(n, rep(0, p), sigma=Sigma)

	# eclairs decomposition
	Sigma.eclairs = eclairs( Y, compute="corr" )

	a = tr(Sigma.eclairs)
	b = sum(eigen(getCor(Sigma.eclairs))$values)
	checkEqualsNumeric(a,b)
}




test_kappa = function(){

	library(Rfast)
	set.seed(1)
	n = 800 
	p = 12

	# create correlation matrix
	Sigma = autocorr.mat(p, .9)

	# draw data from correlation matrix Sigma
	Y = rmvnorm(n, rep(0, p), sigma=Sigma)

	# eclairs decomposition
	Sigma.eclairs = eclairs( Y, compute="corr" )

	a = kappa(Sigma.eclairs)
	b = kappa(getCor(Sigma.eclairs), exact=TRUE)
	checkEqualsNumeric(a,b)
}









test_quadForm = function(){

	library(Rfast)
	set.seed(1)
	n = 800 
	p = 12

	# create correlation matrix
	Sigma = autocorr.mat(p, .9)

	# draw data from correlation matrix Sigma
	Y = rmvnorm(n, rep(0, p), sigma=Sigma*4.5)

	# eclairs decomposition
	Sigma.eclairs = eclairs( Y )

	# return scalar
	a = quadForm( Sigma.eclairs, Y[1,])
	b = Y[1,] %*% solve(getCov(Sigma.eclairs)) %*% ( Y[1,])
	checkEqualsNumeric(a,b)

	# return matrix
	a = quadForm( Sigma.eclairs, Y[1:2,])
	b = Y[1:2,] %*% solve(getCov(Sigma.eclairs)) %*% t( Y[1:2,])
	checkEqualsNumeric(a,b)
}








