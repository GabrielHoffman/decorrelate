# Gabriel Hoffman
# July 12, 2021
#
# Examining risk function

library(clusterGeneration)
library(Rfast)


tr = function(x) sum(diag(x))

minvsqrt = function(S){
	dcmp = eigen(S)
	with(dcmp, vectors %*% diag(1/sqrt(values)) %*% t(vectors))
}


msqrt = function(S){
	dcmp = eigen(S)
	with(dcmp, vectors %*% diag(sqrt(values)) %*% t(vectors))
}

# lowest eigen values
omega = function(C) min(eigen(C)$values)


# largest eigen values
alpha = function(C) max(eigen(C)$values)

n = 100
p = 3
set.seed(1)

Sigma = genPositiveDefMat(p)$Sigma 
min(eigen(Sigma)$values)

X = matrnorm(n,p) / sqrt(n*p)
tr(crossprod(X))

Y = X %*% chol(Sigma)

min(eigen(Sigma)$values) > 1 - tr(sqrt(eigen(Sigma)$values))/tr(crossprod(X))

cor(Y)

lambda = 0.1
SigmaHat = (1-lambda) * cor(Y) + diag(lambda,p)

cor(Y %*% minvsqrt(SigmaHat))



risk1 = function(lambda){

	SigmaHat = (1-lambda) * cor(Y) + diag(lambda,p)
	tr(Sigma %*% solve(SigmaHat)) + 
	tr(X %*% msqrt(Sigma) %*% solve(SigmaHat) %*% t(msqrt(Sigma)) %*% t(X)) - 
	2*tr(X %*% msqrt(Sigma) %*% minvsqrt(SigmaHat) %*% t(X))
}

x = seq(1e-4, 1-1e-4, length.out=1000)
y = sapply(x, risk1)
which.min(y)

plot(x,y, type='l')



risk2 = function(lambda){

	nu = 1
	Y_scale = scale(Y) / sqrt(n-1)
	SigmaHat = (1-lambda) * crossprod(Y_scale) + diag(lambda,p)
	d = svd(Y_scale)$d
	s = (1-lambda)*d^2 + lambda
	# eigen(SigmaHat)$values
	U = eigen(SigmaHat)$vectors

	# tr(Sigma %*% solve(SigmaHat)) 
	# tr(Sigma %*% U %*% diag(1/s) %*% t(U)) 
	tr(t(U) %*% Sigma %*% U %*% diag(1/s)) +

	# tr(X %*% msqrt(Sigma) %*% solve(SigmaHat) %*% t(msqrt(Sigma)) %*% t(X)) 
	tr( t(U) %*% t(msqrt(Sigma)) %*% t(X) %*% X %*% msqrt(Sigma) %*% U %*% diag(1/s)) - 

	# 2*tr(X %*% msqrt(Sigma) %*% minvsqrt(SigmaHat) %*% t(X))

	# 2*tr(X %*% msqrt(Sigma) %*% U %*% diag(s^-.5) %*% t(U)%*% t(X))
	2*tr(t(U)%*%t(X) %*% X %*% msqrt(Sigma) %*% U %*% diag(s^-.5) )
}

y = sapply(x, risk2)
which.min(y)

plot(x,y, type='l')

# derivatives
#############

library(Deriv)
f = function(x){
	((1-x)*dSq + x * n)^-.5
}

Deriv(f)
Deriv(f, nderiv=2)

g = Deriv(f, nderiv=2)
dSq = 3

plot(x,g)

x = .2
g(x)

s = (1-x)*dSq + x * n
3/4 * (n-dSq)^2 / s^2.5


.e3 <- s
.e4 <- sqrt(.e3)
0.5 * ((0.5 * (.e3/.e4) + .e4) * (n - dSq)^2/(.e3 * .e4)^2)


0.5 * ((0.5 * sqrt(s) + sqrt(s)) * (n - dSq)^2/s^3)


0.5 * (1.5*sqrt(s) * (n - dSq)^2/s^3)


3/4 *  (n - dSq)^2/s^2.5




drisk = function(lambda){

	nu = 1
	Y_scale = scale(Y) / sqrt(n-1)
	SigmaHat = (1-lambda) * crossprod(Y_scale) + diag(lambda,p)
	d = svd(Y_scale)$d
	s = (1-lambda)*d^2 + nu*lambda
	# eigen(SigmaHat)$values
	U = eigen(SigmaHat)$vectors
 
	# tr(t(U) %*% Sigma %*% U %*% diag(-(nu- d^2)/s^2)) +
	# tr( t(U) %*% t(msqrt(Sigma)) %*% t(X) %*% X %*% msqrt(Sigma) %*% U %*% diag(-(nu- d^2)/s^2)) - 
	# 2*tr(t(U)%*%t(X) %*% X %*% msqrt(Sigma) %*% U %*% diag(-0.5*(nu- d^2)/s^1.5) )

	# tr(t(U) %*% msqrt(Sigma) %*% (diag(1,p) + crossprod(X)) %*% msqrt(Sigma) %*% U %*% diag(-(nu- d^2)/s^2)) - 
	# 2*tr(t(U)%*%t(X) %*% X %*% msqrt(Sigma) %*% U %*% diag(-0.5*(nu- d^2)/s^1.5) )

	- tr(t(U) %*% msqrt(Sigma) %*% (diag(1,p) + crossprod(X)) %*% msqrt(Sigma) %*% U %*% diag((nu- d^2)/s^2)) +
	tr(t(U)%*%t(X) %*% X %*% msqrt(Sigma) %*% U %*% diag((nu- d^2)/s^1.5) )
}

y = sapply(x, drisk)
which.min( (y)^2)

plot(x,y, type='l')
abline(h=0, lty="dashed")

# always positive when 

tr(t(U)%*%t(X) %*% X %*% msqrt(Sigma) %*% U %*% diag((nu- d^2)/s^1.5) ) > tr(t(U) %*% msqrt(Sigma) %*% (diag(1,p) + crossprod(X)) %*% msqrt(Sigma) %*% U %*% diag((nu- d^2)/s^2))




# larger
tr(t(U)%*%t(X) %*% X %*% msqrt(Sigma) %*% U %*% diag((nu- d^2)/s^1.5) )
tr(t(X) %*% X %*% msqrt(Sigma)) * max((nu- d^2)/s^1.5)
tr(t(X) %*% X) * alpha(msqrt(Sigma)) * max((nu- d^2)/s^1.5)

# smaller
tr(t(U) %*% msqrt(Sigma) %*% (diag(1,p) + crossprod(X)) %*% msqrt(Sigma) %*% U %*% diag((nu- d^2)/s^2))
tr(Sigma %*% (diag(1,p) + crossprod(X))) * min((nu- d^2)/s^2)
(tr(Sigma) + tr(Sigma %*% crossprod(X))) * min((nu- d^2)/s^2)
(tr(Sigma) + alpha(Sigma) * tr(crossprod(X))) * min((nu- d^2)/s^2)





















