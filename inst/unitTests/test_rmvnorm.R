

test_rmvnorm = function(){

  library(Rfast) 
  library(RUnit) 

  X <- matrnorm(10000, 10)
  mu <- colmeans(X)
  Sigma <- cov(X)
  ecl <- eclairs(X, lambda=1e-8)

  a = dmvnorm(X, mu, Sigma, log=TRUE) 
  b = dmvnorm_eclairs( X, mu, ecl, log=TRUE)
  checkEqualsNumeric(a,b)

  a = dmvnorm(X, mu, Sigma, log=FALSE) 
  b = dmvnorm_eclairs( X, mu, ecl, log=FALSE)
  checkEqualsNumeric(a,b)

  a = dmvt(X, mu, Sigma, nu=5, log=TRUE) 
  b = dmvt_eclairs( X, mu, ecl, nu=5,  log=TRUE)
  checkEqualsNumeric(a,b)

  a = dmvt(X, mu, Sigma, nu=5, log=FALSE) 
  b = dmvt_eclairs( X, mu, ecl, nu=5,  log=FALSE)
  checkEqualsNumeric(a,b)
}