
# # library(Rfast)
# # library(decorrelate)
# # X = matrnorm(10000, 3) %*% diag(seq(3)+3) 
# # Y = matrnorm(10000, 4) %*% diag(seq(4)+3)
# # # X = whiten(X, lambda=0)
# # # Y = whiten(Y, lambda=0)
# # k = 4
# # k.x = 3
# # k.y = 4
# # compute = "cov"


# fastcca2 = function(X, Y, k = min(dim(X), dim(Y)), k.x = min(dim(X)), k.y = min(dim(Y)), lambda.x = NULL, lambda.y = NULL, compute=c("covariance", "correlation")){

#   compute <- match.arg(compute)

#   if (!is.matrix(X)) {
#     X <- as.matrix(X)
#   }
#   if (!is.matrix(Y)) {
#     Y <- as.matrix(Y)
#   }

#   if (anyNA(X) | anyNA(Y)) {
#     stop("No NA values are allowed in data")
#   }

#   # mean-center columns
#   X <- scale(X, scale = FALSE)
#   Y <- scale(Y, scale = FALSE)

#   # get dimensions
#   n1 <- nrow(X)
#   n2 <- nrow(Y)
#   p1 <- ncol(X)
#   p2 <- ncol(Y)

#   if (n1 != n2) {
#     stop("X and Y must have same number of rows")
#   }

#   # ensure ranks are properly bounded
#   k <- min(dim(X), dim(Y), k)
#   k.x <- min(dim(X), k.x)
#   k.y <- min(dim(Y), k.y)

#   # decomposition of each dataset
#   ecl.x = eclairs(X, lambda = 0, compute="cor")#, k = k.x, lambda = lambda.x)
#   ecl.y = eclairs(Y, lambda = 0, compute="cor")#, k = k.y, lambda = lambda.y)


#   ecl.x = eclairs(X, lambda = 0, compute="cov")
#   ecl.y = eclairs(Y, lambda = 0, compute="cov")

#   C = with(ecl.x, tcrossprod(U,V)) %*% 
#   	  with(ecl.y, tcrossprod(V,U))
#   dcmp = svd( C )

#   A = crossprod(decorrelate::whiten(X, lambda=0), 
#   				decorrelate::whiten(Y, lambda=0)) / (n1-1)
#   C/A

#   # simplest, from CCA::rcc()
#   res = geigen(cov(X,Y), cov(X), cov(Y))
#   dcmp2 = svd(A)
#   decorrelate(t(dcmp2$u), ecl.x)
#   decorrelate(t(dcmp2$v), ecl.y)

#   res2 = geigen(A, diag(3), diag(4))
  
#   diag(cor(res$Lmat, res2$Lmat))
#   diag(cor(res$Mmat, res2$Mmat))

#   # need to incorporate ecl.x$sigma and ecl.y$sigma


#   # NOTE:
#   # 

#   # set diagonals to be positive
#   dcmp$v <- eachrow(dcmp$v, sign(diag(dcmp$v)), "*")
#   dcmp$u <- eachrow(dcmp$u, sign(diag(dcmp$u)), "*")

#   # with(ecl.x, U %*% diag(1/sqrt(dSq)) %*% t(V))


#   qx <- qr(X)
#   qy <- qr(Y)
#   dx <- qx$rank
#   dy <- qy$rank
#   x.coef = backsolve(qx$qr, dcmp$u)
#   y.coef = backsolve(qy$qr, dcmp$v)

#   fit = cancor(X, Y)

#   dcmp$d^2 - fit$cor^2
#   x.coef / fit$xcoef
#   y.coef / fit$ycoef[,1:3]


#   C2 = crossprod(ecl.x$V, ecl.y$V)
#   dcmp2 = svd( C2 )

#   ecl.x$V %*% dcmp2$u





#   C.std = crossprod(decorrelate::whiten(X, lambda=0), 
#   					decorrelate::whiten(Y, lambda=0)) / n1



#   # cross-product matrix
#   # \tilde D_x ^{-\frac{1}{2}}  D_x U_x^T U_y D_y \tilde D_y ^{-\frac{1}{2}} 
#   d_til_x = with(ecl.x, sqrt(dSq / ((1-lambda)*dSq + lambda * nu)))
#   d_til_y = with(ecl.y, sqrt(dSq / ((1-lambda)*dSq + lambda * nu)))

#   C = diag(d_til_x^-0.5) %*% diag(sqrt(ecl.x$dSq)) %*% crossprod(ecl.x$V, ecl.y$V) %*% diag(sqrt(ecl.y$dSq)) %*% diag(d_til_y^-0.5) 
#   dcmp = svd( t(C)  )


#   qx <- qr(X)
#   qy <- qr(Y)
#   dx <- qx$rank
#   dy <- qy$rank

#   x.coef = backsolve(qx$qr[1L:dx, 1L:dx, drop = FALSE], dcmp$v)
#   y.coef = backsolve(qy$qr[1L:dy, 1L:dy, drop = FALSE], dcmp$u)

#   fit = cancor(X, Y)

#   dcmp$d^2 - fit$cor^2
#   t(x.coef) - fit$xcoef
#   t(y.coef) - fit$ycoef

#   cxx = cov(X)
#   cyy <- cov(Y)
#   cxy <- cov(X, Y)
#   cyx <- t(cxy)
#   solve(cyy, cyx) %*% solve(cxx, cxy)
#   crossprod(C.cc)


#   qx <- qr(X)
#   qy <- qr(Y)
#   nr = nrow(X)
#   C.cc = qr.qty(qx, qr.qy(qy, diag(1, nr, dy)))[1L:dx, , drop = FALSE]

#   fit = cancor(X, Y)


#   C.std = crossprod(whiten(X), whiten(Y)) / n1
#   dcmp.std = svd(C.std)
#   fit = yacca::cca(X, Y, reg.param=0)

#   # SVD of cross-product matrix
#   if( k > min(dim(C))/3 ){
#   	dcmp = svd( C )
#   }else{
#   	dcmp = irlba(C, nu = k, nv = k)
#   }

#   # Backtransform coefficients


# }







