
# 1) predictions are the same when lambda = 0
# 2) coefs are different but 
#     fit1$x.coefs %*% ginv(fit1$y.coefs) are the same
# 3) d is different
# y.vars and x.vars are different? get U instead?

q()
R

library(decorrelate)
library(Rfast)
library(MASS)
library(RUnit)

# source("R/multiwaycca.R")
multiwaycca = decorrelate:::multiwaycca

dat1 <- LifeCycleSavings[, 1:2,drop=FALSE] 
dat2 <- LifeCycleSavings[, 3:5,drop=FALSE]
dat3 <- LifeCycleSavings[, 5,drop=FALSE]

lambda = 0
fit1 = fastcca(dat1, dat2, lambda.x=lambda, lambda.y=lambda)
fit2 <- multiwaycca(list( dat1, dat2), 2, lambdas = c(lambda, lambda))

fit2$coef[[1]] %*% ginv(fit2$coef[[2]])
fit1$x.coefs %*% ginv(fit1$y.coefs)







y.pred1 = predict(fit1, X = dat1)
y.pred2 = predict(fit2, dat1, given="dataset_1", target="dataset_2")

checkEqualsNumeric(y.pred1, y.pred2)


ecl1 = eclairs(as.matrix(dat1), lambda=lambda)
ecl2 = eclairs(as.matrix(dat2), lambda=lambda)

fit1 = fastcca(ecl1, dat2, lambda.x=lambda, lambda.y=lambda)
fit2 <- multiwaycca(list(ecl1, dat2), 2, lambdas = c(lambda, lambda))

y.pred1 = predict(fit1, X = dat1)
y.pred2 = predict(fit2, dat1, given="dataset_1", target="dataset_2")

checkEqualsNumeric(y.pred1, y.pred2)





fit2$coef[[1]] %*% ginv(fit2$coef[[2]])
fit1$x.coefs %*% ginv(fit1$y.coefs)



cor(fit1$x.vars, fit2$vars[[1]])
cor(fit1$y.vars, fit2$vars[[2]])


fit$x.coefs
fit$y.coefs

fit2$coef





fit2



fit$x.coefs / fit2$coef[[1]]
fit$y.coefs / fit2$coef[[2]]



head(fit$x.vars)
head(fit2$vars[[1]])










datasets = list( dat1, dat2)
ks = c(NULL, NULL)
lambdas = c(0, 0)

Xlist = datasets






multiCCA_concat <- function(Xlist, k = 2) {
  # Xlist: list of datasets with same number of rows
  # k: number of canonical components
  
  m <- length(Xlist) # number of datasets
  n <- nrow(Xlist[[1]])
  
  # Center each dataset
  Xc <- lapply(Xlist, function(X) scale(X, center = TRUE, scale = TRUE))
  
  # SVD + whitening
  svd_list <- lapply(Xc, function(X) {
    sv <- svd(X)
    # Keep nonzero singular values
    keep <- which(sv$d > 1e-12)
    list(
      u = sv$u[, keep, drop = FALSE],
      d = sv$d[keep],
      v = sv$v[, keep, drop = FALSE]
    )
  })
  
  # Whitened data: U = X V D^{-1}
  Z <- lapply(seq_len(m), function(i) {
    Xc[[i]] %*% svd_list[[i]]$v %*% diag(1 / svd_list[[i]]$d)
  })
  
  # Concatenate whitened data horizontally
  Zcat <- do.call(cbind, Z)
  
  # SVD of concatenated whitened data
  svd_cat <- svd(Zcat)
  
  # Take top-k shared components
  U_k <- svd_cat$u[, 1:k, drop = FALSE]
  
  # Map back to original space to get canonical weights
  Wlist <- lapply(seq_len(m), function(i) {
    svd_list[[i]]$v %*% diag(1 / svd_list[[i]]$d) %*% 
      t(Z[[i]]) %*% U_k
  })
  
  # Canonical variates (scores) for each dataset
  Clist <- lapply(seq_len(m), function(i) Xc[[i]] %*% Wlist[[i]])
  
  list(weights = Wlist,
       scores = Clist,
       shared = U_k,
       singular_values = svd_cat$d[1:k])
}

# ===== Example usage =====
set.seed(123)
n <- 50
X1 <- matrix(rnorm(n * 5), n, 5)
X2 <- matrix(rnorm(n * 6), n, 6)
X3 <- matrix(rnorm(n * 4), n, 4)

res <- multiCCA_concat(list(X1, X2, X3), k = 2)

# Canonical variates correlations
cor(res$scores[[1]], res$scores[[2]])
cor(res$scores[[1]], res$scores[[3]])
cor(res$scores[[2]], res$scores[[3]])

