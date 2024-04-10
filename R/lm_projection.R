# Gabriel Hoffman


# Consider a regression with response y, covariates X, and additional variable
# A[,j] that is of primary interest but must be evaluated for many values of j.
# Linear regresion is cubic iin the number of covariates But the constant
# covariates can be used to project the reponse and the new covariates into a
# new lower dimensional space.  With 10 covariates, this can be 20X faster Plus
# this takes advantage of sparseMatrix A, while lm does not







lm.projection <- function(y, X) {
  # code for full method\t Py = y - X %*% (tcrossprod(solve(crossprod(X)),X)
  # %*% y) Pa = a - X %*% (tcrossprod(solve(crossprod(X)),X) %*% a) beta =
  # solve(crossprod(Pa)) %*% crossprod(Pa, Py) sse = sum((Py - Pa %*%
  # beta)^2) sqrt(solve(crossprod(Pa)) * sse/(n-ncol(X)-ncol(a)))

  if (is.matrix(y) && (ncol(y) > 1)) {
    stop("Only one response variable is allowed")
  }
  # M <- tcrossprod(solve(crossprod(X)), X)
  M <- solve(crossprod(X), t(X))

  Py <- y - X %*% (M %*% y)
  # Pa = a - X %*% (M %*% a)

  list(M = M, Py = Py, X = X)
}


# @examples # Hypothesis test of Sepal.Width from full model 
# fit = lm(Petal.Width ~ Sepal.Length + Sepal.Width, iris)
# coef(summary(fit))['Sepal.Width',]

# # # Hypothesis test of Sepal.Width using pre-fit model 
# obj = decorrelate:::lm.projection(iris$Petal.Width, model.matrix(~Sepal.Length, iris))
# decorrelate:::lm.test(obj, iris$Sepal.Width)
lm.test <- function(obj, A, two.sided = TRUE) {
  if (!is.matrix(A)) {
    A <- matrix(A, ncol = 1)
  }
  n <- nrow(A)

  Pa <- as.matrix(A - obj$X %*% (obj$M %*% A))

  D <- solve(crossprod(Pa))
  beta <- D %*% crossprod(Pa, obj$Py)

  sse <- sum((obj$Py - Pa %*% beta)^2)

  rdf <- n - ncol(obj$X) - ncol(A)
  stderr <- sqrt(diag(D) * sse / rdf)

  if (two.sided) {
    pvalues <- 2 * pt(abs(beta / stderr), rdf, lower.tail = FALSE)
  } else {
    # one sided test\t\t
    pvalues <- pt(beta / stderr, rdf, lower.tail = FALSE)
  }

  res <- data.frame(
    Estimate = beta, `Std. Error` = stderr, `t value` = beta / stderr,
    `Pr(>|t|)` = pvalues, check.names = FALSE
  )
  rownames(res) <- paste0("A", colnames(A))
  as.matrix(res)
}




