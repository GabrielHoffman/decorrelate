% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/whiten.R
\name{whiten}
\alias{whiten}
\title{Decorrelation projection + eclairs}
\usage{
whiten(X, k = ncol(X), lambda = NULL)
}
\arguments{
\item{X}{matrix to be transformed so *columns* are independent}

\item{k}{the rank of the low rank component}

\item{lambda}{specify lambda and override value estimated by \code{eclairs()}}
}
\value{
data rotated and scaled according to the regularized sample covariance of the input data
}
\description{
Efficient decorrelation projection using eclairs decomposition
}
\examples{
library(Rfast)

n <- 800 # number of samples
p <- 200 # number of features

# create correlation matrix
Sigma <- autocorr.mat(p, .9)

# draw data from correlation matrix Sigma
Y <- rmvnorm(n, rep(0, p), sigma = Sigma * 5.1, seed = 1)

# eclairs decomposition
ecl <- eclairs(Y)

# whitened Y
Y.transform <- decorrelate(Y, ecl)

# Combine eclairs and decorrelate into one step
Y.transform2 <- whiten(Y)

}
