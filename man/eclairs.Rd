% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/eclairs.R
\name{eclairs}
\alias{eclairs}
\title{Estimate covariance/correlation with low rank and shrinkage}
\usage{
eclairs(
  X,
  k,
  lambda = NULL,
  compute = c("covariance", "correlation"),
  n.samples = nrow(X)
)
}
\arguments{
\item{X}{data matrix with n samples as rows and p features as columns}

\item{k}{the rank of the low rank component}

\item{lambda}{shrinkage parameter. If not specified, it is estimated from the data.}

\item{compute}{evaluate either the \code{"covariance"} or \code{"correlation"} of \code{X}}

\item{n.samples}{number of samples data is from.  Usually \code{nrow(X)}, but can be other values in special cases.}
}
\value{
\link{eclairs} object storing:
\describe{
 \item{U: }{orthonormal matrix with k columns representing the low rank component}
 \item{dSq: }{eigen-values so that \eqn{U diag(d^2) U^T} is the low rank component}
 \item{lambda: }{shrinkage parameter \eqn{\lambda} for the scaled diagonal component}
 \item{sigma: }{standard deviations of input columns}
 \item{nu: }{diagonal value, \eqn{\nu}, of target matrix in shrinkage}
 \item{n: }{number of samples (i.e. rows) in the original data}
 \item{p: }{number of features (i.e. columns) in the original data}
 \item{k: }{rank of low rank component}
 \item{rownames: }{sample names from the original matrix}
 \item{colnames: }{features names from the original matrix}
 \item{method: }{method used for decomposition}
 \item{call: }{the function call}
}
}
\description{
Estimate the covariance/correlation between columns as the weighted sum of a low rank matrix and a scaled identity matrix.  The weight acts to shrink the sample correlation matrix towards the identity matrix or the sample covariance matrix towards a scaled identity matrix with constant variance.  An estimate of this form is useful because it is fast, and enables fast operations downstream.  The method is based on the Gaussian Inverse Wishart Empirical Bayes (GIW-EB) model.
}
\details{
Compute \eqn{U}, \eqn{d^2} to approximate the correlation matrix between columns of data matrix X by \eqn{U diag(d^2 (1-\lambda)) U^T + I\nu\lambda}.  When computing the covariance matrix, scale by the standard deviation of each feature.
}
\examples{
library(Rfast)

n <- 800 # number of samples
p <- 200 # number of features

# create correlation matrix
Sigma <- autocorr.mat(p, .9)

# draw data from correlation matrix Sigma
Y <- rmvnorm(n, rep(0, p), sigma = Sigma * 5.1, seed = 1)
rownames(Y) <- paste0("sample_", seq(n))
colnames(Y) <- paste0("gene_", seq(p))

# eclairs decomposition: covariance
ecl <- eclairs(Y, compute = "covariance")

ecl

# eclairs decomposition: correlation
ecl <- eclairs(Y, compute = "correlation")

ecl

}
