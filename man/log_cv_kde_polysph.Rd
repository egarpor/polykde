% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RcppExports.R
\name{log_cv_kde_polysph}
\alias{log_cv_kde_polysph}
\title{Cross-validation for the polyspherical kernel density estimator}
\usage{
log_cv_kde_polysph(X, d, h, weights = as.numeric(c()), wrt_unif = FALSE,
  normalized = TRUE, intrinsic = FALSE, norm_X = FALSE, kernel = 1L,
  kernel_type = 1L, k = 10)
}
\arguments{
\item{X}{a matrix of size \code{c(n, sum(d) + r)} with the sample.}

\item{d}{vector of size \code{r} with dimensions.}

\item{h}{vector of size \code{r} with bandwidths.}

\item{weights}{weights for each observation. If provided, a vector of size
\code{n} with the weights for multiplying each kernel. If not provided,
set internally to \code{rep(1 / n, n)}, which gives the standard estimator.}

\item{wrt_unif}{flag to return a density with respect to the uniform
measure. If \code{FALSE} (default), the density is with respect to the
Lebesgue measure.}

\item{normalized}{flag to compute the normalizing constant of the kernel
and include it in the kernel density estimator. Defaults to \code{TRUE}.}

\item{intrinsic}{use the intrinsic distance, instead of the
extrinsic-chordal distance, in the kernel? Defaults to \code{FALSE}.}

\item{norm_X}{ensure a normalization of the data? Defaults to \code{FALSE}.}

\item{kernel}{kernel employed: \code{1} for von Mises--Fisher (default);
\code{2} for Epanechnikov; \code{3} for softplus.}

\item{kernel_type}{type of kernel employed: \code{1} for product kernel
(default); \code{2} for spherically symmetric kernel.}

\item{k}{softplus kernel parameter. Defaults to \code{10.0}.}
}
\value{
A column vector of size \code{c(n, 1)} with the evaluation of the
logarithm of the cross-validated kernel density estimator.
}
\description{
Computes the logarithm of the cross-validated kernel density
estimator: \eqn{\log \hat{f}_{-i}(\boldsymbol{X}_i;\boldsymbol{h})},
\eqn{i = 1, \ldots, n.}
}
\examples{
# Simple check on S^1 x S^2
n <- 5
d <- c(1, 2)
h <- c(0.2, 0.2)
X <- r_unif_polysph(n = n, d = d)
log_cv_kde_polysph(X = X, d = d, h = h)
kde_polysph(x = X[1, , drop = FALSE], X = X[-1, ], d = d, h = h, log = TRUE)
}
