% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/samplers.R
\name{r_kde_polysph}
\alias{r_kde_polysph}
\title{Sample from polyspherical kernel density estimator}
\usage{
r_kde_polysph(n, X, d, h, kernel = 1, kernel_type = 1, k = 10,
  intrinsic = FALSE, norm_X = FALSE)
}
\arguments{
\item{n}{sample size.}

\item{X}{a matrix of size \code{c(n, sum(d) + r)} with the sample.}

\item{d}{vector of size \code{r} with dimensions.}

\item{h}{vector of size \code{r} with bandwidths.}

\item{kernel}{kernel employed: \code{1} for von Mises--Fisher (default);
\code{2} for Epanechnikov; \code{3} for softplus.}

\item{kernel_type}{type of kernel employed: \code{1} for product kernel
(default); \code{2} for spherically symmetric kernel.}

\item{k}{softplus kernel parameter. Defaults to \code{10.0}.}

\item{intrinsic}{use the intrinsic distance, instead of the
extrinsic-chordal distance, in the kernel? Defaults to \code{FALSE}.}

\item{norm_X}{ensure a normalization of the data?}
}
\value{
A matrix of size \code{c(n, sum(d) + r)} with the sample.
}
\description{
Simulates from the distribution defined by a polyspherical
kernel density estimator on \eqn{\mathcal{S}^{d_1} \times \ldots \times
\mathcal{S}^{d_r}}.
}
\details{
The function uses \code{\link{r_kern_polysph}} to sample from the
considered kernel.
}
\examples{
# Simulated data on (S^1)^2
n <- 50
samp <- r_path_s1r(n = n, r = 2, k = c(1, 2), angles = TRUE)
plot(samp, xlim = c(-pi, pi), ylim = c(-pi, pi), col = rainbow(n),
     axes = FALSE, xlab = "", ylab = "", pch = 16, cex = 0.75)
points(torus_to_angles(r_kde_polysph(n = 10 * n, X = angles_to_torus(samp),
                                     d = c(1, 1), h = c(0.1, 0.1))),
       col = "black", pch = 16, cex = 0.2)
sdetorus::torusAxis()

# Simulated data on S^2
n <- 50
samp <- r_path_s2r(n = n, r = 1, sigma = 0.1, kappa = 5,
                   spiral = TRUE)[, , 1]
sc3d <- scatterplot3d::scatterplot3d(
  samp, xlim = c(-1, 1), ylim = c(-1, 1), zlim = c(-1, 1),
  xlab = "", ylab = "", zlab = "", color = rainbow(n), pch = 16
)
xyz <- r_kde_polysph(n = 10 * n, X = samp, d = 2, h = 0.1)
sc3d$points3d(xyz[, 1], xyz[, 2], xyz[, 3], col = "black", pch = 16,
              cex = 0.2)
}
