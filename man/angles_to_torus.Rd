% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/angles.R
\name{angles_to_torus}
\alias{angles_to_torus}
\alias{torus_to_angles}
\title{Conversion between the angular and Cartesian coordinates of the torus}
\usage{
angles_to_torus(theta)

torus_to_angles(x)
}
\arguments{
\item{theta}{matrix of size \code{c(n, d)} with the angles.}

\item{x}{matrix of size \code{c(n, 2 * d)} with the Cartesian coordinates on
\eqn{(\mathcal{S}^1)^d}. Assumed to be of unit norm by pairs of coordinates
in the rows.}
}
\value{
\itemize{
\item{\code{angles_to_torus}: the matrix \code{x}.}
\item{\code{torus_to_angles}: the matrix \code{theta}.}
}
}
\description{
Transforms the angles \eqn{(\theta_1,\ldots,\theta_d)} in
\eqn{[-\pi,\pi)^d} into the Cartesian coordinates
\deqn{(\cos(x_1), \sin(x_1),\ldots,\cos(x_d), \sin(x_d))}
of the torus \eqn{(\mathcal{S}^1)^d}, and vice versa.
}
\examples{
# Check changes of coordinates
torus_to_angles(angles_to_torus(c(0, pi / 3, pi / 2)))
torus_to_angles(angles_to_torus(rbind(c(0, pi / 3, pi / 2), c(0, 1, -2))))
angles_to_torus(torus_to_angles(c(0, 1, 1, 0)))
angles_to_torus(torus_to_angles(rbind(c(0, 1, 1, 0), c(0, 1, 0, 1))))
}
