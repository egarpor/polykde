% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils.R
\name{softplus}
\alias{softplus}
\title{Stable computation of the softplus function}
\usage{
softplus(t)
}
\arguments{
\item{t}{vector or matrix.}
}
\value{
The softplus function evaluated at \code{t}.
}
\description{
Computes the softplus function \eqn{\log(1+e^{t})} in a
numerically stable way for large absolute values of \eqn{t}.
}
\examples{
curve(softplus(10 * (1 - (1 - x) / 0.1)), from = -1, to = 1)
}
