

#' @title \code{polykde}: Polyspherical Kernel Density Estimation
#'
#' @description Tools for kernel density estimation on the polysphere,
#' hypersphere, and circle. Includes functions for density estimation,
#' regression estimation, ridge estimation, bandwidth selection, kernels,
#' samplers, and homogeneity tests. Companion package to García-Portugués and
#' Meilán-Vila (2024) <\doi{10.48550/arXiv.2411.04166}>.
#'
#' @author Eduardo García-Portugués.
#' @references
#' García-Portugués, E. and Meilán-Vila, A. (2024). Kernel density estimation
#' with polyspherical data and its applications. \emph{arXiv:2411.04166}.
#' \doi{10.48550/arXiv.2411.04166}.
#' @name polykde-package
#' @import graphics Rcpp RcppProgress stats
#' @importFrom utils setTxtProgressBar txtProgressBar
#' @useDynLib polykde, .registration = TRUE
#' @aliases polykde polykde-package
"_PACKAGE"
