
#' @title \code{polykde}: Polyspherical Kernel Density Estimation
#'
#' @description Tools for kernel density estimation on the polysphere,
#' hypersphere, and circle. Includes functions for density estimation,
#' regression estimation, ridge estimation, bandwidth selection, kernels,
#' samplers, and homogeneity tests. Companion package to García-Portugués and
#' Meilán-Vila (2024) <\doi{10.48550/arXiv.2411.04166}> and
#' and García-Portugués and Meilán-Vila (2023)
#' <\doi{10.1007/978-3-031-32729-2_4}>.
#'
#' @author Eduardo García-Portugués.
#' @references
#' García-Portugués, E. and Meilán-Vila, A. (2023). Hippocampus shape analysis
#' via skeletal models and kernel smoothing. In Larriba, Y. (Ed.),
#' \emph{Statistical Methods at the Forefront of Biomedical Advances},
#' pp. 63--82. Springer, Cham. \doi{10.1007/978-3-031-32729-2_4}.
#'
#' García-Portugués, E. and Meilán-Vila, A. (2024). Kernel density estimation
#' with polyspherical data and its applications. \emph{arXiv:2411.04166}.
#' \doi{10.48550/arXiv.2411.04166}.
#' @name polykde-package
#' @import graphics Rcpp RcppProgress stats
#' @importFrom utils setTxtProgressBar txtProgressBar
#' @useDynLib polykde, .registration = TRUE
#' @aliases polykde polykde-package
"_PACKAGE"
