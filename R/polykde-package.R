
#' @title \code{polykde}: Polyspherical Kernel Density Estimation
#'
#' @description Kernel density estimation on the polysphere, (hyper)sphere, and
#' circle. Includes functions for density estimation, regression estimation,
#' ridge estimation, bandwidth selection, kernels, samplers, and homogeneity
#' tests. Companion package to García-Portugués and Meilán-Vila (2025)
#' <\doi{10.1080/01621459.2025.2521898}> and García-Portugués and Meilán-Vila
#' (2023) <\doi{10.1007/978-3-031-32729-2_4}>.
#'
#' @author Eduardo García-Portugués.
#' @references
#' García-Portugués, E. and Meilán-Vila, A. (2025). Kernel density estimation
#' with polyspherical data and its applications. \emph{Journal of the American
#' Statistical Association}, to appear. \doi{10.1080/01621459.2025.2521898}.
#'
#' García-Portugués, E. and Meilán-Vila, A. (2023). Hippocampus shape analysis
#' via skeletal models and kernel smoothing. In Larriba, Y. (Ed.),
#' \emph{Statistical Methods at the Forefront of Biomedical Advances},
#' pp. 63--82. Springer, Cham. \doi{10.1007/978-3-031-32729-2_4}.
#' @name polykde-package
#' @import graphics Rcpp RcppProgress stats
#' @importFrom utils setTxtProgressBar txtProgressBar
#' @useDynLib polykde, .registration = TRUE
#' @aliases polykde polykde-package
"_PACKAGE"
