

#' @title Density of the product of von Mises--Fisher distributions on the
#' polysphere
#'
#' @description Computes the density of the product of von Mises--Fisher
#' densities on the polysphere.
#'
#' @inheritParams kde_polysph
#' @inheritParams r_vmf_polysph
#' @return A vector of size \code{nx} with the evaluated density.
#' @examples
#' # Simple check of integration on S^1 x S^2
#' d <- c(1, 2)
#' mu <- c(0, 1, 0, 1, 0)
#' kappa <- c(1, 1)
#' x <- r_vmf_polysph(n = 1e4, d = d, mu = mu, kappa = kappa)
#' mean(1 / d_vmf_polysph(x = x, d = d, mu = mu, kappa = kappa)) /
#'   prod(rotasym::w_p(p = d + 1))
#' @export
d_vmf_polysph <- function(x, d, mu, kappa, log = FALSE) {

  kde_polysph(x = x, X = rbind(mu), d = d, h = 1 / sqrt(kappa), kernel = 1,
              wrt_unif = FALSE, norm_x = TRUE, norm_X = TRUE, log = log)

}


#' @title Density of the uniform distribution on the polysphere
#'
#' @description Computes the density of the uniform distribution on the
#' polysphere.
#'
#' @inheritParams kde_polysph
#' @return A vector of size \code{nx} with the evaluated density.
#' @examples
#' # Simple check of integration on S^1 x S^2
#' d <- c(1, 2)
#' x <- r_unif_polysph(n = 1e4, d = d)
#' mean(1 / d_unif_polysph(x = x, d = d)) / prod(rotasym::w_p(p = d + 1))
#' @export
d_unif_polysph <- function(x, d, log = FALSE) {

  if (is.null(dim(x))) {

    x <- rbind(x)

  }
  stopifnot(ncol(x) == sum(d + 1))
  log_dens <- rep(-sum(rotasym::w_p(p = d + 1, log = TRUE)), nrow(x))
  if (!log) {

    log_dens <- exp(log_dens)

  }
  return(log_dens)

}
