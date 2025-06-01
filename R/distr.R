

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


#' @title Density of the product of von Mises--Fisher distributions on the
#' polysphere
#'
#' @description Computes the density of an \eqn{m}-mixture of product of von
#' Mises--Fisher densities on the polysphere.
#'
#' @inheritParams kde_polysph
#' @inheritParams r_mvmf_polysph
#' @return A vector of size \code{nx} with the evaluated density.
#' @examples
#' # Simple check of integration on S^1 x S^2
#' d <- c(1, 2)
#' mu <- rbind(c(0, 1, 0, 1, 0), c(1, 0, 1, 0, 0))
#' kappa <- rbind(c(5, 2), c(1, 2))
#' prop <- c(0.7, 0.3)
#' x <- r_mvmf_polysph(n = 1e4, d = d, mu = mu, kappa = kappa, prop = prop)
#' mean(1 / d_mvmf_polysph(x = x, d = d, mu = mu, kappa = kappa, prop = prop)) /
#'   prod(rotasym::w_p(p = d + 1))
#' @export
d_mvmf_polysph <- function(x, d, mu, kappa, prop, log = FALSE) {

  # Check mixture inputs
  r <- length(d)
  m <- length(prop)
  mu <- rbind(mu)
  kappa <- cbind(kappa)
  if (nrow(kappa) != m || ncol(kappa) != r) {

    stop("kappa size is not c(m, r).")

  }
  if (nrow(mu) != m || ncol(mu) != sum(d + 1)) {

    stop("mu size is not c(m, sum(d + 1)).")

  }
  if (abs(sum(prop) - 1) > 1e-15) {

    stop("prop does not add to one.")

  }

  # Density computation
  dens <- matrix(0, nrow = nrow(x), ncol = 1)
  for (k in seq_len(m)) {

    dens <- dens +
      prop[k] * kde_polysph(x = x, X = mu[k, , drop = FALSE], d = d,
                            h = 1 / sqrt(kappa[k, ]), kernel = 1,
                            wrt_unif = FALSE, norm_x = TRUE, norm_X = TRUE,
                            log = FALSE)

  }
  if (log) {

    dens <- log(dens)

  }
  return(dens)

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


#' @title Fast evaluation of the von Mises--Fisher normalizing constant
#'
#' @description Computes the normalizing constant of the von Mises--Fisher
#' distribution on the sphere \eqn{\mathcal{S}^{p-1}} as in
#' \code{\link[rotasym]{c_vMF}} but using a fast spline approximation for the
#' logarithm of (a limited number of) Bessel functions.
#'
#' @param p positive integer giving the dimension of the ambient space that
#' contains the sphere \eqn{\mathcal{S}^{p-1}}. Can be a vector only for
#' \code{spline = FALSE}.
#' @param kappa concentration parameter of the von Mises--Fisher distribution.
#' Can be a vector.
#' @inheritParams log_besselI_scaled
#' @return A vector of the size of \code{kappa} with the normalizing constants.
#' @examples
#' polykde:::fast_log_c_vMF(p = 3, kappa = 1:10, spline = FALSE)
#' polykde:::fast_log_c_vMF(p = 3, kappa = 1:10, spline = TRUE)
#' @noRd
fast_log_c_vMF <- function(p, kappa, spline = FALSE) {

  log_bessel <- log_besselI_scaled(nu = 0.5 * (p - 2), x = kappa,
                                   spline = spline)
  log_c_vMF <- (0.5 * (p - 2)) * log(kappa) - (0.5 * p) * log(2 * pi) -
    kappa - log_bessel
  log_c_vMF[kappa == 0] <- -rotasym::w_p(p = p, log = TRUE)
  return(log_c_vMF)

}
