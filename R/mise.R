
#' @title Computation of the exact MISE for the polyspherical kernel density
#' estimator with von Mises--Fisher kernel for mixtures of von Mises--Fisher
#' distributions
#'
#' @description The function \code{exact_mise_vmf_polysph} computes the
#' exact Mean Integrated Squared Error (MISE) of the kernel density estimator
#' on the polysphere \eqn{\mathcal{S}^{d_1} \times \cdots \times
#' \mathcal{S}^{d_r}} with respect to a density \eqn{f_m(\boldsymbol{x}) =
#' \sum_{j=1}^m p_j
#' f_{\mathrm{vMF}}(\boldsymbol{x}; \boldsymbol{\mu}_j, \kappa_j)} of an
#' \eqn{m}-mixture of von Mises--Fisher distributions. The MISE is
#' \deqn{\mathrm{MISE}_m[\hat{f}(\cdot;\boldsymbol{h})]=\mathbb{E}\left[
#' \int_{\mathcal{S}^{d_1} \times \cdots \times \mathcal{S}^{d_r}}
#' \left(\hat{f}(\boldsymbol{x};\boldsymbol{h})-f_m(\boldsymbol{x})\right)^2
#' \,\mathrm{d}\boldsymbol{x}\right]}
#' and can be computed exactly from Propositions 4 and 5 in García-Portugués
#' et al. (2013), using importance-sampling Monte Carlo to evaluate the matrices
#' \eqn{\boldsymbol{\Psi_1}} and \eqn{\boldsymbol{\Psi_2}}.
#'
#' @param h,log_h matrix of size \code{c(k, r)} with \code{k} vectors of
#' (log)bandwidths.
#' @param n sample size. Either a scalar or a vector of size \code{k}.
#' @param mu a matrix of size \code{c(m, sum(d + 1))} with the means of
#' each mixture components in the rows.
#' @param kappa a matrix of size \code{c(m, r)} with the concentrations of
#' each mixture components in the rows.
#' @param prop a vector of size \code{m} with the proportions of the mixture
#' components.
#' @inheritParams kde_polysph
#' @param M_psi number of importance-sampling Monte Carlo samples. Defaults to
#' \code{1e4}.
#' @param seed_psi seed for approximating the matrices \eqn{\boldsymbol{\Psi_1}}
#' and \eqn{\boldsymbol{\Psi_2}}. Defaults to \code{NULL} (no seed is fixed).
#' @inheritParams log_besselI_scaled
#' @return For \code{exact_mise_vmf*}, a list with the following components:
#' \item{mise}{vector of size \code{k} with the evaluated MISEs.}
#' \item{Psi_0}{matrix \eqn{\boldsymbol{\Psi_0}}.}
#' \item{Psi_1}{matrix \eqn{\boldsymbol{\Psi_1}}.}
#' \item{Psi_2}{matrix \eqn{\boldsymbol{\Psi_2}}.}
#' \item{log_Dq_h}{vector of size \code{k} with normalizing constant.}
#' For \code{log1p_mise_exact}, the \code{log(mise)} is returned.
#' @references
#' García-Portugués, E., Crujeiras, R. M., and González-Manteiga, W. (2013).
#' Kernel density estimation for directional-linear data. \emph{Journal of
#' Multivariate Analysis}, 121:152--175. \doi{10.1016/j.jmva.2013.06.009}.
#' @examples
#' h <- seq(0.01, 1, l = 100)
#' log_mise_h <- log(polykde:::exact_mise_vmf(
#'   h = h, n = 100, mu = rbind(c(0, 1), c(1, 0)), kappa = c(5, 2),
#'   prop = c(0.7, 0.3), d = 1, spline = TRUE, seed_psi = 1)$mise)
#' plot(h, log_mise_h)
#' log_mise_h <- log(polykde:::exact_mise_vmf_polysph(
#'   h = cbind(h, h, h), n = 100,
#'   mu = rbind(c(0, 0, 1, 0, 1, 0, 1), c(0, 0, -1, 0, -1, 0, 1)),
#'   kappa = rbind(1:3, 3:1), prop = c(0.7, 0.3), d = c(2, 1, 1),
#'   spline = TRUE, seed_psi = 1)$mise)
#' plot(h, log_mise_h)
#' @name exact_mise


#' @noRd
#' @rdname exact_mise
exact_mise_vmf_polysph <- function(h, n, mu, kappa, prop, d, M_psi = 1e4,
                                   seed_psi = NULL, spline = FALSE) {

  # Check mixture inputs
  r <- length(d)
  mu <- rbind(mu)
  kappa <- cbind(kappa)
  m <- length(prop)
  h <- rbind(h)
  stopifnot(nrow(kappa) == m)
  stopifnot(ncol(kappa) == r)
  stopifnot(ncol(mu) == sum(d + 1))
  stopifnot(ncol(h) == r)
  if (length(n) == 1) {

    n <- rep(n, nrow(h))

  } else {

    if (length(n) != nrow(h)) {

      stop("n must be a scalar or have the same number of rows as h.")

    }

  }

  # Loop on spheres
  ind_dj <- comp_ind_dj(d)
  log_Dq_h <- 0
  Psi_0 <- matrix(1, nrow = m, ncol = m)
  Psi_1 <- Psi_2 <- array(1, dim = c(nrow(h), m, m))
  for (j in seq_len(r)) {

    mise_j <- exact_mise_vmf(h = h[, j], n = n,
                             mu = mu[, (ind_dj[j] + 1):ind_dj[j + 1]],
                             kappa = kappa[, j], prop = prop, d = d[j],
                             M_psi = M_psi, seed_psi = ifelse(is.null(seed_psi),
                                                              NULL,
                                                              seed_psi + j - 1),
                             spline = spline)

    # Exploiting structure in Proposition 5 of "Kernel density estimation for
    # directional-linear data" (https://doi.org/10.1016/j.jmva.2013.06.009)
    log_Dq_h <- log_Dq_h + mise_j$log_Dq_h
    Psi_0 <- Psi_0 * mise_j$Psi_0
    Psi_1 <- Psi_1 * mise_j$Psi_1
    Psi_2 <- Psi_2 * mise_j$Psi_2

  }

  # MISE
  mise <- numeric(max(nrow(h), length(n)))
  for (i in seq_len(nrow(h))) {

    mise[i] <- exp(-(log_Dq_h[i] + log(n[i]))) +
      drop(t(prop) %*% ((1 - 1 / n[i]) * Psi_2[i, , ] -
                          2 * Psi_1[i, , ] + Psi_0) %*% prop)

  }
  return(list("mise" = mise, "Psi_0" = Psi_0, "Psi_1" = Psi_1, "Psi_2" = Psi_2,
              "log_Dq_h" = unname(log_Dq_h)))

}


#' @noRd
#' @rdname exact_mise
exact_mise_vmf <- function(h, n, mu, kappa, prop, d, M_psi = 1e4,
                           seed_psi = NULL, spline = FALSE) {

  # Check mixture inputs
  m <- length(prop)
  mu <- rbind(mu)
  stopifnot(abs(sum(prop) - 1) < 1e-15)
  stopifnot(ncol(mu) == d + 1)
  stopifnot(nrow(mu) == m)
  stopifnot(length(kappa) == m)
  stopifnot(length(d) == 1)
  if (length(n) == 1) {

    n <- rep(n, length(h))

  } else {

    if (length(n) != length(h)) {

      stop("n must be a scalar or have the same length as h.")

    }

  }

  # mu_i * kappa_i
  mu_kappa <- mu * kappa

  # Psi_0
  # ||kappa_i mu_i + kappa_j mu_j||^2 = kappa_i^2 + kappa_j^2
  #                                    + 2 * kappa_i * kappa_j * mu_i'mu_j
  log_C_kappa <- fast_log_c_vMF(p = d + 1, kappa = kappa, spline = spline)
  kappa_mu_i_kappa_mu_j <-
    sqrt(outer(kappa^2, kappa^2, "+") + 2 * tcrossprod(mu_kappa))
  log_C_kappa_mu_i_kappa_mu_j <-
    fast_log_c_vMF(p = d + 1, kappa = kappa_mu_i_kappa_mu_j, spline = spline)
  Psi_0 <- exp(outer(log_C_kappa, log_C_kappa, "+") -
                 log_C_kappa_mu_i_kappa_mu_j)

  # Set seeds for the Monte Carlos in the Psi_1 and Psi_2
  if (!is.null(seed_psi)) {

    # old_seed <- .Random.seed
    # on.exit({.Random.seed <<- old_seed})
    set.seed(seed_psi, kind = "Mersenne-Twister")

  }

  # Psi_1 and Psi_2 using weighted MC
  Psi_1 <- Psi_2 <- array(NA, dim  = c(length(h), m, m))
  log_C_h <- fast_log_c_vMF(p = d + 1, kappa = 1 / h^2, spline = spline)
  for (j in seq_len(m)) {

    # Sample vMF(mu[j], kappa[j])
    vmf_samp_j <- rotasym::r_vMF(n = M_psi, mu = mu[j, ], kappa = kappa[j])

    for (i in seq_len(j)) {

      # Sample vMF(mu[i], kappa[i])
      vmf_samp_i <- rotasym::r_vMF(n = M_psi, mu = mu[i, ], kappa = kappa[i])

      # Combined sample for estimation based on mixture
      vmf_samp_ij <- rbind(vmf_samp_i, vmf_samp_j)

      # Evaluate mixture
      log_dens_ij <- drop(log(0.5 * (
        rotasym::d_vMF(x = vmf_samp_ij, mu = mu[i, ], kappa = kappa[i]) +
          rotasym::d_vMF(x = vmf_samp_ij, mu = mu[j, ], kappa = kappa[j]))))

      # Compute C(||x/h^2 + kappa_i mu_i||) and C(||x/h^2 + kappa_j mu_j||)
      h_x_kappa_mu_i <- h_x_kappa_mu_j <-
        matrix(0, nrow = length(h), ncol = 2 * M_psi)
      for (k in seq_along(h)) {

        h_x_kappa_mu_i[k, ] <- sqrt(colSums((t(vmf_samp_ij / h[k]^2) +
                                               kappa[i] * mu[i, ])^2))
        h_x_kappa_mu_j[k, ] <- sqrt(colSums((t(vmf_samp_ij / h[k]^2) +
                                               kappa[j] * mu[j, ])^2))

      }
      log_C_h_x_kappa_mu_i <- fast_log_c_vMF(p = d + 1, kappa = h_x_kappa_mu_i,
                                             spline = spline)
      log_C_h_x_kappa_mu_j <- fast_log_c_vMF(p = d + 1, kappa = h_x_kappa_mu_j,
                                             spline = spline)

      # Common terms to Psi_1 and Psi_2
      pre_Psi_12 <- log_C_h + log_C_kappa[i] + log_C_kappa[j] -
        log_C_h_x_kappa_mu_i
      pre_Psi_12 <- t(t(pre_Psi_12) - log_dens_ij)
      kappa[j] * drop(vmf_samp_ij %*% mu[j, ])

      # Psi_1
      Psi_1[, i, j] <-
        rowMeans(exp(t(t(pre_Psi_12) +
                         kappa[j] * drop(vmf_samp_ij %*% mu[j, ]))))

      # Psi_2
      Psi_2[, i, j] <-
        rowMeans(exp(pre_Psi_12 + log_C_h - log_C_h_x_kappa_mu_j))

      # Symmetrize Psi_1 and Psi_2
      Psi_1[, j, i] <- Psi_1[, i, j]
      Psi_2[, j, i] <- Psi_2[, i, j]

    }

  }

  # Dq constant (careful: it is the inverse of what is reported in the paper!)
  log_C_h <- fast_log_c_vMF(p = d + 1, kappa = 1 / h^2, spline = spline)
  log_Dq_h <- -2 * log_C_h +
    fast_log_c_vMF(p = d + 1, kappa = 2 / h^2, spline = spline)

  # MISE in Proposition 4 of "Kernel density estimation for
  # directional-linear data" (https://doi.org/10.1016/j.jmva.2013.06.009)
  mise <- numeric(max(length(h), length(n)))
  for (i in seq_along(h)) {

    mise[i] <- exp(-(log_Dq_h[i] + log(n[i]))) +
      drop(t(prop) %*% ((1 - 1 / n[i]) * Psi_2[i, , ] -
                          2 * Psi_1[i, , ] + Psi_0) %*% prop)

  }
  return(list("mise" = mise, "Psi_0" = Psi_0, "Psi_1" = Psi_1, "Psi_2" = Psi_2,
              "log_Dq_h" = unname(log_Dq_h)))

}


#' @noRd
#' @rdname exact_mise
log1p_mise_exact <- function(log_h, n, d, mu, kappa, prop, M_psi = 1e4,
                             seed_psi = NULL, spline = FALSE) {

  log1p(exact_mise_vmf_polysph(h = exp(log_h), n = n, mu = mu, kappa = kappa,
                               prop = prop, d = d, M_psi = M_psi,
                               seed_psi = seed_psi, spline = spline)$mise)

}


#' @title Exact MISE bandwidth selection for polyspherical kernel density
#' estimator for a mixture of von Mises--Fisher distributions
#'
#' @description Computes the bandwidths
#' \deqn{\boldsymbol{h}_{\mathrm{MISE}}=\arg\min_{\boldsymbol{h}>0}
#' \mathrm{MISE}_m[\hat{f}(\cdot;\boldsymbol{h})]}
#' of the polyspherical kernel density estimator
#' \eqn{\mathcal{S}^{d_1} \times \cdots \times \mathcal{S}^{d_r}} for an
#' \eqn{m}-mixture of von Mises--Fisher densities
#' \eqn{f_m(\boldsymbol{x}) = \sum_{j=1}^m p_j
#' f_{\mathrm{vMF}}(\boldsymbol{x}; \boldsymbol{\mu}_j, \kappa_j)}.
#'
#' @inheritParams kde_polysph
#' @inheritParams exact_mise
#' @param bw0 initial bandwidth vector for minimizing the MISE loss. Can be
#' also a matrix of initial bandwidth vectors.
#' @inheritParams bw_rot_polysph
#' @inheritParams log_besselI_scaled
#' @param ... further arguments passed to \code{\link{nlm}}.
#' @return A list with entries \code{bw} (optimal bandwidth) and \code{opt},
#' the latter containing the output of \code{\link[stats]{nlm}}.
#' @examples
#' polykde:::bw_mise_polysph(n = 100, d = 2, bw0 = 1,
#'                           mu = rbind(c(1, 0, 0)), kappa = rbind(1),
#'                           prop = 1, seed_psi = 1)
#' polykde:::bw_mise_polysph(n = 1000, d = 2, bw0 = 1,
#'                           mu = rbind(c(1, 0, 0)), kappa = rbind(1),
#'                           prop = 1, seed_psi = 1)
#' @noRd
bw_mise_polysph <- function(n, d, bw0 = NULL, upscale = FALSE, deriv = 0,
                            mu = NULL, kappa = NULL, prop = NULL,
                            M_psi = 1e4, seed_psi = NULL, spline = FALSE, ...) {

  # Get r
  r <- length(d)

  # Initial search?
  if (is.null(bw0)) {

    stop("bw0 must be provided.")
    # bw0 <- bw_mrot_polysph(X = X, d = d, kernel = 1, upscale = FALSE
    #                       deriv = 0, mu = mu, kappa = kappa, prop = prop)
    # bw0 <- bw0$bw

  } else {

    bw0 <- rbind(bw0)
    if (r != ncol(bw0)) {

      stop("bw0 and d are incompatible.")

    }
    if (nrow(bw0) > 1) {

      # Name arguments in apply() to avoid mismatches between MARGIN and M_psi
      mise_bw0 <- apply(X = log(bw0), MARGIN = 1, FUN = log1p_mise_exact,
                        n = n, d = d, mu = mu, kappa = kappa, prop = prop,
                        M_psi = M_psi, seed_psi = seed_psi, spline = spline)
      ind_best <- which.min(mise_bw0)
      bw0 <- bw0[ind_best, ]

    }

  }

  # Search h_MISE
  opt <- nlm(f = log1p_mise_exact, p = log(bw0), n = n, d = d,
             mu = mu, kappa = kappa, prop = prop, M_psi = M_psi,
             seed_psi = seed_psi, spline = spline, ...)
  bw <- exp(opt$estimate)

  # Upscale?
  if (upscale > 0) {

    n_up <- n^(1 / (d * r + 4)) * n^(-1 / (d * r + 2 * deriv + 4))
    bw <- bw * n_up

  }
  return(list("bw" = unname(bw), "opt" = opt))

}
