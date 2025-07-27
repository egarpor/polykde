
#' @title Computation of the exact MISE for the polyspherical kernel density
#' estimator with von Mises--Fisher kernel for mixtures of von Mises--Fisher
#' distributions
#'
#' @description The function \code{mise_vmf_polysph} computes the
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
#' @inheritParams r_mvmf_polysph
#' @inheritParams kde_polysph
#' @param M_psi number of importance-sampling Monte Carlo samples. Defaults to
#' \code{1e4}.
#' @param seed_psi seed for approximating the matrices \eqn{\boldsymbol{\Psi_1}}
#' and \eqn{\boldsymbol{\Psi_2}}. Defaults to \code{NULL} (no seed is fixed).
#' @inheritParams log_besselI_scaled
#' @return For \code{mise_vmf*}, a list with the following components:
#' \item{mise}{vector of size \code{k} with the evaluated MISEs.}
#' \item{Psi_0}{matrix \eqn{\boldsymbol{\Psi_0}}.}
#' \item{Psi_1}{matrix \eqn{\boldsymbol{\Psi_1}}.}
#' \item{Psi_2}{matrix \eqn{\boldsymbol{\Psi_2}}.}
#' \item{log_Dq_h}{vector of size \code{k} with normalizing constant.}
#' For \code{log1p_mise}, the \code{log(mise)} is returned.
#' @references
#' García-Portugués, E., Crujeiras, R. M., and González-Manteiga, W. (2013).
#' Kernel density estimation for directional-linear data. \emph{Journal of
#' Multivariate Analysis}, 121:152--175. \doi{10.1016/j.jmva.2013.06.009}.
#' @examples
#' h <- seq(0.01, 1, l = 100)
#' log_mise_h <- log(polykde:::mise_vmf(
#'   h = h, n = 100, mu = rbind(c(0, 1), c(1, 0)), kappa = c(5, 2),
#'   prop = c(0.7, 0.3), d = 1, spline = TRUE, seed_psi = 1)$mise)
#' plot(h, log_mise_h)
#' log_mise_h <- log(polykde:::mise_vmf_polysph(
#'   h = cbind(h, h, h), n = 100,
#'   mu = rbind(c(0, 0, 1, 0, 1, 0, 1), c(0, 0, -1, 0, -1, 0, 1)),
#'   kappa = rbind(1:3, 3:1), prop = c(0.7, 0.3), d = c(2, 1, 1),
#'   spline = TRUE, seed_psi = 1)$mise)
#' plot(h, log_mise_h)
#' @name mise


#' @noRd
#' @rdname mise
mise_vmf_polysph <- function(h, n, mu, kappa, prop, d, M_psi = 1e4,
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
  stopifnot(abs(sum(prop) - 1) < 1e-15)
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

    mise_j <- mise_vmf(h = h[, j], n = n,
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
#' @rdname mise
mise_vmf <- function(h, n, mu, kappa, prop, d, M_psi = 1e4, seed_psi = NULL,
                     spline = FALSE) {

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

  # Psi_1 and Psi_2 using importance-sampling Monte Carlo
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
#' @rdname mise
log1p_mise <- function(log_h, n, d, mu, kappa, prop, M_psi = 1e4,
                       seed_psi = NULL, spline = FALSE) {

  log1p(mise_vmf_polysph(h = exp(log_h), n = n, mu = mu, kappa = kappa,
                         prop = prop, d = d, M_psi = M_psi,
                         seed_psi = seed_psi, spline = spline)$mise)

}


#' @title Exact MISE bandwidth selection for polyspherical kernel density
#' estimator for a mixture of von Mises--Fisher distributions
#'
#' @description Computes the bandwidths
#' \deqn{\boldsymbol{h}_{\mathrm{MISE}}=\arg\min_{\boldsymbol{h}>0}
#' \mathrm{MISE}_m[\hat{f}(\cdot;\boldsymbol{h})]}
#' of the kernel density estimator on the polysphere
#' \eqn{\mathcal{S}^{d_1} \times \cdots \times \mathcal{S}^{d_r}} for an
#' \eqn{m}-mixture of von Mises--Fisher densities
#' \eqn{f_m(\boldsymbol{x}) = \sum_{j=1}^m p_j
#' f_{\mathrm{vMF}}(\boldsymbol{x}; \boldsymbol{\mu}_j, \kappa_j)}.
#'
#' @inheritParams kde_polysph
#' @inheritParams mise_vmf_polysph
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
bw_mise_polysph <- function(n, d, bw0 = NULL, mu, kappa, prop, M_psi = 1e4,
                            seed_psi = NULL, spline = FALSE, ...) {

  # Get r
  r <- length(d)

  # log1p_mise with fixed arguments
  log1p_mise_args <- function(log_h) {

    log1p_mise(log_h = log_h, n = n, d = d, mu = mu, kappa = kappa, prop = prop,
               M_psi = M_psi, seed_psi = seed_psi, spline = spline)

  }

  # Initial search
  if (is.null(bw0)) {

    stop("bw0 must be provided.")

  } else {

    bw0 <- rbind(bw0)
    if (r != ncol(bw0)) {

      stop("bw0 and d are incompatible.")

    }
    if (nrow(bw0) > 1) {

      # Name arguments in apply() to avoid mismatches between MARGIN and M_psi
      mise_bw0 <- apply(X = log(bw0), MARGIN = 1, FUN = log1p_mise_args)
      ind_best <- which.min(mise_bw0)
      bw0 <- bw0[ind_best, ]

    }

  }

  # Search h_MISE
  opt <- nlm(f = log1p_mise_args, p = log(bw0), ...)
  bw <- exp(opt$estimate)
  return(list("bw" = unname(bw), "opt" = opt))

}


#' @title Computation of the exact ISE for the polyspherical kernel density
#' estimator with von Mises--Fisher kernel for mixtures of von Mises--Fisher
#' distributions
#'
#' @description Computes the exact Integrated Squared Error (ISE) of the kernel
#' density estimator on the polysphere
#' \eqn{\mathcal{S}^{d_1} \times \cdots \times \mathcal{S}^{d_r}} with respect
#' to a density \eqn{f_m(\boldsymbol{x}) = \sum_{j=1}^m p_j
#' f_{\mathrm{vMF}}(\boldsymbol{x}; \boldsymbol{\mu}_j, \kappa_j)} of an
#' \eqn{m}-mixture of von Mises--Fisher distributions. The ISE is
#' \deqn{\mathrm{ISE}_m[\hat{f}(\cdot;\boldsymbol{h})]
#' =\|\hat{f}(\cdot;\boldsymbol{h})-f_m\|_2^2
#' =\int_{\mathcal{S}^{d_1} \times \cdots \times \mathcal{S}^{d_r}}
#' \left(\hat{f}(\boldsymbol{x};\boldsymbol{h})-f_m(\boldsymbol{x})\right)^2
#' \,\mathrm{d}\boldsymbol{x}.}
#'
#' @param h,log_h matrix of size \code{c(k, r)} with \code{k} vectors of
#' (log)bandwidths.
#' @inheritParams mise
#' @inheritParams r_mvmf_polysph
#' @param x_mvmf Monte Carlo sample of the mixture to conduct importance
#' sampling. Computed internally if \code{NULL} (default).
#' @param f_mvmf density evaluation on \code{x_mvmf}. Computed internally if
#' \code{NULL} (default).
#' @param seed_psi seed for approximating the ISE if \code{exact = FALSE}.
#' Defaults to \code{NULL} (no seed is fixed).
#' @param exact use the exact von Mises--Fisher formula for the ISE? It
#' escalates quadratically on the sample size. Defaults to \code{FALSE}.
#' @param p power to compute the \eqn{\|\cdot\|_p} norm if \code{exact = FALSE}.
#' Defaults to \code{2}.
#' @return A list with the following fields:
#' \item{ise}{vector of size \code{k} with the evaluated ISEs.}
#' \item{Psi_01}{vector of size \code{k} with the first term of the ISE. If
#' \code{exact = FALSE}, it is \code{NULL}.}
#' \item{Psi_02}{scalar with the second term of the ISE. If
#' \code{exact = FALSE}, it is \code{NULL}.}
#' \item{Psi_03}{vector of size \code{k} with the third term of the ISE. If
#' \code{exact = FALSE}, it is \code{NULL}.}
#' For \code{log1p_ise}, the \code{log(ise)} is returned.
#' @examples
#' M <- 1e4
#' n <- 200
#' mu <-  rbind(c(1, 0, 0), c(-1, 0, 0))
#' kappa <- cbind(c(1, 1))
#' prop <- c(0.5, 0.5)
#' x_mvmf <- r_mvmf_polysph(n = M, d = 2, mu = mu, kappa = kappa, prop = prop)
#' f_mvmf <- d_mvmf_polysph(x = x_mvmf, d = 2, mu = mu, kappa = kappa,
#'                          prop = prop)
#' X <- r_mvmf_polysph(n = n, d = 2, mu = mu, kappa = kappa, prop = prop)
#' h <- 10^seq(-2, 1, l = 100)
#' plot(h, polykde:::ise_vmf_polysph(X = X, d = 2, h = h, x_mvmf = x_mvmf,
#'                                         f_mvmf = f_mvmf)$ise))
#' abline(v = polykde:::bw_ise_polysph(
#'   X = X, d = 2, bw0 = 1, x_mvmf = x_mvmf, f_mvmf = f_mvmf)$bw, col = 2)
#' @noRd
ise_vmf_polysph <- function(X, d, h, mu, kappa, prop, M_psi = 1e4,
                            x_mvmf = NULL, f_mvmf = NULL, seed_psi = NULL,
                            spline = FALSE, exact = FALSE, p = 2) {

  # Check mixture inputs
  h <- cbind(h)
  r <- length(d)
  stopifnot(ncol(h) == r)
  stopifnot(ncol(X) == sum(d + 1))
  if (exact || is.null(x_mvmf) || is.null(f_mvmf)) {

    m <- length(prop)
    mu <- rbind(mu)
    stopifnot(abs(sum(prop) - 1) < 1e-15)
    stopifnot(nrow(mu) == m)
    stopifnot(ncol(mu) == sum(d + 1))
    stopifnot(nrow(kappa) == m)
    stopifnot(ncol(kappa) == r)

  }

  # Exact computation for mixture of vMFs or importance-sampling Monte Carlo?
  if (exact) {

    # Only sphere
    if (r > 1) {

      stop("exact = TRUE not implemented yet for  r > 1.")

    }

    # Common terms: mu_i * kappa_i and 1 / h^2
    n <- nrow(X)
    kappa <- drop(kappa)
    mu_kappa <- mu * kappa
    sqrt_cross_X <- sqrt(2 * pmax(1 + tcrossprod(X), 0))
    cross_mu_X <- 2 * tcrossprod(mu_kappa, X)

    # Second Psi_0 term -- \sum_{i,j=1}^r p_i p_j \Psi_0(mu_i, mu_j)

    # ||kappa_i mu_i + kappa_j mu_j||^2 = kappa_i^2 + kappa_j^2
    #                                    + 2 * kappa_i * kappa_j * mu_i'mu_j
    log_C_kappa_prop <- log(prop) +
      fast_log_c_vMF(p = d + 1, kappa = kappa, spline = spline)
    kappa_mu_i_kappa_mu_j <- sqrt(outer(kappa^2, kappa^2, "+") +
                                    2 * tcrossprod(mu_kappa))
    log_C_kappa_mu_i_kappa_mu_j <-
      fast_log_c_vMF(p = d + 1, kappa = kappa_mu_i_kappa_mu_j, spline = spline)
    Psi_02 <- sum(exp(outer(log_C_kappa_prop, log_C_kappa_prop, "+") -
                        log_C_kappa_mu_i_kappa_mu_j))

    # First Psi_0 term -- n^{-2} \sum_{i,j=1}^n \Psi_0(X_i, X_j)

    # Loop on h
    psi_0_h <- function(h_k) {

      # ||h^{-2} X_i + h^{-2} X_j||^2 = 2 * h^{-4} * (1 - X_i'X_j)
      h2_X_i_X_j <- sqrt_cross_X / h_k^2
      log_C_h_X_i_h_X_j <- fast_log_c_vMF(p = d + 1, kappa = h2_X_i_X_j,
                                          spline = spline)
      log_C_h2 <- fast_log_c_vMF(p = d + 1, kappa = 1 / h_k^2, spline = spline)
      return(sum(exp(2 * log_C_h2 - log_C_h_X_i_h_X_j)))

    }
    Psi_01 <- sapply(h, psi_0_h) / n^2

    # Third Psi_0 term -- n^{-1} \sum_{i=1}^n \sum_{j=1}^r p_j \Psi_0(X_i, mu_j)

    # Loop on h
    psi_0_h_kappa <- function(h_k) {

      # ||h^{-2} X_i + kappa_j mu_j||^2 = h^{-4} + kappa_j^2
      #                                  + 2 * h^{-2} * kappa_j * X_i'mu_j
      h2_X_i_kappa_j_mu_j <- sqrt(1 / h_k^4 + kappa^2 + cross_mu_X / h_k^2)
      log_C_h_X_i_kappa_mu_j <-
        fast_log_c_vMF(p = d + 1, kappa = h2_X_i_kappa_j_mu_j, spline = spline)
      log_C_h2_kappa_prop <- fast_log_c_vMF(p = d + 1, kappa = 1 / h_k^2,
                                            spline = spline) +
        log(prop) + fast_log_c_vMF(p = d + 1, kappa = kappa, spline = spline)
      return(sum(exp(log_C_h2_kappa_prop - log_C_h_X_i_kappa_mu_j)))

    }
    Psi_03 <- sapply(h, psi_0_h_kappa) / n

    # ISE
    ise <- sqrt(Psi_01 + Psi_02 - 2 * Psi_03)

  } else {

    # Set seeds for the Monte Carlo
    if (!is.null(seed_psi)) {

      # old_seed <- .Random.seed
      # on.exit({.Random.seed <<- old_seed})
      set.seed(seed_psi, kind = "Mersenne-Twister")

    }

    # Importance-sampling Monte Carlo of L^p norm
    if (is.null(x_mvmf) || is.null(f_mvmf)) {

      x_mvmf <- r_mvmf_polysph(n = M_psi, d = d, mu = mu, kappa = kappa,
                               prop = prop)
      f_mvmf <- d_mvmf_polysph(x = x_mvmf, d = d, mu = mu, kappa = kappa,
                               prop = prop)

    }
    ise <- rep(NA, nrow(h))
    for (j in seq_len(nrow(h))) {

      f_hat <- kde_polysph(x = x_mvmf, X = X, d = d, h = h[j, ], kernel = 1)
      ise[j] <- mean(abs(f_hat - f_mvmf)^p / f_mvmf)^(1 / p)

    }
    Psi_01 <- Psi_02 <- Psi_03 <- NULL

  }
  return(list("ise" = ise,
              "Psi_01" = Psi_01, "Psi_02" = Psi_02, "Psi_03" = Psi_03))

}


#' @noRd
#' @rdname ise
log1p_ise <- function(log_h, X, d, mu, kappa, prop, M_psi = 1e4,
                      x_mvmf = NULL, f_mvmf = NULL, seed_psi = NULL,
                      spline = FALSE, exact = FALSE, p = 2) {

  log1p(ise_vmf_polysph(X = X, d = d, h = exp(log_h), mu = mu, kappa = kappa,
                        prop = prop, M_psi = M_psi, x_mvmf = x_mvmf,
                        f_mvmf = f_mvmf, seed_psi = seed_psi,
                        spline = spline, exact = exact, p = p)$ise)

}


#' @title Exact ISE bandwidth selection for polyspherical kernel density
#' estimator for a mixture of von Mises--Fisher distributions
#'
#' @description Computes the bandwidths
#' \deqn{\boldsymbol{h}_{\mathrm{ISE}}=\arg\min_{\boldsymbol{h}>0}
#' \mathrm{ISE}_m[\hat{f}(\cdot;\boldsymbol{h})]}
#' of the kernel density estimator on the polysphere
#' \eqn{\mathcal{S}^{d_1} \times \cdots \times \mathcal{S}^{d_r}} for an
#' \eqn{m}-mixture of von Mises--Fisher densities
#' \eqn{f_m(\boldsymbol{x}) = \sum_{j=1}^m p_j
#' f_{\mathrm{vMF}}(\boldsymbol{x}; \boldsymbol{\mu}_j, \kappa_j)}.
#'
#' @inheritParams ise_vmf
#' @inheritParams mise
#' @inheritParams r_mvmf_polysph
#' @param bw0 initial bandwidth vector for minimizing the ISE loss. Can be
#' also a matrix of initial bandwidth vectors.
#' @return A list with entries \code{bw} (optimal bandwidth) and \code{opt},
#' the latter containing the output of \code{\link[stats]{nlm}}.
#' @examples
#' n <- 200
#' mu <-  rbind(c(1, 0, 0), c(-1, 0, 0))
#' kappa <- cbind(c(1, 1))
#' prop <- c(0.5, 0.5)
#' samp <- r_mvmf_polysph(n = n, d = 2, mu = mu, kappa = kappa, prop = prop)
#' polykde:::bw_mise_polysph(n = n, d = 2, bw0 = 5, mu = mu, kappa = kappa,
#'                           prop = prop, seed_psi = 1)
#' polykde:::bw_ise_polysph(X = samp, d = 2, bw0 = 5, mu = mu, kappa = kappa,
#'                          prop = prop, seed_psi = 1)
#' @noRd
bw_ise_polysph <- function(X, d, bw0 = NULL, mu, kappa, prop, M_psi = 1e4,
                           x_mvmf = NULL, f_mvmf = NULL, seed_psi = NULL,
                           spline = FALSE, exact = FALSE, p = 2, ...) {

  # Set seeds for the Monte Carlo
  if (!is.null(seed_psi)) {

    # old_seed <- .Random.seed
    # on.exit({.Random.seed <<- old_seed})
    set.seed(seed_psi, kind = "Mersenne-Twister")

  }

  # Loss function
  if (is.null(x_mvmf) || is.null(f_mvmf)) {

    x_mvmf <- r_mvmf_polysph(n = M_psi, d = d, mu = mu, kappa = kappa,
                             prop = prop)
    f_mvmf <- d_mvmf_polysph(x = x_mvmf, d = d, mu = mu, kappa = kappa,
                             prop = prop)

  }

  # Get r
  r <- length(d)

  # log1p_ise with fixed arguments
  log1p_ise_args <- function(log_h) {

    log1p_ise(log_h = log_h, X = X, d = d, mu = mu, kappa = kappa,
              prop = prop, x_mvmf = x_mvmf, f_mvmf = f_mvmf, spline = spline,
              exact = exact, p = p)

  }

  # Initial search
  if (is.null(bw0)) {

    stop("bw0 must be provided.")

  } else {

    bw0 <- rbind(bw0)
    if (r != ncol(bw0)) {

      stop("bw0 and d are incompatible.")

    }
    if (nrow(bw0) > 1) {

      # Name arguments in apply() to avoid mismatches between MARGIN and M_psi
      ise_bw0 <- apply(X = log(bw0), MARGIN = 1, FUN = log1p_ise_args)
      ind_best <- which.min(ise_bw0)
      bw0 <- bw0[ind_best, ]

    }

  }

  # Search h_ISE
  opt <- nlm(f = log1p_ise_args, p = log(bw0), ...)
  bw <- exp(opt$estimate)
  return(list("bw" = unname(bw), "opt" = opt))

}
