
#' @title Kernels on the hypersphere and their derivatives
#'
#' @description An isotropic kernel \eqn{L} on \eqn{\mathcal{S}^d} and its
#' normalizing constant are such that \eqn{\int_{\mathcal{S}^d} c(h, d, L)
#' L\left(\frac{1 - \boldsymbol{x}'\boldsymbol{y}}{h^2}\right) d\boldsymbol{x}
#' = 1} (extrinsic-chordal distance) or \eqn{\int_{\mathcal{S}^d} c(h,d,L)
#' L\left(\frac{\cos^{-1}(\boldsymbol{x}'\boldsymbol{y})^2}{2h^2}\right)
#' d\boldsymbol{x} = 1} (intrinsic distance).
#' @inheritParams kde_polysph
#' @param t kernel argument. A positive scalar.
#' @param squared square the kernel? Only for \code{deriv == 0}. Defaults to
#' \code{FALSE}.
#' @param deriv kernel derivative. Must be \code{0}, \code{1}, or \code{2}.
#' Defaults to \code{0}.
#' @param inc_sfp include \code{softplus(k)} in the constant? Defaults to
#' \code{TRUE}.
#' @param log compute the logarithm of the constant?
#' @inheritParams L
#' @param intrinsic consider the intrinsic distance? Defaults to \code{FALSE}.
#' @param y center of the kernel.
#' @details The gradient and Hessian are computed for the functions
#' \eqn{\boldsymbol{x} \mapsto
#' L\left(\frac{1 - \boldsymbol{x}'\boldsymbol{y}}{h^2}\right)}.
#' @return For \code{L}, a vector with the kernel evaluated at \code{t}. For
#' \code{grad_L}, a vector with the gradient evaluated at \code{x}. For
#' \code{hess_L}, a matrix with the Hessian evaluated at \code{x}.
#' @examples
#' # Constants in terms of h
#' h_grid <- seq(0.01, 4, l = 100)
#' r <- 2
#' d <- 2
#' dr <- rep(d, r)
#' c_vmf <- sapply(h_grid, function(hi)
#'   log(c_kern(h = rep(hi, r), d = dr, kernel = 1, kernel_type = 2)))
#' c_epa <- sapply(h_grid, function(hi)
#'   log(c_kern(h = rep(hi, r), d = dr, kernel = 2, kernel_type = 2)))
#' c_sfp <- sapply(h_grid, function(hi)
#'   log(c_kern(h = rep(hi, r), d = dr, kernel = 3, k = 1, kernel_type = 2)))
#' plot(h_grid, c_epa, type = "l", ylab = "Constant", xlab = "h", col = 2)
#' lines(h_grid, c_sfp, col = 3)
#' lines(h_grid, c_vmf, col = 1)
#' abline(v = sqrt(2), lty = 2, col = 2)
#'
#' # Kernel and its derivatives
#' h <- 0.5
#' x <- c(sqrt(2), -sqrt(2), 0) / 2
#' y <- c(-sqrt(2), sqrt(3), sqrt(3)) / 3
#' L(t = (1 - sum(x * y)) / h^2)
#' grad_L(x = x, y = y, h = h)
#' hess_L(x = x, y = y, h = h)
#' @name kernel


#' @rdname kernel
#' @export
L <- function(t, kernel = "1", squared = FALSE, deriv = 0, k = 10,
              inc_sfp = TRUE) {

  # Stop if invalid kernel
  if (!(kernel %in% c(1:3, "1", "2", "3"))) {

    stop("\"kernel\" must be 1, 2, or 3.")

  }

  # Which derivative?
  if (deriv == 0) {

    L_0 <- switch(kernel,
                  "1" = exp(-t) * (0 <= t),
                  "2" = (1 - t) * (0 <= t & t < 1),
                  "3" = softplus(k * (1 - t)) /
                    ifelse(inc_sfp, softplus(k), 1) * (0 <= t)
    )

    if (squared) L_0 <- L_0 * L_0
    return(L_0)

  } else if (deriv == 1) {

    L_1 <- switch(kernel,
                  "1" = -exp(-t) * (0 <= t),
                  "2" = (-1) * (0 <= t & t < 1),
                  "3" = -k / (ifelse(inc_sfp, softplus(k), 1) *
                                (1 + exp(-k * (1 - t)))) * (0 <= t)
    )
    return(L_1)

  } else if (deriv == 2) {

    L_2 <- switch(kernel,
                  "1" = exp(-t) * (0 <= t),
                  "2" = rep(0, length(t)) * (0 < t),
                  "3" = k^2 / ifelse(inc_sfp, softplus(k), 1) *
                    exp(-k * (1 - t)) / (1 + exp(-k * (1 - t)))^2 * (0 <= t)
    )
    return(L_2)

  } else {

    stop("deriv must be 0, 1, or 2.")

  }

}


#' @rdname kernel
#' @export
c_kern <- function(h, d, kernel = "1", kernel_type = "1", k = 10, log = FALSE,
                   inc_sfp = TRUE, intrinsic = FALSE) {

  if (intrinsic) {

    if (kernel_type == "1") {

      const <- rotasym::w_p(p = d) * sapply(seq_along(d), function(i)
        integrate(f = function(t)
          L(t = acos(t)^2 / (2 * h[i]^2), kernel = kernel, k = k,
            inc_sfp = inc_sfp) * (1 - t^2)^(d[i] / 2 - 1),
          lower = -1, upper = 1)$value)
      return(switch(log + 1, 1 / const, -log(const)))

    } else if (kernel_type == "2") {

      if (kernel == "1") {

        const <- rotasym::w_p(p = d) * sapply(seq_along(d), function(i)
          integrate(f = function(t)
            L(t = acos(t)^2 / (2 * h[i]^2), kernel = "1") *
              (1 - t^2)^(d[i] / 2 - 1), lower = -1, upper = 1)$value)
        const <- prod(const)
        return(switch(log + 1, 1 / const, -log(const)))

      } else {

        # Asymptotic form
        intrinsic <- FALSE

      }

    } else {

      stop("\"kernel_type\" must be 1 or 2.")

    }

  }

  # A separate if is needed because in the previous one intrinsic <- FALSE
  if (!intrinsic) {

    if (kernel == "1") {

      # Case d = 2 vs. general case
      if (all(d == 2)) {

        const <- -log(1 / h^2) + log(2 * pi) + log1p(-exp(-2 / h^2))

      } else {

        const <- (0.5 * (d + 1)) * log(2 * pi) + (d - 1) * log(h) +
          log(besselI(x = 1 / h^2, nu = 0.5 * (d - 1), expon.scaled = TRUE))

      }

      # Product or spherically symmetric
      if (kernel_type == "1") {

        return(switch(log + 1, exp(-const), -const))

      } else if (kernel_type == "2") {

        return(switch(log + 1, exp(-sum(const)), -sum(const)))

      } else {

        stop("\"kernel_type\" must be 1 or 2.")

      }

    } else if (kernel == "2") {

      # Product or spherically symmetric
      if (kernel_type == "1") {

        m_h <- pmax(-1, 1 - h^2)
        F_d <- sapply(seq_along(d), function(i)
          drop(sphunif::p_proj_unif(x = m_h[i], p = d[i] + 1)))
        const <- rotasym::w_p(p = d + 1) * (1 - h^(-2)) * (1 - F_d) +
          rotasym::w_p(p = d) * (1 - m_h^2)^(d / 2) / (d * h^2)
        return(switch(log + 1, 1 / const, -log(const)))

      } else if (kernel_type == "2") {

        # Exact form
        if (all(h == h[1]) && all(d == 2)) {

          r <- length(d)
          if (h[1] < sqrt(2)) {

            const <- r * log(2 * pi * h[1]^2) - lfactorial(r + 1)
            return(switch(log + 1, exp(-const), -log(const)))

          } else {

            l <- seq(ceiling(r - h[1]^2 / 2), r, by = 1)
            const <- sum(
              exp(r * log(2 * pi * h[1]^2) + lchoose(r, l) -
                    lfactorial(r + 1)) *
                (-1)^(r + l) * (1 - 2 * (r - l) / h[1]^2)^(r + l))
            return(switch(log + 1, 1 / const, -log(const)))

          }

          # Asymptotic form
        } else {

          d_tilde <- sum(d)
          const <- sum(d * log(h)) + (d_tilde / 2) * log(2 * pi) +
            log(4 / (d_tilde * (d_tilde + 2))) - lgamma(d_tilde / 2)
          return(switch(log + 1, exp(-const), -const))

        }

      } else {

        stop("\"kernel_type\" must be 1 or 2.")

      }

    } else if (kernel == "3") {

      # Product or spherically symmetric
      if (kernel_type == "1") {

        # N <- 1280
        # t_k <- sphunif::Gauss_Legen_nodes(N = N)
        # w_k <- sphunif::Gauss_Legen_weights(N = N)
        # const <- rotasym::w_p(p = d) / softplus(k) *
        #   sapply(seq_along(d), function(i)
        #     sum(w_k * softplus(k * (1 - (1 - t_k) / h[i]^2)) *
        #           (1 - t_k^2)^(d[i] / 2 - 1)))
        const <- rotasym::w_p(p = d) * sapply(seq_along(d), function(i)
          integrate(f = function(t)
            L(t = (1 - t) / h[i]^2, kernel = 3, k = k, inc_sfp = inc_sfp) *
              (1 - t^2)^(d[i] / 2 - 1), lower = -1, upper = 1)$value)
        return(switch(log + 1, 1 / const, -log(const)))

      } else if (kernel_type == "2") {

        # Exact form
        if (all(h == h[1]) && all(d == 2)) {

          r <- length(d)
          l <- seq(0, r, by = 1)
          mu <- k * (1 - r / h[1]^2) + k * (2 * l - r) / h[1]^2
          const <- sum(
            exp(r * log(2 * pi * h[1]^2) - r * log(k) -
                  ifelse(inc_sfp, log(softplus(k)), 0) +
                  log(-polylog_minus_exp_mu(mu = mu, s = r + 1)) +
                  lchoose(r, l)) * (-1)^(r + l))
          return(switch(log + 1, 1 / const, -log(const)))

          # Asymptotic form
        } else {

          d_tilde <- sum(d)
          const <- sum(d * log(h)) + (d_tilde / 2) * log(2 * pi / k) +
            log(-polylog_minus_exp_mu(mu = k, s = d_tilde / 2 + 1)) -
            ifelse(inc_sfp, log(softplus(k)), 0)
          return(switch(log + 1, exp(-const), -const))

        }

      } else {

        stop("\"kernel_type\" must be 1 or 2.")

      }

    } else {

      stop("\"kernel\" must be 1, 2, or 3.")

    }

  }

}


#' @rdname kernel
#' @export
grad_L <- function(x, y, h, kernel = 1, k = 10) {

  -drop(L(t = (1 - x %*% y) / h^2, kernel = kernel, deriv = 1, k = k)) *
    y / h^2

}


#' @rdname kernel
#' @export
hess_L <- function(x, y, h, kernel = 1, k = 10) {

  drop(L(t = (1 - x %*% y) / h^2, kernel = kernel, deriv = 2, k = k)) *
    tcrossprod(y) / h^4

}


#' @title Sample from the angular kernel density
#'
#' @description Simulation from the angular density function of an isotropic
#' kernel on \eqn{\mathcal{S}^d}.
#'
#' @inheritParams r_unif_polysph
#' @inheritParams kde_polysph
#' @return A vector of length \code{n}.
#' @examples
#' hist(r_g_kern(n = 1e3, d = 2, h = 1, kernel = "1"), breaks = 30,
#'      probability = TRUE, main = "", xlim = c(-1, 1))
#' hist(r_g_kern(n = 1e3, d = 2, h = 1, kernel = "2"), breaks = 30,
#'      probability = TRUE, main = "", xlim = c(-1, 1))
#' hist(r_g_kern(n = 1e3, d = 2, h = 1, kernel = "3"), breaks = 30,
#'      probability = TRUE, main = "", xlim = c(-1, 1))
#' @export
r_g_kern <- function(n, d, h, kernel = "1", k = 10) {

  stopifnot(length(d) == 1)
  if (kernel == "1") {

    return(rotasym::r_g_vMF(n = n, p = d + 1, kappa = 1 / h^2))

  } else if (kernel == "2") {

    # Inversion method with explicit solution for d = 2
    u <- runif(n)
    m_h <- max(-1, 1 - h^2)
    if (d == 2) {

      e <- 1 - h^2 + sqrt((h^2 - 1)^2 - m_h + (m_h + u * (1 - m_h)) *
                            (2 * h^2 - 1 + m_h))

    } else {

      F_epa <- function(t) {
        (rotasym::w_p(p = d) * c_kern(h = h, d = d, kernel = "2") / h^2) * (
          (h^2 - 1) * rotasym::w_p(p = d + 1) / rotasym::w_p(p = d) *
            (drop(sphunif::p_proj_unif(x = t, p = d + 1)) -
               drop(sphunif::p_proj_unif(x = m_h, p = d + 1))) -
            ((1 - t^2)^(d / 2) - (1 - m_h^2)^(d / 2)) / d
        )
      }
      e <- sapply(u, function(v) uniroot(f = function(t) F_epa(t) - v,
                                         lower = m_h, upper = 1,
                                         tol = 1e-6)$root)

    }
    return(e)

  } else if (kernel == "3") {

    if (d == 2) {

      # Inversion method with explicit solution for d = 2
      u <- runif(n)
      a <- polylog_minus_exp_mu(mu = k * (1 - 2 / h^2), s = 2)
      den <- a - polylog_minus_exp_mu(mu = k, s = 2)
      F_sfp <- function(t) {
        (a - polylog_minus_exp_mu(mu = k * (1 - (1 - t) / h^2), s = 2)) / den
      }
      e <- sapply(u, function(v) uniroot(f = function(t) F_sfp(t) - v,
                                         lower = -1, upper = 1,
                                         tol = 1e-6)$root)
      # TODO: this can be made faster by performing series inversion on the
      # polylogarithm function https://math.stackexchange.com/a/2840919

    } else {

      # Acceptance-rejection from vMF kernel

      # Get M
      t_grid <- seq(-1, 1, l = 100)
      t_grid <- (1 - t_grid) / h^2
      const_f <- c_kern(h = h, d = d, kernel = kernel, kernel_type = "1", k = k)
      const_g <- c_kern(h = h, d = d, kernel = "1", kernel_type = "1")
      dens_f <- L(t = t_grid, kernel = kernel, k = k)
      dens_g <- L(t = t_grid, kernel = "1")
      dens_ratio <- dens_g / dens_f
      M <- 1 / min(dens_ratio[is.finite(dens_ratio)])
      ratio <- const_f / (const_g * M)

      # Acceptance-rejection
      V <- numeric(n)
      i <- 1
      while (i <= n) {

        # vMF sample
        y <- r_g_kern(n = 1, d = d, h = h, kernel = "1")
        t <- (1 - y) / h^2

        # Rejection step
        u <- runif(1)
        if (u < (ratio * L(t = t, kernel = kernel, k = k) /
                 L(t = t, kernel = "1"))) {

          V[i] <- y
          i <- i + 1

        }

      }
      return(V)

    }

  } else {

    stop("\"kernel\" must be 1, 2, or 3.")

  }

}


#' @title Polyspherical kernel moments and efficiencies
#'
#' @description Computes moments of kernels on \eqn{\mathcal{S}^{d_1} \times
#' \cdots \times \mathcal{S}^{d_r}} and efficiencies of kernels on
#' \eqn{(\mathcal{S}^d)^r}.
#'
#' @param d either a vector with the dimensions or a scalar with the
#' common dimension of each hypersphere (a scalar).
#' @param r number of polyspheres of the same dimension.
#' @inheritParams kde_polysph
#' @param kernel_ref reference kernel to which compare the efficiency. Uses the
#' same codification as the \code{kernel}. Defaults to \code{"2"}.
#' @param kernel_type type of kernel. Must be either \code{"prod"} (product
#' kernel, default) or \code{"sph"} (spherically symmetric kernel).
#' @param kernel_ref_type type of the reference kernel. Must be either
#' \code{"prod"} (product kernel) or \code{"sph"} (spherically symmetric kernel,
#' default).
#' @param ... further parameters passed to \code{\link{integrate}}, such as
#' \code{upper}, \code{abs.tol}, \code{rel.tol}, etc.
#' @param bias,squared flags parameters for computing the numerator constants
#' in the bias and variance constants.
#' @return For \code{b_d}, a vector with the first kernel moment on
#' each hypersphere (common if \code{kernel_type = "sph"}). For \code{v_d},
#' a vector with the second kernel moment if \code{kernel_type = "prod"}, or a
#' scalar if \code{kernel_type = "sph"}. For \code{eff_kern}, a scalar with the
#' kernel efficiency.
#' @examples
#' # Kernel moments
#' b_d(kernel = 2, d = c(2, 3), kernel_type = "prod")
#' v_d(kernel = 2, d = c(2, 3), kernel_type = "prod")
#' b_d(kernel = 2, d = c(2, 3), kernel_type = "sph")
#' v_d(kernel = 2, d = c(2, 3), kernel_type = "sph")
#'
#' # Kernel efficiencies
#' eff_kern(d = 2, r = 1, kernel = "1")
#' eff_kern(d = 2, r = 1, kernel = "2")
#' eff_kern(d = 2, r = 1, k = 10, kernel = "3")
#' @export
eff_kern <- function(d, r, k = 10, kernel, kernel_type = c("prod", "sph")[1],
                     kernel_ref = "2", kernel_ref_type = c("prod", "sph")[2],
                     ...) {

  # Check d and r
  stopifnot(length(d) == 1)
  stopifnot(length(r) == 1)

  # Make kernel_type and kernel_ref_type character
  if (is.numeric(kernel_type)) {

    kernel_type <- switch(kernel_type, "1" = "prod", "2" = "sph",
                          stop("\"kernel_type\" must be 1 or 2."))

  }
  if (is.numeric(kernel_ref_type)) {

    kernel_ref_type <- switch(kernel_ref_type, "1" = "prod", "2" = "sph",
                              stop("\"kernel_type\" must be 1 or 2."))

  }

  # AMISE constant of the kernel, without raising to the power 4 / (d * r + 4)
  C_d_r <- function(kern, kernel_type) {

    if (kernel_type == "prod") {

      C <- v_d(kernel = kern, d = d, kernel_type = "prod",
               k = k, ...)^r *
        b_d(kernel = kern, d = d, kernel_type = "prod",
            k = k, ...)^(d * r / 2)

    } else if (kernel_type == "sph") {

      C <- v_d(kernel = kern, d = rep(d, r), kernel_type = "sph",
               k = k, ...) *
        b_d(kernel = kern, d = rep(d, r), kernel_type = "sph",
            k = k, ...)[1]^(d * r / 2)

    } else {

      stop("\"kernel_type\" must be either \"prod\" or \"sph\".")

    }
    return(C)

  }

  # Efficiency
  # eff <- (C_d_r(kern = kernel_ref, kernel_type = kernel_ref_type) /
  #           C_d_r(kern = kernel, kernel_type = kernel_type))^((d * r + 4) / 4)
  eff <- C_d_r(kern = kernel_ref, kernel_type = kernel_ref_type) /
    C_d_r(kern = kernel, kernel_type = kernel_type)
  return(eff)

}


#' @rdname eff_kern
#' @export
b_d <- function(kernel, d, k = 10, kernel_type = c("prod", "sph")[1], ...) {

  # Make kernel_type and kernel_ref_type character
  if (is.numeric(kernel_type)) {

    kernel_type <- switch(kernel_type, "1" = "prod", "2" = "sph",
                          stop("\"kernel_type\" must be 1 or 2."))

  }

  if (is.function(kernel)) {

    if (kernel_type == "prod") {

      bd <- sapply(d, function(di) {
        num <- integrate(function(t) kernel(t) * t^(di / 2),
                         lower = 0, upper = Inf, ...)$value
        den <- integrate(function(t) kernel(t) * t^(di / 2 - 1),
                         lower = 0, upper = Inf, ...)$value
        return(num / (di * den))
      })

    } else if (kernel_type == "sph") {

      d_tilde <- sum(d)
      num <- integrate(function(t) kernel(t) * t^(d_tilde / 2),
                       lower = 0, upper = Inf, ...)$value
      den <- integrate(function(t) kernel(t) * t^(d_tilde / 2 - 1),
                       lower = 0, upper = Inf, ...)$value
      bd <- num / (d_tilde * den)

    } else {

      stop("\"kernel_type\" must be either \"prod\" or \"sph\".")

    }

  } else {

    if (kernel == "1") {

      bd <- rep(1 / 2, ifelse(kernel_type == "prod", length(d), 1))

    } else if (kernel == "2") {

      if (kernel_type == "prod") {

        bd <- 1 / (d + 4)

      } else if (kernel_type == "sph") {

        d_tilde <- sum(d)
        bd <- 1 / (d_tilde + 4)

      } else {

        stop("\"kernel_type\" must be either \"prod\" or \"sph\".")

      }

    } else if (kernel == "3") {

      if (kernel_type == "prod") {

        Li_num <- sapply(d, function(di)
          polylog_minus_exp_mu(s = di / 2 + 2, mu = k))
        Li_den <- sapply(d, function(di)
          polylog_minus_exp_mu(s = di / 2 + 1, mu = k))
        bd <- Li_num / (2 * k * Li_den)

      } else if (kernel_type == "sph") {

        d_tilde <- sum(d)
        Li_num <- polylog_minus_exp_mu(s = d_tilde / 2 + 2, mu = k)
        Li_den <- polylog_minus_exp_mu(s = d_tilde / 2 + 1, mu = k)
        bd <- Li_num / (2 * k * Li_den)

      } else {

        stop("\"kernel_type\" must be either \"prod\" or \"sph\".")

      }

    } else {

      stop("\"kernel\" must be 1 (vMF), 2 (Epa), 3 (softplus), or a function.")

    }

  }

  # Ensure returned value is a vector
  if (kernel_type == "sph") {

    bd <- rep(bd, length(d))

  }
  return(bd)

}


#' @rdname eff_kern
#' @export
v_d <- function(kernel, d, k = 10, kernel_type = c("prod", "sph")[1], ...) {

  # Make kernel_type and kernel_ref_type character
  if (is.numeric(kernel_type)) {

    kernel_type <- switch(kernel_type, "1" = "prod", "2" = "sph",
                          stop("\"kernel_type\" must be 1 or 2."))

  }

  if (is.function(kernel)) {

    if (kernel_type == "prod") {

      vd <- sapply(d, function(di) {
        num <- integrate(function(t) kernel(t)^2 * t^(di / 2 - 1),
                         lower = 0, upper = Inf, ...)$value
        den <- integrate(function(t) kernel(t) * t^(di / 2 - 1),
                         lower = 0, upper = Inf, ...)$value
        return(exp(-(di / 2 - 1) * log(2) - rotasym::w_p(p = di, log = TRUE) +
                     log(num) - 2 * log(den)))
      })

    } else if (kernel_type == "sph") {

      d_tilde <- sum(d)
      num <- integrate(function(t) kernel(t)^2 * t^(d_tilde / 2 - 1),
                       lower = 0, upper = Inf, ...)$value
      den <- integrate(function(t) kernel(t) * t^(d_tilde / 2 - 1),
                       lower = 0, upper = Inf, ...)$value
      vd <- exp(-(d_tilde / 2 - 1) * log(2) -
                  rotasym::w_p(p = d_tilde, log = TRUE) +
                  log(num) - 2 * log(den))

    } else {

      stop("\"kernel_type\" must be either \"prod\" or \"sph\".")

    }

  } else {

    if (kernel == "1") {

      if (kernel_type == "prod") {

        vd <- 1 / (2 * sqrt(pi))^d

      } else if (kernel_type == "sph") {

        d_tilde <- sum(d)
        vd <- 1 / (2 * sqrt(pi))^d_tilde

      } else {

        stop("\"kernel_type\" must be either \"prod\" or \"sph\".")

      }

    } else if (kernel == "2") {

      if (kernel_type == "prod") {

        vd <- exp(log(4) + lgamma(d / 2 + 2) - (d / 2) * log(2 * pi) -
                    log(d + 4))

      } else if (kernel_type == "sph") {

        d_tilde <- sum(d)
        vd <- exp(log(4) + lgamma(d_tilde / 2 + 2) -
                    (d_tilde / 2) * log(2 * pi) - log(d_tilde + 4))

      } else {

        stop("\"kernel_type\" must be either \"prod\" or \"sph\".")

      }

    } else if (kernel == "3") {

      if (kernel_type == "prod") {

        vd <- exp(
          d * log(k) + log(J_d_k(d = d, k = k, ...)) -
            ((d / 2) * log(2 * pi) + lgamma(d / 2) +
               2 * log(abs(polylog_minus_exp_mu(s = d / 2 + 1, mu = k)))))

      } else if (kernel_type == "sph") {

        d_tilde <- sum(d)
        vd <- exp(
          d_tilde * log(k) + log(J_d_k(d = d_tilde, k = k, ...)) -
            ((d_tilde / 2) * log(2 * pi) + lgamma(d_tilde / 2) +
               2 * log(abs(polylog_minus_exp_mu(s = d_tilde / 2 + 1, mu = k)))))

      } else {

        stop("\"kernel_type\" must be either \"prod\" or \"sph\".")

      }

    } else {

      stop("\"kernel\" must be 1 (vMF), 2 (Epa), 3 (softplus), or a function.")

    }

  }
  return(vd)

}


#' @rdname eff_kern
#' @keywords internal
lambda_h <- function(d, h = NULL, kernel = "1", bias = FALSE, squared = FALSE,
                     k = 10, ...) {

  if (is.null(h)) {

    lambda <- rotasym::w_p(p = d) * 2^(d / 2 - 1) * sapply(d, function(di) {
      integrate(function(t) L(t = t, kernel = kernel, squared = squared,
                              k = k) * t^(di / 2 - !bias),
                lower = 0, upper = Inf, ...)$value
    })

  } else {

    lambda <- rotasym::w_p(p = d) * sapply(d, function(di) {
      integrate(function(t) L(t = t, kernel = kernel, squared = squared,
                              k = k) * t^(di / 2 - !bias) *
                  (2 - h^2 * t)^(di / 2 - 1),
                lower = 0, upper = 2 / h^2, ...)$value
    })

  }
  return(lambda)

}


#' @rdname eff_kern
#' @keywords internal
lambda_vmf_h <- function(d, h = NULL, bias = FALSE, squared = FALSE) {

  if (bias) {

    lambda <- ((d - 1) / 2) * log(2) + 0.5 * log(pi) +
                 log(besselI(x = 1 / h^2, nu = (d - 1) / 2,
                             expon.scaled = TRUE) -
                       besselI(x = 1 / h^2, nu = (d + 1) / 2,
                               expon.scaled = TRUE)) +
                 lgamma(d / 2) - 3 * log(h)

  } else if (squared) {

    lambda <- 0.5 * log(pi) + log(besselI(x = 2 / h^2, nu = (d - 1) / 2,
                                  expon.scaled = TRUE)) +
      lgamma(d / 2) - log(h)

  } else {

    lambda <- ((d - 1) / 2) * log(2) + 0.5 * log(pi) +
      log(besselI(x = 1 / h^2, nu = (d - 1) / 2, expon.scaled = TRUE)) +
      lgamma(d / 2) - log(h)

  }
  return(exp(rotasym::w_p(p = d, log = TRUE) + lambda))

}


#' @title Fixed-point algorithm for the weighted Epanechnikov kernel
#'
#' @description Solves the equation
#' \deqn{\boldsymbol{\beta} = \mathrm{diag}(\boldsymbol{d}^{\odot (-1)})
#' \boldsymbol{R}^{-1} \boldsymbol{\beta}^{\odot(-1)}}
#' to determine the weights of the weighted Epanechnikov kernel.
#'
#' @inheritParams eff_kern
#' @param R curvature matrix as outputted by \code{\link{curv_vmf_polysph}}.
#' @param tol tolerance for the fixed-point algorithm. Defaults to \code{1e-10}.
#' @return Vector \eqn{\boldsymbol{\beta}} of weights.
#' @examples
#' r <- 2
#' d <- seq_len(r)
#' R <- curv_vmf_polysph(kappa = 10 * d, d = d)
#' polykde:::beta0_R(d = d, R = R)
#' @keywords internal
beta0_R <- function(d, R, tol = 1e-10) {

  # # Simple fixed-point algorithm
  # b <- rep(1, length(d))
  # b_old <- rep(0, length(d))
  # inv_d_R <- diag(1 / d) %*% R
  # while (sum(abs(b - b_old)) > tol) {
  #
  #   b_old <- b
  #   b <- inv_d_R %*% (1 / b_old)
  #
  # }
  # return(b)

  # Fixed-point algorithm
  r <- nrow(R)
  inv_d_R <- diag(1 / d, nrow = r, ncol = r) %*% R
  f <- function(beta) inv_d_R %*% (1 / beta)
  fp <- FixedPoint::FixedPoint(Function = f, Inputs = rep(1, r),
                               Method = "Anderson")
  return(fp$FixedPoint)

}
