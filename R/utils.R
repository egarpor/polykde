
#' @title Stable computation of the softplus function
#'
#' @description Computes the softplus function \eqn{\log(1+e^{t})} in a
#' numerically stable way for large absolute values of \eqn{t}.
#'
#' @param t vector or matrix.
#' @return The softplus function evaluated at \code{t}.
#' @examples
#' curve(softplus(10 * (1 - (1 - x) / 0.1)), from = -1, to = 1)
#' @export
softplus <- function(t) {

  # Evaluate the function in a stable for either positive or negative t's
  res <- numeric(length(t))
  ind_pos <- t >= 0
  res[ind_pos] <- t[ind_pos] + log1p(exp(-t[ind_pos]))
  res[!ind_pos] <- log1p(exp(t[!ind_pos]))
  return(res)

}


#' @rdname dist_polysph
#' @export
dist_polysph_matrix <- function(x, ind_dj, norm_x = FALSE, norm_y = FALSE,
                                std = TRUE) {

  n <- nrow(x)
  dij <- matrix(0, nrow = n, ncol = n)
  for (i in seq_len(n)[-1]) {
    for (j in seq_len(i - 1)) {

      dij[i, j] <- drop(dist_polysph(x = x[i, , drop = FALSE],
                                     y = x[j, , drop = FALSE],
                                     ind_dj = ind_dj, norm_x = norm_x,
                                     norm_y = norm_y, std = std))

    }
  }
  dij <- s(dij, add = TRUE)
  colnames(dij) <- rownames(dij) <- seq_len(n)
  return(dij)

}


#' @title Index of spheres on a polysphere
#'
#' @description Given Cartesian coordinates of polyspherical data, computes the
#' \code{0}-based indexes at which the Cartesian coordinates for each sphere
#' start and end.
#'
#' @inheritParams kde_polysph
#' @return A vector of size \code{sum(d) + 1}.
#' @examples
#' # Example on (S^1)^3
#' d <- c(1, 1, 1)
#' comp_ind_dj(d = d)
#' comp_ind_dj(d = d) + 1
#'
#' # Example on S^1 x S^2
#' d <- c(1, 2)
#' comp_ind_dj(d = d)
#' comp_ind_dj(d = d) + 1
#' @export
comp_ind_dj <- function(d) {

  ind_dj <- rep(0, length(d) + 1)
  ind_dj[-1] <- d + 1
  return(cumsum(ind_dj))

}


#' @title Merge a list of lists
#'
#' @description Binds the entries of a list of lists by rows or columns.
#'
#' @param lists lists whose entries are to be binded. Must have a common
#' structure (same fields with the same dimensions).
#' @param bind bind operator, either \code{"rbind"} or \code{"cbind"}.
#' @return A list with the same structure as the input lists, but with the
#' entries binded.
#' @examples
#' polykde:::bind_lists(lists = list(list(1:3, 4:6), list(7:9, 10:12)))
#' @noRd
bind_lists <- function(lists, bind = "rbind") {

  stopifnot(bind %in% c("rbind", "cbind"))
  n <- length(lists)
  merged <- lists[[which(!sapply(lists, is.null))[1]]]
  for (j in seq_along(merged)) {

    bind_j <- ifelse(length(dim(merged[[j]])) == 3,
                     function(...)
                       abind::abind(..., along = ifelse(bind == "rbind", 1, 2)),
                     bind)
    merged[[j]] <- do.call(what = bind_j,
                           args = lapply(seq_len(n),
                                         function(i) lists[[i]][[j]]))

  }
  return(merged)

}


#' @title Polylogarithm function with negative argument
#'
#' @description Computation of the polylogarithm \eqn{\mathrm{Li}_s(-e^\mu)}.
#'
#' @param mu vector with exponents of the negative argument.
#' @param s vector with indexes of the polylogarithm.
#' @inheritParams eff_kern
#' @param upper upper limit of integration. Defaults to \code{Inf}.
#' @details If \code{s} is an integer, 1/2, 3/2, or 5/2, then routines from
#' the \href{https://www.gnu.org/software/gsl/}{GSL library} to compute
#' Fermi--Dirac integrals are called. Otherwise, numerical integration is used.
#' @return A vector of size \code{length(mu)} or \code{length(s)} with the
#' values of the polylogarithm.
#' @examples
#' polylog_minus_exp_mu(mu = 1:5, s = 1)
#' polylog_minus_exp_mu(mu = 1, s = 1:5)
#' polylog_minus_exp_mu(mu = 1:5, s = 1:5)
#' @export
polylog_minus_exp_mu <- function(mu, s, upper = Inf, ...) {

  # gsl::fermi_dirac_int() computes the complete Fermi-Dirac integral
  # with an **integer** index of j:
  # F_j(x) = (1 / \Gamma(j + 1))
  #          \int_0^\infty (t^j /(\exp(t - x) + 1)) dt.
  # -Li_s(-exp(mu)) = (1 / \Gamma(s))
  #                   \int_0^\infty (t^(s - 1) /(\exp(t - mu) + 1)) dt
  #                 = F_{s - 1}(mu)
  # gsl::fermi_dirac_int() automatically rounds j to an integer when
  # returning F_j

  # Vectorization in s and mu
  l <- max(c(length(mu), length(s)))
  if (length(mu) == 1 && length(s) > 1) {

    mu <- rep(mu, length.out = l)

  } else if (length(s) == 1 && length(mu) > 1) {

    s <- rep(s, length.out = l)

  } else if (length(mu) != length(s)) {

    stop("mu and s must have equal length if both are not of length 1.")

  }
  poly <- numeric(l)

  # Indices
  ind_int <- which(abs(s - round(s)) < sqrt(.Machine$double.eps))
  ind_mhalf <- which(s == 1 / 2)
  ind_half <- which(s == 3 / 2)
  ind_3half <- which(s == 5 / 2)
  ind_other <- setdiff(seq_len(l), c(ind_int, ind_mhalf, ind_half, ind_3half))

  # F_{s - 1} for integer s
  if (length(ind_int) > 0) {

    poly[ind_int] <- -gsl::fermi_dirac_int(j = s[ind_int] - 1, x = mu[ind_int])

  }

  # F_{-1 / 2}
  if (length(ind_mhalf) > 0) {

    poly[ind_mhalf] <- -gsl::fermi_dirac_mhalf(x = mu[ind_mhalf])

  }

  # F_{1 / 2}
  if (length(ind_half) > 0) {

    poly[ind_half] <- -gsl::fermi_dirac_half(x = mu[ind_half])

  }

  # F_{3 / 2}
  if (length(ind_3half) > 0) {

    poly[ind_3half] <- -gsl::fermi_dirac_3half(x = mu[ind_3half])

  }

  # F_{s - 1} for non-integer s
  if (length(ind_other) > 0) {

    poly[ind_other] <- -sapply(ind_other, function(i)
      integrate(function(t) {
        exp((s[i] - 1) * log(t) - (lgamma(s[i]) + log1p(exp(t - mu[i]))))
        }, lower = 0, upper = upper, ...)$value)

  }

  # Polylogarithm
  return(poly)

}


#' @title Computes the integral \eqn{J_{d, k}}
#'
#' @description Computes the integral \eqn{J_{d, k}=\int_0^\infty
#' e^{2 \log(\log(1 + e^{k(1 - t)}))} t^{d / 2 - 1} dt}.
#'
#' @inheritParams eff_kern
#' @inheritParams polylog_minus_exp_mu
#' @return A vector of size \code{length(d)} with the values of the integral.
#' @examples
#' polykde:::J_d_k(d = 1:5, k = 10)
#' @noRd
J_d_k <- function(d, k = 10, upper = Inf, ...) {

  sapply(d, function(di)
    integrate(function(t)
      exp(2 * log(log1p(exp(k * (1 - t)))) + (di / 2 - 1) * log(t)),
      lower = 0, upper = upper, ...)$value)

}


#' @title Fast and stable evaluation of \eqn{\log(e^{-x} \mathcal{I}_{\nu}(x))}
#'
#' @description Computes a fast and stable approximation of the logarithm of
#' the scaled modified Bessel function of the first kind. Spline interpolation
#' is used for the tabulated orders \eqn{\nu = 0, 0.5, 1, \ldots, 24.5}; any
#' other order falls back to a direct evaluation.
#'
#' @param nu a scalar with the order of the Bessel function.
#' @param x vector with evaluation points for the Bessel function.
#' @param spline if \code{TRUE}, uses a fast spline interpolation for the values
#' of the Bessel function when \code{nu} is a single value in
#' \code{seq(0, 24.5, by = 0.5)}; for any other \code{nu} (off-grid or a vector
#' of orders), it falls back to the standard \code{\link[base]{besselI}}
#' computation. If \code{FALSE}, always uses \code{\link[base]{besselI}}.
#' Defaults to \code{FALSE}.
#' @details When \code{spline = TRUE} and \code{nu} is a single value in
#' \code{seq(0, 24.5, by = 0.5)}, the (scaled, log) Bessel function is
#' approximated by cubic spline interpolation in \code{u = log(x)} over \code{x}
#' in \code{[1e-5, 1e4]}; below \code{1e-5} an exact small-argument series
#' \eqn{\nu \log(x/2) - \log\Gamma(\nu + 1) - x} is used. For any other
#' \code{nu} (off-grid or vectorized) or \code{spline = FALSE}, the standard
#' \code{\link[base]{besselI}} is used.
#'
#' In both cases, for \code{x} larger than \code{1e4} a large-argument
#' asymptotic approximation with \code{\link[Bessel]{besselIasym}} is used. For
#' \code{nu} larger than \code{100} (only reachable when \code{spline = FALSE}),
#' a large-order asymptotic approximation with
#' \code{\link[Bessel]{besselI.nuAsym}} (DLMF 10.41) is used and takes
#' precedence over the large-argument one.
#' @return A vector of size \code{length(x)} with the evaluated function.
#' @examples
#' curve(polykde:::log_besselI_scaled(nu = 0.5, x = x, spline = TRUE),
#'       from = 0.01, to = 10, n = 200)
#' curve(polykde:::log_besselI_scaled(nu = 0.5, x = x, spline = FALSE),
#'       n = 200, add = TRUE, col = 2)
#' @noRd
log_besselI_scaled <- function(nu, x, spline = FALSE) {

  # Check length of nu and see if vectorized
  if (!(length(x) == 1 || length(nu) == 1 || length(nu) == length(x))) {

    stop("nu and x must have equal lengths or one being a single number.")

  } else {

    single_nu <- length(nu) == 1

  }

  # Remove NAs in x and propagate to nu
  x_na <- is.na(x)
  x <- x[!x_na]
  if (anyNA(nu)) {

    stop("nu must not contain NAs.")

  }
  if (all(x_na)) {

    return(as.numeric(rep(NA, max(length(x_na), max(length(nu))))))

  }
  if (!single_nu) {

    nu <- nu[!x_na]

  }

  # The spline tables only cover a single nu in bessel_nus; for any other nu
  # (off-grid or vectorized), set spline = FALSE.
  if (spline && !(single_nu && any(abs(bessel_nus - nu) < 1e-8))) {

    spline <- FALSE

  }

  # Spline interpolation and asymptotic approximations
  if (spline) {

    # Three regimes (spline requires a single nu, so x is never broadcast here):
    # analytic small-argument series below bessel_x_lo, cubic spline in
    # u = log(x) on [bessel_x_lo, bessel_x_hi), and large-argument asymptotics
    # at/above bessel_x_hi
    res <- rep(NA_real_, length(x))
    ind_small <- x < bessel_x_lo
    ind_asymp <- x >= bessel_x_hi
    ind_spl <- !ind_small & !ind_asymp

    # Cubic spline in u = log(x)
    if (any(ind_spl)) {

      key <- sprintf("%03d", round(10 * nu))
      res[ind_spl] <- get(paste0("log_besselI_scaled_spline_", key))(
        log(x[ind_spl]))

    }

    # Analytic small-argument series (handles x = 0; the nu = 0 branch avoids
    # 0 * log(0) = NaN)
    if (any(ind_small)) {

      xs <- x[ind_small]
      res[ind_small] <- if (nu == 0) -xs else
        nu * log(xs / 2) - lgamma(nu + 1) - xs

    }

    # Large-argument asymptotic approximation
    if (any(ind_asymp)) {

      res[ind_asymp] <- Bessel::besselIasym(x = x[ind_asymp],
                                            nu = nu, k.max = 10,
                                            expon.scaled = TRUE,
                                            log = TRUE)

    }

  } else {

    # Broadcast nu and x to a common length
    n <- max(length(x), length(nu))
    if (length(x) == 1) {

      x <- rep(x, n)

    }
    if (length(nu) == 1) {

      nu <- rep(nu, n)

    }

    # nu-asymptotics (large order) takes priority over x-asymptotics (large
    # argument), as the latter is invalid unless x >> nu^2
    ind_nu_asymp <- nu > 100
    ind_x_asymp <- !ind_nu_asymp & (x >= 10000)
    ind_std <- !ind_nu_asymp & !ind_x_asymp

    # Combine the three regimes
    res <- rep(NA_real_, n)
    if (any(ind_nu_asymp)) {

      # https://dlmf.nist.gov/10.41#E3
      res[ind_nu_asymp] <- Bessel::besselI.nuAsym(x = x[ind_nu_asymp],
                                                  nu = nu[ind_nu_asymp],
                                                  k.max = 5,
                                                  expon.scaled = TRUE,
                                                  log = TRUE)

    }
    if (any(ind_x_asymp)) {

      res[ind_x_asymp] <- Bessel::besselIasym(x = x[ind_x_asymp],
                                              nu = nu[ind_x_asymp], k.max = 10,
                                              expon.scaled = TRUE, log = TRUE)

    }
    if (any(ind_std)) {

      res[ind_std] <- log(besselI(x = x[ind_std], nu = nu[ind_std],
                                  expon.scaled = TRUE))

    }

  }

  # Add NAs to x
  if (single_nu) {

    res_na <- rep(NA, length(x_na))
    res_na[!x_na] <- res

  } else {

    res_na <- rep(NA, length(nu))
    res_na[!x_na] <- res

  }
  return(res_na)

}


#' @title Stable evaluation of \eqn{\log(\sum_{i=1}^n e^{x_i})}
#'
#' @description Computes \eqn{\log(\sum_{i=1}^n e^{x_i})} using the "LogSumExp
#' trick".
#'
#' @param logs vector with \eqn{x_i} values.
#' @param avg replace the sum by the average? Defaults to \code{FALSE}.
#' @return A vector of size \code{length(x)} with the evaluated function.
#' @examples
#' logs <- c(1e5, 1)
#' log(sum(exp(logs)))
#' polykde:::log_sum_exp(logs)
#' @noRd
log_sum_exp <- function(logs, avg = FALSE) {

  logs_M <- max(logs)
  return(logs_M + log(sum(exp(logs - logs_M))) - avg * log(length(logs)))

}


#' @title Stable evaluation of \eqn{\log(1 - \exp(-x))}
#'
#' @description Computes \eqn{\log(1 - \exp(-x))} accurately for \eqn{x\geq 0}.
#'
#' @param x vector with \eqn{x} values.
#' @return A vector of size \code{length(x)} with the evaluated function.
#' @examples
#' x <- c(1e-30, 1e-15, 1e-10, 10, 100)
#' polykde:::log1m_exp(x)
#' log(1 - exp(-x))
#' @noRd
log1m_exp <- function(x) {

  # Smaller than log(2)
  res <- rep(NA, length(x))
  ind_lt <- !is.na(x) & x <= log(2)
  res[ind_lt] <- log(-expm1(-x[ind_lt]))

  # Greater than log(2)
  ind_gt <- !is.na(x) & x > log(2)
  res[ind_gt] <- log1p(-exp(-x[ind_gt]))
  return(res)

}


#' @title Stable evaluation of \eqn{\exp(x) - \exp(y)}
#'
#' @description Computes \eqn{\exp(x) - \exp(y)} through log-scale and keeping
#' track of the sign.
#'
#' @param log_p vector with \eqn{x} values.
#' @param log_n vector with \eqn{y} values.
#' @param tol tolerance for considering the log-values equal.
#' @return A list with entries \code{log_abs} (vector with
#' \eqn{\log(|\exp(x) - \exp(y)|)} and \code{sgn} (vector with the signs of
#' \eqn{\exp(x) - \exp(y)}).
#' @examples
#' log_p <- c(10, 5, 1)
#' log_n <- rev(log_p)
#' log_diff <- polykde:::log_diff_exp(log_p = log_p, log_n = log_n)
#' log_diff$sgn * exp(log_diff$log_abs)
#' exp(log_p) - exp(log_n)
#' @noRd
log_diff_exp <- function(log_p, log_n, tol = 1e-15) {

  # Check lengths
  stopifnot(length(log_p) == length(log_n))

  # Signs (the near-zero tolerance branch must win over the strict sign checks)
  log_diff <- log_p - log_n
  sgn <- log_abs <- rep(NA, length(log_p))
  near_zero <- abs(log_diff) <= tol
  sgn[!near_zero & log_diff > 0] <- 1
  sgn[!near_zero & log_diff < 0] <- -1
  sgn[near_zero] <- 0

  # Positive case
  ind_sgn_pos <- !is.na(sgn) & sgn > 0
  log_abs[ind_sgn_pos] <- log_p[ind_sgn_pos] +
    log1m_exp(log_diff[ind_sgn_pos])

  # Negative case
  ind_sgn_neg <- !is.na(sgn) & sgn < 0
  log_abs[ind_sgn_neg] <- log_n[ind_sgn_neg] +
    log1m_exp(-log_diff[ind_sgn_neg])

  # Near-zero case: the difference is null at the working tolerance, so
  # log|exp(x) - exp(y)| = -Inf. Set it directly, as log1m_exp() would be
  # evaluated at a negative argument (giving NaN) whenever log_diff < 0
  ind_sgn_zero <- !is.na(sgn) & sgn == 0
  log_abs[ind_sgn_zero] <- -Inf
  return(list(log_abs = log_abs, sgn = sgn))

}


#' @title Stable evaluation of
#' \eqn{\operatorname{arcsinh}(sign(x) * \exp(\log(|x|)))}
#'
#' @description Computes
#' \eqn{\operatorname{arcsinh}(sign(x) * \exp(\log(|x|)))}. Useful to
#' logarithmically scale while preserving the sign of \eqn{x}.
#'
#' @param log_abs vector with \eqref{log(|x|)} values.
#' @param sgn vector with signs of \eqn{x}.
#' @return A vector of size \code{length(x)} with the evaluated function.
#' @examples
#' polykde:::asinh_log(log_abs = 10, sgn = 1)
#' asinh(exp(10))
#' polykde:::asinh_log(log_abs = 10, sgn = -1)
#' asinh(-exp(10))
#' @noRd
asinh_log <- function(log_abs, sgn) {

  # We use that for s in {-1, 1},
  # asinh(s * |x|) = s * log(|x| + sqrt(x^2 + 1))
  # and hence
  # asinh(s * exp(log(|x|))) = s * log(exp(log(|x|)) +
  #                             sqrt(exp(2 * log(|x|)) + 1))
  #                          = s * log(exp(log(|x|)) *
  #                             (1 + sqrt(1 + exp(-2 * log(|x|)))))
  #                          = s * (log(|x|) +
  #                             log1p(sqrt(1 + exp(-2 * log(|x|)))))
  # for log(|x|) >= 0.
  # Indeed, if s = 1:
  # asinh(s * |x|) = log(s * |x| + sqrt(x^2 + 1))
  #                = log(|x| + sqrt(x^2 + 1))
  # If s = -1:
  # asinh(s * |x|) = log(s * |x| + sqrt(x^2 + 1))
  #                = log(-|x| + sqrt(x^2 + 1))
  #                = log(1 / (|x| + sqrt(x^2 + 1)))
  #                = -log(|x| + sqrt(x^2 + 1)),
  # since 1 / (y + sqrt(y^2 + 1)) = -y + sqrt(y^2 + 1)
  # For log(|x|) < 0, use
  # asinh(s * exp(log(|x|))) = s * (log(|x|) +
  #                             log(1 + sqrt(1 + exp(-2 * log(|x|)))))
  #                          = s * log(exp(log(|x|)) +
  #                              sqrt(exp(2 * log(|x|)) + 1))
  # ash <- sgn * (log_abs + log1p(sqrt(1 + exp(-2 * log_abs))))

  # Check lengths
  stopifnot(length(log_abs) == length(sgn))

  # Positive case
  ash <- rep(NA, length(log_abs))
  log_abs_pos <- !is.na(log_abs) & log_abs >= 0
  ash[log_abs_pos] <-
    sgn[log_abs_pos] * (log_abs[log_abs_pos] +
                          log1p(sqrt(1 + exp(-2 * log_abs[log_abs_pos]))))

  # Negative case
  log_abs_neg <- !is.na(log_abs) & log_abs < 0
  ash[log_abs_neg] <-
    sgn[log_abs_neg] * log(exp(log_abs[log_abs_neg]) +
                              sqrt(1 + exp(2 * log_abs[log_abs_neg])))
  ash[sgn == 0] <- 0
  return(ash)

}
