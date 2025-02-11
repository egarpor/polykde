
#' @title Stable computation of the softplus function
#'
#' @description Computes the softplus function \eqn{\log(1+e^{t})} in a
#' numerically stable way for large absolute values of \eqn{t}.
#'
#' @param t numeric vector or matrix.
#' @return The softplus function evaluated at \code{t}.
#' @examples
#' h <- 0.1
#' k <- 1
#' curve(polykde:::softplus(k * (1 - (1 - x) / h^2)), from = -1, to = 1)
#' @keywords internal
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
  for (i in 2:n) {
    for (j in seq_len(i - 1)) {

      dij[i, j] <- drop(dist_polysph(x = x[i, , drop = FALSE],
                                     y = x[j, , drop = FALSE],
                                     ind_dj = ind_dj, norm_x = norm_x,
                                     norm_y = norm_y, std = std))

    }
  }
  dij <- s(dij, add = TRUE)
  colnames(dij) <- rownames(dij) <- 1:n
  return(dij)

}


#' @title Index of hyperspheres on a polysphere
#'
#' @description Given Cartesian coordinates of polyspherical data, computes
#' the \code{0}-based indexes at which the Cartesian coordinates for each
#' hypersphere start and end.
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
#' lists <- list(list(1:3, 4:6), list(7:9, 10:12))
#' polykde:::bind_lists(lists = lists, bind = "rbind")
#' @keywords internal
bind_lists <- function(lists, bind = "rbind") {

  stopifnot(bind %in% c("rbind", "cbind"))
  n <- length(lists)
  merged <- lists[[which(!sapply(lists, is.null))[1]]]
  for (j in seq_len(length(merged))) {

    bind_j <- ifelse(length(dim(merged[[j]])) == 3,
                     function(...)
                       abind::abind(..., along = ifelse(bind == "rbind", 1, 2)),
                     bind)
    merged[[j]] <- do.call(what = bind_j,
                           args = lapply(1:n, function(i) lists[[i]][[j]]))

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
#' @param upper upper limit of integration.
#' @details If \code{s} is an integer, 1/2, 3/2, or 5/2, then routines from
#' the \href{https://www.gnu.org/software/gsl/}{GSL library} to compute
#' Fermi--Dirac integrals are called. Otherwise, numerical integration is used.
#' @return A vector of size \code{length(mu)} or \code{length(s)} with the
#' values of the polylogarithm.
#' @examples
#' polykde:::polylog_minus_exp_mu(mu = 1:5, s = 1)
#' polykde:::polylog_minus_exp_mu(mu = 1, s = 1:5)
#' polykde:::polylog_minus_exp_mu(mu = 1:5, s = 1:5)
#' @keywords internal
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
        # t^(s[i] - 1) / (gamma(s[i]) * (exp(t - mu[i]) + 1))
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
#' @keywords internal
J_d_k <- function(d, k = 10, upper = Inf, ...) {

  sapply(d, function(di)
    integrate(function(t)
      exp(2 * log(log1p(exp(k * (1 - t)))) + (di / 2 - 1) * log(t)),
      lower = 0, upper = upper, ...)$value)

}
