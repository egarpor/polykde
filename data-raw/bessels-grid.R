
# Log-Bessel spline tables consumed by log_besselI_scaled(..., spline = TRUE)
#
# The (scaled, log) Bessel function log(exp(-x) * I_nu(x)) is smooth in
# u = log(x), so it is interpolated on a uniform grid in u. Below x_lo an
# analytic small-argument series is used and above x_hi a large-argument
# asymptotic (Bessel::besselIasym) is used, so the stored spline only needs to
# cover [x_lo, x_hi].

# Grid parameters
bessel_nus <- seq(0, 24.5, by = 0.5)
bessel_x_lo <- 1e-5
bessel_x_hi <- 1e4
bessel_N <- 1400

# Uniform knots in u = log(x)
u_bessel <- seq(log(bessel_x_lo), log(bessel_x_hi), length.out = bessel_N)
x_bessel <- exp(u_bessel)

# Analytic leading term of log(exp(-x) * I_nu(x)) as x -> 0
log_besselI_series <- function(x, nu) {
  if (nu == 0) -x else nu * log(x / 2) - lgamma(nu + 1) - x
}

# Evaluate the (scaled, log) Bessel on the grid, one column per order
logI_bessel <- matrix(NA_real_, nrow = bessel_N, ncol = length(bessel_nus))
for (j in seq_along(bessel_nus)) {

  nu <- bessel_nus[j]
  cat("nu =", nu, "\n")
  y <- log(besselI(x = x_bessel, nu = nu, expon.scaled = TRUE))

  # Fill any underflow (-Inf) with the analytic series for small x. With
  # x_lo = 1e-5 and nu <= 24.5 this never triggers, but it keeps the
  # builder robust if x_lo or the order range are pushed further.
  bad <- !is.finite(y) & x_bessel < 1e-2
  if (any(bad)) {

    y[bad] <- log_besselI_series(x_bessel[bad], nu)

  }
  stopifnot(all(is.finite(y)))
  logI_bessel[, j] <- y

}

# Save tables
save(list = c("u_bessel", "bessel_nus", "bessel_x_lo", "bessel_x_hi",
              "logI_bessel"),
     file = "../R/sysdata.rda", compress = "xz")
