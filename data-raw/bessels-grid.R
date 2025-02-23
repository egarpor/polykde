
# Mesh for interpolation
x_bessel <- c(10^seq(-10, -5, by = 1), seq(1e-4, 1, by = 1e-4),
              seq(1 + 1e-2, 10, by = 1e-3), seq(10 + 1e-1, 100, by = 1e-2),
              seq(100 + 1e0, 1e3, by = 1e0), seq(1000 + 1e1, 5e4, by = 2e1))

# Evaluated log-Bessels
log_besselI_scaled_00_grid <-
  log(besselI(x = x_bessel, nu = 0, expon.scaled = TRUE))
log_besselI_scaled_05_grid <-
  log(besselI(x = x_bessel, nu = 0.5, expon.scaled = TRUE))
log_besselI_scaled_10_grid <-
  log(besselI(x = x_bessel, nu = 1, expon.scaled = TRUE))
log_besselI_scaled_15_grid <-
  log(besselI(x = x_bessel, nu = 1.5, expon.scaled = TRUE))
log_besselI_scaled_20_grid <-
  log(besselI(x = x_bessel, nu = 2, expon.scaled = TRUE))
log_besselI_scaled_25_grid <-
  log(besselI(x = x_bessel, nu = 2.5, expon.scaled = TRUE))
log_besselI_scaled_30_grid <-
  log(besselI(x = x_bessel, nu = 3, expon.scaled = TRUE))
log_besselI_scaled_35_grid <-
  log(besselI(x = x_bessel, nu = 3.5, expon.scaled = TRUE))
log_besselI_scaled_40_grid <-
  log(besselI(x = x_bessel, nu = 4, expon.scaled = TRUE))
log_besselI_scaled_45_grid <-
  log(besselI(x = x_bessel, nu = 4.5, expon.scaled = TRUE))

# Change -Inf to -1e6
log_besselI_scaled_00_grid[is.infinite(log_besselI_scaled_00_grid)] <- -1e6
log_besselI_scaled_05_grid[is.infinite(log_besselI_scaled_05_grid)] <- -1e6
log_besselI_scaled_10_grid[is.infinite(log_besselI_scaled_10_grid)] <- -1e6
log_besselI_scaled_15_grid[is.infinite(log_besselI_scaled_15_grid)] <- -1e6
log_besselI_scaled_20_grid[is.infinite(log_besselI_scaled_20_grid)] <- -1e6
log_besselI_scaled_25_grid[is.infinite(log_besselI_scaled_25_grid)] <- -1e6
log_besselI_scaled_30_grid[is.infinite(log_besselI_scaled_30_grid)] <- -1e6
log_besselI_scaled_35_grid[is.infinite(log_besselI_scaled_35_grid)] <- -1e6
log_besselI_scaled_40_grid[is.infinite(log_besselI_scaled_40_grid)] <- -1e6
log_besselI_scaled_45_grid[is.infinite(log_besselI_scaled_45_grid)] <- -1e6

# Save tables
save(list = c("x_bessel",
              "log_besselI_scaled_00_grid", "log_besselI_scaled_05_grid",
              "log_besselI_scaled_10_grid", "log_besselI_scaled_15_grid",
              "log_besselI_scaled_20_grid", "log_besselI_scaled_25_grid",
              "log_besselI_scaled_30_grid", "log_besselI_scaled_35_grid",
              "log_besselI_scaled_40_grid", "log_besselI_scaled_45_grid"),
     file = "../R/sysdata.rda", compress = "xz")

# Accuracy interpolations
x <- seq(1e-8, 1e4, l = 1e3)
summary(log_besselI_scaled(nu = 0, x = x, spline = TRUE) -
          log_besselI_scaled(nu = 0, x = x, spline = FALSE))
summary(log_besselI_scaled(nu = 0.5, x = x, spline = TRUE) -
          log_besselI_scaled(nu = 0.5, x = x, spline = FALSE))
summary(log_besselI_scaled(nu = 1, x = x, spline = TRUE) -
          log_besselI_scaled(nu = 1, x = x, spline = FALSE))
summary(log_besselI_scaled(nu = 2, x = x, spline = TRUE) -
          log_besselI_scaled(nu = 2, x = x, spline = FALSE))
summary(log_besselI_scaled(nu = 3.5, x = x, spline = TRUE) -
          log_besselI_scaled(nu = 3.5, x = x, spline = FALSE))
summary(log_besselI_scaled(nu = 4, x = x, spline = TRUE) -
          log_besselI_scaled(nu = 4, x = x, spline = FALSE))
