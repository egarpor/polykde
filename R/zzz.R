
.onLoad <- function(libname = find.package("polykde"), pkgname = "polykde") {

  # Build splines in u = log(x) from (u_bessel, bessel_nus, logI_bessel) in
  # sysdata.rda
  env <- asNamespace(pkgname)
  for (j in seq_along(bessel_nus)) {

    key <- sprintf("%03d", round(10 * bessel_nus[j]))
    assign(x = paste0("log_besselI_scaled_spline_", key),
           value = splinefun(x = u_bessel, y = logI_bessel[, j],
                             method = "fmm"),
           pos = env)

  }

  # CRAN NOTE avoidance
  if (getRversion() >= "2.15.1") {

    utils::globalVariables(c("u_bessel", "bessel_nus", "bessel_x_lo",
                             "bessel_x_hi", "logI_bessel"))

  }
  invisible()

}
