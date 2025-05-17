
.onLoad <- function(libname = find.package("polykde"), pkgname = "polykde") {

  # Assign global variables
  env <- asNamespace(pkgname)
  nus <- seq(0, 60, by = 5)
  nus_char <- sprintf("%02d", nus)
  for (nu_i in nus_char) {

    assign(x = paste0("log_besselI_scaled_spline_", nu_i),
           value =
             splinefun(x = x_bessel,
                       y = get(paste0("log_besselI_scaled_", nu_i, "_grid"),
                               envir = env)),
           pos = env)

  }

  # CRAN NOTE avoidance
  if (getRversion() >= "2.15.1") {

    utils::globalVariables(c("x_bessel",
                             paste0("log_besselI_scaled_", nus_char, "_grid")))

  }
  invisible()

}
