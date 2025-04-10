
# Mesh for interpolation
x_bessel <- c(10^seq(-10, -5, by = 1),
              seq(1e-4, 1, by = 1e-4),
              seq(1 + 1e-2, 10, by = 1e-3),
              seq(10 + 1e-1, 100, by = 1e-2),
              seq(100 + 1e0, 1e3, by = 1e0),
              seq(1000 + 1e1, 5e4, by = 2e1))

# Evaluate log-Bessels
nus <- seq(0, 60, by = 5)
nus_char <- sprintf("%02d", nus)
for (nu_i in nus_char) {

  # Call Bessel function
  cat("nu =", nu_i, "\n")
  bessel_nu_i <- log(besselI(x = x_bessel, nu = as.integer(nu_i) / 10,
                             expon.scaled = TRUE))

  # Change -Inf to -1e6
  bessel_nu_i[is.infinite(bessel_nu_i)] <- -1e6

  # Save object
  assign(x = paste0("log_besselI_scaled_", nu_i, "_grid"),
         value = bessel_nu_i)

}

# Save tables
save(list = c("x_bessel", paste0("log_besselI_scaled_", nus_char, "_grid")),
     file = "../R/sysdata.rda", compress = "xz")
