
## von Mises--Fisher density

test_that("Check d_vmf_polysph by inverse weighted sampling", {

  M <- 1e3
  d <- c(1, 2)
  mu <- c(0, 1, 0, 1, 0)
  kappa <- c(1, 1)
  x <- r_vmf_polysph(n = M, d = d, mu = mu, kappa = kappa)
  dens <- d_vmf_polysph(x = x, d = d, mu = mu, kappa = kappa)
  expect_equal(mean(1 / dens) / prod(rotasym::w_p(p = d + 1)), 1,
               tolerance = 1e-1)

})

test_that("fast_log_c_vMF() works properly with spline = TRUE / FALSE", {

  expect_equal(polykde:::fast_log_c_vMF(p = 5, kappa = 0:10, spline = TRUE),
               polykde:::fast_log_c_vMF(p = 5, kappa = 0:10, spline = FALSE))

})

## Uniform density

test_that("Check d_unif_polysph by inverse weighted sampling", {

  M <- 1e3
  d <- c(1, 2)
  x <- r_unif_polysph(n = M, d = d)
  dens <- d_unif_polysph(x = x, d = d)
  expect_equal(mean(1 / dens) / prod(rotasym::w_p(p = d + 1)), 1,
               tolerance = 1e-1)

})

test_that("Check d_unif_polysph with log", {

  M <- 1e3
  d <- c(1, 2)
  x <- c(0, 1, 0, 1, 0)
  expect_equal(d_unif_polysph(x = x, d = d, log = TRUE),
               log(d_unif_polysph(x = x, d = d)))

})
