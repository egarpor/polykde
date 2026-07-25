
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

  expect_equal(fast_log_c_vMF(p = 5, kappa = 0:10, spline = TRUE),
               fast_log_c_vMF(p = 5, kappa = 0:10, spline = FALSE))
})

test_that("fast_log_c_vMF(spline = TRUE) works for high p (off-grid nu)", {
  expect_no_error(fast_log_c_vMF(p = 60, kappa = 1:5, spline = TRUE))
  expect_equal(fast_log_c_vMF(p = 60, kappa = 1:5, spline = TRUE),
               fast_log_c_vMF(p = 60, kappa = 1:5, spline = FALSE))
})

## Mixture of von Mises--Fisher densities

test_that("Check d_mvmf_polysph by inverse weighted sampling", {

  M <- 1e3
  d <- c(1, 2)
  mu <- rbind(c(0, 1, 0, 1, 0), c(0, -1, 0, -1, 0))
  kappa <- rbind(c(1, 1), c(2, 1))
  prop <- c(0.7, 0.3)
  x <- r_mvmf_polysph(n = M, d = d, mu = mu, kappa = kappa, prop = prop)
  dens <- d_mvmf_polysph(x = x, d = d, mu = mu, kappa = kappa, prop = prop)
  expect_equal(mean(1 / dens) / prod(rotasym::w_p(p = d + 1)), 1,
               tolerance = 1e-1)

})

test_that("Mixture functions accept a length-r kappa vector when m = 1", {

  d <- c(1, 2)
  mu <- c(0, 1, 0, 1, 0)
  kappa <- c(5, 2)
  x <- r_mvmf_polysph(n = 5, d = d, mu = mu, kappa = kappa, prop = 1)
  expect_equal(dim(x), c(5L, sum(d + 1)))
  dens <- d_mvmf_polysph(x = x, d = d, mu = mu, kappa = kappa, prop = 1)
  expect_length(dens, 5)
  expect_true(all(is.finite(dens)))

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
