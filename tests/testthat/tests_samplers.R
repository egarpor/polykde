
# Randomize testing
r <- 2 # sample(1:3, size = 1)
d <- rep(2, 2) # rpois(r, 2) + 1
h <- runif(r, 0.25, 0.50)
mu <- r_unif_polysph(n = 1, d = d)
X <- r_kern_polysph(n = 10, d = d, mu = mu, h = rep(0.2, r), kernel = 1)
n <- 300

test_that("r_unif_polysph", {

  ind <- cumsum(c(1, d + 1))
  for (j in seq_len(r)) {
    ind_j <- ind[j]:(ind[j + 1] - 1)
    pval <- sphunif::unif_test(data = r_unif_polysph(n = n, d = d)[, ind_j],
                               type = "Bingham")$p.value
    expect_gt(pval, 0.01)
  }

})

test_that("Check r_kern_polysph", {

  ind <- cumsum(c(1, d + 1))
  for (kernel in 1:2) {
    data <- r_kern_polysph(n = n, d = d, mu = mu, h = h, kernel = kernel)
    for (j in seq_len(r)) {
      ind_j <- ind[j]:(ind[j + 1] - 1)
      pval <- sphunif::unif_test(data = data[, ind_j],
                                 type = "Bingham")$p.value
      expect_lt(pval, 0.10)
    }
  }

})

test_that("Check r_kern_polysph by inverse weighted sampling", {

  h <- c(2, 2)
  for (kernel in 1:3) {

    data <- r_kern_polysph(n = n, d = d, mu = mu, h = h, kernel = kernel,
                           kernel_type = 1)
    dens <- kde_polysph(x = data, X = mu, d = d, h = h, kernel = kernel,
                        kernel_type = 1)
    inv_dens <- 1 / dens
    expect_equal(mean(inv_dens[is.finite(inv_dens)]) /
                   prod(rotasym::w_p(p = d + 1)), 1,
                 tolerance = 1e-1)

  }

})

test_that("Check r_kern_polysph by inverse weighted sampling", {

  h <- c(2, 2)
  for (kernel in 1:3) {

    data <- r_kern_polysph(n = n, d = d, mu = mu, h = h, kernel = kernel,
                           kernel_type = 2)
    dens <- kde_polysph(x = data, X = mu, d = d, h = h, kernel = kernel,
                        kernel_type = 2)
    inv_dens <- 1 / dens
    expect_equal(mean(inv_dens[is.finite(inv_dens)]) /
                   prod(rotasym::w_p(p = d + 1)), 1,
                 tolerance = 1e-1)

  }

})

test_that("Check coherency between kde_polysph(), c_kern(), and L()", {

  for (kernel in 1:3) {

    data <- r_kern_polysph(n = n, d = d, mu = mu, h = h, kernel = 1)
    dens1 <- kde_polysph(x = data, X = mu, d = d, h = h, kernel = kernel,
                         kernel_type = 1)
    dens2 <- prod(c_kern(h = h, d = d, kernel = kernel, kernel_type = 1)) *
      L((1 - data[, 1:3] %*% mu[1:3]) / h[1]^2, kernel = kernel) *
      L((1 - data[, 4:6] %*% mu[4:6]) / h[2]^2, kernel = kernel)
    expect_equal(dens1, dens2, tolerance = 1e-2)

    dens1 <- kde_polysph(x = data, X = mu, d = d, h = h, kernel = kernel,
                         kernel_type = 2)
    dens2 <- c_kern(h = h, d = d, kernel = kernel, kernel_type = 2) *
      L((1 - data[, 1:3] %*% mu[1:3]) / h[1]^2 +
          (1 - data[, 4:6] %*% mu[4:6]) / h[2]^2, kernel = kernel)
    expect_equal(dens1, dens2, tolerance = 1e-2)

  }

})
