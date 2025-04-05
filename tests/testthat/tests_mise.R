
test_that("exact_mise_vmf() works properly", {

  # Parameters
  M <- 1e4
  n <- 2
  d <- 1
  m <- 1
  mu <- rbind(c(1, 0))
  kappa <- 5
  prop <- 1
  h <- 1

  # Sample X's and evaluate density
  set.seed(42)
  N1 <- 1e3
  X <- lapply(seq_along(prop), function(m)
    r_vmf_polysph(n = round(N1 * prop[m]), d = d, mu = mu[m, ],
                  kappa = kappa))
  X <- do.call(rbind, X)
  f_X <- drop(kde_polysph(x = X, X = mu, d = d, h = 1 / sqrt(kappa),
                          weights = prop))

  # Sample Y's to evaluate the kde
  N2 <- 1e3
  Y <- lapply(seq_len(N2), function(nj) {
    Y_j <- lapply(seq_along(prop), function(m) {
      r_vmf_polysph(n = round(n * prop[m]), d = d, mu = mu[m, ], kappa = kappa)
    })
    do.call(rbind, Y_j)
  })

  # Simulate sample and compute kde
  kde_f_2 <- sapply(seq_len(N2), function(k) {
    kde <- tryCatch(kde_polysph(x = X, X = Y[[k]][1:n, ], d = d, h = h),
                    error = function(e) NA)
    (kde - f_X)^2
  })
  kde_f_2 <- rowMeans(kde_f_2, na.rm = TRUE)

  expect_equal(exact_mise_vmf(h = h, n = n, mu = mu, kappa = kappa, prop = prop,
                              d = d, seed_psi = 42, spline = TRUE)$mise,
               mean(kde_f_2 / f_X),
               tolerance = 1e-2)

})

test_that("exact_mise_vmf_polysph() and exact_mise_vmf() equal on the sphere", {

  h <- 0.5
  expect_equal(exact_mise_vmf(h = h, n = 100, mu = rbind(c(0, 1), c(1, 0)),
                              kappa = c(5, 2), prop = c(0.7, 0.3), d = 1,
                              seed_psi = 1, spline = TRUE),
               exact_mise_vmf_polysph(h = h, n = 100,
                                      mu = rbind(c(0, 1), c(1, 0)),
                                      kappa = c(5, 2), prop = c(0.7, 0.3),
                                      d = 1, seed_psi = 1, spline = TRUE))

})

test_that("bw_mise_polysph() minimizes the MISE on the sphere", {

  r <- 1
  m <- rpois(1, 3) + 1
  d <- rpois(r, 3) + 1
  mu <- r_unif_polysph(n = m, d = d)
  kappa <- matrix(abs(rnorm(m * r, sd = 2)), nrow = m, ncol = r)
  prop <- runif(m)
  prop <- prop / sum(prop)

  n <- 10
  bw0 <- cbind(10^seq(log10(0.1), log10(5), l = 10))
  log1p_mise_bw0 <- sapply(bw0, function(h) {
    log1p_mise_exact(log_h = log(h), n = n, mu = mu, kappa = kappa, prop = prop,
                     d = d, seed_psi = 1, spline = TRUE)
  })
  bw_mise <- bw_mise_polysph(n = n, d = d, bw0 = bw0, mu = mu, kappa = kappa,
                             prop = prop, seed_psi = 1, spline = TRUE)
  plot(bw0, log1p_mise_bw0, type = "o", xlab = "h", ylab = "log1p_mise",
       xlim = range(c(bw0, bw_mise$bw)))
  points(bw_mise$bw, bw_mise$opt$minimum, col = "red", pch = 19)
  expect_lte(bw_mise$opt$minimum, min(log1p_mise_bw0))

})
