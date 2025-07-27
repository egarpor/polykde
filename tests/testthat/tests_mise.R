
test_that("mise_vmf() does the job on computing exactly the MISE", {

  # Parameters
  M <- 1e4
  n <- 3
  d <- rpois(1, lambda = 1) + 1
  m <- 2
  kappa <- rep(5, m)
  mu <- r_unif_polysph(n = m, d = d)
  prop <- rep(0.5, m)
  h <- 1

  # Sample X's and evaluate density
  set.seed(42)
  N1 <- 1e3
  X <- r_mvmf_polysph(n = N1, d = d, mu = mu, kappa = kappa, prop = prop)
  f_X <- drop(d_mvmf_polysph(x = X, d = d, mu = mu, kappa = kappa, prop = prop))

  # Sample Y's to evaluate the kde
  N2 <- 1e3
  Y <- lapply(seq_len(N2), function(nj) {
    r_mvmf_polysph(n = n, d = d, mu = mu, kappa = kappa, prop = prop)
  })

  # Simulate sample and compute kde
  kde_f_2 <- sapply(seq_len(N2), function(k) {
    kde <- tryCatch(kde_polysph(x = X, X = Y[[k]][1:n, ], d = d, h = h),
                    error = function(e) NA)
    (kde - f_X)^2
  })
  kde_f_2 <- rowMeans(kde_f_2, na.rm = TRUE)

  # Exact vs. Monte Carlo
  expect_equal(drop(mise_vmf(h = h, n = n, mu = mu, kappa = kappa,
                             prop = prop, d = d, seed_psi = 42,
                             spline = TRUE)$mise),
               mean(kde_f_2 / f_X),
               tolerance = 1e-2)

})

test_that("mise_vmf() is properly vectorized in h and n", {

  h <- c(0.1, 0.2)
  n <- c(10, 20)
  mu <- rbind(c(0, 1), c(1, 0))
  kappa <- c(5, 2)
  prop <- c(0.7, 0.3)
  expect_equal(mise_vmf(h = h, n = n, mu = mu, kappa = kappa,
                        prop = prop, d = 1, seed_psi = 1, M_psi = 10,
                        spline = TRUE)$mise,
               c(mise_vmf(h = h[1], n = n[1], mu = mu, kappa = kappa,
                          prop = prop, d = 1, seed_psi = 1, M_psi = 10,
                          spline = TRUE)$mise,
                 mise_vmf(h = h[2], n = n[2], mu = mu, kappa = kappa,
                          prop = prop, d = 1, seed_psi = 1, M_psi = 10,
                          spline = TRUE)$mise))
  expect_equal(mise_vmf(h = h, n = n[1], mu = mu, kappa = kappa,
                        prop = prop, d = 1, seed_psi = 1, M_psi = 10,
                        spline = TRUE)$mise,
               c(mise_vmf(h = h[1], n = n[1], mu = mu, kappa = kappa,
                          prop = prop, d = 1, seed_psi = 1, M_psi = 10,
                          spline = TRUE)$mise,
                 mise_vmf(h = h[2], n = n[1], mu = mu, kappa = kappa,
                          prop = prop, d = 1, seed_psi = 1, M_psi = 10,
                          spline = TRUE)$mise))

})

test_that("mise_vmf_polysph() and mise_vmf() equal on the sphere", {

  h <- 0.5
  expect_equal(mise_vmf(h = h, n = 100, mu = rbind(c(0, 1), c(1, 0)),
                        kappa = c(5, 2), prop = c(0.7, 0.3), d = 1,
                        seed_psi = 1, spline = TRUE),
               mise_vmf_polysph(h = h, n = 100,
                                mu = rbind(c(0, 1), c(1, 0)),
                                kappa = c(5, 2), prop = c(0.7, 0.3),
                                d = 1, seed_psi = 1, spline = TRUE))

})

test_that("bw_mise_polysph() minimizes the MISE on the sphere", {

  # Parameters
  r <- 1
  m <- rpois(1, 3) + 1
  d <- rpois(r, 3) + 1
  mu <- r_unif_polysph(n = m, d = d)
  kappa <- matrix(abs(rnorm(m * r, sd = 2)), nrow = m, ncol = r)
  prop <- runif(m)
  prop <- prop / sum(prop)
  n <- 5

  # Minimization of ISE
  bw0 <- cbind(10^seq(log10(0.2), log10(5), l = 10))
  log1p_mise_bw0 <- sapply(bw0, function(h) {
    log1p_mise(log_h = log(h), n = n, mu = mu, kappa = kappa, prop = prop,
               d = d, seed_psi = 1, spline = TRUE)
  })
  bw_mise <- bw_mise_polysph(n = n, d = d, bw0 = bw0, mu = mu, kappa = kappa,
                             prop = prop, seed_psi = 1, spline = TRUE)
  plot(bw0, log1p_mise_bw0, type = "o", xlab = "h", ylab = "log1p_mise",
       xlim = range(c(bw0, bw_mise$bw)))
  points(bw_mise$bw, bw_mise$opt$minimum, col = "red", pch = 19)
  expect_lte(bw_mise$opt$minimum, min(log1p_mise_bw0))

})

test_that("ise_vmf_polysph() does the job on computing exactly the ISE", {

  # Parameters and sample
  n <- 3
  d <- rpois(1, lambda = 1) + 1
  m <- 2
  kappa <- rep(5, m)
  mu <- r_unif_polysph(n = m, d = d)
  prop <- rep(0.5, m)
  h <- seq(0.1, 1, l = 10)
  X <- r_mvmf_polysph(n = n, d = d, mu = mu, kappa = kappa, prop = prop)

  # Exact vs. Monte Carlo
  expect_equal(
    ise_vmf_polysph(X = X, d = d, h = h, mu = mu, kappa = kappa,
                    prop = prop, spline = TRUE, exact = TRUE)$ise,
    ise_vmf_polysph(X = X, d = d, h = h, mu = mu, kappa = kappa,
                    prop = prop, seed_psi = 42, M_psi = 1e4,
                    spline = TRUE, exact = FALSE)$ise,
    tolerance = 5e-2
  )

})

test_that("bw_ise_polysph() minimizes the ISE on the sphere", {

  # Parameters and sample
  r <- 1
  m <- rpois(1, 3) + 1
  d <- rpois(r, 3) + 1
  mu <- r_unif_polysph(n = m, d = d)
  kappa <- matrix(abs(rnorm(m * r, sd = 2)), nrow = m, ncol = r)
  prop <- runif(m)
  prop <- prop / sum(prop)
  n <- 5
  X <- r_mvmf_polysph(n = n, d = d, mu = mu, kappa = kappa, prop = prop)

  # Minimization of ISE
  bw0 <- cbind(10^seq(log10(0.2), log10(5), l = 20))
  log1p_ise_bw0 <- sapply(bw0, function(h) {
    log1p_ise(log_h = log(h), X = X, mu = mu, kappa = kappa, prop = prop,
              d = d, spline = TRUE, exact = TRUE)
  })
  bw_ise <- bw_ise_polysph(X = X, d = d, bw0 = bw0, mu = mu, kappa = kappa,
                           prop = prop, seed_psi = 1, spline = TRUE,
                           exact = TRUE)
  plot(bw0, log1p_ise_bw0, type = "o", xlab = "h", ylab = "log1p_mise",
       xlim = range(c(bw0, bw_ise$bw)))
  points(bw_ise$bw, bw_ise$opt$minimum, col = "red", pch = 19)
  expect_lte(bw_ise$opt$minimum, min(log1p_ise_bw0))

})
