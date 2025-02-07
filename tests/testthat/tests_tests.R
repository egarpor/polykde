
## Well-computed JSD for k = 2

r <- 3
d <- rep(2, r)
h1 <- rep(0.25, r)
h2 <- h1
n1 <- 100
n2 <- 50
n <- n1 + n2
p <- c(n1, n2) / n
mu1 <- r_unif_polysph(n = 1, d = d)
mu2 <- r_unif_polysph(n = 1, d = d)
kappa1 <- rep(1, r)
kappa2 <- rep(3, r)
X1 <- r_vmf_polysph(n = n1, mu = mu1, kappa = kappa1, d = d)
X2 <- r_vmf_polysph(n = n2, mu = mu2, kappa = kappa2, d = d)
X <- rbind(X1, X2)
labels <- rep(c(TRUE, FALSE), times = c(n1, n2))
M <- 1e4

# Monte Carlo divergence
samp_1 <- r_kde_polysph(n = M, X = X1, d = d, h = h1)
samp_2 <- r_kde_polysph(n = M, X = X2, d = d, h = h2)
samp_0 <- rbind(samp_1[1:round(p[1] * M), ], samp_2[1:round(p[2] * M), ])
H_1 <- -mean(kde_polysph(x = samp_1, X = X1, d = d, h = h1, log = TRUE))
H_2 <- -mean(kde_polysph(x = samp_2, X = X2, d = d, h = h2, log = TRUE))
H_0 <- -mean(log(p[1] * kde_polysph(x = samp_0, X = X1, d = d, h = h1) +
                   p[2] * kde_polysph(x = samp_0, X = X2, d = d, h = h2)))

# Monte Carlo with sample (naive)
H_1_n <- -mean(kde_polysph(x = X1, X = X1, d = d, h = h1, log = TRUE))
H_2_n <- -mean(kde_polysph(x = X2, X = X2, d = d, h = h2, log = TRUE))
H_0_n <- -mean(c(log(p[1] * kde_polysph(x = X, X = X1, d = d, h = h1) +
                        p[2] * kde_polysph(x = X, X = X2, d = d, h = h2))))

# Cross-validation (with two approaches for f0_cv)
H_1_cv <- -mean(log_cv_kde_polysph(X = X1, d = d, h = h1))
H_2_cv <- -mean(log_cv_kde_polysph(X = X2, d = d, h = h2))
H_0_cv <- -mean(c(log(p[1] * exp(log_cv_kde_polysph(X = X1, d = d, h = h1)) +
                        p[2] * kde_polysph(x = X1, X = X2, d = d, h = h2)),
                  log(p[1] * kde_polysph(x = X2, X = X1, d = d, h = h1) +
                        p[2] * exp(log_cv_kde_polysph(X = X2, d = d, h = h2)))))
H_0_cv_bis <- -mean(log(c(
  1 / (1 - p[1] / n1) *
    (((n1 - 1) / n1) * p[1] * exp(log_cv_kde_polysph(X = X1, d = d, h = h1)) +
       p[2] * kde_polysph(x = X1, X = X2, d = d, h = h2)),
  1 / (1 - p[2] / n2) *
    (p[1] * kde_polysph(x = X2, X = X1, d = d, h = h1) +
       ((n2 - 1) / n2) * p[2] * exp(log_cv_kde_polysph(X = X2, d = d, h = h2)))
  )))

# More pure cross-validation
log_K01 <- drop(kde_polysph(x = X1[1, , drop = FALSE],
                            X = X1[1, , drop = FALSE],
                            d = d, h = h1, log = TRUE))
log_K02 <- drop(kde_polysph(x = X2[1, , drop = FALSE],
                            X = X2[1, , drop = FALSE],
                            d = d, h = h2, log = TRUE))
H_1_cv2 <- -mean(log1p(exp(log((n1 - 1) / n1) +
                             log_cv_kde_polysph(X = X1, d = d, h = h1) -
                             (log_K01 - log(n1))))) -
  (log_K01 - log(n1))
H_2_cv2 <- -mean(log1p(exp(log((n2 - 1) / n2) +
                             log_cv_kde_polysph(X = X2, d = d, h = h2) -
                             (log_K02 - log(n1))))) -
  (log_K02 - log(n2))
H_0_cv2_1 <- -sum(log1p(exp(log(
  (p[1] * (n1 - 1) / n1) * exp(log_cv_kde_polysph(X = X1, d = d, h = h1)) +
    p[2] * kde_polysph(x = X1, X = X2, d = d, h = h2)) -
    (log_K01 + log(p[1] / n1))))) -
  (log_K01 + log(p[1] / n1)) * n1
H_0_cv2_2 <- -sum(log1p(exp(log(
  (p[2] * (n2 - 1) / n2) * exp(log_cv_kde_polysph(X = X2, d = d, h = h2)) +
    p[1] * kde_polysph(x = X2, X = X1, d = d, h = h1)) -
    (log_K02 + log(p[2] / n2))))) -
  (log_K02 + log(p[2] / n2)) * n2
H_0_cv2 <- (H_0_cv2_1 + H_0_cv2_2) / (n1 + n2)

# Connection of cross-validation with naive Monte Carlo
H_1_cv3 <- -mean(log1p(exp(log((n1 - 1) / n1) +
                             log_cv_kde_polysph(X = X1, d = d, h = h1) -
                             (log_K01 - log(n1)))))
H_2_cv3 <- -mean(log1p(exp(log((n2 - 1) / n2) +
                             log_cv_kde_polysph(X = X2, d = d, h = h2) -
                             (log_K02 - log(n1)))))
H_0_cv3_1 <- -sum(log1p(exp(log(
  (p[1] * (n1 - 1) / n1) * exp(log_cv_kde_polysph(X = X1, d = d, h = h1)) +
    p[2] * kde_polysph(x = X1, X = X2, d = d, h = h2)) -
    (log_K01 + log(p[1] / n1)))))
H_0_cv3_2 <- -sum(log1p(exp(log(
  (p[2] * (n2 - 1) / n2) * exp(log_cv_kde_polysph(X = X2, d = d, h = h2)) +
    p[1] * kde_polysph(x = X2, X = X1, d = d, h = h1)) -
    (log_K02 + log(p[2] / n2)))))
H_0_cv3 <- (H_0_cv3_1 + H_0_cv3_2) / (n1 + n2)

# Checks
H_0 - (p[1] * H_1 + p[2] * H_2)
H_0_cv - (p[1] * H_1_cv + p[2] * H_2_cv)
H_0_cv_bis - (p[1] * H_1_cv + p[2] * H_2_cv)
H_0_n - (p[1] * H_1_n + p[2] * H_2_n)
H_0_cv2 - (p[1] * H_1_cv2 + p[2] * H_2_cv2)
H_0_cv3 - (p[1] * H_1_cv3 + p[2] * H_2_cv3) - sum(log(p) * p)
-sum(log(p) * p)

test_that("Jensen--Shannon distance with Monte Carlo and k = 2", {
  expect_equal(unname(hom_test_poly(X = X, d = d, labels = labels, type = "jsd",
                                    h = h1, B = 1, M = M,
                                    cv_jsd = 123)$statistic),
               H_0 - (p[1] * H_1 + p[2] * H_2),
               tolerance = 5e-2)
})

test_that("Jensen--Shannon distance with cv_jsd = 1 and k = 2", {
  skip("Unstable")
  expect_equal(unname(hom_test_poly(X = X, d = d, labels = labels, type = "jsd",
                                    h = h1, B = 1, M = M,
                                    cv_jsd = 1)$statistic),
               H_0_cv - (p[1] * H_1_cv + p[2] * H_2_cv))
})

test_that("Jensen--Shannon distance with cv_jsd = 2 and k = 2", {
  expect_equal(unname(hom_test_poly(X = X, d = d, labels = labels, type = "jsd",
                                    h = h1, B = 1, M = M,
                                    cv_jsd = 2)$statistic),
               H_0_n - (p[1] * H_1_n + p[2] * H_2_n))
})

## Well-computed JSD for k = 3

r <- 3
d <- rep(2, r)
h1 <- rep(0.25, r)
h2 <- h1
h3 <- h1
n1 <- 100
n2 <- 50
n3 <- 25
n <- n1 + n2 + n3
p <- c(n1, n2, n3) / n
mu1 <- r_unif_polysph(n = 1, d = d)
mu2 <- r_unif_polysph(n = 1, d = d)
mu3 <- r_unif_polysph(n = 1, d = d)
kappa1 <- rep(1, r)
kappa2 <- rep(3, r)
kappa3 <- rep(2, r)
X1 <- r_vmf_polysph(n = n1, mu = mu1, kappa = kappa1, d = d)
X2 <- r_vmf_polysph(n = n2, mu = mu2, kappa = kappa2, d = d)
X3 <- r_vmf_polysph(n = n3, mu = mu3, kappa = kappa3, d = d)
X <- rbind(X1, X2, X3)
labels <- rep(1:3, times = c(n1, n2, n3))
M <- 1e4

# Monte Carlo divergence
samp_1 <- r_kde_polysph(n = M, X = X1, d = d, h = h1)
samp_2 <- r_kde_polysph(n = M, X = X2, d = d, h = h2)
samp_3 <- r_kde_polysph(n = M, X = X3, d = d, h = h3)
samp_0 <- rbind(samp_1[1:round(p[1] * M), ],
                samp_2[1:round(p[2] * M), ],
                samp_3[1:round(p[3] * M), ])
H_1 <- -mean(kde_polysph(x = samp_1, X = X1, d = d, h = h1, log = TRUE))
H_2 <- -mean(kde_polysph(x = samp_2, X = X2, d = d, h = h2, log = TRUE))
H_3 <- -mean(kde_polysph(x = samp_3, X = X3, d = d, h = h3, log = TRUE))
H_0 <- -mean(log(p[1] * kde_polysph(x = samp_0, X = X1, d = d, h = h1) +
                   p[2] * kde_polysph(x = samp_0, X = X2, d = d, h = h2) +
                   p[3] * kde_polysph(x = samp_0, X = X3, d = d, h = h3)))

# Monte Carlo with sample (naive)
H_1_n <- -mean(kde_polysph(x = X1, X = X1, d = d, h = h1, log = TRUE))
H_2_n <- -mean(kde_polysph(x = X2, X = X2, d = d, h = h2, log = TRUE))
H_3_n <- -mean(kde_polysph(x = X3, X = X3, d = d, h = h3, log = TRUE))
H_0_n <- -mean(c(log(p[1] * kde_polysph(x = X, X = X1, d = d, h = h1) +
                       p[2] * kde_polysph(x = X, X = X2, d = d, h = h2) +
                       p[3] * kde_polysph(x = X, X = X3, d = d, h = h3))))

# Cross-validation (with two approaches for f0_cv)
H_1_cv <- -mean(log_cv_kde_polysph(X = X1, d = d, h = h1))
H_2_cv <- -mean(log_cv_kde_polysph(X = X2, d = d, h = h2))
H_3_cv <- -mean(log_cv_kde_polysph(X = X3, d = d, h = h3))
H_0_cv <- -mean(c(log(p[1] * exp(log_cv_kde_polysph(X = X1, d = d, h = h1)) +
                        p[2] * kde_polysph(x = X1, X = X2, d = d, h = h2) +
                        p[3] * kde_polysph(x = X1, X = X3, d = d, h = h3)),
                  log(p[1] * kde_polysph(x = X2, X = X1, d = d, h = h1) +
                        p[2] * exp(log_cv_kde_polysph(X = X2, d = d, h = h2)) +
                        p[3] * kde_polysph(x = X2, X = X3, d = d, h = h3)),
                  log(p[1] * kde_polysph(x = X3, X = X1, d = d, h = h1) +
                        p[2] * kde_polysph(x = X3, X = X2, d = d, h = h2) +
                        p[3] * exp(log_cv_kde_polysph(X = X3, d = d, h = h3)))))

test_that("Jensen--Shannon distance with Monte Carlo and k = 3", {
  skip("Unstable")
  expect_equal(unname(hom_test_poly(X = X, d = d, labels = labels, type = "jsd",
                                    h = h1, B = 1, M = M,
                                    cv_jsd = 123)$statistic),
               H_0 - (p[1] * H_1 + p[2] * H_2 + p[3] * H_3),
               tolerance = 1e-2)
})

test_that("Jensen--Shannon distance with cv_jsd = 1 and k = 3", {
  skip("Unstable")
  expect_equal(unname(hom_test_poly(X = X, d = d, labels = labels, type = "jsd",
                                    h = h1, B = 1, M = M,
                                    cv_jsd = 1)$statistic),
               H_0_cv - (p[1] * H_1_cv + p[2] * H_2_cv + p[3] * H_3_cv))
})

test_that("Jensen--Shannon distance with cv_jsd = 2 and k = 3", {
  expect_equal(unname(hom_test_poly(X = X, d = d, labels = labels, type = "jsd",
                                    h = h1, B = 1, M = M,
                                    cv_jsd = 2)$statistic),
               H_0_n - (p[1] * H_1_n + p[2] * H_2_n + p[3] * H_3_n))
})
