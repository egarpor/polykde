
# Randomize testing
r <- 5
d <- rep(2, r)
ind_dj <- comp_ind_dj(d = d)
n <- 100
X <- rnorm(n)
Y <- t(sapply(1:n, function(i) {
  mu <- c(cos(2 * pi * X[i]), sin(2 * pi * X[i]), 0.1)
  r_kern_polysph(n = 1, d = d, mu = rbind(rep(mu, r)), h = rep(1, r),
                 norm_mu = TRUE)
}))
Y_avg <- proj_polysph(x = rbind(colMeans(Y)), ind_dj = ind_dj)

# # Timing
# microbenchmark::microbenchmark(
#   bw_cv_kre_polysph(X = X, Y = Y, d = d, h = h, p = 0, fast = FALSE, plot_cv = FALSE),
#   bw_cv_kre_polysph(X = X, Y = Y, d = d, h = h, p = 0, fast = TRUE, plot_cv = FALSE),
#   times = 1e3)
# r = 10; d = 2; n = 100
# Unit: milliseconds
# expr       min        lq
# bw_cv_kre_polysph(X = X, Y = Y, d = d, h = h, p = 0, fast = FALSE, plot_cv = FALSE) 313.29404 322.22574
# bw_cv_kre_polysph(X = X, Y = Y, d = d, h = h, p = 0, fast = TRUE, plot_cv = FALSE)  20.66351  22.71197
# mean    median        uq      max neval cld
# 332.12150 325.80480 331.54650 610.0708  1000   b
# 31.03197  23.25623  30.83821 144.7515  1000  a
# r = 100; d = 2; n = 200
# Unit: milliseconds
# expr       min        lq
# bw_cv_kre_polysph(X = X, Y = Y, d = d, h = h, p = 0, fast = FALSE, plot_cv = FALSE) 961.26206 988.12673
# bw_cv_kre_polysph(X = X, Y = Y, d = d, h = h, p = 0, fast = TRUE, plot_cv = FALSE)  74.81003  81.48791
# mean    median        uq       max neval cld
# 1009.0932 995.36549 1006.7707 1584.8514  1000   b
# 114.0201  85.91021  170.0534  361.3551  1000  a

test_that("Equivalent bandwidth selection with fast and slow CV in NW", {
  expect_equal(bw_cv_kre_polysph(X = X, Y = Y, d = d, p = 0, fast = FALSE)$cvm,
               bw_cv_kre_polysph(X = X, Y = Y, d = d, p = 0, fast = TRUE)$cvm)
})

# test_that("Equivalent bandwidth selection with fast and slow CV in LL", {
#   expect_equal(bw_cv_kre_polysph(X = X, Y = Y, d = d, p = 1, fast = FALSE)$cvm,
#                bw_cv_kre_polysph(X = X, Y = Y, d = d, p = 1, fast = TRUE)$cvm)
# })

test_that("NW interpolates for very small bandwidth", {
  expect_lt(suppressWarnings(
    max(dist_polysph(x = kre_polysph(x = X, X = X, Y = Y,
                                     d = d, h = 1e-5, p = 0),
                     y = Y, ind_dj = ind_dj, std = FALSE))),
    1e-2)
})

test_that("NW becomes the unconditional mean for very large bandwidth", {
  expect_lt(max(dist_polysph(x = kre_polysph(x = X, X = X, Y = Y,
                                             d = d, h = 100, p = 0),
                             y = Y_avg, ind_dj = ind_dj, std = FALSE)),
            1e-2)
})

test_that("LL interpolates for very small bandwidth", {
  expect_lt(suppressWarnings(
    max(dist_polysph(x = kre_polysph(x = X, X = X, Y = Y,
                                     d = d, h = 1e-5, p = 1),
                     y = Y, ind_dj = ind_dj, std = FALSE))),
    1e-2)
})

test_that("LL becomes the unconditional mean for very large bandwidth", {
  skip("Unstable")
  expect_lt(max(dist_polysph(x = kre_polysph(x = X, X = X, Y = Y,
                                             d = d, h = 100, p = 1),
                             y = Y_avg, ind_dj = ind_dj, std = FALSE)),
            1e-2)
})

# Visualize interpolation curve -- LL extrapolates, NW stops
x_grid <- seq(0, 1, l = 200)
x_grid2 <- seq(-1, 2, l = 200)
n <- 100
X <- seq(0, 1, l = n)
Y <- r_path_s2r(n = n, r = 1, Theta = cbind(c(1, 0, 0)), sigma = 0.05,
                spiral = TRUE)[, , 1]

sc3 <- scatterplot3d::scatterplot3d(Y, pch = 16, xlim = c(-1, 1),
                                    ylim = c(-1, 1), zlim = c(-1, 1),
                                    xlab = "", ylab = "", zlab = "",
                                    box = FALSE, axis = FALSE, grid = FALSE)
h0 <- 0.01
sc3$points3d(kre_polysph(x = x_grid, X = X, Y = Y, d = 2, h = h0, p = 0),
             pch = 16, type = "l", col = 2, lwd = 2)
sc3$points3d(kre_polysph(x = x_grid2, X = X, Y = Y, d = 2, h = h0, p = 0),
             pch = 16, type = "l", col = 2)
h1 <- 0.1
sc3$points3d(kre_polysph(x = x_grid, X = X, Y = Y, d = 2, h = h1, p = 1),
             pch = 16, type = "l", col = 3, lwd = 2)
sc3$points3d(kre_polysph(x = x_grid2, X = X, Y = Y, d = 2, h = h1, p = 1),
             pch = 16, type = "l", col = 3)
L <- 100; tth <- seq(0, 2 * pi, l = L); pph <- seq(0, pi, l = L)
for (ph in seq(0, pi, l = 21))
  sc3$points3d(DirStats::to_sph(th = tth, ph = rep(ph, L)), type = "l",
               col = gray(0.5, alpha = 0.25))
for (th in seq(-pi, pi, l = 21))
  sc3$points3d(DirStats::to_sph(th = rep(th, L), ph = pph), type = "l",
               col = gray(0.5, alpha = 0.25))
