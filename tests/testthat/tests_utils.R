
## proj_polysph()

r <- 3
d <- c(2, 3, 1)
n <- 2
X <- r_unif_polysph(n = n, d = d)
ind_dj <- comp_ind_dj(d)

test_that("Normalization in dist_polysph()", {
  expect_equal(proj_polysph(x = X, ind_dj = ind_dj), X)
  expect_equal(proj_polysph(x = 1:n * X, ind_dj = ind_dj), X)
})

## dist_polysph()

r <- 2
d <- c(2, 3)
n <- 10
X1 <- r_unif_polysph(n = n, d = d)
X2 <- r_unif_polysph(n = n, d = d)
ind_dj <- comp_ind_dj(d)

test_that("Normalization in dist_polysph()", {
  expect_equal(
    dist_polysph(x = X1, y = X2, ind_dj = ind_dj,
                 norm_x = FALSE, norm_y = FALSE),
    dist_polysph(x = X1, y = X2, ind_dj = ind_dj,
                 norm_x = TRUE, norm_y = TRUE))
})

test_that("Simple case in dist_polysph()", {
  expect_equal(
    drop(dist_polysph(x = X1, y = X2, ind_dj = ind_dj, std = FALSE)),
    sqrt(
      acos(rowSums(X1[, 1:3] * X2[, 1:3]))^2 +
        acos(rowSums(X1[, 4:7] * X2[, 4:7]))^2
    ))
})

test_that("Distance between the same points is zero in dist_polysph()", {
  expect_equal(drop(dist_polysph(x = X1, y = X1, ind_dj = ind_dj,
                                 norm_x = TRUE, norm_y = TRUE)),
               rep(0, n))
  expect_equal(drop(dist_polysph(x = X2, y = X2, ind_dj = ind_dj,
                                 norm_x = TRUE, norm_y = TRUE)),
               rep(0, n))
  expect_equal(drop(dist_polysph(x = X1, y = X1, ind_dj = ind_dj,
                                 norm_x = TRUE, norm_y = TRUE)),
               drop(dist_polysph(x = X1, y = X1, ind_dj = ind_dj,
                                 norm_x = FALSE, norm_y = FALSE)))
  expect_equal(drop(dist_polysph(x = X2, y = X2, ind_dj = ind_dj,
                                 norm_x = TRUE, norm_y = TRUE)),
               drop(dist_polysph(x = X2, y = X2, ind_dj = ind_dj,
                                 norm_x = FALSE, norm_y = FALSE)))
})

## dist_polysph_cross()

test_that("dist_polysph_cross() equals dist_polysph_matrix() for x = y", {
  expect_equal(
    c(dist_polysph_matrix(x = X1, ind_dj = ind_dj, std = FALSE)),
    c(dist_polysph_cross(x = X1, y = X1, ind_dj = ind_dj, std = FALSE)),
    tolerance = 1e-6)
})

test_that("dist_polysph_cross() equals dist_polysph()", {
  expect_equal(
    cbind(dist_polysph(x = X1[1:3, ], y = X2[1, , drop = FALSE],
                       ind_dj = ind_dj),
          dist_polysph(x = X1[1:3, ], y = X2[2, , drop = FALSE],
                       ind_dj = ind_dj)),
    dist_polysph_cross(x = X1[1:3, ], y = X2[1:2, ], ind_dj = ind_dj),
    tolerance = 1e-6)
  expect_equal(
    dist_polysph(x = X1[1, , drop = FALSE], y = X2[1, , drop = FALSE],
                 ind_dj = ind_dj),
    dist_polysph_cross(x = X1[1, , drop = FALSE], y = X2[1, , drop = FALSE],
                       ind_dj = ind_dj),
    tolerance = 1e-6)
})

## diamond_crossprod()

# Randomize testing
r <- 2
d <- c(2, 3)
n <- 3
X <- r_unif_polysph(n = n, d = d)
ind_dj <- comp_ind_dj(d)

test_that("Simple case in diamond_crossprod()", {
  for (i in seq_len(n)) {
    expect_equal(
      diamond_crossprod(X = X, ind_dj = ind_dj)[i, , ],
      rbind(cbind(tcrossprod(X[i, 1:3]), tcrossprod(X[i, 1:3], X[i, 4:7])),
            cbind(tcrossprod(X[i, 4:7], X[i, 1:3]), tcrossprod(X[i, 4:7])))
      )
  }
})

## s()

test_that("Simmetrization using s()", {
  expect_true(isSymmetric(s(r_unif_polysph(n = 4, d = 3))))
})

## proj_P()

# Randomize testing
r <- rpois(1, lambda = 3) + 1
d <- rpois(r, lambda = 2) + 1
x <- r_unif_polysph(n = 1, d = d)
ind_dj <- comp_ind_dj(d)

proj_P <- function(x, d) {

  # Projections
  ind <- cumsum(c(1, d + 1))
  P <- matrix(0, nrow = ind[length(d) + 1] - 1, ncol = ind[length(d) + 1] - 1)
  for (j in seq_along(d)) {

    ind_j <- ind[j]:(ind[j + 1] - 1)
    P[ind_j, ind_j] <- -tcrossprod(x[ind_j])

  }
  diag(P) <- diag(P) + 1
  return(P)

}

test_that("AP does the job", {
  expect_equal(proj_P(x = x, d = d), AP(x = x, v = x, ind_dj = ind_dj)$P)
  expect_equal(diag(1, nrow = ncol(x), ncol = ncol(x)) - proj_P(x = x, d = d),
               AP(x = x, v = x, ind_dj = ind_dj, orth = TRUE)$P)
})

## polylog_minus_exp_mu()

k <- c(0.1, 0.5, 1, 10, 100)
test_that("polylog_minus_exp_mu() is smooth on integer indexes", {
  for (s in 1:5) {
    expect_equal(polylog_minus_exp_mu(s = s, mu = k),
                 polylog_minus_exp_mu(s = s + 1e-6, mu = k),
                 tolerance = 1e-5)
  }
})
test_that("polylog_minus_exp_mu() is smooth on half arguments", {
  for (s in c(0.5, 1.5, 2.5)) {
    expect_equal(polylog_minus_exp_mu(s = s, mu = k),
                 polylog_minus_exp_mu(s = s + 1e-6, mu = k),
                 tolerance = 1e-5)
  }
})
