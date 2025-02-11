
test_that("angles_to_torus() and torus_to_angles() are inverses", {
  d <- 3
  x <- r_unif_polysph(n = 10, d = rep(1, d))
  theta <- torus_to_angles(x)
  expect_equal(x, angles_to_torus(torus_to_angles(x)))
  expect_equal(theta, torus_to_angles(angles_to_torus(theta)))
  expect_equal(x[1, , drop = FALSE],
               angles_to_torus(torus_to_angles(x[1, ])))
  expect_equal(theta[1, , drop = FALSE],
               torus_to_angles(angles_to_torus(theta[1, ])))
})

test_that("Error if p is odd in angles_to_torus()", {
  expect_error(torus_to_angles(matrix(0, nrow = 1, ncol = 3)))
})

test_that("angles_to_sph() and sph_to_angles() are inverses, d = 1", {
  d <- 1
  x <- r_unif_polysph(n = 10, d = d)
  theta <- sph_to_angles(x)
  expect_equal(x, angles_to_sph(sph_to_angles(x)))
  expect_equal(theta, sph_to_angles(angles_to_sph(theta)))
  expect_equal(x[1, , drop = FALSE],
               angles_to_sph(sph_to_angles(x[1, ])))
  expect_equal(theta[1, , drop = FALSE],
               sph_to_angles(angles_to_sph(theta[1, ])))
})

test_that("angles_to_sph() and sph_to_angles() are inverses, d = 3", {
  d <- 3
  x <- r_unif_polysph(n = 10, d = d)
  theta <- sph_to_angles(x)
  expect_equal(x, angles_to_sph(sph_to_angles(x)))
  expect_equal(theta, sph_to_angles(angles_to_sph(theta)))
  expect_equal(x[1, , drop = FALSE],
               angles_to_sph(sph_to_angles(x[1, ])))
  expect_equal(theta[1, , drop = FALSE],
               sph_to_angles(angles_to_sph(theta[1, ])))
})

test_that("polysph_to_angles() and angles_to_polysph() are inverses", {
  d <- 1:5
  x <- r_unif_polysph(n = 10, d = d)
  theta <- polysph_to_angles(x, d = d)
  expect_equal(x, angles_to_polysph(polysph_to_angles(x, d = d), d = d))
  expect_equal(theta, polysph_to_angles(angles_to_polysph(theta, d = d), d = d))
  expect_equal(x[1, , drop = FALSE],
               angles_to_polysph(polysph_to_angles(x[1, ], d = d), d = d))
  expect_equal(theta[1, , drop = FALSE],
               polysph_to_angles(angles_to_polysph(theta[1, ], d = d), d = d))
})
