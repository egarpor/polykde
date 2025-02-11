
test_that("interp_srep() works forward and backward", {

  expect_equal(
    interp_polysph(x = c(1, 0), y = c(0, 1),
                   ind_dj = comp_ind_dj(d = 1), N = 10),
    interp_polysph(x = c(0, 1), y = c(1, 0),
                   ind_dj = comp_ind_dj(d = 1), N = 10)[10:1, ])

})

test_that("view_srep() does not return error with static = TRUE", {

  base <- r_unif_polysph(n = 50, d = 2)
  dirs <- base
  radii <- runif(nrow(base), min = 0.5, max = 1)
  bdry <- base + radii * dirs
  expect_no_error(view_srep(base = base, dirs = dirs, bdry = bdry,
                            radii = radii, static = TRUE, texts = 1:50))
  expect_no_error(view_srep(base = base, dirs = dirs, radii = radii,
                            static = TRUE))
  expect_no_error(view_srep(base = base, dirs = dirs, bdry = bdry,
                            static = TRUE))
  expect_no_error(view_srep(base = base, dirs = dirs, bdry = bdry,
                            static = TRUE, show_base = FALSE,
                            show_seg = FALSE, show_bdry_pt = FALSE,
                            show_bdry = FALSE))

})
