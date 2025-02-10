

#' @rdname euler_ridge
#' @export
parallel_euler_ridge <- function(x, X, d, h, h_euler, N = 1e3, eps = 1e-5,
                                 keep_paths = FALSE, cores = 1, ...) {

  # Starting values
  nx <- nrow(x)

  # Parallel backend
  old_dopar <- doFuture::registerDoFuture()
  old_plan <- future::plan(future::multisession(), workers = cores)
  options(future.rng.onMisuse = "ignore")
  on.exit({

    with(old_dopar, foreach::setDoPar(fun = fun, data = data, info = info))
    future::plan(old_plan)
    options(future.rng.onMisuse = NULL)

  })
  `%op%` <- foreach::`%dopar%`

  # Measure progress?
  if (requireNamespace("progressr", quietly = TRUE)) {

    prog <- progressr::progressor(along = 1:nx)

  }

  # Empty Euler
  p <- ncol(x)
  empty_euler <- list("ridge_y" = rbind(rep(NA, p)), "log_dens_y" = rbind(NA),
                      "paths" = array(NA, dim = c(1, p, N + 1)),
                      "start_x" = NA, "iter" = rbind(NA), "conv" = rbind(NA),
                      "error" = NA)

  # Eulers
  k <- 0
  eu_list <- foreach::foreach(k = 1:nx, .combine = c, .inorder = TRUE,
                              .multicombine = TRUE, .maxcombine = 100,
                              .packages = "polykde",
                              .export = "empty_euler") %op% {

    # Euler algorithm, with error handling
    eu <- try(euler_ridge(x = x[k, , drop = FALSE], X = X, d = d, h = h,
                          h_euler = h_euler, N = N, eps = eps,
                          keep_paths = keep_paths, ...))
    if (inherits(eu, what = "try-error")) {

      empty_euler$start_x <- x[k, , drop = FALSE]
      empty_euler$error <- eu
      eu <- empty_euler

    } else {

      eu$error <- NA

    }

    # Signal progress
    if (requireNamespace("progressr", quietly = TRUE)) {

      prog()

    }

    # Return Euler
    list(eu)

  }

  # Return a single list
  euler <- bind_lists(lists = eu_list, bind = "rbind")
  euler$h <- eu_list[[1]]$h
  euler$d <- eu_list[[1]]$d
  return(euler)

}


#' @rdname euler_ridge
#' @export
block_euler_ridge <- function(x, X, d, h, h_euler, ind_blocks, N = 1e3,
                              eps = 1e-5, keep_paths = FALSE, cores = 1, ...) {

  # Determine begin and end of the S^d's
  r <- length(d)
  ind_dj <- comp_ind_dj(d = d)
  ini_dj <- ind_dj[-(r + 1)] + 1
  end_dj <- ind_dj[-1]

  # Compute
  blocks <- sort(unique(ind_blocks))
  n_blocks <- length(blocks)
  e_block_j <- list(n_blocks)

  # Stop if all d's are not equal
  if (!all(d == d[1])) {

    stop("All d's must be equal for splitting in blocks.")

  }

  # Sanity checks for ind_blocks
  if (length(ind_blocks) != r) {

    stop("Length of ind_blocks is not r.")

  }
  if (!all(1:n_blocks %in% blocks)) {

    stop("The unique elements of ind_blocks are not 1:n_blocks.")

  }
  for (blk in blocks) {

    j <- which(ind_blocks == blk) # Which Sj are in block?
    r_j <- length(j)
    d_j <- d[j]
    ind_j <- unlist(apply(cbind(ini_dj[j], end_dj[j]), 1, function(x) x[1]:x[2],
                          simplify = FALSE))
    X_j <- X[, ind_j]
    if (any(abs(rowSums(X_j^2) / r_j - 1) > 1e-10)) {

      stop(paste0("ind_blocks defines blocks with non-unit norm for ",
                  "block ", blk, "!"))

    }

  }

  # Block ridges
  keep_blocks <- list()
  for (blk in blocks) {

    # Progress
    message(blk, " / ", n_blocks)

    # Block definition
    j <- which(ind_blocks == blk) # Which Sj are in block?
    r_j <- length(j)
    d_j <- d[j]
    ind_j <- unlist(apply(cbind(ini_dj[j], end_dj[j]), 1, function(x) x[1]:x[2],
                          simplify = FALSE))
    keep_blocks[[blk]] <- j

    # Data in the block
    x_j <- x[, ind_j]
    X_j <- X[, ind_j]
    h_j <- h[j]
    h_euler_j <- h_euler[j]

    # Run block fit
    e_block_j[[blk]] <- parallel_euler_ridge(x = x_j, X = X_j, d = d_j,
                                             h = h_j, h_euler = h_euler_j,
                                             N = N, eps = eps,
                                             keep_paths = keep_paths, ...,
                                             cores = cores)

    # Avoid issues on merging due to unequal sphere dimensions
    e_block_j[[blk]]$d <- t(e_block_j[[blk]]$d)
    e_block_j[[blk]]$h <- t(e_block_j[[blk]]$h)

  }

  # Merge components
  e_block <- bind_lists(lists = e_block_j, bind = "cbind")
  e_block$d <- t(e_block$d)
  e_block$h <- t(e_block$h)

  # Use ranks to achieve a fast reordering of the blocks
  # reord_blocks <- unlist(lapply(blocks, function(blk)
  #  which(ind_blocks == blk)))
  # ind_reord_blocks <- unlist(lapply(1:r, function(j)
  #  which(j == reord_blocks)))
  ind_reord_blocks <- rank(ind_blocks, ties.method = "first")

  # Create map of begin and end spheres in the block-ordered data
  d_reord <- d[order(ind_blocks)]
  ind_dj_reord <- comp_ind_dj(d = d_reord)
  ini_dj_reord <- ind_dj_reord[-(r + 1)] + 1
  end_dj_reord <- ind_dj_reord[-1]

  # Index to revert the block-ordering back to the original
  ind_X_blocks <- unlist(apply(cbind(ini_dj_reord, end_dj_reord),
                               1, function(x) x[1]:x[2],
                               simplify = FALSE)[ind_reord_blocks])

  # Revert block-orderings
  e_block$ridge_y <- e_block$ridge_y[, ind_X_blocks]
  e_block$lamb_norm_y <- e_block$lamb_norm_y[, ind_X_blocks]
  if (keep_paths) e_block$paths <- e_block$paths[, ind_X_blocks, ]
  e_block$start_x <- e_block$start_x[, ind_X_blocks]
  e_block$d <- c(e_block$d)[ind_reord_blocks]
  e_block$h <- c(e_block$h)[ind_reord_blocks]

  # Final checks
  if (any(e_block$start_x != x)) {

    warning("Unequal x and $start_x: wrong recollection of blocks!")

  }
  if (any(e_block$h != h)) {

    warning("Unequal h and $h: wrong recollection of blocks!")

  }
  if (any(e_block$d != d)) {

    warning("Unequal d and $d: wrong recollection of blocks!")

  }
  return(e_block)

}


#' @title Clean ridge points coming from spurious fits
#'
#' @description Remove points from the ridge that are spurious. The cleaning is
#' done by removing end points in the Euler algorithm that did not converge,
#' do not have a negative second eigenvalue, or are in low-density regions.
#'
#' @param e outcome from \code{\link{euler_ridge}} or
#' \code{\link{parallel_euler_ridge}}.
#' @param p_out proportion of outliers to remove. Defaults to \code{NULL} (no
#' cleaning).
#' @inheritParams euler_ridge
#' @return A list with the same structure as that returned by
#' \code{\link{euler_ridge}}, but with the spurious points. The removed points
#' are informed in the \code{removed} field.
#' @examples
#' ## Test on S^2 with some spurious end points
#'
#' # Sample
#' r <- 1
#' d <- 2
#' n <- 50
#' ind_dj <- comp_ind_dj(d = d)
#' set.seed(987202226)
#' X <- r_path_s2r(n = n, r = r, spiral = FALSE, Theta = cbind(c(1, 0, 0)),
#'                 sigma = 0.2)[, , 1]
#' col_X <- rep(gray(0), n)
#' col_X_alp <- rep(gray(0, alpha = 0.25), n)
#'
#' # Euler
#' h_rid <- 0.5
#' h_eu <- h_rid^2
#' N <- 30
#' eps <- 1e-6
#' X0 <- r_unif_polysph(n = n, d = d)
#' Y <- euler_ridge(x = X0, X = X, d = d, h = h_rid, h_euler = h_eu,
#'                  N = N, eps = eps, keep_paths = TRUE)
#' Y_removed <- clean_euler_ridge(e = Y, X = X)$removed
#' col_X[Y_removed] <- 2
#' col_X_alp[Y_removed] <- 2
#'
#' # Visualization
#' i <- N # Between 1 and N
#' sc3 <- scatterplot3d::scatterplot3d(Y$paths[, , 1], color = col_X_alp,
#'                                     pch = 19, xlim = c(-1, 1),
#'                                     ylim = c(-1, 1), zlim = c(-1, 1),
#'                                     xlab = "x", ylab = "y", zlab = "z")
#' sc3$points3d(rbind(Y$paths[, , i]), col = col_X, pch = 16, cex = 0.75)
#' invisible(sapply(seq_len(nrow(Y$paths)), function(k) {
#'   sc3$points3d(t(Y$paths[k, , ]), col = col_X_alp[k], type = "l")
#' }))
#' @export
clean_euler_ridge <- function(e, X, p_out = NULL) {

  # Check all end points have a negative second non-null eigenvalue
  lambda_2 <- apply(e$lamb_norm_y, 1, function(lamb) {
    lamb <- lamb[abs(lamb) > 1e-10]
    lamb[2]
  })
  lambda_2_neg <- !is.na(lambda_2) & (lambda_2 < 0)
  message("End points with non-negative lambda_2: ", sum(!lambda_2_neg))

  # Check convergences
  conv <- apply(e$conv, 1, prod) == 1
  message("End points that are non-convergent: ", sum(!conv))

  # Remove low-density (with respect to the original sample) outliers
  if (!is.null(p_out)) {

    log_kde_x <- log_cv_kde_polysph(X = X, d = e$d, h = e$h)
    log_kde_y <- kde_polysph(x = e$ridge_y, X = X, d = e$d, h = e$h, log = TRUE)
    out <- log_kde_y < quantile(log_kde_x, probs = p_out)
    message("End points marked as outliers: ", sum(out))

  } else {

    out <- FALSE

  }

  # Keep only convergent points with negative second eigenvalue
  keep <- !out & lambda_2_neg & conv
  message("Removed end points: ", sum(!keep))
  for (j in c(1:7)[-4]) e[[j]] <- e[[j]][keep, , drop = FALSE]
  e$paths <- e$paths[keep, , , drop = FALSE]

  # Add information on removed points
  e$removed <- !keep
  return(e)

}


#' @title Index a ridge curve, creating the smoothed and indexed ridge
#'
#' @description Indexing of an unordered collection of points defining the
#' estimated density ridge curve. The indexing is done by a multidimensional
#' scaling map to the real line, while the smoothing is done by local polynomial
#' regression for polyspherical-on-scalar regression.
#'
#' @param endpoints end points of the ridge algorithm to be indexed.
#' @inheritParams euler_ridge
#' @param l_index length of the grid index used for finding projections.
#' Defaults to \code{1e3}.
#' @param f_index factor with the range of the grid for finding ridge
#' projections. Defaults to \code{2}, which is twice the range of MDS indexes.
#' @param probs_scores probabilities for indexing the ridge on the
#' \code{probs_scores}-quantiles of the scores. Defaults to
#' \code{seq(0, 1, l = 101)}.
#' @param verbose show diagnostic plots?
#' @param type_bwd type of cross-validation bandwidth for Nadaraya--Watson,
#' either \code{"min"} (minimizer of the cross-validation loss) or \code{"1se"}
#' (the "one standard error rule", smoother). Defaults to \code{"1se"}.
#' @inheritParams kre_polysph
#' @param ... further parameters passed to \code{\link{bw_cv_kre_polysph}}.
#' @return A list with the following fields:
#' \item{scores_X}{TODO.}
#' \item{projs_X}{TODO.}
#' \item{ord_X}{TODO.}
#' \item{scores_grid}{TODO.}
#' \item{ridge_grid}{TODO.}
#' \item{mds_index}{TODO.}
#' \item{ridge_fun}{TODO.}
#' \item{h}{TODO.}
#' \item{probs_scores}{TODO.}
#' @examples
#' ## Test on S^2
#'
#' # Sample
#' r <- 1
#' d <- 2
#' n <- 50
#' ind_dj <- comp_ind_dj(d = d)
#' set.seed(987204452)
#' X <- r_path_s2r(n = n, r = r, spiral = FALSE, Theta = cbind(c(1, 0, 0)),
#'                 sigma = 0.2)[, , 1]
#' col_X <- rep(gray(0), n)
#' col_X_alp <- rep(gray(0, alpha = 0.25), n)
#'
#' # Euler
#' h_rid <- 0.5
#' h_eu <- h_rid^2
#' N <- 30
#' eps <- 1e-6
#' X0 <- r_unif_polysph(n = n, d = d)
#' Y <- euler_ridge(x = X0, X = X, d = d, h = h_rid, h_euler = h_eu,
#'                  N = N, eps = eps, keep_paths = TRUE)
#' ind_rid <- index_ridge(endpoints = Y$ridge_y, X = X, d = d)
#' # TODO
#'
#' # Visualization
#' i <- N # Between 1 and N
#' sc3 <- scatterplot3d::scatterplot3d(Y$paths[, , 1], color = col_X_alp,
#'                                     pch = 19, xlim = c(-1, 1),
#'                                     ylim = c(-1, 1), zlim = c(-1, 1),
#'                                     xlab = "x", ylab = "y", zlab = "z")
#' sc3$points3d(rbind(Y$paths[, , i]), col = col_X, pch = 16, cex = 0.75)
#' invisible(sapply(seq_len(nrow(Y$paths)), function(k) {
#'   sc3$points3d(t(Y$paths[k, , ]), col = col_X_alp[k], type = "l")
#' }))
#' sc3$points3d(ind_rid$ridge_fun(seq(0, 1, l = 10)), col = "red", pch = 16, cex = 1.5, type = "l")
#' @export
index_ridge <- function(endpoints, X, d, l_index = 1e3, f_index = 2,
                        probs_scores = seq(0, 1, l = 101), verbose = FALSE,
                        type_bwd = c("1se", "min")[1], p = 0, ...) {

  # Distance matrix between end points
  ind_dj <- comp_ind_dj(d = d)
  dij <- dist_polysph_matrix(x = endpoints, ind_dj = ind_dj, norm_x = TRUE,
                             norm_y = TRUE, std = FALSE)

  # Indexing via MDS
  mds_index <- smacof::mds(delta = dij, ndim = 1, itmax = 1e3, eps = 1e-6)
  mds_index <- drop(mds_index$conf)
  if (verbose) {

    # MDS-indexes
    plot(density(mds_index), main = "MDS indexes")
    rug(mds_index)

  }

  # Smoothing of the ridge using the learned index
  h_nw <- bw_cv_kre_polysph(X = mds_index, Y = endpoints, d = d,
                            p = p, plot_cv = verbose, ...)
  h_nw <- switch(type_bwd, "min" = h_nw$h_min, "1se" = h_nw$h_1se,
                 stop("type_bwd must be \"min\" or \"1se\""))
  index_grid <- seq(f_index * min(mds_index), f_index * max(mds_index),
                    l = l_index)
  ridge_grid <- kre_polysph(x = index_grid, X = mds_index, Y = endpoints,
                            d = d, h = h_nw, p = p)

  # Scores and projections onto the smoothed ridge
  index_scores_ridge_grid <- sapply(seq_len(nrow(X)), function(i) {
    which.min(dist_polysph(x = ridge_grid, y = X[i, , drop = FALSE],
                           ind_dj = ind_dj, norm_x = TRUE,
                           norm_y = TRUE, std = FALSE))
  })
  scores_X <- index_grid[index_scores_ridge_grid]
  projs_X <- ridge_grid[index_scores_ridge_grid, ]

  # Recenter scores
  mds_index <- mds_index - mean(scores_X)
  scores_X <- scores_X - mean(scores_X)
  if (verbose) {

    # MDS-indexes
    plot(density(scores_X), main = "Scores")
    rug(scores_X)

  }

  # Ridge path indexed by the quantiles of the sample scores
  scores_grid <- quantile(scores_X, probs = probs_scores)
  ridge_scores_grid <- kre_polysph(x = scores_grid, X = mds_index,
                                   Y = endpoints, d = d, h = h_nw, p = p)

  # Ridge handle
  r_hat <- function(t) {

    kre_polysph(x = t, X = mds_index, Y = endpoints, d = d, h = h_nw, p = p)

  }

  # List
  return(list("scores_X" = scores_X, "projs_X" = projs_X,
              "ord_X" = order(index_scores_ridge_grid),
              "scores_grid" = scores_grid, "ridge_grid" = ridge_scores_grid,
              "mds_index" = mds_index, "ridge_fun" = r_hat, "h" = h_nw,
              "probs_scores" = probs_scores))

}
