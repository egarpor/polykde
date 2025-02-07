

#' @title Parallel Euler algorithm
#'
#' @inheritParams euler_ridge
#' @param ... further parameters passed to \code{\link{euler_ridge}}.
#' @param cores cores to use. Defaults to \code{1}.
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


#' @title Blockwise Euler algorithm
#'
#' @param ind_blocks indexes of the blocks, a vector or length \code{r}.
#' @inheritParams euler_ridge
#' @inheritParams parallel_euler_ridge
#' @param ... further parameters passed to \code{\link{euler_ridge}}.
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

  # Revert bloc-orderings
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
#' @param e outcome from \code{\link{euler_ridge}} or
#' \code{\link{parallel_euler_ridge}}.
#' @param p_out proportion of outliers to remove. Defaults to \code{NULL} (no
#' cleaning).
#' @inheritParams euler_ridge
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
#' @param endpoints TODO.
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
