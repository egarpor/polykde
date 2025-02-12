
#' @title Interpolation on the polysphere
#'
#' @description Creates a sequence of points on the polysphere linearly
#' interpolating between two points extrinsically.
#'
#' @param x a vector of size \code{sum(d) + r} with the begin point.
#' @param y a vector of size \code{sum(d) + r} with the end point.
#' @inheritParams proj_polysph
#' @param N number of points in the sequence. Defaults to \code{10}.
#' @return A matrix of size \code{c(N, sum(d) + r)} with the interpolation.
#' @examples
#' interp_polysph(x = c(1, 0), y = c(0, 1), ind_dj = comp_ind_dj(d = 1))
#' @export
interp_polysph <- function(x, y, ind_dj, N = 10) {

  # Checks
  stopifnot(N > 1)
  stopifnot(is.vector(x) && is.vector(y))

  # Interpolation
  t <- seq(0, 1, l = N)
  xy <- t %o% x + (1 - t) %o% y

  # Projection
  return(proj_polysph(x = xy, ind_dj = ind_dj))

}


#' @title s-rep viewer
#'
#' @description Plots a skeletal representation (s-rep) object based on its
#' three-dimensional base, spokes, and boundary.
#'
#' @param base base points, a matrix of size \code{c(nx, 3)}.
#' @param dirs directions of spokes, a matrix of size \code{c(nx, 3)} with
#' unit vectors.
#' @param bdry boundary points, a matrix of size \code{c(nx, 3)}.
#' @param radii radii of spokes, a vector of size \code{nx}.
#' @param show_base,show_base_pt show base and base grid? Default to
#' \code{TRUE}.
#' @param show_bdry,show_bdry_pt show boundary and boundary grid? Default to
#' \code{TRUE}.
#' @param show_seg show segments? Defaults to \code{TRUE}.
#' @param col_base,col_bndy,col_seg colors for the base, boundary, and segments.
#' Default to \code{"red"}, \code{"blue"}, and \code{"green"}, respectively.
#' @param static use static (\code{\link[scatterplot3d]{scatterplot3d}}) or
#' interactive (\code{\link[rgl]{plot3d}}) plot? Default to \code{TRUE}.
#' @param texts add text labels? If given, it should be a vector of size
#' \code{nx} with the labels. Defaults to \code{NULL}.
#' @param cex_base,cex_bdry size of the base and boundary points.
#' @param cex_texts size of the text labels. Defaults to \code{1}.
#' @param lwd_seg width of the segments.
#' @param alpha_base,alpha_bdry transparencies for base and boundary. Default to
#' \code{0.1} and \code{0.15}, respectively.
#' @param r_texts magnification of the radius to separate the text labels.
#' Defaults to \code{1.25}.
#' @param alpha_ashape3d_base,alpha_ashape3d_bdry alpha parameters for
#' \code{\link[alphashape3d]{ashape3d}}. Default to \code{NULL}.
#' @param lit lit parameter passed to \code{\link[rgl]{material3d}}.
#' @param ... further arguments to be passed to \code{\link[rgl]{plot3d}} or
#' \code{\link[scatterplot3d]{scatterplot3d}}.
#' @return Creates a static or interactive plot.
#' @examples
#' base <- r_unif_polysph(n = 50, d = 2)
#' dirs <- base
#' radii <- runif(nrow(base), min = 0.5, max = 1)
#' bdry <- base + radii * dirs
#' view_srep(base = base, dirs = dirs, bdry = bdry, radii = radii,
#'           texts = 1:50, xlim = c(-2, 2), ylim = c(-2, 2), zlim = c(-2, 2))
#' @export
view_srep <- function(base, dirs, bdry, radii, show_base = TRUE,
                      show_base_pt = TRUE, show_bdry = TRUE,
                      show_bdry_pt = TRUE, show_seg = TRUE, col_base = "red",
                      col_bndy = "blue", col_seg = "green", static = TRUE,
                      texts = NULL, cex_base = ifelse(static, 0.5, 6),
                      cex_bdry = ifelse(static, 1, 8),
                      lwd_seg = ifelse(static, 1, 2), cex_texts = 1,
                      alpha_base = 0.1, alpha_bdry = 0.15, r_texts = 1.25,
                      alpha_ashape3d_base = NULL, alpha_ashape3d_bdry = NULL,
                      lit = FALSE, ...) {

  # Checks
  stopifnot(is.matrix(base) & is.matrix(dirs))
  stopifnot(all(dim(base) == dim(dirs)))

  # bdry from base + dirs
  if (missing(bdry)) {

    stopifnot(!missing(base) & !missing(radii) & !missing(dirs))
    stopifnot(is.vector(radii) & (length(radii) == nrow(base)))
    bdry <- base + radii * dirs

  }

  # radii from base + bdry
  if (missing(radii)) {

    stopifnot(!missing(base) & !missing(bdry))
    radii <- sqrt(rowSums((bdry - base)^2))

  }

  # Show base points
  if (static) {

    sc3 <- scatterplot3d::scatterplot3d(base, color = col_base,
                                        tick.marks = FALSE, xlab = "",
                                        ylab = "", zlab = "", pch = 16,
                                        type = ifelse(show_base_pt, "p", "n"),
                                        cex.symbols = cex_base, ...)

  } else {

    rgl::plot3d(base, col = col_base, axes = FALSE, box = FALSE,
                xlab = "", ylab = "", zlab = "", size = cex_base,
                type = ifelse(show_base_pt, "p", "n"),
                aspect = apply(bdry, 2, function(x) diff(range(x))),
                lit = lit, ...)

    # Show base surface
    if (show_base) {

      # Remove almost equal points --- otherwise ashape3d() will crash R
      D <- as.matrix(dist(base))
      D[lower.tri(D, diag = TRUE)] <- NA
      duplicated_base <- which(D < 1e-6, arr.ind = TRUE)[, 1]

      # Plot surface
      if (is.null(alpha_ashape3d_base)) {

        stop("alpha_ashape3d_base must be provided")

      }
      ash_base <- alphashape3d::ashape3d(base[-duplicated_base, ],
                                         alpha = alpha_ashape3d_base)
      m_base <- rgl::as.mesh3d(ash_base, triangles = FALSE)
      rgl::shade3d(m_base, col = col_base, alpha = alpha_base, lit = lit,
                   meshColor = "faces", smooth = FALSE)
      rgl::addNormals(x = m_base)

    }

  }

  # Show segments
  if (show_seg) {

    if (length(col_seg) == 1) {

      col_seg <- rep(col_seg, nrow(base))

    }
    for (i in seq_len(nrow(base))) {

      if (static) {

        sc3$points3d(rbind(base[i, ], bdry[i, ]), col = col_seg[i],
                     lwd = lwd_seg, type = "l")

      } else {

        rgl::segments3d(rbind(base[i, ], bdry[i, ]), col = col_seg[i],
                        lwd = lwd_seg, lit = lit)

      }

    }

  }

  # Show boundary surface
  if (show_bdry_pt) {

    if (static) {

      sc3$points3d(bdry, col = col_bndy, pch = 16, cex = cex_bdry)

    } else {

      rgl::points3d(bdry, col = col_bndy, size = cex_bdry, lit = lit)

    }

  }
  if (show_bdry && !static) {

    # Remove almost equal points --- otherwise ashape3d() will crash R
    D <- as.matrix(dist(bdry))
    D[lower.tri(D, diag = TRUE)] <- NA
    duplicated_bdry <- which(D < 1e-6, arr.ind = TRUE)[, 1]

    # Plot surface
    if (is.null(alpha_ashape3d_base)) {

      stop("alpha_ashape3d_base must be provided")

    }
    ash_bndy <- alphashape3d::ashape3d(bdry[-duplicated_bdry, ],
                                       alpha = alpha_ashape3d_bdry)
    m_bndy <- rgl::as.mesh3d(ash_bndy)
    rgl::shade3d(m_bndy, alpha = alpha_bdry, col = col_bndy, lit = lit)
    rgl::addNormals(x = m_bndy)

  }

  # Add texts
  if (!is.null(texts)) {

    stopifnot(!missing(radii) & !missing(dirs))
    if (static) {

      text(sc3$xyz.convert(base + r_texts * radii * dirs), labels = texts,
           col = col_bndy, pos = 2, cex = cex_texts)

    } else {

      rgl::text3d(base + r_texts * radii * dirs, texts = texts,
                  col = col_bndy, pos = 0, cex = cex_texts, family = "sans")

    }

  }

}


#' #' @title Movie weighted by distances
#' #'
#' #' @description TODO
#' #'
#' #' @param X TODO
#' #' @inheritParams proj_polysph
#' #' @param ord TODO
#' #' @param L number of frames. If \code{NULL}, \code{nx} frames are taken.
#' #' @return TODO
#' #' @examples
#' #' X <- polykde:::interp_polysph(x = c(1, 0), y = c(0, 1),
#' #'                               ind_dj = comp_ind_dj(d = 1), N = 5)
#' #' polykde:::interp_srep(X = X, ind_dj = comp_ind_dj(d = 1), ord = 1:5, L = 10)
#' #' polykde:::interp_srep(X = X, ind_dj = comp_ind_dj(d = 1), ord = 5:1, L = 10)
#' #' @keywords internal
#' interp_srep <- function(X, ind_dj, ord, L = NULL) {
#'
#'   # Ordered
#'   X_ord <- X[ord, ]
#'
#'   # How many frames per index
#'   if (is.null(L)) {
#'
#'     L <- nrow(X_ord) - 1
#'     frames <- rep(1, L)
#'     ind_frames <- seq_len(L)
#'
#'   } else {
#'
#'     dists <- dist_polysph(x = X_ord[-1, ], y = X_ord[-nrow(X_ord), ],
#'                           ind_dj = ind_dj)
#'     frames <- round(c(L * dists / sum(dists)))
#'     ind_frames <- which(frames > 0)
#'     frames <- frames[ind_frames]
#'
#'   }
#'
#'   # Interpolation
#'   int <- X_ord[ind_frames[1], ]
#'   for (i in seq_along(ind_frames)) {
#'
#'     int_i <- interp_polysph(x = X_ord[ind_frames[i], ],
#'                             y = X_ord[ind_frames[i] + 1, ],
#'                             ind_dj = ind_dj, N = frames[i] + 1)[-1, ]
#'     int <- rbind(int, int_i)
#'
#'   }
#'   return(int)
#'
#' }
#'
#'
#' #' @title View s-reps across time
#' #'
#' #' @description TODO
#' #'
#' #' @param avg_base average base points, a matrix of size \code{c(nx, 3)}.
#' #' @param avg_dirs average directions of spokes, a matrix of size
#' #' \code{c(nx, 3)} with unit vectors.
#' #' @param avg_radii average radii of spokes, a vector of size \code{nx}.
#' #' @param moving_spokes TODO
#' #' @inheritParams interp_srep
#' #' @param col_pal color palette for
#' #' @inheritParams view_srep
#' #' @param ... further arguments passed to \code{\link{view_srep}}.
#' #' @return Creates a static or interactive plot.
#' #' @examples
#' #' base <- r_unif_polysph(n = 50, d = 2)
#' #' dirs <- base
#' #' radii <- runif(nrow(base), min = 0.5, max = 1)
#' #' bdry <- base + radii * dirs
#' #' polykde:::view_srep_movie(avg_base = base, avg_radii = radii,
#' #'                           avg_dirs = dirs, moving_spokes = 1:25,
#' #'                           X = dirs)
#' #' @keywords internal
#' view_srep_movie <- function(avg_base, avg_radii, avg_dirs, moving_spokes,
#'                             X, ord, col_pal = grDevices::rainbow, L = NULL,
#'                             static = TRUE, ...) {
#'
#'   # Obtain ind_dj_nice
#'   r_nice <- sum(moving_spokes)
#'   d_nice <- rep(2, r_nice)
#'   ind_dj_nice <- rep(0, r_nice + 1)
#'   ind_dj_nice[-1] <- d_nice + 1
#'   ind_dj_nice <- cumsum(ind_dj_nice)
#'
#'   # Interpolation of directions
#'   rid_dirs <- interp_srep(X = X, ind_dj = ind_dj_nice,
#'                           ord = ord, L = L)
#'   L <- nrow(rid_dirs)
#'   cols <- col_pal(L)
#'   rid_dirs_k <- avg_dirs
#'
#'   # Plot s-reps
#'   if (static) {
#'
#'     angle <- 40
#'     manipulate::manipulate({
#'
#'       # Update directions
#'       rid_dirs_k[moving_spokes, ] <- matrix(rid_dirs[k, ], ncol = 3,
#'                                             byrow = TRUE)
#'
#'       # Plot s-rep
#'       view_srep(base = avg_base, dirs = rid_dirs_k, radii = avg_radii,
#'                 show_base = FALSE, show_base_pt = TRUE, show_bdry = TRUE,
#'                 show_bdry_pt = TRUE, show_seg = TRUE,
#'                 col_bndy = ifelse(moving_spokes, cols[k], "gray"),
#'                 col_base = 1, col_seg = ifelse(moving_spokes, cols[k], "gray"),
#'                 static = static, angle = angle, ...)
#'       # theta = theta, phi = phi)
#'
#'     }, k = manipulate::slider(min = 1, max = L, initial = 1, step = 1),
#'     angle = manipulate::slider(min = 0, max = 360, initial = 40, step = 1))
#'     # theta = manipulate::slider(min = -180, max = 180,
#'     #                            initial = -135, step = 1),
#'     # phi = manipulate::slider(min = -90, max = 90, initial = 0, step = 1))
#'
#'   } else {
#'
#'     rgl::open3d()
#'     rgl::par3d(windowRect = c(80, 125, 1280, 826), zoom = 0.78)
#'     rgl::mfrow3d(nr = ceiling(sqrt(L)), nc = ceiling(sqrt(L)),
#'                  sharedMouse = TRUE)
#'     for (k in seq_len(L)) {
#'
#'       # Update directions
#'       rid_dirs_k[moving_spokes, ] <- matrix(rid_dirs[k, ], ncol = 3,
#'                                             byrow = TRUE)
#'
#'       # Plot s-rep
#'       view_srep(base = avg_base, dirs = rid_dirs_k, radii = avg_radii,
#'                 show_base = FALSE, show_base_pt = FALSE, show_bdry = TRUE,
#'                 show_bdry_pt = FALSE, show_seg = TRUE,
#'                 col_seg = moving_spokes + 2, static = static)
#'
#'     }
#'
#'   }
#'
#' }
