
// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::depends(RcppProgress)]]
#include <RcppArmadillo.h>
#include <Rmath.h>
#include <R.h>
#include <progress.hpp>
#include <progress_bar.hpp>

// Declaration for functions
arma::mat proj_polysph(arma::mat x, arma::uvec ind_dj);
Rcpp::List proj_grad_kde_polysph(arma::mat x, arma::mat X, arma::uvec d,
                                 arma::vec h, Rcpp::NumericVector weights,
                                 bool wrt_unif, bool normalized,
                                 bool norm_x, bool norm_X, arma::uword kernel,
                                 arma::uword kernel_type, double k,
                                 bool proj_alt, bool fix_u1, bool sparse);
arma::vec dist_polysph(arma::mat x, arma::mat y, arma::uvec ind_dj,
                       bool norm_x, bool norm_y, bool std);
arma::vec kde_polysph(arma::mat x, arma::mat X, arma::uvec d, arma::vec h,
                      Rcpp::NumericVector weights, bool log, bool wrt_unif,
                      bool normalized, bool intrinsic, bool norm_x, bool norm_X,
                      arma::uword kernel, arma::uword kernel_type, double k);

//' @title Euler algorithms for polyspherical density ridge estimation
//'
//' @description Functions to perform density ridge estimation on the
//' polysphere \eqn{\mathcal{S}^{d_1} \times \cdots \times \mathcal{S}^{d_r}}
//' through the Euler algorithm in standard, parallel, or block mode.
//'
//' @param x a matrix of size \code{c(nx, sum(d) + r)} with the starting points
//' for the Euler algorithm.
//' @inheritParams kde_polysph
//' @inheritParams proj_grad_kde_polysph
//' @param h_euler vector of size \code{r} with the advance steps in the Euler
//' method. Set internally as \code{h} if not provided.
//' @param N maximum number of Euler iterations. Defaults to \code{1e3}.
//' @param eps convergence tolerance. Defaults to \code{1e-5}.
//' @param keep_paths keep the Euler paths to the ridge? Defaults to
//' \code{FALSE}.
//' @param show_prog display a progress bar for \code{x}? Defaults to
//' \code{TRUE}.
//' @param show_prog_j display a progress bar for \code{N}? Defaults to
//' \code{FALSE}.
//' @param ind_blocks indexes of the blocks, a vector or length \code{r}.
//' @param ... further parameters passed to \code{\link{euler_ridge}}.
//' @param cores cores to use. Defaults to \code{1}.
//' @details \code{euler_ridge} is the main function to perform density ridge
//' estimation through the Euler algorithm from the starting values \code{x}
//' to initiate the ridge path. The function \code{euler_ridge_parallel}
//' parallelizes on the starting values \code{x}. The function
//' \code{euler_ridge_block} runs the Euler algorithm marginally in blocks
//' of hyperspheres, instead of jointly in the whole polysphere. This function
//' requires that all the dimensions are the same.
//' @return The three functions return a list with the following fields:
//' \item{ridge_y}{a matrix of size \code{c(nx, sum(d) + r)} with the end
//' points of Euler algorithm defining the estimated ridge.}
//' \item{lamb_norm_y}{a matrix of size \code{c(nx, sum(d) + r)} with the
//' Hessian eigenvalues (largest to smallest) evaluated at end points.}
//' \item{log_dens_y}{a column vector of size \code{c(nx, 1)} with the
//' logarithm of the density at end points.}
//' \item{paths}{an array of size \code{c(nx, sum(d) + r, N + 1)} containing
//' the Euler paths.}
//' \item{start_x}{a matrix of size \code{c(nx, sum(d) + r)} with the starting
//' points for the Euler algorithm.}
//' \item{iter}{a column vector of size \code{c(nx, 1)} counting the iterations
//' required for each point.}
//' \item{conv}{a column vector of size \code{c(nx, 1)} with convergence flags.}
//' \item{d}{vector \code{d}.}
//' \item{h}{bandwidth used for the kernel density estimator.}
//' \item{error}{a column vector of size \code{c(nx, 1)} indicating if errors
//' were found for each path.}
//' @examples
//' ## Test on S^2 with a small circle trend
//'
//' # Sample
//' r <- 1
//' d <- 2
//' n <- 50
//' ind_dj <- comp_ind_dj(d = d)
//' set.seed(987204452)
//' X <- r_path_s2r(n = n, r = r, spiral = FALSE, Theta = cbind(c(1, 0, 0)),
//'                 sigma = 0.35)[, , 1]
//' col_X_alp <- viridis::viridis(n, alpha = 0.25)
//' col_X <- viridis::viridis(n)
//'
//' # Euler
//' h_rid <- 0.5
//' h_eu <- h_rid^2
//' N <- 30
//' eps <- 1e-6
//' Y <- euler_ridge(x = X, X = X, d = d, h = h_rid, h_euler = h_eu,
//'                  N = N, eps = eps, keep_paths = TRUE)
//' Y
//'
//' # Visualization
//' i <- N # Between 1 and N
//' sc3 <- scatterplot3d::scatterplot3d(Y$paths[, , 1], color = col_X_alp,
//'                                     pch = 19, xlim = c(-1, 1),
//'                                     ylim = c(-1, 1), zlim = c(-1, 1),
//'                                     xlab = "x", ylab = "y", zlab = "z")
//' sc3$points3d(rbind(Y$paths[, , i]), col = col_X, pch = 16, cex = 0.75)
//' invisible(sapply(seq_len(nrow(Y$paths)), function(k) {
//'   sc3$points3d(t(Y$paths[k, , ]), col = col_X_alp[k], type = "l")
//' }))
//' @export
// [[Rcpp::export]]
Rcpp::List euler_ridge(arma::mat x, arma::mat X, arma::uvec d, arma::vec h,
                       Rcpp::NumericVector h_euler =
                         Rcpp::NumericVector::create(),
                       Rcpp::NumericVector weights =
                         Rcpp::NumericVector::create(), bool wrt_unif = false,
                       bool normalized = true, bool norm_x = false,
                       bool norm_X = false, arma::uword kernel = 1,
                       arma::uword kernel_type = 1, double k = 10.0,
                       arma::uword N = 1e3, double eps = 1e-5,
                       bool keep_paths = false, bool proj_alt = true,
                       bool fix_u1 = true, bool sparse = false,
                       bool show_prog = true, bool show_prog_j = false) {

  // Number of starting points
  arma::uword m = x.n_rows;

  // How many S^dj's
  arma::uword r = d.n_elem;

  // Dimensions
  arma::uword p = X.n_cols;
  if (p != x.n_cols) {

    Rcpp::stop("Dimension of X mismatches with x.");

  }
  if (r != h.n_elem) {

    Rcpp::stop("Size of h mismatches with d.");

  }
  if (p != arma::accu(d + 1)) {

    Rcpp::stop("Dimension of X mismatches with d.");

  }

  // Indexes with begin and end of each S^dj
  arma::uvec ind_dj = arma::conv_to<arma::uvec>::from(arma::zeros(r + 1));
  ind_dj.tail(r) = d + 1;
  ind_dj = arma::cumsum(ind_dj);

  // Transform h_euler to arma::vec
  arma::vec h_eu = Rcpp::as<arma::vec>(Rcpp::wrap(h_euler));

  // h_euler given?
  if (h_eu.n_elem == 0) {

    // Fill with h
    h_eu = h;

  }

  // Compute h_eu_diamond to use in h ◇ eta
  arma::rowvec h_eu_diamond = arma::zeros(p).t();
  for (arma::uword dj = 0; dj < r; dj++) {

    h_eu_diamond.subvec(ind_dj(dj), ind_dj(dj + 1) - 1).fill(h_eu(dj));

  }

  // Iterations and convergences
  arma::uvec iter(m);
  iter.fill(N);
  arma::uvec conv(m);
  conv.fill(0);

  // Progress bar
  Progress prog(m, show_prog);
  Progress prog_j(N, show_prog_j);

  // Save paths?
  arma::cube paths(1, 1, 1);
  paths.fill(arma::datum::nan);
  if (keep_paths) {

    paths.set_size(m, p, N + 1);
    paths.fill(arma::datum::nan);
    paths.slice(0) = x;

  }

  // Kde-normalized eigenvalues
  arma::uword n_lamb = sparse ? 1 : p;
  arma::mat lamb_norm = arma::zeros(m, n_lamb);

  // Loop on grid points
  arma::mat y = x;
  for (arma::uword i = 0; i < m; i++) {

    // i-th grid point, to be updated
    arma::rowvec yi_old = y.row(i);
    arma::rowvec yi_new = yi_old;

    // Euler loop
    for (arma::uword j = 0; j < N; j++) {

      // Projected normalized gradient and eigenvalues
      Rcpp::List proj = proj_grad_kde_polysph(yi_old, X, d, h, weights,
                                              wrt_unif, normalized, norm_x,
                                              norm_X, kernel, kernel_type, k,
                                              proj_alt, fix_u1, sparse);
      arma::rowvec eta = proj["eta"];
      arma::rowvec lamb = proj["lamb_norm"];

      // Advance
      yi_new = yi_old + h_eu_diamond % eta;

      // Project to polysphere
      yi_new = proj_polysph(yi_new, ind_dj);

      // Save paths?
      if (keep_paths) {

        paths(arma::span(i), arma::span::all, arma::span(j + 1)) = yi_new;

      }

      // Save last value
      y.row(i) = yi_new;

      // Convergence? Use standardized distance. Normalize arguments for better
      // accuracy.
      if (arma::as_scalar(dist_polysph(yi_old, yi_new, ind_dj,
                                       true, true, true)) < eps) {

        // Save iteration and convergence flag
        iter(i) = j + 1;
        conv(i) = 1;

        // Save eigenvalues
        lamb_norm.row(i) = lamb;

        // Set final values in the path
        if (keep_paths && j < (N - 1)) {

          paths(arma::span(i), arma::span::all, arma::span(j + 2, N)) =
            arma::repmat(yi_new.t(), 1, N - j - 1);

        }

        // Stop Euler
        break;

      } else {

        // Update
        yi_old = yi_new;

      }

      // Update
      prog_j.increment();

    }

    // Update and check for interruptions
    prog.increment();
    Rcpp::checkUserInterrupt();

  }

  // Compute log-density final points
  arma::vec log_dens = kde_polysph(y, X, d, h, weights, true, wrt_unif,
                                   normalized, false, norm_x, norm_X, kernel,
                                   kernel_type, k);

  // Return a Rcpp list
  return Rcpp::List::create(Rcpp::Named("ridge_y") = y,
                            Rcpp::Named("lamb_norm_y") = lamb_norm,
                            Rcpp::Named("log_dens_y") = log_dens,
                            Rcpp::Named("paths") = paths,
                            Rcpp::Named("start_x") = x,
                            Rcpp::Named("iter") = iter,
                            Rcpp::Named("conv") = conv,
                            Rcpp::Named("d") = d,
                            Rcpp::Named("h") = h);

}
