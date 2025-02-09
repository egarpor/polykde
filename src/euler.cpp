
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

//' @title Euler step for density ridge estimation
//'
//' @description TODO
//'
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
//' @return TODO
//' @examples
//' # TODO
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
