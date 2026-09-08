
// [[Rcpp::depends("RcppArmadillo")]]
#include <RcppArmadillo.h>
#include <Rmath.h>
#include <R.h>

// Declaration for functions
arma::mat proj_polysph(arma::mat x, arma::uvec ind_dj);
arma::vec kde_polysph(arma::mat x, arma::mat X, arma::uvec d, arma::vec h,
                      Rcpp::NumericVector weights, bool log, bool wrt_unif,
                      bool normalized, bool intrinsic, bool norm_x, bool norm_X,
                      arma::uword kernel, arma::uword kernel_type, double k,
                      bool loo);
arma::mat sfp(arma::mat t);

// Constants
const double log_two = std::log(2.0);
const double log_two_M_PI = std::log(2.0 * M_PI);
const double log_M_PI = std::log(M_PI);


//' @title Polyspherical kernel density estimator
//'
//' @description Computes the kernel density estimator for data on the
//' polysphere \eqn{\mathbb{S}^{d_1} \times \cdots \times \mathbb{S}^{d_r}}.
//' Given a sample \eqn{\boldsymbol{X}_1,\ldots,\boldsymbol{X}_n}, this
//' estimator is
//' \deqn{\hat{f}(\boldsymbol{x};\boldsymbol{h})=\sum_{i=1}^n
//' L_{\boldsymbol{h}}(\boldsymbol{x},\boldsymbol{X}_i)}
//' for a kernel \eqn{L} and a vector of bandwidths \eqn{\boldsymbol{h}}.
//'
//' @param x a matrix of size \code{c(nx, sum(d) + r)} with the evaluation
//' points.
//' @param X a matrix of size \code{c(n, sum(d) + r)} with the sample.
//' @param d vector of size \code{r} with dimensions.
//' @param h vector of size \code{r} with bandwidths.
//' @param weights weights for each observation. If provided, a vector of size
//' \code{n} with the weights for multiplying each kernel. If not provided, set
//' internally to \code{rep(1 / n, n)}, which gives the standard estimator.
//' @param log compute the logarithm of the density? Defaults to \code{FALSE}.
//' @param wrt_unif flag to return a density with respect to the uniform
//' measure. If \code{FALSE} (default), the density is with respect to the
//' Lebesgue measure.
//' @param normalized flag to compute the normalizing constant of the kernel and
//' include it in the kernel density estimator. Defaults to \code{TRUE}.
//' @param intrinsic use the intrinsic distance, instead of the
//' extrinsic-chordal distance, in the kernel? Defaults to \code{FALSE}.
//' @param norm_x,norm_X ensure a normalization of the data? Defaults to
//' \code{FALSE}.
//' @param kernel kernel employed: \code{1} for von Mises--Fisher (default);
//' \code{2} for Epanechnikov; \code{3} for softplus.
//' @param kernel_type type of kernel employed: \code{1} for product kernel
//' (default); \code{2} for spherically symmetric kernel.
//' @param k softplus kernel parameter. Defaults to \code{10.0}.
//' @param loo leave-one-out: skip the \eqn{i}-th observation when evaluating
//' at \eqn{\boldsymbol{X}_i}? Requires \code{x = X} and \eqn{n \ge 2}.
//' Defaults to \code{FALSE}.
//' @return A column matrix of size \code{c(nx, 1)} with the evaluation of the
//' kernel density estimator.
//' @references
//' García-Portugués, E. and Meilán-Vila, A. (2026). Kernel density estimation
//' with polyspherical data and its applications. \emph{Journal of the American
//' Statistical Association}, 121(553):427--439.
//' \doi{10.1080/01621459.2025.2521898}.
//' @seealso \code{\link{log_cv_kde_polysph}},
//' \code{\link{grad_hess_kde_polysph}}, \code{\link{r_kde_polysph}},
//' \code{\link{bw_rot_polysph}}, \code{\link{bw_cv_polysph}},
//' \code{\link{bw_mrot_polysph}}.
//' @examples
//' # Simple check on S^1 x S^2
//' n <- 1e3
//' d <- c(1, 2)
//' mu <- c(0, 1, 0, 0, 1)
//' kappa <- c(5, 5)
//' h <- c(0.2, 0.2)
//' X <- r_vmf_polysph(n = n, d = d, mu = mu, kappa = kappa)
//' kde_polysph(x = rbind(mu), X = X, d = d, h = h)
//' d_vmf_polysph(x = rbind(mu), d = d, mu = mu, kappa = kappa)
//' @export
// [[Rcpp::export]]
arma::vec kde_polysph(arma::mat x, arma::mat X, arma::uvec d, arma::vec h,
                      Rcpp::NumericVector weights =
                        Rcpp::NumericVector::create(), bool log = false,
                      bool wrt_unif = false, bool normalized = true,
                      bool intrinsic = false, bool norm_x = false,
                      bool norm_X = false, arma::uword kernel = 1,
                      arma::uword kernel_type = 1, double k = 10.0,
                      bool loo = false) {

  // Sample size
  arma::uword n = X.n_rows;

  // How many S^dj's
  arma::uword r = d.n_elem;

  // Dimensions
  arma::uword p = X.n_cols;
  if (p != x.n_cols) {

    Rcpp::stop("Dimensions of X and x mismatch.");

  }
  if (loo && !arma::approx_equal(x, X, "absdiff", 1e-10)) {

    Rcpp::stop("loo = TRUE requires x = X.");

  }
  if (loo && n < 2) {

    Rcpp::stop("loo = TRUE requires n >= 2.");

  }
  if (p != arma::accu(d + 1)) {

    Rcpp::stop("Dimension of X mismatches with d.");

  }
  if (r != h.n_elem) {

    Rcpp::stop("Size of h mismatches with d.");

  }

  // 1 / h^2
  arma::vec inv_h2 = 1.0 / arma::square(h);

  // d as double to avoid casting issues
  arma::vec dd = arma::conv_to<arma::vec>::from(d);

  // Transform NumericVector to arma::vec
  arma::vec log_weights = Rcpp::as<arma::vec>(weights);

  // Are weights given?
  if (log_weights.n_elem == 0) {

    // Fill with -log(n) sample size
    log_weights.set_size(n);
    log_weights.fill(-std::log(static_cast<double>(n)));

  } else {

    // Normalize weights to be positive and sum 1, and take its logarithm
    log_weights = log_weights.clamp(0, arma::datum::inf);
    log_weights /= arma::accu(log_weights);
    log_weights = arma::log(log_weights);

  }

  // Log-uniform area measure
  arma::rowvec log_w_d = arma::zeros<arma::rowvec>(r);
  if (wrt_unif) {

    log_w_d = (log_two + 0.5 * (dd + 1.0) * log_M_PI -
      arma::lgamma(0.5 * (dd + 1.0))).t();

  }
  double log_w_d_sum = arma::accu(log_w_d);

  // Kernel normalizing constants
  arma::rowvec log_const_h = arma::zeros<arma::rowvec>(r);
  double log_const_sph_h = 0;
  if (normalized) {

    if (!intrinsic) {

      if (kernel == 1) { // vMF

        // Special case (S^2)^r: use closed-form for the log-exp-scaled Bessels
        // to avoid overflows
        // I_{1/2}(z) = sqrt(2 / (pi * z)) * sinh(z)
        //            = 1 / sqrt(2 * pi * z) * (exp(z) - exp(-z))
        //            = 1 / sqrt(2 * pi * z) * exp(z) * (1 - exp(-2 * z))
        // log(I_{1/2}(z) * exp(-z)) = log(1 / sqrt(2 * pi * z))
        //                             + log(1 - exp(-2 * z))
        //                           = -0.5 * log(2 * pi) - 0.5 * log(z)
        //                             + log1p(-exp(-2 * z))
        if (arma::all(d == 2)) {

          // Log-normalizing constant
          arma::vec log_exps_bessel_no_cts =
            arma::log1p(-arma::exp(-2.0 * inv_h2));
          log_const_h = (arma::log(inv_h2) - log_two_M_PI -
            log_exps_bessel_no_cts).t();

        // General case: use log-exp-scaled Bessels to avoid overflows as
        // possible
        } else {

          // Bessels
          arma::vec log_exps_bessel = arma::zeros(r);
          for (arma::uword dj = 0; dj < r; dj++) {

            // Codification: bessel_i(x, nu, expon.scaled + 1)
            log_exps_bessel(dj) =
              std::log(R::bessel_i(inv_h2(dj), 0.5 * (d(dj) - 1.0), 2));

          }

          // Log-normalizing constant
          log_const_h = (0.5 * (dd - 1.0) % arma::log(inv_h2) -
            0.5 * (dd + 1.0) * log_two_M_PI - log_exps_bessel).t();

        }

        // Normalizing constant for spherically symmetric version: sum of
        // normalizing constants because of the product kernel
        log_const_sph_h = arma::accu(log_const_h);

        // Non-finite vMF normalizing constant?
        if (!log_const_h.is_finite()) {

          Rcpp::stop("Non-finite vMF normalizing constant: bandwidth too small?");

        }

      } else if (kernel == 2) { // Epanechnikov product

        Rcpp::Function c("c_kern");
        Rcpp::NumericVector const_epa = c(h, d, 2, kernel_type, 1,
                                          true, false, false);
        log_const_h = const_epa;
        log_const_sph_h = const_epa(0);

      } else if (kernel == 3) { // Softplus product

        // Remove division by sfp(k) in the kernel because then the kernels are
        // computed without this normalizing constant to save computations
        Rcpp::Function c("c_kern");
        Rcpp::NumericVector const_sfp = c(h, d, 3, kernel_type, k,
                                          true, false, false);
        log_const_h = const_sfp;
        log_const_sph_h = const_sfp(0);

      } else {

        Rcpp::stop("kernel must be 1, 2, or 3.");

      }

    } else {

      Rcpp::Function c("c_kern");
      Rcpp::NumericVector const_kernel = c(h, d, kernel, kernel_type, k,
                                           true, false, true);
      log_const_h = const_kernel;
      log_const_sph_h = const_kernel(0);

    }

  }

  // Indexes with begin and end of each S^dj
  arma::uvec ind_dj = arma::conv_to<arma::uvec>::from(arma::zeros(r + 1));
  ind_dj.tail(r) = d + 1;
  ind_dj = arma::cumsum(ind_dj);

  // Normalization?
  if (norm_X) {

    X = proj_polysph(X, ind_dj);

  }
  if (norm_x) {

    x = proj_polysph(x, ind_dj);

  }

  // Loop on evaluation points for computing the log-kde
  arma::uword nx = x.n_rows;
  arma::vec log_kde = arma::zeros(nx);
  arma::mat x_prime_X = arma::zeros(n, r);
  arma::vec log_wd_const_weights = log_weights +
    (log_w_d_sum + log_const_sph_h);
  for (arma::uword ki = 0; ki < nx; ki++) {

    // Compute argument of the kernels for extrinsic or intrinsic distances
    if (!intrinsic) {

      // Compute (1 - x'X) / h^2 for each S^dj
      for (arma::uword dj = 0; dj < r; dj++) {

        arma::mat x_prime_Xj = X.cols(ind_dj(dj), ind_dj(dj + 1) - 1);
        x_prime_Xj.each_row() %=
          x(ki, arma::span(ind_dj(dj), ind_dj(dj + 1) - 1));
        x_prime_X.col(dj) = (1.0 - arma::sum(x_prime_Xj, 1)) * inv_h2(dj);

      }

    } else {

      // Compute arccos(x'X)^2 / (2 * h^2) for each S^dj
      for (arma::uword dj = 0; dj < r; dj++) {

        arma::mat x_prime_Xj = X.cols(ind_dj(dj), ind_dj(dj + 1) - 1);
        x_prime_Xj.each_row() %=
          x(ki, arma::span(ind_dj(dj), ind_dj(dj + 1) - 1));
        x_prime_X.col(dj) = arma::square(arma::acos(
          arma::clamp(arma::sum(x_prime_Xj, 1), -1.0, 1.0))) *
            (0.5 * inv_h2(dj));

      }

    }

    // Product or spherically symmetric kernel?
    arma::vec log_L = arma::zeros(n);
    if (kernel_type == 1) {

      // vMF kernel
      if (kernel == 1) {

        // Multiply by area measure and normalizing constant (later exp(-x) is
        // taken, so the constants multiply)
        x_prime_X.each_row() -= log_w_d + log_const_h;

        // Product of normalized kernels across spheres
        log_L = arma::sum(-x_prime_X, 1) + log_weights;

      // Epanechnikov kernel
      } else if (kernel == 2) {

        // Apply [1 - (1 - x'X) / h^2] * 1_{0 <= (1 - x'X) / h^2 <= 1} and log()
        // This sets (1 - x'X) / h^2 > 1 to 1, which afterwards become
        // log(1 - 1) = log(0) = -Inf (as intended)
        x_prime_X.clamp(-arma::datum::inf, 1.0);
        x_prime_X = arma::log1p(-x_prime_X);

        // Product of normalized kernels across spheres
        x_prime_X.each_row() += log_w_d + log_const_h;
        log_L = arma::sum(x_prime_X, 1) + log_weights;

      // Softplus kernel
      } else if (kernel == 3) {

        // Apply softplus kernel and log()
        x_prime_X = arma::log(sfp(k * (1.0 - x_prime_X)));

        // Product of normalized kernels across spheres
        x_prime_X.each_row() += log_w_d + log_const_h;
        log_L = arma::sum(x_prime_X, 1) + log_weights;

      }

    } else if (kernel_type == 2) {

      // Sum along all S^dj's to fed the kernels
      arma::vec x_prime_X_sum = sum(x_prime_X, 1);

      // vMF kernel
      if (kernel == 1) {

        // Log-kernel times log-normalizing constants and log-weights
        log_L = -x_prime_X_sum + log_wd_const_weights;

      // Epanechnikov kernel
      } else if (kernel == 2) {

        // Apply [1 - sum_j((1 - x_j'X) / h_j^2)] *
        //   1_{0 <= sum_j((1 - x_j'X) / h_j^2) <= 1} and log()
        // This sets sum_j((1 - x_j'X) / h_j^2) > 1 to 1, which afterwards
        // become log(1 - 1) = log(0) = -Inf (as intended)
        x_prime_X_sum.clamp(-arma::datum::inf, 1.0);
        x_prime_X_sum = arma::log1p(-x_prime_X_sum);

        // Log-kernel times log-normalizing constants and log-weights
        log_L = x_prime_X_sum + log_wd_const_weights;

      // Softplus kernel
      } else if (kernel == 3) {

        // Apply softplus kernel and log()
        x_prime_X_sum = arma::log(sfp(k * (1.0 - x_prime_X_sum)));

        // Log-kernel times log-normalizing constants and log-weights
        log_L = x_prime_X_sum + log_wd_const_weights;

      }

    }

    // Leave-one-out: drop the ki-th observation from the sum
    if (loo) {

      log_L(ki) = -arma::datum::inf;

    }

    // Sum kernels with LogSumExp trick to avoid overflows in large
    // exponentials -- but be careful, as if max_log_L = -inf, then
    // log_L - max_log_L = -inf + inf = nan! (and the subtraction of max_log_L
    // really does not make sense as all the elements are -inf)
    double max_log_L = arma::max(log_L);
    if (max_log_L > -arma::datum::inf) {

      log_L -= max_log_L;

    }

    // Log-kde computed using the LogSumExp trick
    log_kde(ki) = max_log_L + std::log(arma::accu(arma::exp(log_L)));

  }

  // Leave-one-out: renormalize the weights of the remaining observations, so
  // that they sum one. A weight equal to one leaves no observations behind, so
  // the -Inf - (-Inf) = NaN arising there is a null density.
  if (loo) {

    log_kde -= arma::log1p(-arma::exp(log_weights));
    log_kde.replace(arma::datum::nan, -arma::datum::inf);

  }

  // Return log-kde or kde
  if (!log) {

    log_kde = arma::exp(log_kde);

  }
  return log_kde;

}


//' @title Cross-validation for the polyspherical kernel density estimator
//'
//' @description Computes the logarithm of the cross-validated kernel density
//' estimator: \eqn{\log \hat{f}_{-i}(\boldsymbol{X}_i;\boldsymbol{h})},
//' \eqn{i = 1, \ldots, n.}
//'
//' @inheritParams kde_polysph
//' @param norm_X ensure a normalization of the data? Defaults to \code{FALSE}.
//' @return A column matrix of size \code{c(n, 1)} with the evaluation of the
//' logarithm of the cross-validated kernel density estimator.
//' @examples
//' # Simple check on S^1 x S^2
//' n <- 5
//' d <- c(1, 2)
//' h <- c(0.2, 0.2)
//' X <- r_unif_polysph(n = n, d = d)
//' log_cv_kde_polysph(X = X, d = d, h = h)
//' kde_polysph(x = X[1, , drop = FALSE], X = X[-1, ], d = d, h = h, log = TRUE)
//' @export
// [[Rcpp::export]]
arma::vec log_cv_kde_polysph(arma::mat X, arma::uvec d, arma::vec h,
                             Rcpp::NumericVector weights =
                               Rcpp::NumericVector::create(),
                             bool wrt_unif = false, bool normalized = true,
                             bool intrinsic = false, bool norm_X = false,
                             arma::uword kernel = 1,
                             arma::uword kernel_type = 1, double k = 10.0) {

  // Sample size
  arma::uword n = X.n_rows;

  // How many S^dj's
  arma::uword r = d.n_elem;

  // Dimensions
  arma::uword p = X.n_cols;
  if (p != arma::accu(d + 1)) {

    Rcpp::stop("Dimension of X mismatches with d.");

  }
  if (r != h.n_elem) {

    Rcpp::stop("Size of h mismatches with d.");

  }
  if (n < 2) {

    Rcpp::stop("X must have at least 2 rows (n >= 2) for cross-validation.");

  }

  // Indexes with begin and end of each S^dj
  arma::uvec ind_dj = arma::conv_to<arma::uvec>::from(arma::zeros(r + 1));
  ind_dj.tail(r) = d + 1;
  ind_dj = arma::cumsum(ind_dj);

  // Normalization?
  if (norm_X) {

    X = proj_polysph(X, ind_dj);

  }

  // Leave-one-out kde evaluated at the sample, except for
  // vMF-product-extrinsic kernel that admits a more efficient computation.
  if (kernel != 1 || kernel_type != 1 || intrinsic) {

    return kde_polysph(X, X, d, h, weights, true, wrt_unif, normalized,
                       intrinsic, false, false, kernel, kernel_type, k, true);

  } else {

    // For the extrinsic product vMF kernel, the log-kernel is
    // log L_h(X_i, X_j) = C + sum_dj (X_{i,dj}'X_{j,dj} - 1) / h_dj^2. The
    // log-kernel matrix is then, up to constants, the Gram matrix X (X S)',
    // with S scaling the coordinates of S^dj by 1 / h_dj^2, which can be
    // blockwise computed. C is computed with a single kde call.
    arma::vec inv_h2 = 1.0 / arma::square(h);
    double C = arma::as_scalar(kde_polysph(X.row(0), X.row(0), d, h,
                                           Rcpp::NumericVector::create(), true,
                                           wrt_unif, normalized, false, false,
                                           false, kernel, kernel_type, k,
                                           false)) - arma::accu(inv_h2);

    // Normalized log-weights, as in kde_polysph()
    arma::vec w = arma::ones(n);
    if (weights.size() > 0) {

      w = arma::clamp(Rcpp::as<arma::vec>(weights), 0.0, arma::datum::inf);

    }
    arma::vec log_w = arma::log(w / arma::accu(w));

    // Sample with the coordinates of S^dj scaled by 1 / h_dj^2
    arma::mat XS = X;
    for (arma::uword dj = 0; dj < r; dj++) {

      XS.cols(ind_dj(dj), ind_dj(dj + 1) - 1) *= inv_h2(dj);

    }

    // Row-wise LogSumExp of the blocks, excluding the diagonal
    arma::uword b = std::max<arma::uword>(1, std::min(n, 250000 / n));
    arma::vec log_cv = arma::zeros(n);
    for (arma::uword i0 = 0; i0 < n; i0 += b) {

      arma::uword i1 = std::min(i0 + b, n) - 1;
      arma::mat G = XS.rows(i0, i1) * X.t();
      G.each_row() += log_w.t();
      G.diag(i0).fill(-arma::datum::inf);

      // Rows without finite entries (all the other weights are zero) get a zero
      // shift, so that they give -Inf instead of NaN
      arma::vec m = arma::max(G, 1);
      m.replace(-arma::datum::inf, 0.0);
      G.each_col() -= m;
      log_cv.subvec(i0, i1) = m + arma::log(arma::sum(arma::exp(G), 1));

    }

    // Add the constants and renormalize the weights of the remaining
    // observations, so that they sum one, as in kde_polysph().
    log_cv += C - arma::log1p(-arma::exp(log_w));
    log_cv.replace(arma::datum::nan, -arma::datum::inf);
    return log_cv;

  }

}
