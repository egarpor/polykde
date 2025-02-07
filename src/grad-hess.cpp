
// [[Rcpp::depends("RcppArmadillo")]]
#include <RcppArmadillo.h>
#include <Rmath.h>
#include <R.h>

// Declaration for functions
arma::mat proj_polysph(arma::mat x, arma::uvec ind_dj);
arma::cube diamond_crossprod(arma::mat X, arma::uvec ind_dj);
arma::mat s(arma::mat A, bool add = false);
Rcpp::List grad_hess_kde_polysph(arma::mat x, arma::mat X, arma::uvec d,
                                 arma::vec h, Rcpp::NumericVector weights,
                                 bool projected, bool proj_alt,
                                 bool norm_grad_hess, bool log, bool wrt_unif,
                                 bool normalized, bool norm_x, bool norm_X,
                                 arma::uword kernel, arma::uword kernel_type,
                                 double k);
Rcpp::List AP(arma::rowvec x, arma::rowvec v, arma::uvec ind_dj, bool orth);

// Constants
const double log_two = std::log(2.0);
const double log_two_M_PI = std::log(2.0 * M_PI);
const double log_M_PI = std::log(M_PI);


//' @title Gradient and Hessian of the polyspherical kde
//'
//' @inheritParams kde_polysph
//' @param projected compute the \emph{projected} gradient and Hessian that
//' accounts for the radial projection? Defaults to \code{TRUE}.
//' @param proj_alt alternative projection. Defaults to \code{TRUE}.
//' @param norm_grad_hess normalize the gradient and Hessian dividing by the
//' kde? Defaults to \code{FALSE}.
//' @export
// [[Rcpp::export]]
Rcpp::List grad_hess_kde_polysph(arma::mat x, arma::mat X, arma::uvec d,
                                 arma::vec h, Rcpp::NumericVector weights =
                                   Rcpp::NumericVector::create(),
                                 bool projected = true, bool proj_alt = true,
                                 bool norm_grad_hess = false, bool log = false,
                                 bool wrt_unif = false, bool normalized = true,
                                 bool norm_x = false, bool norm_X = false,
                                 arma::uword kernel = 1,
                                 arma::uword kernel_type = 1, double k = 10.0) {

  /*
   * Code replicated from kde_polysph(). NEW indicates the main differences.
   */

  // Stop if kernel_type is not product (NEW)
  if (kernel_type != 1) {

    Rcpp::stop("Only kernel_type = 1 is supported.");

  }

  // Sample size
  arma::uword n = X.n_rows;

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

  // 1 / h^2
  arma::vec inv_h2 = 1.0 / arma::square(h);

  // d as double to avoid casting issues
  arma::vec dd = arma::conv_to<arma::vec>::from(d);

  // Transform NumericVector to arma::vec
  arma::vec log_weights = Rcpp::as<arma::vec>(Rcpp::wrap(weights));

  // Are weights given?
  if (log_weights.n_elem == 0) {

    // Fill with -log(n) sample size
    log_weights.set_size(n);
    log_weights.fill(-std::log(n));

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

  // Kernel normalizing constants
  arma::rowvec log_const_h = arma::zeros<arma::rowvec>(r);
  if (normalized) {

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
        arma::vec log_exps_bessel_no_cts = arma::log1p(-arma::exp(-2.0 * inv_h2));
        log_const_h = (arma::log(inv_h2) - log_two_M_PI -
          log_exps_bessel_no_cts).t();

      // General case: use log-exp-scaled Bessels to avoid overflows as possible
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

      // Non-finite vMF normalizing constant?
      if (!log_const_h.is_finite()) {

        Rcpp::stop("Non-finite vMF normalizing constant: bandwidth too small?");

      }

    } else if (kernel == 2 && kernel_type == 1) { // Epanechnikov product

      Rcpp::Function c("c_kern");
      Rcpp::NumericVector const_epa = c(h, d, 2, 1, true, false);
      log_const_h = const_epa;

    } else if (kernel == 3 && kernel_type == 1) { // Softplus product

      Rcpp::Function c("c_kern");
      Rcpp::NumericVector const_sfp = c(h, d, 3, k, true, false);
      log_const_h = const_sfp;
      // Remove division by sfp(k) in the kernel because then the kernels are
      // computed without this normalizing constant to save computations

    } else {

      Rcpp::stop("Unsupported combination of (\"kernel\", \"kernel_type\"). Only these are available: (1, 1), (2, 1), (3, 1), (1, 2). Consider computing the estimator unnormalized?");

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

  // -2 * log(h) for diamond multiplication in gradient (NEW)
  arma::rowvec log_inv_h2_diamond_grad = arma::zeros(p).t();
  for (arma::uword dj = 0; dj < r; dj++) {

    log_inv_h2_diamond_grad.subvec(ind_dj(dj), ind_dj(dj + 1) - 1).fill(h(dj));

  }
  log_inv_h2_diamond_grad = -2.0 * arma::log(log_inv_h2_diamond_grad);

  // -2 * log(h) for diamond multiplication in Hessian (NEW)
  arma::mat log_inv_h2_diamond_hess = arma::zeros(p, p);
  for (arma::uword di = 0; di < r; di++) {
    for (arma::uword dj = 0; dj < r; dj++) {

      log_inv_h2_diamond_hess(
        arma::span(ind_dj(di), ind_dj(di + 1) - 1),
        arma::span(ind_dj(dj), ind_dj(dj + 1) - 1)).fill(h(di) * h(dj));

    }
  }
  log_inv_h2_diamond_hess = -2.0 * arma::log(log_inv_h2_diamond_hess);

  // Create Xi ◇ Xi' (NEW)
  arma::cube Xi_diamond = diamond_crossprod(X, ind_dj);

  // Loop on evaluation points for computing the log-kde, gradient, and Hessian
  arma::uword nx = x.n_rows;
  arma::vec log_kde = arma::zeros(nx);
  arma::mat grad_kde = arma::zeros(nx, p);
  arma::cube hess_kde = arma::zeros(nx, p, p);
  for (arma::uword ki = 0; ki < nx; ki++) {

    // Compute (1 - x'X) / h^2 for each S^dj
    arma::mat x_prime_X = arma::zeros(n, r);
    for (arma::uword dj = 0; dj < r; dj++) {

      arma::mat x_prime_Xj = X.cols(ind_dj(dj), ind_dj(dj + 1) - 1);
      x_prime_Xj.each_row() %=
        x(ki, arma::span(ind_dj(dj), ind_dj(dj + 1) - 1));
      x_prime_X.col(dj) = (1.0 - arma::sum(x_prime_Xj, 1)) * inv_h2(dj);

    }

    // vMF kernel
    arma::vec log_L = arma::zeros(n);
    if (kernel == 1) {

      // Multiply by area measure and normalizing constant (later exp(-x) is
      // taken, so the constants multiply)
      x_prime_X.each_row() -= log_w_d + log_const_h;

      // Product of normalized kernels across spheres
      log_L = arma::sum(-x_prime_X, 1) + log_weights;

    }

    // Get maximum kernel logarithm and scale kernels (NEW)
    double max_log_L = arma::max(log_L);
    if (max_log_L > -arma::datum::inf) {

      log_L -= max_log_L;

    }

    // Log-kde without maximum log, computed using the LogSumExp trick (NEW)
    log_kde(ki) = std::log(arma::accu(arma::exp(log_L)));

    // Gradient: multiply by weights and sum (NEW)
    grad_kde.row(ki) = arma::exp(log_L).t() * X;

    // // Gradient: multiply by maximum kernel and apply bandwidths (NEW)
    // grad_kde.row(ki) %= arma::exp(log_inv_h2_diamond_grad + max_log_L);

    // Normalized gradient: apply bandwidths and take ratio with kde. Recall
    // the max_log_L cancels out (NEW).
    grad_kde.row(ki) %= arma::exp(log_inv_h2_diamond_grad - log_kde(ki));

    // Hessian: multiply by weights and sum (NEW)
    arma::cube Xi_diamond_ki = Xi_diamond;
    for (arma::uword i = 0; i < n; i++) {

      Xi_diamond_ki.row(i) *= std::exp(log_L(i));

    }
    hess_kde.row(ki) = arma::sum(Xi_diamond_ki, 0);

    // // Hessian: multiply by maximum kernel and apply bandwidths (NEW)
    // hess_kde.row(ki) %= arma::exp(log_inv_h2_diamond_hess + max_log_L);

    // Normalized Hessian: apply bandwidths and take ratio with kde. Recall
    // the max_log_L cancels out (NEW).
    hess_kde.row(ki) %= arma::exp(log_inv_h2_diamond_hess - log_kde(ki));

    // Add maximum log to log-kde (NEW)
    log_kde(ki) += max_log_L;

    // Not kde-normalized gradient and Hessian? (NEW)
    if (!norm_grad_hess) {

      // Multiply by kde to undo normalization (NEW)
      double e = std::exp(log_kde(ki));
      grad_kde.row(ki) *= e;
      hess_kde.row(ki) *= e;

    }

  }

  // Project the gradient and Hessian (NEW)
  if (projected) {

    // Gradient
    arma::mat proj_grad_kde = grad_kde;

    // Loop on components
    for (arma::uword dj = 0; dj < r; dj++) {

      // j-terms
      arma::uword ini_dj = ind_dj(dj);
      arma::uword end_dj = ind_dj(dj + 1) - 1;
      arma::mat proj_grad_kde_dj = grad_kde.cols(ini_dj, end_dj);
      arma::mat Idj = arma::eye(d(dj) + 1, d(dj) + 1);

      // Loop on evaluation points
      for (arma::uword i = 0; i < nx; i++) {

        arma::rowvec x_i_dj = x.row(i).subvec(ini_dj, end_dj);
        proj_grad_kde_dj.row(i) = proj_grad_kde_dj.row(i) *
          (Idj - x_i_dj.t() * x_i_dj);

      }

      // Save projected gradient
      proj_grad_kde.cols(ini_dj, end_dj) = proj_grad_kde_dj;

    }

    // Hessian
    arma::cube proj_hess_kde = hess_kde;

    // Project Hessian to the tangent
    if (proj_alt) {

      // Loop on evaluation points
      for (arma::uword i = 0; i < nx; i++) {

        Rcpp::List AP_i = AP(x.row(i), grad_kde.row(i), ind_dj, false);
        arma::mat A_i = AP_i["A"];
        arma::mat P_i = AP_i["P"];
        arma::mat proj_hess_kde_i = proj_hess_kde.row(i);
        proj_hess_kde.row(i) = P_i * (proj_hess_kde_i - A_i) * P_i;

      }

    // Hessian from radial projection
    } else {

      // Loop on the j-component
      for (arma::uword dj = 0; dj < r; dj++) {

        // j-terms
        arma::uword ini_dj = ind_dj(dj);
        arma::uword end_dj = ind_dj(dj + 1) - 1;
        arma::mat Idj = arma::eye(d(dj) + 1, d(dj) + 1);

        // Loop on the k-component
        for (arma::uword dk = 0; dk < r; dk++) {

          // k-terms
          arma::uword ini_dk = ind_dj(dk);
          arma::uword end_dk = ind_dj(dk + 1) - 1;
          arma::cube proj_hess_kde_djk = hess_kde(arma::span::all,
                                                  arma::span(ini_dj, end_dj),
                                                  arma::span(ini_dk, end_dk));

          // Distinguish diagonal and off-diagonal cases
          if (dj == dk) {

            // Loop on evaluation points
            for (arma::uword i = 0; i < nx; i++) {

              arma::rowvec x_i_dj = x.row(i).subvec(ini_dj, end_dj);
              arma::rowvec grad_kde_i_dj =
                grad_kde.row(i).subvec(ini_dj, end_dj);
              arma::mat proj_hess_kde_djk_i = proj_hess_kde_djk.row(i);
              arma::mat I_proj_dj = Idj - x_i_dj.t() * x_i_dj;

              proj_hess_kde_djk.row(i) =
                I_proj_dj * proj_hess_kde_djk_i * I_proj_dj;
              //proj_hess_kde_djk.row(i) -= s(I_proj_dj *
              //  (Idj * arma::as_scalar(grad_kde_i_dj * x_i_dj.t()) +
              //  2 * grad_kde_i_dj.t() * x_i_dj), false);
              proj_hess_kde_djk.row(i) -= s(x_i_dj.t() * grad_kde_i_dj, true) +
                (Idj - 3 * x_i_dj.t() * x_i_dj) *
                arma::as_scalar(grad_kde_i_dj * x_i_dj.t());

            }

          } else {

            arma::mat Idk = arma::eye(d(dk) + 1, d(dk) + 1);

            // Loop on evaluation points
            for (arma::uword i = 0; i < nx; i++) {

              arma::rowvec x_i_dj = x.row(i).subvec(ini_dj, end_dj);
              arma::rowvec x_i_dk = x.row(i).subvec(ini_dk, end_dk);
              arma::mat proj_hess_kde_djk_i = proj_hess_kde_djk.row(i);

              // Caution, the order (j, k) is the opposite to (k, j) in the math
              // of the Hippocampus chapter
              proj_hess_kde_djk.row(i) = (Idj - x_i_dj.t() * x_i_dj) *
                proj_hess_kde_djk_i * (Idk - x_i_dk.t() * x_i_dk);

            }

          }

          // Save projected Hessian
          proj_hess_kde(arma::span::all,
                        arma::span(ini_dj, end_dj),
                        arma::span(ini_dk, end_dk)) = proj_hess_kde_djk;

        }

      }

    }

    // Overwrite unprojected gradient and Hessian
    grad_kde = proj_grad_kde;
    hess_kde = proj_hess_kde;

  }

  // Return log-kde or kde
  if (!log) {

    log_kde = arma::exp(log_kde);

  }

  // Return a Rcpp list (NEW)
  return Rcpp::List::create(Rcpp::Named("dens") = log_kde,
                            Rcpp::Named("grad") = grad_kde,
                            Rcpp::Named("hess") = hess_kde);

}


//' @title Projected gradient of the polyspherical kde
//'
//' @inheritParams kde_polysph
//' @inheritParams grad_hess_kde_polysph
//' @param fix_u1 ensure the \eqn{u_1} vector is different from \eqn{x}?
//' Prevents the Euler algorithm to "surf the ridge". Defaults to \code{TRUE}.
//' @param sparse use a sparse eigendecomposition of the Hessian? Defaults to
//' \code{FALSE}.
//' @export
// [[Rcpp::export]]
Rcpp::List proj_grad_kde_polysph(arma::mat x, arma::mat X, arma::uvec d,
                                 arma::vec h, Rcpp::NumericVector weights =
                                   Rcpp::NumericVector::create(),
                                 bool wrt_unif = false, bool normalized = true,
                                 bool norm_x = false, bool norm_X = false,
                                 arma::uword kernel = 1,
                                 arma::uword kernel_type = 1,
                                 double k = 10.0, bool proj_alt = true,
                                 bool fix_u1 = true, bool sparse = false) {

  // How many S^dj's
  arma::uword r = d.n_elem;

  // Indexes with begin and end of each S^dj
  arma::uvec ind_dj = arma::conv_to<arma::uvec>::from(arma::zeros(r + 1));
  ind_dj.tail(r) = d + 1;
  ind_dj = arma::cumsum(ind_dj);

  // Compute normalized gradient and Hessian
  Rcpp::List grad_hess = grad_hess_kde_polysph(x, X, d, h, weights, true,
                                               proj_alt, true, false, wrt_unif,
                                               normalized, norm_x, norm_X,
                                               kernel, kernel_type, k);
  arma::mat grad = grad_hess["grad"];
  arma::cube hess = grad_hess["hess"];

  // Diagonal matrix
  arma::uword p = X.n_cols;
  arma::mat I = arma::eye(p, p);

  // Eigenvalues and eigenvectors
  arma::uword n_lamb = sparse ? (r + r/2) : p;
  n_lamb = std::min(n_lamb, p);
  arma::uword nx = x.n_rows;
  arma::mat lambda = arma::zeros(nx, n_lamb);
  arma::mat u1 = arma::zeros(nx, p);

  // Loop on evaluation points
  arma::mat eta = arma::zeros(nx, p);
  arma::vec eigval = arma::zeros(n_lamb);
  arma::mat eigvec = arma::zeros(p, n_lamb);
  for (arma::uword ki = 0; ki < nx; ki++) {

    // Sparse or dense eigen decomposition
    arma::mat hess_norm = hess.row(ki);
    if (sparse) {

      // Call to arma::eigs_sym() for computing the r + 2 eigenvalues with
      // largest algebraic value ("la") and its associated eigenvectors. We
      // need more than r eigenvalues since it may be that the first r are
      // numerical zeros. The returned eigenvalues are stored in ascending
      // order.
      arma::sp_mat H = arma::sp_mat(hess_norm);
      bool eig = arma::eigs_sym(eigval, eigvec, H, n_lamb, "la");

      // Check errors
      if (!eig) {

        Rcpp::stop("Sparse eigendecomposition eigs_sym() failed.");

      }

      // Check that there are at least r + 1 eigenvalues. It might be that less
      // than the required eigenvalues are returned!
      if (eigval.n_elem < (r + 1)) {

        Rcpp::stop("Sparse eigendecomposition eigs_sym() returned only %d eigenvectors (required: %d)",
                   eigval.n_elem, r + 1);

      }

      // Filter the null
      arma::uvec where_nonzero_eigvals = arma::find(arma::abs(eigval) > 1e-10);

      // Save eigenvalues
      lambda.row(ki) = arma::reverse(eigval).t();

      // Ensure the first eigenvector does not have null eigenvalue. See
      // reasons below. Save the first eigenvector with non-zero eigenvalue
      // as u1.
      if (fix_u1) {

        arma::uword where_u1 = arma::as_scalar(where_nonzero_eigvals.tail(1));
        u1.row(ki) = eigvec.col(where_u1).t();

      } else {

        u1.row(ki) = eigvec.tail_cols(1).t();

      }

    } else {

      // Decomposition -- eigenvalues are in ascending (algebraic) order
      // The divide-and-conquer "dc" is considerably faster for large matrices
      bool eig = arma::eig_sym(eigval, eigvec, hess_norm, "dc");

      // Check errors
      if (!eig) {

        Rcpp::stop("Dense eigendecomposition eig_sym() failed.");

      }

      // Check that there are r null eigenvalues. It is mathematically impossible
      // that there are less than r null eigenvalues, but it may happen due to
      // a loss of accuracy with small bandwidths.
      arma::uvec where_nonzero_eigvals = arma::find(arma::abs(eigval) > 1e-10);
      arma::uword n_zero_eigvals = p - where_nonzero_eigvals.n_elem;
      if (n_zero_eigvals != r) {

        Rcpp::warning("%d null eigenvalues in projected Hessian (expected: %d): accuracy loss due to small bandwidth?",
                      n_zero_eigvals, r);

      }

      // Save eigenvalues in reverse order (largest to smallest)
      lambda.row(ki) = arma::reverse(eigval).t();

      // Ensure the first eigenvector does not have null eigenvalue. That would
      // correspond to x if r = 1 and it would imply projecting *onto* the
      // direction of maximal curvature u1, causing a displacement along the
      // ridge. It might happen if the other eigenvalues are negative or closer
      // to zero for some reason. So save the first eigenvector with non-zero
      // eigenvalue as u1.
      if (fix_u1) {

        arma::uword where_u1 = arma::as_scalar(where_nonzero_eigvals.tail(1));
        u1.row(ki) = eigvec.col(where_u1).t();

      } else {

        u1.row(ki) = eigvec.col(n_lamb - 1).t();

      }

    }

    // Projection -- both approaches are actually the same since grad'x = 0!
    // if (proj_tang) {
    //
    //   arma::mat IP = AP(x.row(ki), x.row(ki), ind_dj, false)["P"];
    //   eta.row(ki) = grad.row(ki) * (IP - u1.row(ki).t() * u1.row(ki));
    //
    // } else {
    //
    //   eta.row(ki) = grad.row(ki) * (I - u1.row(ki).t() * u1.row(ki));
    //
    // }
    eta.row(ki) = grad.row(ki) * (I - u1.row(ki).t() * u1.row(ki));

  }

  // Return a Rcpp list
  return Rcpp::List::create(Rcpp::Named("eta") = eta,
                            Rcpp::Named("u1") = u1,
                            Rcpp::Named("lamb_norm") = lambda);

}
