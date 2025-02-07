
// [[Rcpp::depends("RcppArmadillo")]]
#include <RcppArmadillo.h>
#include <Rmath.h>
#include <R.h>

// Declaration for functions
arma::mat proj_polysph(arma::mat x, arma::uvec ind_dj);
arma::vec dist_polysph(arma::mat x, arma::mat y, arma::uvec ind_dj,
                       bool norm_x, bool norm_y, bool std);


//' @title Stable computation of the softplus function
//'
//' @inheritParams softplus
//' @examples
//' curve(log(sfp(rbind(5 * (1 - x)))), from = -10, to = 10)
//' @export
// [[Rcpp::export]]
arma::mat sfp(arma::mat t) {

  // return t + arma::log1p(arma::exp(-t));
  // The function is going to be evaluated for negative t's
  return arma::log1p(arma::exp(t));

}

//' @title Polyspherical projection
//'
//' @inheritParams kde_polysph
//' @param ind_dj 0-based index separating the spheres. Computed using
//' \code{\link{comp_ind_dj}}.
//' @export
// [[Rcpp::export]]
arma::mat proj_polysph(arma::mat x, arma::uvec ind_dj) {

  // How many S^dj's
  arma::uword r = ind_dj.n_elem - 1;

  // Dimensions
  arma::uword p = x.n_cols;
  if (p != ind_dj(r)) {

    Rcpp::stop("Dimension of x mismatches with ind_dj.");

  }

  // Normalization process
  for (arma::uword dj = 0; dj < r; dj++) {

    arma::vec norm = arma::sqrt(arma::sum(
     arma::square(x.cols(ind_dj(dj), ind_dj(dj + 1) - 1)), 1));
    x.cols(ind_dj(dj), ind_dj(dj + 1) - 1).each_col() /= norm;

    // // SLOWER
    // arma::uword ini_dj = ind_dj(dj);
    // arma::uword end_dj = ind_dj(dj + 1) - 1;
    // x.cols(ini_dj, end_dj) = arma::normalise(x.cols(ini_dj, end_dj), 2, 1);

  }
  return x;

}


//' @title Polyspherical distance
//'
//' @param x a matrix of size \code{c(n, sum(d) + r)}.
//' @param y either a matrix of the same dimension of \code{x} or a vector of
//' length \code{sum(d) + r}.
//' @inheritParams proj_polysph
//' @param norm_x,norm_y ensure a normalization of the data?
//' @param std standardize distance to \eqn{[0,1]}? Uses that the maximum
//' distance is \eqn{\sqrt{r}\pi}. Defaults to \code{TRUE}.
//' @export
// [[Rcpp::export]]
arma::vec dist_polysph(arma::mat x, arma::mat y, arma::uvec ind_dj,
                       bool norm_x = false, bool norm_y = false,
                       bool std = true) {

  // How many S^dj's
  arma::uword r = ind_dj.n_elem - 1;

  // Repeat rows of y there is only one row
  if (y.n_rows == 1) {

    y = arma::repmat(y, x.n_rows, 1);

  }

  // Dimensions
  arma::uword p = x.n_cols;
  if (p != y.n_cols || x.n_rows != y.n_rows) {

    Rcpp::stop("Dimension of x and y mismatch.");

  }
  if (p != ind_dj(r)) {

    Rcpp::stop("Dimension of x mismatches with ind_dj.");

  }

  // Normalization?
  if (norm_x) {

    x = proj_polysph(x, ind_dj);

  }
  if (norm_y) {

    y = proj_polysph(y, ind_dj);

  }

  // Epsilon + one
  double e1 = 1.0 + 1e-15;

  // Product
  arma::mat xy = x % y;

  // Distances on each sphere
  arma::mat dist_j = arma::zeros(x.n_rows, p);
  for (arma::uword dj = 0; dj < r; dj++) {

    arma::vec prod = arma::sum(xy.cols(ind_dj(dj), ind_dj(dj + 1) - 1), 1);
    double max_abs = arma::max(arma::abs(prod));
    if (max_abs > 1.0) {

      prod.clamp(-1.0, 1.0);
      if (max_abs > e1) {

        Rcpp::warning("max_i |x_j'y_j| = %g for j = %d. Clamping.",
                      max_abs, dj + 1);

      }

    }
    dist_j.col(dj) = arma::acos(prod);

  }

  // Distance
  arma::vec dist = arma::sqrt(arma::sum(arma::square(dist_j), 1));

  // Standardize?
  if (std) {

    dist /= arma::datum::pi * std::sqrt(r);

  }
  return dist;

}


//' @title Cross polyspherical distance
//'
//' @inheritParams dist_polysph
//' @param y either a matrix of the same dimension of \code{x} or a vector of
//' length \code{sum(d) + r}.
//' @export
// [[Rcpp::export]]
arma::mat dist_polysph_cross(arma::mat x, arma::mat y, arma::uvec ind_dj,
                             bool norm_x = false, bool norm_y = false,
                             bool std = true) {

  // How many S^dj's
  arma::uword r = ind_dj.n_elem - 1;

  // Dimensions
  arma::uword p = x.n_cols;
  if (p != y.n_cols) {

    Rcpp::stop("Dimension of x and y mismatch.");

  }
  if (p != ind_dj(r)) {

    Rcpp::stop("Dimension of x mismatches with ind_dj.");

  }

  // Normalization?
  if (norm_x) {

    x = proj_polysph(x, ind_dj);

  }
  if (norm_y) {

    y = proj_polysph(y, ind_dj);

  }

  // Distances to each y
  arma::mat dist = arma::zeros(x.n_rows, y.n_rows);
  for (arma::uword i = 0; i < y.n_rows; i++) {

    dist.col(i) = dist_polysph(x, y.row(i), ind_dj, false, false, false);

  }

  // Standardize?
  if (std) {

    dist /= arma::datum::pi * std::sqrt(r);

  }
  return dist;

}


//' @title Compute cube \eqn{X_i \diamond X_i'}
//'
//' @inheritParams kde_polysph
//' @inheritParams proj_polysph
//' @export
// [[Rcpp::export]]
arma::cube diamond_crossprod(arma::mat X, arma::uvec ind_dj) {

  // How many S^dj's
  arma::uword r = ind_dj.n_elem - 1;

  // Dimensions
  arma::uword n = X.n_rows;
  arma::uword p = X.n_cols;
  if (p != ind_dj(r)) {

    Rcpp::stop("Dimension of X mismatches with ind_dj.");

  }

  // Loop on the cube rows
  arma::cube X_diamond = arma::zeros(n, p, p);
  for (arma::uword i = 0; i < n; i++) {

    // i-th observation
    arma::vec X_ki = X.row(i).t();

    // Loop on the j-component // TODO: symmetryze for faster computation
    for (arma::uword dj = 0; dj < r; dj++) {

      // j-terms
      arma::uword ini_dj = ind_dj(dj);
      arma::uword end_dj = ind_dj(dj + 1) - 1;
      arma::vec X_ki_dj = X_ki.subvec(ini_dj, end_dj);

      // Loop on the k-component
      for (arma::uword dk = 0; dk < r; dk++) {

        // k-terms
        arma::uword ini_dk = ind_dj(dk);
        arma::uword end_dk = ind_dj(dk + 1) - 1;
        arma::vec X_ki_dk = X_ki.subvec(ini_dk, end_dk);

        // Fill (dj, dk)-block
        X_diamond(arma::span(i), arma::span(ini_dj, end_dj),
                  arma::span(ini_dk, end_dk)) = X_ki_dj * X_ki_dk.t();

      }
    }
  }

  return X_diamond;

}


//' @title Symmetrize a matrix A with (A+A')/2
//'
//' @param A matrix.
//' @param add return simply the addition A + A'? Defaults to \code{FALSE}
//' @export
// [[Rcpp::export]]
arma::mat s(arma::mat A, bool add = false) {

  if (add) {

    return A + A.t();

  } else {

    return 0.5 * (A + A.t());

  }

}


//' @title Projection matrices P and A
//'
//' @description The \eqn{jj}-block of \eqn{P} is \eqn{I_{d_j} - x_j x_j'}. The
//' \eqn{jj}-block of \eqn{A} is \eqn{(x_j' * v_j) * I_{d_j}}.
//'
//' @param x,v row vectors of size \code{sum(d) + r}.
//' @inheritParams proj_polysph
//' @param orth return the orthogonal complement of \eqn{P}, \eqn{I - P}?
//' @export
// [[Rcpp::export]]
Rcpp::List AP(arma::rowvec x, arma::rowvec v, arma::uvec ind_dj,
              bool orth = false) {

  // How many S^dj's
  arma::uword r = ind_dj.n_elem - 1;

  // Dimensions
  arma::uword p = x.n_elem;
  if (p != ind_dj(r)) {

    Rcpp::stop("Dimension of x mismatches with ind_dj.");

  }

  // Fill diagonal blocks of the matrices
  arma::mat P = arma::zeros(p, p);
  arma::mat A = P;
  for (arma::uword dj = 0; dj < r; dj++) {

    arma::uword ini_dj = ind_dj(dj);
    arma::uword end_dj = ind_dj(dj + 1) - 1;
    arma::rowvec xj = x.subvec(ini_dj, end_dj);
    arma::rowvec vj = v.subvec(ini_dj, end_dj);

    P(arma::span(ini_dj, end_dj), arma::span(ini_dj, end_dj)) = xj.t() * xj;
    A(arma::span(ini_dj, end_dj), arma::span(ini_dj, end_dj)).diag().fill(
        arma::as_scalar(xj * vj.t()));

  }

  // Orthogonal complement?
  if (!orth) {

    P *= -1.0;
    P.diag() += 1;

  }

  // Return a Rcpp list
  return Rcpp::List::create(Rcpp::Named("P") = P, Rcpp::Named("A") = A);

}

