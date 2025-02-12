
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
//' @description Computes the softplus function \eqn{\log(1+e^{t})} in a
//' numerically stable way for large absolute values of \eqn{t}.
//'
//' @inheritParams softplus
//' @return The softplus function evaluated at \code{t}.
//' @examples
//' curve(log(polykde:::sfp(rbind(5 * (1 - x)))), from = -10, to = 10)
//' @keywords internal
// [[Rcpp::export]]
arma::mat sfp(arma::mat t) {

  // return t + arma::log1p(arma::exp(-t));
  // The function is going to be evaluated for negative t's
  return arma::log1p(arma::exp(t));

}


//' @title Projection onto the polysphere
//'
//' @description Projects points on \eqn{\mathbb{R}^{d_1 + \cdots + d_r + r}}
//' onto the polysphere \eqn{\mathcal{S}^{d_1} \times \cdots \times
//' \mathcal{S}^{d_r}} by normalizing each block of \eqn{d_j} coordinates.
//'
//' @param x a matrix of size \code{c(n, sum(d) + r)}.
//' @param ind_dj \code{0}-based index separating the blocks of spheres that
//' is computed with \code{\link{comp_ind_dj}}.
//' @return A matrix of size \code{c(n, sum(d) + r)} with the projected points.
//' @examples
//' # Example on (S^1)^2
//' d <- c(1, 1)
//' x <- rbind(c(2, 0, 1, 1))
//' polykde:::proj_polysph(x, ind_dj = comp_ind_dj(d))
//' @keywords internal
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
//' @description Computation of the distance between points \eqn{\boldsymbol{x}}
//' and \eqn{\boldsymbol{y}} on the polysphere
//' \eqn{\mathcal{S}^{d_1} \times \cdots \times \mathcal{S}^{d_r}}:
//' \deqn{\sqrt{\sum_{j=1}^r
//' d_{\mathcal{S}^{d_j}}(\boldsymbol{x}_j, \boldsymbol{y}_j)^2},}
//' where \eqn{d_{\mathcal{S}^{d_j}}(\boldsymbol{x}_j, \boldsymbol{y}_j)=
//' \cos^{-1}(\boldsymbol{x}_j' \boldsymbol{y}_j)}.
//'
//' @inheritParams proj_polysph
//' @param y either a matrix of size \code{c(m, sum(d) + r)} or a vector of
//' length \code{sum(d) + r}.
//' @inheritParams proj_polysph
//' @param norm_x,norm_y ensure a normalization of the data? Default to
//' \code{FALSE}.
//' @param std standardize distance to \eqn{[0,1]}? Uses that the maximum
//' distance is \eqn{\sqrt{r}\pi}. Defaults to \code{TRUE}.
//' @return
//' \itemize{
//' \item{\code{dist_polysph}: a vector of size \code{n} with the distances
//' between \code{x} and \code{y}.}
//' \item{\code{dist_polysph_matrix}: a matrix of size \code{c(n, n)} with the
//' pairwise distances of \code{x}.}
//' \item{\code{dist_polysph_cross}: a matrix of distances of size
//' \code{c(n, m)} with the cross distances between \code{x} and \code{y}.}
//' }
//' @examples
//' # Example on S^2 x S^3 x S^1
//' d <- c(2, 3, 1)
//' ind_dj <- comp_ind_dj(d)
//' n <- 3
//' x <- r_unif_polysph(n = n, d = d)
//' y <- r_unif_polysph(n = n, d = d)
//'
//' # Distances of x to y
//' dist_polysph(x = x, y = y, ind_dj = ind_dj, std = FALSE)
//' dist_polysph(x = x, y = y[1, , drop = FALSE], ind_dj = ind_dj, std = FALSE)
//'
//' # Pairwise distance matrix of x
//' dist_polysph_matrix(x = x, ind_dj = ind_dj, std = FALSE)
//'
//' # Cross distances between x and y
//' dist_polysph_cross(x = x, y = y, ind_dj = ind_dj, std = FALSE)
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

  // One + epsilon
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


//' @rdname dist_polysph
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


//' @title Diamond cross-product
//'
//' @description Given a matrix \eqn{\boldsymbol{X}} whose \eqn{n} rows are on a
//' polysphere \eqn{\mathcal{S}^{d_1} \times \cdots \times \mathcal{S}^{d_r}},
//' the function computes the cube whose rows are
//' \eqn{\boldsymbol{X}_i \diamond \boldsymbol{X}_i'}, \eqn{i = 1, \ldots, n},
//' and \eqn{\diamond} is a block-by-block product.
//'
//' @inheritParams kde_polysph
//' @inheritParams proj_polysph
//' @return An array of size \code{c(nrow(X), ncol(X), ncol(X))}.
//' @examples
//' d <- c(1, 2)
//' X <- r_unif_polysph(n = 2, d = d)
//' polykde:::diamond_crossprod(X = X, ind_dj = comp_ind_dj(d))
//' @keywords internal
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


//' @title Symmetrize a matrix
//'
//' @description Symmetrizes a matrix \eqn{\boldsymbol{A}} by returning
//' \eqn{(\boldsymbol{A} + \boldsymbol{A}') / 2}.
//'
//' @param A a matrix.
//' @param add return simply the addition
//' \eqn{\boldsymbol{A} + \boldsymbol{A}'}? Defaults to \code{FALSE}
//' @return A symmetric matrix with the same dimensions as \code{A}.
//' @examples
//' A <- matrix(rnorm(4), nrow = 2, ncol = 2)
//' polykde:::s(A)
//' @keywords internal
// [[Rcpp::export]]
arma::mat s(arma::mat A, bool add = false) {

  if (add) {

    return A + A.t();

  } else {

    return 0.5 * (A + A.t());

  }

}


//' @title Projection matrices \eqn{\boldsymbol{P}} and \eqn{\boldsymbol{A}}
//'
//' @description Computation of the projection matrices \eqn{\boldsymbol{P}}
//' and \eqn{\boldsymbol{A}}. The \eqn{jj}-block of \eqn{\boldsymbol{P}} is
//' \eqn{\boldsymbol{I}_{d_j} - \boldsymbol{x}_j \boldsymbol{x}_j'}. The
//' \eqn{jj}-block of \eqn{\boldsymbol{A}} is \eqn{(\boldsymbol{x}_j'
//' \boldsymbol{v}_j) \boldsymbol{I}_{d_j}}, \eqn{j=1,\ldots,r}.
//'
//' @param x,v row vectors of size \code{sum(d) + r}.
//' @inheritParams proj_polysph
//' @param orth return the orthogonal complement of \eqn{\boldsymbol{P}},
//' \eqn{\boldsymbol{I} - \boldsymbol{P}}? Defaults to \code{FALSE}.
//' @return A list with the matrices \eqn{\boldsymbol{P}} and
//' \eqn{\boldsymbol{A}}. Both matrices have size
//' \code{c(sum(d) + r, sum(d) + r)}.
//' @examples
//' d <- c(1, 2)
//' x <- r_unif_polysph(n = 1, d = d)
//' v <- r_unif_polysph(n = 1, d = d)
//' polykde:::AP(x = x, v = v, ind_dj = comp_ind_dj(d))
//' @keywords internal
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
