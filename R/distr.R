

#' @title TODO
#'
#' @description TODO
#'
#' @param x TODO
#' @param d TODO
#' @param mu TODO
#' @param kappa TODO
#' @param log TODO
#' @return TODO
#' @examples
#' # TODO
#' @export
d_vmf_polysph <- function(x, d,mu, kappa, log = FALSE) {

  kde_polysph(x = x, X = rbind(mu), d = d, h = 1 / sqrt(kappa),
              kernel = 1, wrt_unif = FALSE, norm_x = TRUE, norm_X = TRUE,
              log = log)

}


#' @title TODO
#'
#' @description TODO
#'
#' @inheritParams d_vmf_polysph
#' @return TODO
#' @examples
#' # TODO
#' @export
d_unif_polysph <- function(x, d, log = FALSE) {

  if (is.null(dim(x))) {
    x <- rbind(x)
  }
  stopifnot(ncol(x) == sum(d + 1))
  log_dens <- rep(-sum(rotasym::w_p(p = d + 1, log = TRUE)), nrow(x))
  if (!log) {

    log_dens <- exp(log_dens)

  }
  return(log_dens)

}
