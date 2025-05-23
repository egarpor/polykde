// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// euler_ridge
Rcpp::List euler_ridge(arma::mat x, arma::mat X, arma::uvec d, arma::vec h, Rcpp::NumericVector h_euler, Rcpp::NumericVector weights, bool wrt_unif, bool normalized, bool norm_x, bool norm_X, arma::uword kernel, arma::uword kernel_type, double k, arma::uword N, double eps, bool keep_paths, bool proj_alt, bool fix_u1, bool sparse, bool show_prog, bool show_prog_j);
RcppExport SEXP _polykde_euler_ridge(SEXP xSEXP, SEXP XSEXP, SEXP dSEXP, SEXP hSEXP, SEXP h_eulerSEXP, SEXP weightsSEXP, SEXP wrt_unifSEXP, SEXP normalizedSEXP, SEXP norm_xSEXP, SEXP norm_XSEXP, SEXP kernelSEXP, SEXP kernel_typeSEXP, SEXP kSEXP, SEXP NSEXP, SEXP epsSEXP, SEXP keep_pathsSEXP, SEXP proj_altSEXP, SEXP fix_u1SEXP, SEXP sparseSEXP, SEXP show_progSEXP, SEXP show_prog_jSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type x(xSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type X(XSEXP);
    Rcpp::traits::input_parameter< arma::uvec >::type d(dSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type h(hSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type h_euler(h_eulerSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type weights(weightsSEXP);
    Rcpp::traits::input_parameter< bool >::type wrt_unif(wrt_unifSEXP);
    Rcpp::traits::input_parameter< bool >::type normalized(normalizedSEXP);
    Rcpp::traits::input_parameter< bool >::type norm_x(norm_xSEXP);
    Rcpp::traits::input_parameter< bool >::type norm_X(norm_XSEXP);
    Rcpp::traits::input_parameter< arma::uword >::type kernel(kernelSEXP);
    Rcpp::traits::input_parameter< arma::uword >::type kernel_type(kernel_typeSEXP);
    Rcpp::traits::input_parameter< double >::type k(kSEXP);
    Rcpp::traits::input_parameter< arma::uword >::type N(NSEXP);
    Rcpp::traits::input_parameter< double >::type eps(epsSEXP);
    Rcpp::traits::input_parameter< bool >::type keep_paths(keep_pathsSEXP);
    Rcpp::traits::input_parameter< bool >::type proj_alt(proj_altSEXP);
    Rcpp::traits::input_parameter< bool >::type fix_u1(fix_u1SEXP);
    Rcpp::traits::input_parameter< bool >::type sparse(sparseSEXP);
    Rcpp::traits::input_parameter< bool >::type show_prog(show_progSEXP);
    Rcpp::traits::input_parameter< bool >::type show_prog_j(show_prog_jSEXP);
    rcpp_result_gen = Rcpp::wrap(euler_ridge(x, X, d, h, h_euler, weights, wrt_unif, normalized, norm_x, norm_X, kernel, kernel_type, k, N, eps, keep_paths, proj_alt, fix_u1, sparse, show_prog, show_prog_j));
    return rcpp_result_gen;
END_RCPP
}
// grad_hess_kde_polysph
Rcpp::List grad_hess_kde_polysph(arma::mat x, arma::mat X, arma::uvec d, arma::vec h, Rcpp::NumericVector weights, bool projected, bool proj_alt, bool norm_grad_hess, bool log, bool wrt_unif, bool normalized, bool norm_x, bool norm_X, arma::uword kernel, arma::uword kernel_type, double k);
RcppExport SEXP _polykde_grad_hess_kde_polysph(SEXP xSEXP, SEXP XSEXP, SEXP dSEXP, SEXP hSEXP, SEXP weightsSEXP, SEXP projectedSEXP, SEXP proj_altSEXP, SEXP norm_grad_hessSEXP, SEXP logSEXP, SEXP wrt_unifSEXP, SEXP normalizedSEXP, SEXP norm_xSEXP, SEXP norm_XSEXP, SEXP kernelSEXP, SEXP kernel_typeSEXP, SEXP kSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type x(xSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type X(XSEXP);
    Rcpp::traits::input_parameter< arma::uvec >::type d(dSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type h(hSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type weights(weightsSEXP);
    Rcpp::traits::input_parameter< bool >::type projected(projectedSEXP);
    Rcpp::traits::input_parameter< bool >::type proj_alt(proj_altSEXP);
    Rcpp::traits::input_parameter< bool >::type norm_grad_hess(norm_grad_hessSEXP);
    Rcpp::traits::input_parameter< bool >::type log(logSEXP);
    Rcpp::traits::input_parameter< bool >::type wrt_unif(wrt_unifSEXP);
    Rcpp::traits::input_parameter< bool >::type normalized(normalizedSEXP);
    Rcpp::traits::input_parameter< bool >::type norm_x(norm_xSEXP);
    Rcpp::traits::input_parameter< bool >::type norm_X(norm_XSEXP);
    Rcpp::traits::input_parameter< arma::uword >::type kernel(kernelSEXP);
    Rcpp::traits::input_parameter< arma::uword >::type kernel_type(kernel_typeSEXP);
    Rcpp::traits::input_parameter< double >::type k(kSEXP);
    rcpp_result_gen = Rcpp::wrap(grad_hess_kde_polysph(x, X, d, h, weights, projected, proj_alt, norm_grad_hess, log, wrt_unif, normalized, norm_x, norm_X, kernel, kernel_type, k));
    return rcpp_result_gen;
END_RCPP
}
// proj_grad_kde_polysph
Rcpp::List proj_grad_kde_polysph(arma::mat x, arma::mat X, arma::uvec d, arma::vec h, Rcpp::NumericVector weights, bool wrt_unif, bool normalized, bool norm_x, bool norm_X, arma::uword kernel, arma::uword kernel_type, double k, bool proj_alt, bool fix_u1, bool sparse);
RcppExport SEXP _polykde_proj_grad_kde_polysph(SEXP xSEXP, SEXP XSEXP, SEXP dSEXP, SEXP hSEXP, SEXP weightsSEXP, SEXP wrt_unifSEXP, SEXP normalizedSEXP, SEXP norm_xSEXP, SEXP norm_XSEXP, SEXP kernelSEXP, SEXP kernel_typeSEXP, SEXP kSEXP, SEXP proj_altSEXP, SEXP fix_u1SEXP, SEXP sparseSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type x(xSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type X(XSEXP);
    Rcpp::traits::input_parameter< arma::uvec >::type d(dSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type h(hSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type weights(weightsSEXP);
    Rcpp::traits::input_parameter< bool >::type wrt_unif(wrt_unifSEXP);
    Rcpp::traits::input_parameter< bool >::type normalized(normalizedSEXP);
    Rcpp::traits::input_parameter< bool >::type norm_x(norm_xSEXP);
    Rcpp::traits::input_parameter< bool >::type norm_X(norm_XSEXP);
    Rcpp::traits::input_parameter< arma::uword >::type kernel(kernelSEXP);
    Rcpp::traits::input_parameter< arma::uword >::type kernel_type(kernel_typeSEXP);
    Rcpp::traits::input_parameter< double >::type k(kSEXP);
    Rcpp::traits::input_parameter< bool >::type proj_alt(proj_altSEXP);
    Rcpp::traits::input_parameter< bool >::type fix_u1(fix_u1SEXP);
    Rcpp::traits::input_parameter< bool >::type sparse(sparseSEXP);
    rcpp_result_gen = Rcpp::wrap(proj_grad_kde_polysph(x, X, d, h, weights, wrt_unif, normalized, norm_x, norm_X, kernel, kernel_type, k, proj_alt, fix_u1, sparse));
    return rcpp_result_gen;
END_RCPP
}
// kde_polysph
arma::vec kde_polysph(arma::mat x, arma::mat X, arma::uvec d, arma::vec h, Rcpp::NumericVector weights, bool log, bool wrt_unif, bool normalized, bool intrinsic, bool norm_x, bool norm_X, arma::uword kernel, arma::uword kernel_type, double k);
RcppExport SEXP _polykde_kde_polysph(SEXP xSEXP, SEXP XSEXP, SEXP dSEXP, SEXP hSEXP, SEXP weightsSEXP, SEXP logSEXP, SEXP wrt_unifSEXP, SEXP normalizedSEXP, SEXP intrinsicSEXP, SEXP norm_xSEXP, SEXP norm_XSEXP, SEXP kernelSEXP, SEXP kernel_typeSEXP, SEXP kSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type x(xSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type X(XSEXP);
    Rcpp::traits::input_parameter< arma::uvec >::type d(dSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type h(hSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type weights(weightsSEXP);
    Rcpp::traits::input_parameter< bool >::type log(logSEXP);
    Rcpp::traits::input_parameter< bool >::type wrt_unif(wrt_unifSEXP);
    Rcpp::traits::input_parameter< bool >::type normalized(normalizedSEXP);
    Rcpp::traits::input_parameter< bool >::type intrinsic(intrinsicSEXP);
    Rcpp::traits::input_parameter< bool >::type norm_x(norm_xSEXP);
    Rcpp::traits::input_parameter< bool >::type norm_X(norm_XSEXP);
    Rcpp::traits::input_parameter< arma::uword >::type kernel(kernelSEXP);
    Rcpp::traits::input_parameter< arma::uword >::type kernel_type(kernel_typeSEXP);
    Rcpp::traits::input_parameter< double >::type k(kSEXP);
    rcpp_result_gen = Rcpp::wrap(kde_polysph(x, X, d, h, weights, log, wrt_unif, normalized, intrinsic, norm_x, norm_X, kernel, kernel_type, k));
    return rcpp_result_gen;
END_RCPP
}
// log_cv_kde_polysph
arma::vec log_cv_kde_polysph(arma::mat X, arma::uvec d, arma::vec h, Rcpp::NumericVector weights, bool wrt_unif, bool normalized, bool intrinsic, bool norm_X, arma::uword kernel, arma::uword kernel_type, double k);
RcppExport SEXP _polykde_log_cv_kde_polysph(SEXP XSEXP, SEXP dSEXP, SEXP hSEXP, SEXP weightsSEXP, SEXP wrt_unifSEXP, SEXP normalizedSEXP, SEXP intrinsicSEXP, SEXP norm_XSEXP, SEXP kernelSEXP, SEXP kernel_typeSEXP, SEXP kSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type X(XSEXP);
    Rcpp::traits::input_parameter< arma::uvec >::type d(dSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type h(hSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type weights(weightsSEXP);
    Rcpp::traits::input_parameter< bool >::type wrt_unif(wrt_unifSEXP);
    Rcpp::traits::input_parameter< bool >::type normalized(normalizedSEXP);
    Rcpp::traits::input_parameter< bool >::type intrinsic(intrinsicSEXP);
    Rcpp::traits::input_parameter< bool >::type norm_X(norm_XSEXP);
    Rcpp::traits::input_parameter< arma::uword >::type kernel(kernelSEXP);
    Rcpp::traits::input_parameter< arma::uword >::type kernel_type(kernel_typeSEXP);
    Rcpp::traits::input_parameter< double >::type k(kSEXP);
    rcpp_result_gen = Rcpp::wrap(log_cv_kde_polysph(X, d, h, weights, wrt_unif, normalized, intrinsic, norm_X, kernel, kernel_type, k));
    return rcpp_result_gen;
END_RCPP
}
// sfp
arma::mat sfp(arma::mat t);
RcppExport SEXP _polykde_sfp(SEXP tSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type t(tSEXP);
    rcpp_result_gen = Rcpp::wrap(sfp(t));
    return rcpp_result_gen;
END_RCPP
}
// proj_polysph
arma::mat proj_polysph(arma::mat x, arma::uvec ind_dj);
RcppExport SEXP _polykde_proj_polysph(SEXP xSEXP, SEXP ind_djSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type x(xSEXP);
    Rcpp::traits::input_parameter< arma::uvec >::type ind_dj(ind_djSEXP);
    rcpp_result_gen = Rcpp::wrap(proj_polysph(x, ind_dj));
    return rcpp_result_gen;
END_RCPP
}
// dist_polysph
arma::vec dist_polysph(arma::mat x, arma::mat y, arma::uvec ind_dj, bool norm_x, bool norm_y, bool std);
RcppExport SEXP _polykde_dist_polysph(SEXP xSEXP, SEXP ySEXP, SEXP ind_djSEXP, SEXP norm_xSEXP, SEXP norm_ySEXP, SEXP stdSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type x(xSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type y(ySEXP);
    Rcpp::traits::input_parameter< arma::uvec >::type ind_dj(ind_djSEXP);
    Rcpp::traits::input_parameter< bool >::type norm_x(norm_xSEXP);
    Rcpp::traits::input_parameter< bool >::type norm_y(norm_ySEXP);
    Rcpp::traits::input_parameter< bool >::type std(stdSEXP);
    rcpp_result_gen = Rcpp::wrap(dist_polysph(x, y, ind_dj, norm_x, norm_y, std));
    return rcpp_result_gen;
END_RCPP
}
// dist_polysph_cross
arma::mat dist_polysph_cross(arma::mat x, arma::mat y, arma::uvec ind_dj, bool norm_x, bool norm_y, bool std);
RcppExport SEXP _polykde_dist_polysph_cross(SEXP xSEXP, SEXP ySEXP, SEXP ind_djSEXP, SEXP norm_xSEXP, SEXP norm_ySEXP, SEXP stdSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type x(xSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type y(ySEXP);
    Rcpp::traits::input_parameter< arma::uvec >::type ind_dj(ind_djSEXP);
    Rcpp::traits::input_parameter< bool >::type norm_x(norm_xSEXP);
    Rcpp::traits::input_parameter< bool >::type norm_y(norm_ySEXP);
    Rcpp::traits::input_parameter< bool >::type std(stdSEXP);
    rcpp_result_gen = Rcpp::wrap(dist_polysph_cross(x, y, ind_dj, norm_x, norm_y, std));
    return rcpp_result_gen;
END_RCPP
}
// diamond_crossprod
arma::cube diamond_crossprod(arma::mat X, arma::uvec ind_dj);
RcppExport SEXP _polykde_diamond_crossprod(SEXP XSEXP, SEXP ind_djSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type X(XSEXP);
    Rcpp::traits::input_parameter< arma::uvec >::type ind_dj(ind_djSEXP);
    rcpp_result_gen = Rcpp::wrap(diamond_crossprod(X, ind_dj));
    return rcpp_result_gen;
END_RCPP
}
// diamond_rcrossprod
arma::cube diamond_rcrossprod(arma::mat X, arma::uvec ind_dj);
RcppExport SEXP _polykde_diamond_rcrossprod(SEXP XSEXP, SEXP ind_djSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type X(XSEXP);
    Rcpp::traits::input_parameter< arma::uvec >::type ind_dj(ind_djSEXP);
    rcpp_result_gen = Rcpp::wrap(diamond_rcrossprod(X, ind_dj));
    return rcpp_result_gen;
END_RCPP
}
// s
arma::mat s(arma::mat A, bool add);
RcppExport SEXP _polykde_s(SEXP ASEXP, SEXP addSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type A(ASEXP);
    Rcpp::traits::input_parameter< bool >::type add(addSEXP);
    rcpp_result_gen = Rcpp::wrap(s(A, add));
    return rcpp_result_gen;
END_RCPP
}
// AP
Rcpp::List AP(arma::rowvec x, arma::rowvec v, arma::uvec ind_dj, bool orth);
RcppExport SEXP _polykde_AP(SEXP xSEXP, SEXP vSEXP, SEXP ind_djSEXP, SEXP orthSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::rowvec >::type x(xSEXP);
    Rcpp::traits::input_parameter< arma::rowvec >::type v(vSEXP);
    Rcpp::traits::input_parameter< arma::uvec >::type ind_dj(ind_djSEXP);
    Rcpp::traits::input_parameter< bool >::type orth(orthSEXP);
    rcpp_result_gen = Rcpp::wrap(AP(x, v, ind_dj, orth));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_polykde_euler_ridge", (DL_FUNC) &_polykde_euler_ridge, 21},
    {"_polykde_grad_hess_kde_polysph", (DL_FUNC) &_polykde_grad_hess_kde_polysph, 16},
    {"_polykde_proj_grad_kde_polysph", (DL_FUNC) &_polykde_proj_grad_kde_polysph, 15},
    {"_polykde_kde_polysph", (DL_FUNC) &_polykde_kde_polysph, 14},
    {"_polykde_log_cv_kde_polysph", (DL_FUNC) &_polykde_log_cv_kde_polysph, 11},
    {"_polykde_sfp", (DL_FUNC) &_polykde_sfp, 1},
    {"_polykde_proj_polysph", (DL_FUNC) &_polykde_proj_polysph, 2},
    {"_polykde_dist_polysph", (DL_FUNC) &_polykde_dist_polysph, 6},
    {"_polykde_dist_polysph_cross", (DL_FUNC) &_polykde_dist_polysph_cross, 6},
    {"_polykde_diamond_crossprod", (DL_FUNC) &_polykde_diamond_crossprod, 2},
    {"_polykde_diamond_rcrossprod", (DL_FUNC) &_polykde_diamond_rcrossprod, 2},
    {"_polykde_s", (DL_FUNC) &_polykde_s, 2},
    {"_polykde_AP", (DL_FUNC) &_polykde_AP, 4},
    {NULL, NULL, 0}
};

RcppExport void R_init_polykde(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
