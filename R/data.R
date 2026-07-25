

#' @title Skeletal representation of hippocampus shapes
#'
#' @description Contains skeletal representations (s-reps) for the hippocampus
#' shapes of \eqn{177} 6-month-old infants. The s-reps models include \eqn{168}
#' spokes, and hence the directions of the spokes form a sample on the
#' polysphere \eqn{(\mathbb{S}^2)^{168}}. There are \eqn{34} infants that later
#' developed autism.
#'
#' @docType data
#' @format A list with six fields:
#' \describe{
#'   \item{base}{an array of size \code{c(177, 168, 3)} with the coordinates of
#'   the base points.}
#'   \item{bdry}{an array of size \code{c(177, 168, 3)} with the coordinates of
#'   the boundary points.}
#'   \item{dirs}{an array of size \code{c(177, 168, 3)} with the directions of
#'   spokes (unit vectors).}
#'   \item{rads}{a matrix of size \code{c(177, 168)} with the radii of the
#'   spokes.}
#'   \item{ids}{hippocampi identifiers.}
#'   \item{ids_labs}{class labels. \code{TRUE} stands for infants that developed
#'    autism and \code{FALSE} for those that did not.}
#' }
#' @source The fitted s-reps data was shared by Prof. Stephen M. Pizer and Dr.
#' Zhiyuan Liu, see Liu et al. (2021) and the repository
#' \url{https://github.com/ZhiyLiu/shanapy}. The source for the raw hippocampi
#' data is the Infant Brain Imaging Study (IBIS) network.
#' @references
#' García-Portugués, E. and Meilán-Vila, A. (2023). Hippocampus shape analysis
#' via skeletal models and kernel smoothing. In Larriba, Y. (Ed.),
#' \emph{Statistical Methods at the Forefront of Biomedical Advances},
#' pp. 63--82. Springer, Cham. \doi{10.1007/978-3-031-32729-2_4}.
#'
#' Liu, Z., Hong, J., Vicory, J., Damon, J. N. and Pizer, S. M. (2021). Fitting
#' unbranching skeletal structures to objects. \emph{Medical Image Analysis},
#' 70:102020. \doi{10.1016/j.media.2021.102020}
#' @examples
#' \donttest{
#' if (requireNamespace("scatterplot3d", quietly = TRUE)) {
#'   # Load data
#'   data("hippocampus")
#'   attach(hippocampus)
#'
#'   # View ith hippocampus s-rep
#'   i <- 50
#'   view_srep(base = base[i, , ], dirs = dirs[i, , ], bdry = bdry[i, , ])
#' }
#' }
"hippocampus"
