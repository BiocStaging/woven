# AllClasses.R -- S4 class definitions
#
# The S4 "woven" class replaces the previous plain S3 list, so the object
# returned by woven() participates in R's formal class system (validity
# checking, show/summary/plot as real S4 methods) rather than relying on an
# informally-attached class attribute.

#' S4 class for a fitted WOVEN model
#'
#' The object returned by \code{\link{woven}}. Holds the consensus latent
#' scores for every subject, the per-modality projection matrices, and the
#' hyperparameters and bookkeeping needed by \code{\link{woven_scores}},
#' \code{\link{woven_predict}}, \code{\link{woven_metrics}}, and the
#' \code{woven_plot_*} functions.
#'
#' @slot Z numeric matrix, n x K consensus latent scores for all n subjects
#'   (anchors and block-missing); rows with no observed modality are NA.
#' @slot W_list list of V numeric matrices, each p_v x K, the per-modality
#'   projection matrices.
#' @slot Z_anchors list of V numeric matrices, each n_a x K, the
#'   per-modality anchor latent scores before consensus averaging.
#' @slot singular_values numeric vector of length K, the supervised
#'   canonical correlations from the dual eigendecomposition.
#' @slot anchor_idx integer vector, indices (into the original n subjects)
#'   of the fully-observed anchor subjects.
#' @slot Y_anchor labels of the anchor subjects, in \code{anchor_idx} order.
#' @slot Y_levels sorted unique class labels seen during fitting.
#' @slot Y_labels optional original label names (e.g. \code{"CN"}), or
#'   \code{NULL} if \code{Y} was already numeric/integer.
#' @slot K integer, latent dimension.
#' @slot V integer, number of modalities.
#' @slot n integer, total number of subjects.
#' @slot mod_names optional character vector of modality names, or
#'   \code{NULL} if \code{X_list} was unnamed.
#' @slot lambdas numeric vector of length V, per-modality Laplacian
#'   regularization strength used at fit time.
#' @slot gamma_y numeric, supervision strength used at fit time.
#' @slot k_nn integer, k-NN graph size used at fit time.
#' @slot ridge_w numeric, ridge strength used at fit time.
#' @slot screen_top numeric (integer count or \code{Inf}), feature-screening
#'   depth used at fit time.
#' @slot scaled logical, whether modalities were standardized before fitting.
#' @slot scale_center list of V numeric vectors (or \code{NULL}), per-feature
#'   centers used for standardization, reapplied in \code{\link{woven_scores}}.
#' @slot scale_scale list of V numeric vectors (or \code{NULL}), per-feature
#'   scales used for standardization, reapplied in \code{\link{woven_scores}}.
#' @slot fit_mcca list, the raw output of \code{\link{woven_mcca_dual}} for
#'   advanced/diagnostic use.
#'
#' @seealso [woven()]
#' @exportClass woven
setClass("woven",
    representation(
        Z = "matrix",
        W_list = "list",
        Z_anchors = "list",
        singular_values = "numeric",
        anchor_idx = "integer",
        Y_anchor = "ANY",
        Y_levels = "ANY",
        Y_labels = "ANY",
        K = "integer",
        V = "integer",
        n = "integer",
        mod_names = "ANY",
        lambdas = "numeric",
        gamma_y = "numeric",
        k_nn = "integer",
        ridge_w = "numeric",
        screen_top = "numeric",
        scaled = "logical",
        scale_center = "ANY",
        scale_scale = "ANY",
        fit_mcca = "list"
    )
)

setValidity("woven", function(object) {
    errs <- character(0)
    if (nrow(object@Z) != object@n) {
        errs <- c(errs, "nrow(Z) must equal n")
    }
    if (ncol(object@Z) != object@K) {
        errs <- c(errs, "ncol(Z) must equal K")
    }
    if (length(object@W_list) != object@V) {
        errs <- c(errs, "length(W_list) must equal V")
    }
    if (length(object@Z_anchors) != object@V) {
        errs <- c(errs, "length(Z_anchors) must equal V")
    }
    if (length(object@singular_values) != object@K) {
        errs <- c(errs, "length(singular_values) must equal K")
    }
    if (length(object@anchor_idx) != length(object@Y_anchor)) {
        errs <- c(errs, "length(anchor_idx) must equal length(Y_anchor)")
    }
    if (length(errs) == 0L) TRUE else errs
})
