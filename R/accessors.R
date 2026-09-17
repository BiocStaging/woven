# accessors.R -- S4 accessor generics for the "woven" class
#
# Bioconductor convention favors named accessor functions over direct '@'
# slot access in user-facing code (documentation, examples, vignettes).
# These wrap the most commonly used slots of a woven fit; the full slot
# list is documented at ?`woven-class`.

#' Consensus latent scores
#' @param x a \code{\linkS4class{woven}} object
#' @return numeric matrix, n x K
#' @examples
#' data(woven_example)
#' fit <- woven(woven_example$X_complete, Y = woven_example$Y, K = 3L)
#' dim(Z(fit))
#' @export
setGeneric("Z", function(x) standardGeneric("Z"))

#' @rdname Z
#' @export
setMethod("Z", "woven", function(x) x@Z)

#' Per-modality projection matrices
#' @param x a \code{\linkS4class{woven}} object
#' @return named list of V numeric matrices, each p_v x K
#' @examples
#' data(woven_example)
#' fit <- woven(woven_example$X_complete, Y = woven_example$Y, K = 3L)
#' names(WList(fit))
#' @export
setGeneric("WList", function(x) standardGeneric("WList"))

#' @rdname WList
#' @export
setMethod("WList", "woven", function(x) x@W_list)

#' Indices of the anchor (fully-observed) subjects
#' @param x a \code{\linkS4class{woven}} object
#' @return integer vector
#' @examples
#' data(woven_example)
#' fit <- woven(woven_example$X_complete, Y = woven_example$Y, K = 3L)
#' length(anchorIdx(fit))
#' @export
setGeneric("anchorIdx", function(x) standardGeneric("anchorIdx"))

#' @rdname anchorIdx
#' @export
setMethod("anchorIdx", "woven", function(x) x@anchor_idx)

#' Supervised canonical correlations
#' @param x a \code{\linkS4class{woven}} object
#' @return numeric vector of length K
#' @examples
#' data(woven_example)
#' fit <- woven(woven_example$X_complete, Y = woven_example$Y, K = 3L)
#' singularValues(fit)
#' @export
setGeneric("singularValues", function(x) standardGeneric("singularValues"))

#' @rdname singularValues
#' @export
setMethod("singularValues", "woven", function(x) x@singular_values)
