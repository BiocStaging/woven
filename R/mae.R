# mae.R -- MultiAssayExperiment input support for woven()
#
# MultiAssayExperiment (Ramos et al. 2017) is the Bioconductor container for
# exactly the setting woven() targets: several assay platforms measured on
# overlapping but non-identical sets of patients. This file converts a
# MultiAssayExperiment into the named list of subject x feature matrices
# (with all-NA rows for block-missing subjects) that the rest of the package
# operates on, so no internal code needs to know about MultiAssayExperiment.

#' Convert a MultiAssayExperiment to woven's internal X_list format
#'
#' Aligns every experiment in \code{mae} to the shared patient roster
#' (\code{rownames(colData(mae))}) using \code{sampleMap(mae)}, transposing
#' each assay to subject x feature and filling an all-NA row for any subject
#' not present in that experiment. This is the same block-missing
#' representation used throughout the package, so the result can be passed
#' directly as \code{X_list} to any woven function.
#'
#' @param mae a MultiAssayExperiment object
#' @return named list of V numeric matrices, each n x p_v, n = number of
#'   subjects in \code{colData(mae)}, row order matching
#'   \code{rownames(colData(mae))}
#' @keywords internal
.mae_to_xlist <- function(mae) {
    if (!requireNamespace("MultiAssayExperiment", quietly = TRUE)) {
        stop("Package 'MultiAssayExperiment' is required to pass a ",
            "MultiAssayExperiment to woven(). Install it with ",
            "BiocManager::install('MultiAssayExperiment').",
            call. = FALSE
        )
    }
    primary <- rownames(MultiAssayExperiment::colData(mae))
    exps <- MultiAssayExperiment::experiments(mae)
    smap <- as.data.frame(MultiAssayExperiment::sampleMap(mae))

    X_list <- lapply(names(exps), function(nm) {
        ex <- exps[[nm]]
        m <- if (methods::is(ex, "SummarizedExperiment")) {
            SummarizedExperiment::assay(ex)
        } else {
            as.matrix(ex)
        }
        mt <- t(m) # features x samples -> samples x features
        map_i <- smap[smap$assay == nm, , drop = FALSE]
        rownames(mt) <- map_i$primary[match(rownames(mt), map_i$colname)]
        out <- matrix(NA_real_,
            nrow = length(primary), ncol = ncol(mt),
            dimnames = list(primary, colnames(mt))
        )
        common <- intersect(primary, rownames(mt))
        out[common, ] <- mt[common, , drop = FALSE]
        out
    })
    names(X_list) <- names(exps)
    X_list
}
