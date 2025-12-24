#' Gather true fractions
#'
#' \code{gather_true_fractions} Gathers true fractions of all cell types per
#' sample from scRNAseq data
#'
#'
#' @param scRNAseq SingleCellExperiment object of which to gather
#' fractions per sample
#' @param ids character vector with ids of samples
#' @param label_col character for the column name with the cell type labels
#'
#' @return DataFrame with fractions of all cell types per sample
#' @importFrom SingleCellExperiment colData
#' @importFrom methods is
#' @export
#'
#' @examples
#' if (requireNamespace("scRNAseq", quietly = TRUE)) {
#'     library(scRNAseq)
#'     ## Load data
#'     scRNAseq <- scRNAseq::SegerstolpePancreasData()
#'     ## Subset to 1 healthy and 3 type 2 diabetes samples
#'     scRNAseq <- scRNAseq[, scRNAseq$individual %in% c(
#'         "H3",
#'         "T2D1", "T2D2"
#'     )]
#'     ## remove NA cells
#'     scRNAseq <- scRNAseq[, !is.na(scRNAseq$`cell type`)]
#'
#'     ## remove cells with less than 100 in total cohort
#'     celltypes_to_remove <- names(table(scRNAseq$`cell type`)
#'     [(table(scRNAseq$`cell type`) < 100)])
#'     scRNAseq <- scRNAseq[, !scRNAseq$`cell type` %in% celltypes_to_remove]
#'
#'     true_fractions <- gather_true_fractions(scRNAseq,
#'         ids = scRNAseq$individual, label_col = "cell type"
#'     )
#' }
gather_true_fractions <- function(scRNAseq, ids, label_col) {
    if (!is(scRNAseq, "SingleCellExperiment")) {
        stop("scRNAseq is not a SingleCellExperiment object")
    }
    ## Init list to save fractions
    true_fractions <- list()

    ## loop over samples
    for (sample in unique(ids)) {
        temp_scRNAseq <- scRNAseq[, ids == sample]
        temp_true_fractions <- DataFrame(table(
            colData(temp_scRNAseq)[[label_col]]
        ) / ncol(temp_scRNAseq))

        ## add True_fractions of sample to list
        true_fractions[sample] <- list(temp_true_fractions$Freq)
    }

    ## Add all samples to one DataFrame
    true_fractions <- DataFrame(true_fractions,
        row.names = temp_true_fractions$Var1
    )

    ## return DF
    return(true_fractions)
}
