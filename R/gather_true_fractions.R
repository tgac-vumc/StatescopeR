#' Gather true fractions
#'
#' \code{gather_true_fractions} Gathers true fractions of all cell types per
#' sample from scRNAseq data
#'
#'
#' @param data SingleCellExperiment object of which to gather
#' fractions per sample
#'
#' @return DataFrame with fractions of all cell types per sample
#' @importFrom scRNAseq SegerstolpePancreasData
#' @export
#'
#' @examples
#' ## Load data
#' data <- scRNAseq::SegerstolpePancreasData()
#'
#' ## subset to 100 genes for example
#' data <- data[1:100]
#' ## Preprocess data
#' data$donor <- data$individual
#' data$label <- data$`cell type`
#'
#' ## remove NA cells
#' data <- data[, !is.na(data$label)]
#'
#' ## remove duplicates gene names
#' data <- data[!duplicated(rownames(data)), ]
#'
#' ## remove cells with less than 100 in total cohort
#' celltypes_to_remove <- names(table(data$label)[(table(data$label) < 100)])
#' data <- data[, !data$label %in% celltypes_to_remove]
#'
#' true_fractions <- gather_true_fractions(data)
gather_true_fractions <- function(data) {
    ## Init list to save fractions
    true_fractions <- list()

    ## loop over samples
    for (sample in unique(colData(data)$donor)) {
        temp_data <- data[, colData(data)$donor == sample]
        temp_true_fractions <- DataFrame(table(temp_data$label) /
                                                                ncol(temp_data))

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
