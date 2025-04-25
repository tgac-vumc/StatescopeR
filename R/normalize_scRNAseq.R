#' Normalize scRNAseq data
#'
#' \code{normalize_scRNAseq} adds norm and lognorm counts to your
#' SingleCellExperiment object
#'
#' @param SingleCellExperiment object
#'
#' @return SingleCellExperiment object with add normcounts and logcounts assays
#' @import dplyr
#' @import SingleCellExperiment
#' @importFrom methods as
#' @importFrom scRNAseq SegerstolpePancreasData
#' @export
#'
#' @examples
#' ## Load data
#' data <- scRNAseq::SegerstolpePancreasData()
#'
#' ## subset to 100 genes for example
#' data <- data[1:100]
#' data <- normalize_scRNAseq(data)
normalize_scRNAseq <- function(SingleCellExperiment) {
    ## Calculate cp 10k
    norm_counts <- as(as.matrix(as_tibble(as.matrix(counts(SingleCellExperiment))) %>% mutate_all(funs(. / sum(.) * 10000))), "dgCMatrix")
    rownames(norm_counts) <- rownames(SingleCellExperiment) ## add rownames back

    ## assign normalized and lognormalized counts
    normcounts(SingleCellExperiment) <- norm_counts
    logcounts(SingleCellExperiment) <- log1p(normcounts(SingleCellExperiment))

    ## return sc experiment object
    return(SingleCellExperiment)
}
