#' Normalize bulk RNAseq data
#'
#' \code{normalize_bulkRNAseq} normalizes to counts per 10k
#'
#' @param bulk SummarizedExperiment with raw bulk mRNA in the first assay
#'
#' @return SummarizedExperiment with added normalized counts
#' @import dplyr
#' @importFrom scRNAseq SegerstolpePancreasData
#' @importFrom SummarizedExperiment SummarizedExperiment assay assay<-
#' @importFrom SummarizedExperiment assayNames assayNames<-
#' @export
#'
#' @examples
#' #' ## Load data
#' data <- scRNAseq::SegerstolpePancreasData()
#'
#' ## subset to 100 genes for example
#' data <- data[1:100]
#' ## Preprocess data
#' data$donor <- data$individual
#' data$label <- data$`cell type`
#' ## Create and normalized pseudobulk from scRNAseq
#' pseudobulk <- generate_pseudobulk(data)
#' pseudobulk_norm <- normalize_bulkRNAseq(pseudobulk)
normalize_bulkRNAseq <- function(bulk) {
    ## Normalize to cp 10k
    bulk_norm <- DataFrame(as_tibble(assay(bulk)) %>%
        mutate_all(funs(. / sum(.) * 10000)))

    ## Add gene symbols back as rownames
    rownames(bulk_norm) <- rownames(bulk)

    ## add normalized counts to SummarizedExperiment object
    assay(bulk, 2) <- bulk_norm
    assayNames(bulk)[2] <- "normalized_counts" ## set name

    ## return SummarizedExperiment with normalized bulk
    return(bulk)
}
