#' Generate pseudobulk
#'
#' \code{generate_pseudobulk} Creates pseudobulk from raw counts of
#' SingleCellExperiment object
#'
#'
#' @param data SingleCellExperiment object of which to generate pseudobulk
#'
#' @return SummarizedExperiment with pseudobulk from the scRNAseq
#' @import SingleCellExperiment
#' @importFrom S4Vectors DataFrame SimpleList
#' @importFrom SummarizedExperiment SummarizedExperiment
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
#' pseudobulk <- generate_pseudobulk(data)
generate_pseudobulk <- function(data) {
  ## init pseudobulk df
  pseudobulk <- S4Vectors::DataFrame()

  ## loop over samples
  for (sample in unique(data$donor)) {
    ## Subset data on sample
    temp_data <- data[,data$donor == sample]
    ## sum raw counts over cells to create pseudobulk
    temp_pseudobulk <- rowSums(as.array(counts(temp_data)))
    ## add sample to pseudobulk
    pseudobulk[sample] <- temp_pseudobulk
  }

  ## Add gene symbol as rownames
  rownames(pseudobulk) <- rownames(data)

  ## Convert to SummarizedExperiment
  pseudobulk <- SummarizedExperiment(assays = S4Vectors::SimpleList(counts = pseudobulk))
  ## return pseudobulk
  return(pseudobulk)
}

