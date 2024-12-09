#' Normalize bulkRNAseq data
#'
#' \code{normalize_bulkRNAseq} normalizes to counts per 10k
#'
#' @param bulk SummarizedExperiment with raw bulk mRNA in the first assay
#'
#' @return SummarizedExperiment with added normalized counts
#' @export
#'
#' @examples
#' bulk_norm = normalize_bulkRNAseq(bulk)
normalize_bulkRNAseq <- function(bulk) {
  ## Normalize to cp 10k
  bulk_norm = DataFrame(as_tibble(assay(bulk)) %>% mutate_all(funs(./sum(.)*10000)))

  ## Add gene symbols back as rownames
  rownames(bulk_norm) = rownames(bulk)

  ## add normalized counts to SummarizedExperiment object
  assay(bulk, 2) = bulk_norm
  assayNames(bulk)[2] = 'normalized_counts' ## set name

  ## return SummarizedExperiment with normalized bulk
  return(bulk)
}

