#' Normalize scRNAseq data
#'
#' \code{normalize_scRNAseq} adds norm and lognorm counts to your
#' SingleCellExperiment object
#'
#' @param SingleCellExperiment object
#'
#' @return SingleCellExperiment object with add normcounts and logcounts assays
#' @export
#'
#' @examples
#' data = normalize_scRNAseq(data)
normalize_scRNAseq <- function(data) {
  ## Calculate cp 10k
  norm_counts = as(as.matrix(as_tibble(as.matrix(counts(data))) %>% mutate_all(funs(./sum(.)*10000))), "dgCMatrix")
  rownames(norm_counts) = rownames(data) ## add rownames back

  ## assign normalized and lognormalized counts
  normcounts(data) <- norm_counts
  logcounts(data) <- log1p(normcounts(data))

  ## return sc experiment object
  return(data)
}

