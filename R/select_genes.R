#' Select genes
#'
#' \code{select_genes.R} select genes to use for deconvolution
#'
#'
#' @param data SingleCellExperiment object to use for gene selection, should be
#' same as signature dataset
#'
#' @return Vector of genes to use for deconvolution
#' @import scran reticulate
#' @export
#'
#' @examples
#' selected_genes = select_genes(data)
select_genes <- function(data) {
  ## import autogenes
  ag <- import('autogenes')
  print('This takes about 15 minutes with 4k genes and 14 cell types ')
  ## First select hvg
  ## calculate per gene variance
  dec.data = modelGeneVar(data, assay.type = 'logcounts')

  ## select hvg
  hvg_genes = getTopHVGs(dec.data, n = 4000L)

  ## init centroids df
  centroids = data.frame(row.names = hvg_genes)
  ## Calculate centroids for each celltype
  for (ct in unique(data$label)){
    ## subset data on celltype
    temp_data = data[hvg_genes,data$label == ct]
    ## Calculate centroids for all genes
    centroids[ct] = rowMeans(as.array(logcounts(temp_data)))

  }

  ## Select genes with AutoGeneS
  ag$init(t(centroids))
  ag$optimize(ngen= 5000L, seed = 42L, offspring_size = 100L, verbose = FALSE)
  index = ag$select(index=0L)
  selected_genes = rownames(centroids)[index]

  return(selected_genes)

}

