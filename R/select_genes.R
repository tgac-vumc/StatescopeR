#' Select genes using AutoGeneS
#'
#' \code{select_genes.R} select genes using AutoGeneS for deconvolution
#'
#'
#' @param scRNAseq SingleCellExperiment object to use for gene selection, should
#' be same as signature dataset
#' @param fixed_n_features integer number of genes to pick with autogenes,
#' default is NA which lets autogenes itself pick
#' @param n_hvg_genes int which allows the users to choose the number of highly
#' variable genes
#'
#' @return Vector of genes to use for deconvolution
#' @import scran basilisk reticulate
#' @importFrom scRNAseq SegerstolpePancreasData
#' @export
#'
#' @examples
#' # Load scRNAseq
#' scRNAseq <- scRNAseq::SegerstolpePancreasData()
#'
#' ## subset to first 100 genes for example
#' scRNAseq <- scRNAseq[1:100]
#'
#' # remove duplicates gene names
#' scRNAseq <- scRNAseq[!duplicated(rownames(scRNAseq)), ]
#'
#' ## Preprocess scRNAseq
#' scRNAseq$donor <- scRNAseq$individual
#' scRNAseq$label <- scRNAseq$`cell type`
#'
#' ## remove NA cells
#' scRNAseq <- scRNAseq[, !is.na(scRNAseq$label)]
#'
#' # remove cells with less than 100 in total cohort
#' celltypes_to_remove <- names(table(scRNAseq$label)[(table(scRNAseq$label)
#' < 100)])
#' scRNAseq <- scRNAseq[, !scRNAseq$label %in% celltypes_to_remove]
#'
#' ## Normalize to log cp 10k
#' scRNAseq <- normalize_scRNAseq(scRNAseq)
#'
#' ## Select genes by autogenes
#' selected_genes <- select_genes(scRNAseq, 25L, n_hvg_genes = 50L) # 25 genes
select_genes <- function(scRNAseq, fixed_n_features = NA, n_hvg_genes = 3000L) {
    ## First select hvg
    ## calculate per gene variance
    dec.data <- modelGeneVar(scRNAseq, assay.type = "logcounts")

    ## select hvg
    if (nrow(scRNAseq) < n_hvg_genes){
        hvg_genes <- rownames(scRNAseq) ## don't select hvg_genes
    }else hvg_genes <- getTopHVGs(dec.data, n = n_hvg_genes)

    ## init centroids df
    centroids <- data.frame(row.names = hvg_genes)
    ## Calculate centroids for each celltype
    for (ct in unique(scRNAseq$label)) {
        ## subset scRNAseq on celltype
        temp_scRNAseq <- scRNAseq[hvg_genes, scRNAseq$label == ct]
        ## Calculate centroids for all genes
        centroids[ct] <- rowMeans(as.array(logcounts(temp_scRNAseq)))
    }

    ## start basilisk
    proc <- basiliskStart(autogenes, testload = c("autogenes"))

    ## Select genes with AutoGeneS using Basilisk
    selected_genes <- basiliskRun(proc, fun = function( centroids, ngen, seed,
        offspring_size, fixed_n_features) {
            ## import autogenes
            ag <- reticulate::import("autogenes")
            ag$init(t(centroids))
            if (is.na(fixed_n_features)) {
                ag$optimize(ngen = ngen, seed = seed,
                            offspring_size = offspring_size, verbose = FALSE)
            } else {
                ag$optimize(
                    ngen = ngen, nfeatures = fixed_n_features, seed = seed,
                    mode = "fixed", offspring_size = offspring_size,
                    verbose = FALSE)
            }
            index <- ag$select(index = 0L)
            selected_genes <- rownames(centroids)[index]
            selected_genes
        }, centroids = centroids, ngen = 5000L, seed = 42L,
        offspring_size = 100L, fixed_n_features = fixed_n_features)

    ## stop basilisk
    basiliskStop(proc)

    return(selected_genes)
}
