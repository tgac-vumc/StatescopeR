#' Select genes using AutoGeneS
#'
#' \code{select_genes.R} select genes using AutoGeneS for deconvolution
#'
#'
#' @param data SingleCellExperiment object to use for gene selection, should be
#' same as signature dataset
#' @param fixed_n_features integer number of genes to pick with autogenes,
#' default is NA which lets autogenes itself pick
#'
#' @return Vector of genes to use for deconvolution
#' @import scran basilisk reticulate
#' @importFrom scRNAseq SegerstolpePancreasData
#' @export
#'
#' @examples
#' # Load data
#' data <- scRNAseq::SegerstolpePancreasData()
#'
#' ## subset to 100 genes for example
#' data <- data[1:100]
#' data$donor <- data$individual
#' data$label <- data$`cell type`
#' ## remove NA cells
#' data <- data[, !is.na(data$label)]
#'
#' # remove duplicates gene names
#' data <- data[!duplicated(rownames(data)), ]
#'
#' # remove cells with less than 100 in total cohort
#' celltypes_to_remove <- names(table(data$label)[(table(data$label) < 100)])
#' data <- data[, !data$label %in% celltypes_to_remove]
#' data <- normalize_scRNAseq(data)
#' selected_genes <- select_genes(data)
select_genes <- function(data, fixed_n_features = NA) {
    ## First select hvg
    ## calculate per gene variance
    dec.data <- modelGeneVar(data, assay.type = "logcounts")

    ## select hvg
    hvg_genes <- getTopHVGs(dec.data, n = 3000L)

    ## init centroids df
    centroids <- data.frame(row.names = hvg_genes)
    ## Calculate centroids for each celltype
    for (ct in unique(data$label)) {
        ## subset data on celltype
        temp_data <- data[hvg_genes, data$label == ct]
        ## Calculate centroids for all genes
        centroids[ct] <- rowMeans(as.array(logcounts(temp_data)))
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
