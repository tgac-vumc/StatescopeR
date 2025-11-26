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
#' @param labels character vector with cell type labels
#'
#' @return Vector of genes to use for deconvolution
#' @import scran basilisk reticulate
#' @importFrom Matrix rowMeans
#' @importFrom methods is
#' @importFrom SingleCellExperiment logcounts<- logcounts
#' @export
#'
#' @examples
#' if (requireNamespace("scRNAseq", quietly = TRUE)) {
#'     library(scRNAseq)
#'     library(scuttle)
#'     ## Load SegerstolpePancreas data set
#'     scRNAseq <- SegerstolpePancreasData()
#'
#'     ## remove duplicate genes
#'     scRNAseq <- scRNAseq[!duplicated(rownames(scRNAseq)), ]
#'
#'     ## Subset to 1 healthy and 2 type 2 diabetes samples
#'     scRNAseq = scRNAseq[,scRNAseq$individual %in% c('H2',
#'                                                'T2D1', 'T2D2')]
#'     ## remove cells with no cell type label
#'     scRNAseq <- scRNAseq[, !is.na(scRNAseq$`cell type`)]
#'
#'     ## remove rare cell types (<100 cells in total data set)
#'     celltypes_to_remove <-
#'     names(table(scRNAseq$`cell type`)[(table(scRNAseq$`cell type`) < 100)])
#'     scRNAseq <- scRNAseq[, !scRNAseq$`cell type` %in% celltypes_to_remove]
#'
#'     ## remove NA cells
#'     scRNAseq <- scRNAseq[, !is.na(scRNAseq$`cell type`)]
#'
#'     ## Subset to first 1k genes
#'     scRNAseq = scRNAseq[1:1000,]
#'
#'     ## Normalize (cp10k) and logtransform scRNAseq
#'     cpm(scRNAseq) <- scuttle::calculateCPM(scRNAseq)
#'     logcounts(scRNAseq) <- log1p(cpm(scRNAseq)/100)
#'
#'     ## Select genes by autogenes
#'     selected_genes <- select_genes(scRNAseq, 3L, n_hvg_genes = 5L,
#'      labels = scRNAseq$`cell type`) # 3 genes
#' }
select_genes <- function(scRNAseq, fixed_n_features = NA, n_hvg_genes = 3000L,
    labels) {
    ## Check if scRNAseq is actually a SCE
    if (!is(scRNAseq, 'SingleCellExperiment')){
        stop('scRNAseq is not a SingleCellExperiment object')}
    else if (length(labels) != ncol(scRNAseq)){
        stop("labels are not the same length as number of cells in scRNAseq")}
    ## calculate per gene variance
    dec.data <- modelGeneVar(scRNAseq, assay.type = "logcounts")

    ## select hvg
    if (nrow(scRNAseq) < n_hvg_genes) {
        hvg_genes <- rownames(scRNAseq) ## don't select hvg_genes
    } else {hvg_genes <- getTopHVGs(dec.data, n = n_hvg_genes)}

    ## init centroids df
    centroids <- data.frame(row.names = hvg_genes)
    ## Calculate centroids for each celltype
    for (ct in unique(labels)) {
        ## subset scRNAseq on celltype
        temp_scRNAseq <- scRNAseq[hvg_genes, labels == ct]
        ## Calculate centroids for all genes
        centroids[ct] <- Matrix::rowMeans(logcounts(temp_scRNAseq))}

    ## start basilisk
    proc <- basiliskStart(autogenes, testload = c("autogenes"))

    ## Select genes with AutoGeneS using Basilisk
    selected_genes <- basiliskRun(proc, fun = function(centroids, ngen, seed,
    offspring_size, fixed_n_features) {
            ## import autogenes
            ag <- reticulate::import("autogenes")
            ag$init(t(centroids))
            if (is.na(fixed_n_features)) {ag$optimize(ngen = ngen, seed = seed,
                    offspring_size = offspring_size, verbose = FALSE)
                } else {ag$optimize(ngen = ngen, nfeatures = fixed_n_features,
                    seed = seed,mode = "fixed", offspring_size = offspring_size,
                    verbose = FALSE)}
            index <- ag$select(index = 0L)
            selected_genes <- rownames(centroids)[index]
            selected_genes}, centroids = centroids, ngen = 5000L, seed = 42L,
        offspring_size = 100L, fixed_n_features = fixed_n_features)
    ## stop basilisk
    basiliskStop(proc)

    return(selected_genes)
}
