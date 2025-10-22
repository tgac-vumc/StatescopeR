#' Create scRNAseq Signature
#'
#' \code{create_signature} Creates signature from scRNAseq data
#'
#'
#' @param scRNAseq SingleCellExperiment object of which to make signature
#' @param hvg_genes boolean which chooses if mu and omega should be subset to
#' highly variable genes or not
#' @param n_hvg_genes int which allows the users to choose the number of highly
#' variable genes
#' @param labels character vector for the cell type labels
#'
#' @return SimpleList DataFrames for Mu (mean per gene per cell type) and
#' Omega (variance corrected std.dev per gene per cell type)
#' @import scran
#' @importFrom matrixStats rowSds rowVars
#' @importFrom SingleCellExperiment logcounts<- logcounts
#' @importFrom methods is
#' @export
#'
#' @examples
#' if (requireNamespace("scRNAseq", quietly = TRUE)) {
#'     library(scRNAseq)
#'     library(scuttle)
#'     ## Load scRNaseq
#'     scRNAseq <- scRNAseq::SegerstolpePancreasData()
#'
#'     ## remove NA cells
#'     scRNAseq <- scRNAseq[, !is.na(scRNAseq$`cell type`)]
#'
#'     ## Normalize (cp10k) and logtransform scRNAseq
#'     cpm(scRNAseq) <- scuttle::calculateCPM(scRNAseq)
#'     SingleCellExperiment::logcounts(scRNAseq) <- log1p(cpm(scRNAseq) / 100)
#'
#'     ## Create signature
#'     signature <- create_signature(scRNAseq, labels = scRNAseq$`cell type`)
#' }
create_signature <- function(scRNAseq, hvg_genes = FALSE, n_hvg_genes = 3000L,
    labels) {
    if (!is(scRNAseq, 'SingleCellExperiment')){
        stop('scRNAseq is not a SingleCellExperiment object')}
    ## init Mu, Omega & Var
    mu <- DataFrame()
    omega <- DataFrame()
    var <- DataFrame()

    for (celltype in unique(labels)) {
        ## subset scRNAseq on celltype
        temp_scRNAseq <- scRNAseq[, labels == celltype]
        ## Calculate Mu, Omega & Var
        ct_mu <- rowMeans(as.array(logcounts(temp_scRNAseq)))
        ct_omega <- rowSds(as.array(logcounts(temp_scRNAseq)))
        ct_var <- rowVars(as.array(logcounts(temp_scRNAseq)))

        ## Add ct Mu, Omega & var
        mu[celltype] <- ct_mu
        omega[celltype] <- ct_omega
        var[celltype] <- ct_var
    }

    ## Correct Omega by mean-var trend
    new_omega <- omega[, 0] ## make df with same # rows
    for (ct in colnames(omega)) {
        trend <- fitTrendVar(mu[, ct], omega[, ct])$trend
        new_omega[, ct] <- trend(mu[, ct])
    }

    ## subset on hvg_genes if true
    if (hvg_genes) {
        ## calculate per gene variance
        dec.data <- modelGeneVar(scRNAseq, assay.type = "logcounts")

        ## select hvg
        hvg_genes <- getTopHVGs(dec.data, n = n_hvg_genes)

        ## subset mu and omega
        mu <- mu[hvg_genes, ]
        new_omega <- new_omega[hvg_genes, ]
    }

    return(SimpleList(mu = mu, omega = new_omega))
}
