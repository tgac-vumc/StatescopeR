#' Create scRNAseq Signature
#'
#' \code{create_signature} Creates signature from scRNAseq data
#'
#'
#' @param data SingleCellExperiment object of which to make signature
#' @param hvg_genes boolean which chooses if mu and omega should be subset to
#' highly variable genes or not
#'
#' @return SimpleList DataFrames for Mu (mean per gene per cell type) and
#' Omega (variance corrected std.dev per gene per cell type)
#' @import scran
#' @importFrom scRNAseq SegerstolpePancreasData
#' @importFrom matrixStats rowSds rowVars
#' @export
#'
#' @examples
#' ## Load data
#' data <- scRNAseq::SegerstolpePancreasData()
#'
#' ## subset to 100 genes for example
#' data <- data[1:100]
#' data <- normalize_scRNAseq(data)
#' signature <- create_signature(data)
create_signature <- function(data, hvg_genes = FALSE) {
    ## init Mu, Omega & Var
    mu <- DataFrame()
    omega <- DataFrame()
    var <- DataFrame()

    for (celltype in unique(data$label)) {
        ## subset data on celltype
        temp_data <- data[, data$label == celltype]
        ## Calculate Mu, Omega & Var
        ct_mu <- rowMeans(as.array(logcounts(temp_data)))
        ct_omega <- rowSds(as.array(logcounts(temp_data)))
        ct_var <- rowVars(as.array(logcounts(temp_data)))

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
        dec.data <- modelGeneVar(data, assay.type = "logcounts")

        ## select hvg
        hvg_genes <- getTopHVGs(dec.data, n = 3000L)

        ## subset mu and omega
        mu <- mu[hvg_genes, ]
        new_omega <- new_omega[hvg_genes, ]
    }

    return(SimpleList(mu = mu, omega = new_omega))
}
