#' Run Refinement
#'
#' \code{Refinement.R} # Perform Gene Expression Refinement
#'
#' @details
#' This function takes the output from BLADE_deconvolution and refines the
#' cell type specific gene expression by reoptimizing the initial estimates.
#' It reoptimzes by fixing the estimated fractions and weighing the objective
#' function value in a way that it tries to resemble the bulk RNA expression
#' more than initially
#'
#'
#' @param Statescope Statescope obj from BLADE_deconvolution
#' @param signature SimpleList with mu and sigma, respectively mean
#' gene expression and mean-variance corrected variance per cell type
#' @param bulk mRNA to be deconvolved
#' @param cores number of cores to use for paralellization
#'
#' @return updated Statescope S4 object with ct_specific_gep added
#' @import basilisk reticulate
#' @importFrom S4Vectors DataFrame
#' @export
#'
#' @examples
#' ## Load data
#' data <- scRNAseq::SegerstolpePancreasData()
#'
#' ## subset to 100 genes for example
#' data <- data[1:100]
#' ## Preprocess data
#' data$donor <- data$individual
#' data$label <- data$`cell type`
#'
#' ## remove NA cells
#' data <- data[, !is.na(data$label)]
#'
#' ## remove duplicates gene names
#' data <- data[!duplicated(rownames(data)), ]
#'
#' ## remove cells with less than 100 in total cohort
#' celltypes_to_remove <- names(table(data$label)[(table(data$label) < 100)])
#' data <- data[, !data$label %in% celltypes_to_remove]
#'
#' data <- normalize_scRNAseq(data)
#'
#' ## Create and normalized pseudobulk from scRNAseq
#' pseudobulk <- generate_pseudobulk(data)
#'
#' pseudobulk <- normalize_bulkRNAseq(pseudobulk)
#'
#' ## Measure true cell fractions from pseudobulk
#' true_fractions <- gather_true_fractions(data)
#'
#' ## Create signature from scRNAseq for deconvolution
#' signature <- create_signature(data)
#'
#' ## Select genes optimized for deconvolution
#' selected_genes <- select_genes(data)
#'
#' ## Perform Deconvolution with BLADE
#' Statescope <- BLADE_deconvolution(signature, pseudobulk, selected_genes, 1L)
#' Statescope <- Refinement(Statescope, signature, pseudobulk, 1L)
#' ct_specific_gep(Statescope)
Refinement <- function(Statescope, signature, bulk, cores = 1L) {
    ## Prepare Refinement input
    BLADE_obj <- list("final_obj" = BLADE_output(Statescope)[[1]],
        "outs" = BLADE_output(Statescope)[[3]])
    Mu <- as.matrix(signature$mu)
    Omega <- as.matrix(signature$omega)
    genes <- rownames(Mu) # subset bulk on signature genes
    bulk <- as.matrix(assay(bulk, "normalized_counts"))[genes,]

    ## start basilisk & Run Refinement
    setBasiliskShared(TRUE)
    proc <- basiliskStart(deconvolution)
    Statescope <- basiliskRun(proc, fun = function(BLADE_obj, Mu, Omega, bulk,
                                                    cores) {
        ## import BLADE
        reticulate::source_python(system.file("python/BLADE.py",
            package = "StatescopeR"))

        ## Run refinement
        result <- Purify_AllGenes(BLADE_obj, Mu, Omega, bulk, cores)

        ## update BLADE results
        BLADE_output(Statescope) <- result

        ## Gather ct specific gep
        ct_specific_gep <- list()
        for (i in seq_along(colnames(signature$mu))) {
            temp_gep <- t(result[[1]]$Nu[, , i])
            ## Weight by Omega
            omega_weighted_gep <- (temp_gep - colMeans(temp_gep)) *
                result[[1]]$Omega[, i]

            ## Name samples and genes
            omega_weighted_gep_df <- DataFrame(omega_weighted_gep)
            colnames(omega_weighted_gep_df) <- colnames(bulk)
            rownames(omega_weighted_gep_df) <- rownames(bulk)

            ## add to ct_specific_gep list
            ct_specific_gep[colnames(signature$mu)[i]] <- SummarizedExperiment(
                assays = SimpleList(weighted_gep = omega_weighted_gep_df))}

        ## Add cell type specific gene expression to Statescope obj
        ct_specific_gep(Statescope) <- ct_specific_gep
        Statescope
    }, BLADE_obj = BLADE_obj, Mu = Mu, Omega = Omega, bulk = bulk,
                                                                cores = cores)
    ## stop basilisk
    basiliskStop(proc)

    return(Statescope)
}
