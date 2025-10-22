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
#' @param Statescope SummarizedExperiment object from BLADE_deconvolution
#' @param signature SimpleList with mu and sigma, respectively mean
#' gene expression and mean-variance corrected variance per cell type
#' @param bulk mRNA to be refined
#' @param cores number of cores to use for paralellization
#'
#' @return updated SummarizedExperiment object with ct_specific_gep added
#' @import basilisk reticulate
#' @importFrom S4Vectors DataFrame metadata<- metadata SimpleList
#' @importFrom SummarizedExperiment SummarizedExperiment
#' @importFrom methods is
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
#'     ## Subset to 3 healthy and 3 type 2 diabetes samples
#'     scRNAseq = scRNAseq[,scRNAseq$individual %in% c('H2', 'H3', 'H4',
#'                                                'T2D1', 'T2D2', 'T2D3')]
#'     ## remove cells with no cell type label
#'     scRNAseq <- scRNAseq[, !is.na(scRNAseq$`cell type`)]
#'
#'     ## remove very rare cell types (<150 cells in total data set)
#'     celltypes_to_remove <-names(table(scRNAseq$`cell type`)
#'         [(table(scRNAseq$`cell type`) < 150)])
#'     scRNAseq <- scRNAseq[, !scRNAseq$`cell type` %in% celltypes_to_remove]
#'
#'     ## Create pseudobulk and normalize to cp10k
#'     pseudobulk <- aggregateAcrossCells(scRNAseq, ids = scRNAseq$individual)
#'     normcounts(pseudobulk) <- calculateCPM(pseudobulk)/100
#'     pseudobulk = as(pseudobulk, "SummarizedExperiment")
#'     rownames(pseudobulk) = rownames(scRNAseq)
#'     ## Load signature
#'     load(system.file("extdata", "example_signature.RData",
#'         package = "StatescopeR"))
#'
#'     ## Load Deconvolved Statescope object
#'     load(system.file("extdata", "example_Statescope_Deconvolved.RData",
#'         package = "StatescopeR"))
#'
#'     ## Run Refinement
#'     Statescope <- Refinement(Statescope, signature, pseudobulk, 2L)
#'
#'     ## Show cell type specific gene expression profile estimates
#'     S4Vectors::metadata(Statescope)$ct_specific_gep
#'     }
Refinement <- function(Statescope, signature, bulk, cores = 1L) {
    if (!is(Statescope, 'SummarizedExperiment')){ ## Check Statescope input
        stop('Statescope is not a SummarizedExperiment object')}
    ## Prepare Refinement input
    BLADE_obj <- list("final_obj" = metadata(Statescope)$BLADE_output[[1]],
        "outs" = metadata(Statescope)$BLADE_output[[3]])
    Mu <- as.matrix(signature$mu)
    Omega <- as.matrix(signature$omega)
    genes <- rownames(Mu) # subset bulk on signature genes
    bulk_matrix <- as.matrix(assay(bulk[genes, ], "normcounts"))

    ## start basilisk & Run Refinement
    proc <- basiliskStart(deconvolution)
    Statescope <- basiliskRun(proc, fun = function(BLADE_obj, Mu, Omega,
    bulk_matrix, cores) {
            ## import BLADE
            reticulate::source_python(system.file("python/BLADE.py",
                package = "StatescopeR"))
            ## Run refinement
            result <- Purify_AllGenes(BLADE_obj, Mu, Omega, bulk_matrix, cores)

            ## update BLADE results
            S4Vectors::metadata(Statescope)$BLADE_output <- result

            ## Gather ct specific gep
            ct_specific_gep <- list()
            for (i in seq_along(colnames(signature$mu))) {
                temp_gep <- t(result[[1]]$Nu[, , i])
                ## Weight by Omega
                omega_weighted_gep <- (temp_gep - colMeans(temp_gep)) *
                    result[[1]]$Omega[, i]
                ## Name samples and genes
                omega_weighted_gep_df <-
                    S4Vectors::DataFrame(omega_weighted_gep)
                colnames(omega_weighted_gep_df) <- colnames(bulk_matrix)
                rownames(omega_weighted_gep_df) <- rownames(bulk_matrix)
                ## add to ct_specific_gep list
                ct_specific_gep[colnames(signature$mu)[i]] <-
                    SummarizedExperiment::SummarizedExperiment(
                        assays = S4Vectors::SimpleList(weighted_gep =
                                                        omega_weighted_gep_df))}
            ## Add cell type specific gene expression to Statescope obj
            S4Vectors::metadata(Statescope)$ct_specific_gep <- ct_specific_gep
            Statescope}, BLADE_obj = BLADE_obj, Mu = Mu, Omega = Omega,
        bulk_matrix = bulk_matrix, cores = cores)
    ## stop basilisk
    basiliskStop(proc)

    return(Statescope)
}
