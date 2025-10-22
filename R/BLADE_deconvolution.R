#' Run BLADE deconvolution
#'
#' \code{BLADE_deconvolution.R} Runs BLADE to estimate cell fractions from bulk
#' mRNA
#'
#'
#' @param signature SimpleList with mu and sigma, both nGene x nCelltype
#' dataframes, respectively mean gene expression and mean-variance corrected
#' variance per cell type
#' @param bulk nSample x nGene mRNA to be deconvolved
#' @param genes subset of genes to be used for deconvolution
#' @param prior (optional) nSample x nCelltype matrix with prior fraction
#' expectations
#' @param cores number of cores to use for paralellization
#' @param Alpha BLADE Hyperparameter
#' @param Alpha0 BLADE Hyperparameter
#' @param Kappa0 BLADE Hyperparameter
#' @param sY BLADE Hyperparameter
#' @param Nrep Number of BLADE initializations
#' @param Nrepfinal Number of maximum optimization iterations
#'
#' @return SummarizedExperiment object with BLADE output and estimated fractions
#' @import reticulate basilisk
#' @importFrom SummarizedExperiment SummarizedExperiment assay
#' @importFrom S4Vectors DataFrame metadata<-
#' @importFrom methods new
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
#'     ## Subset to 3 healthy and 3 type 2 diabetes samples
#'     scRNAseq = scRNAseq[,scRNAseq$individual %in% c('H2', 'H3', 'H4',
#'                                                'T2D1', 'T2D2', 'T2D3')]
#'     ## remove cells with no cell type label
#'     scRNAseq <- scRNAseq[, !is.na(scRNAseq$`cell type`)]
#'
#'     ## remove very rare cell types (<120 cells in total data set)
#'     celltypes_to_remove <- names(table(scRNAseq$`cell type`)
#'         [(table(scRNAseq$`cell type`) < 150)])
#'     scRNAseq <- scRNAseq[, !scRNAseq$`cell type` %in% celltypes_to_remove]
#'
#'     ## Create pseudobulk and normalize to cp10k
#'     pseudobulk <- aggregateAcrossCells(scRNAseq, ids = scRNAseq$individual)
#'     normcounts(pseudobulk) <- calculateCPM(pseudobulk)/100
#'     pseudobulk = as(pseudobulk, "SummarizedExperiment")
#'     rownames(pseudobulk) = rownames(scRNAseq)
#'
#'     ## Load signature
#'     load(system.file("extdata", "example_signature.RData",
#'         package = "StatescopeR"))
#'
#'     ## Load selected_genes
#'     load(system.file("extdata", "example_selected_genes.RData",
#'         package = "StatescopeR"))
#'
#'     ##  Load prior
#'     load(system.file("extdata", "example_prior.RData",
#'         package = "StatescopeR"))
#'
#'     ## Perform Deconvolution with BLADE
#'     Statescope <- BLADE_deconvolution(
#'         signature, pseudobulk, selected_genes,
#'         prior, 1L,
#'         Nrep = 1L)
#'
#'     ## show estimated fractions
#'     S4Vectors::metadata(Statescope)$fractions
#'     }
BLADE_deconvolution <- function(signature, bulk, genes, prior = NULL,
    cores = 1L, Alpha = 1L, Alpha0 = 1000L, Kappa0 = 1L, sY = 1L, Nrep = 10L,
        Nrepfinal = 1000L) {
    if (any(colSums(as.matrix(assay(bulk, "normcounts"))) != 10000)){
        stop('bulk is not in counts per 10k')}
    ## make matrices from Mu, omega, bulk & prior
    Mu <- as.matrix(signature$mu[genes, ])
    Omega <- as.matrix(signature$omega[genes, ])
    bulk_matrix <- as.matrix(assay(bulk[genes, ], "normcounts"))
    if (is.null(prior)) { # keep null if no prior is given
    } else if (is.list(prior)) { # do nothing if group prior
    } else {prior <- as.matrix(prior)[colnames(bulk), colnames(signature$mu)]}

    ## start basilisk
    proc <- basiliskStart(deconvolution)

    ## Estimate fractions with BLADE using Basilisk
    Statescope <- basiliskRun(proc,fun = function(bulk, Mu, Omega, bulk_matrix,
    prior, Alpha, Alpha0, Kappa0, sY, Nrep, cores, Nrepfinal) {
            ## import BLADE
            reticulate::source_python(system.file("python/BLADE.py",
                package = "StatescopeR"))

            ## Run deconvolution
            result <- Framework_Iterative(Mu, Omega, bulk_matrix,
                Expectation = prior, Alpha = Alpha, Alpha0 = Alpha0,
                Kappa0 = Kappa0, sY = sY, Nrep = Nrep, Njob = cores,
                IterMax = Nrepfinal)

            ## make fractions DF
            fractions <- S4Vectors::DataFrame(t(result[[1]]$ExpF(
                result[[1]]$Beta)), row.names = colnames(Mu))
            colnames(fractions) <- colnames(bulk_matrix)

            ## Make named list from result for use in refinement
            result[[1]] <- list("Alpha" = result[[1]]$Alpha, "Beta" =
                                    result[[1]]$Beta)
            result[[3]] <- NULL # Remove Python object
            ## Save BLADE parameters and fractions to metadata
            S4Vectors::metadata(bulk)$BLADE_output <- result
            S4Vectors::metadata(bulk)$fractions <- fractions
            bulk}, bulk = bulk, Mu = Mu, Omega = Omega,
    bulk_matrix = bulk_matrix,prior = prior, Alpha = Alpha, Alpha0 = Alpha0,
    Kappa0 = Kappa0, sY = sY,Nrep = Nrep, cores = cores, Nrepfinal = Nrepfinal)

    ## stop basilisk
    basiliskStop(proc)

    return(Statescope)
}
