#' Run BLADE deconvolution
#'
#' \code{BLADE_deconvolution.R} Runs BLADE to estimate cell fractions from bulk
#' mRNA
#'
#'
#' @param signature SimpleList with mu and sigma, respectively mean gene
#' expression and mean-variance corrected variance per cell type
#' @param bulk mRNA to be deconvolved
#' @param genes subset of genes to be used for deconvolution
#' @param cores number of cores to use for paralellization
#' @param Alpha BLADE Hyperparameter
#' @param Alpha0 BLADE Hyperparameter
#' @param Kappa0 BLADE Hyperparameter
#' @param sY BLADE Hyperparameter
#' @param Nrep Number of BLADE initializations
#' @param Nrepfinal Number of maximum optimization iterations
#'
#' @return Statescope S4 object
#' @import reticulate basilisk
#' @importFrom scRNAseq SegerstolpePancreasData
#' @importFrom SummarizedExperiment SummarizedExperiment assay
#' @importFrom S4Vectors DataFrame
#' @importFrom methods new
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
#' fractions(Statescope)
#'
BLADE_deconvolution <- function(signature, bulk, genes, cores = 1L, Alpha = 1L,
                                Alpha0 = 1000L, Kappa0 = 1L, sY = 1L,
                                Nrep = 10L, Nrepfinal = 1000L) {
    ## make matrices from Mu, omega & bulk
    Mu <- as.matrix(signature$mu[genes, ])
    Omega <- as.matrix(signature$omega[genes, ])
    bulk <- as.matrix(assay(bulk[genes, ], "normalized_counts"))

    ## start basilisk
    setBasiliskShared(TRUE)
    proc <- basiliskStart(deconvolution)

    ## Estimate fractions with BLADE using Basilisk
    Statescope <- basiliskRun(proc, fun = function(Mu, Omega, bulk, Alpha,
                                                    Alpha0, Kappa0, sY, Nrep,
                                                    cores, Nrepfinal) {
        ## import BLADE
        reticulate::source_python(system.file("python/BLADE.py",
        package = "StatescopeR"))

        ## Run deconvolution
        result <- Framework_Iterative(Mu, Omega,
            bulk, Alpha = Alpha, Alpha0 = Alpha0, Kappa0 = Kappa0, sY = sY,
            Nrep = Nrep, Njob = cores, IterMax = Nrepfinal)

            ## make fractions DF
            fractions <- DataFrame(t(result[[1]]$ExpF(result[[1]]$Beta)),
                row.names = colnames(Mu))
            colnames(fractions) <- colnames(bulk)

            ## Make named list from result for use in refinement
            result[[1]] <- list("Alpha" = result[[1]]$Alpha,
                                "Beta" = result[[1]]$Beta)

            ## Remove Python object
            result[[3]] <- NULL

            ## Save S4 object with Statescope and fractions slots
            Statescope <- new("Statescope", BLADE_output = result,
                                            fractions = fractions)
            Statescope
        }, Mu = Mu, Omega = Omega, bulk = bulk, Alpha = Alpha, Alpha0 = Alpha0,
        Kappa0 = Kappa0, sY = sY, Nrep = Nrep, cores = cores,
        Nrepfinal = Nrepfinal)

    ## stop basilisk
    basiliskStop(proc)

    return(Statescope)
}
