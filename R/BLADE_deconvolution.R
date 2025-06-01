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
#' @return Statescope S4 object
#' @import reticulate basilisk
#' @importFrom scRNAseq SegerstolpePancreasData
#' @importFrom SummarizedExperiment SummarizedExperiment assay
#' @importFrom S4Vectors DataFrame
#' @importFrom methods new
#' @export
#'
#' @examples
#' ## Load scRNAseq
#' scRNAseq <- scRNAseq::SegerstolpePancreasData()
#'
#' ## subset to 100 genes for example
#' scRNAseq <- scRNAseq[1:100]
#' ## Preprocess scRNAseq
#' scRNAseq$donor <- scRNAseq$individual
#' scRNAseq$label <- scRNAseq$`cell type`
#'
#' ## remove NA cells
#' scRNAseq <- scRNAseq[, !is.na(scRNAseq$label)]
#'
#' ## remove duplicates gene names
#' scRNAseq <- scRNAseq[!duplicated(rownames(scRNAseq)), ]
#'
#' ## remove cells with less than 100 in total cohort
#' celltypes_to_remove <-
#'     names(table(scRNAseq$label)[(table(scRNAseq$label) < 100)])
#' scRNAseq <- scRNAseq[, !scRNAseq$label %in% celltypes_to_remove]
#'
#' scRNAseq <- normalize_scRNAseq(scRNAseq)
#'
#' ## Create and normalized pseudobulk from scRNAseq
#' pseudobulk <- generate_pseudobulk(scRNAseq)
#'
#' pseudobulk <- normalize_bulkRNAseq(pseudobulk)
#'
#' ## Create signature from scRNAseq for deconvolution
#' signature <- create_signature(scRNAseq)
#'
#' ## Select genes optimized for deconvolution
#' selected_genes <- select_genes(scRNAseq)
#'
#' ## Optionally create prior expectation
#' prior <- gather_true_fractions(scRNAseq) # Use True sc fractions for this
#' prior[rownames(prior) != "ductal cell", ] <- NA # Keep only ductal cell
#'
#' ## Tranpose it to nSample x nCelltype
#' prior <- t(prior)
#'
#' ## Perform Deconvolution with BLADE
#' Statescope <- BLADE_deconvolution(
#'     signature, pseudobulk, selected_genes,
#'     prior, 1L
#' )
#' fractions(Statescope)
#'
BLADE_deconvolution <- function(signature, bulk, genes, prior = NULL,
    cores = 1L, Alpha = 1L, Alpha0 = 1000L, Kappa0 = 1L, sY = 1L, Nrep = 10L,
    Nrepfinal = 1000L) {
    ## make matrices from Mu, omega, bulk & prior
    Mu <- as.matrix(signature$mu[genes, ])
    Omega <- as.matrix(signature$omega[genes, ])
    bulk <- as.matrix(assay(bulk[genes, ], "normalized_counts"))
    if (is.null(prior)) { # keep null if no prior is given
    } else if (is.list(prior)) { # do nothing if group prior
    } else {prior <- as.matrix(prior)[colnames(bulk), colnames(signature$mu)]}

    ## start basilisk
    proc <- basiliskStart(deconvolution)

    ## Estimate fractions with BLADE using Basilisk
    Statescope <- basiliskRun(proc, fun = function(Mu, Omega, bulk, prior,
    Alpha, Alpha0, Kappa0, sY, Nrep, cores, Nrepfinal) {
            ## import BLADE
            reticulate::source_python(system.file("python/BLADE.py",
                package = "StatescopeR"))

            ## Run deconvolution
            result <- Framework_Iterative(Mu, Omega, bulk, Expectation = prior,
                Alpha = Alpha, Alpha0 = Alpha0, Kappa0 = Kappa0, sY = sY,
                Nrep = Nrep, Njob = cores, IterMax = Nrepfinal)

            ## make fractions DF
            fractions <- DataFrame(t(result[[1]]$ExpF(result[[1]]$Beta)),
                row.names = colnames(Mu))
            colnames(fractions) <- colnames(bulk)

            ## Make named list from result for use in refinement
            result[[1]] <- list("Alpha" = result[[1]]$Alpha,
                                "Beta" = result[[1]]$Beta)
            result[[3]] <- NULL # Remove Python object

            ## Save S4 object with Statescope and fractions slots
            Statescope <- new("Statescope", BLADE_output = result,
                fractions = fractions)
            Statescope
        }, Mu = Mu, Omega = Omega, bulk = bulk, prior = prior, Alpha = Alpha,
        Alpha0 = Alpha0, Kappa0 = Kappa0, sY = sY, Nrep = Nrep, cores = cores,
        Nrepfinal = Nrepfinal)

    ## stop basilisk
    basiliskStop(proc)

    return(Statescope)
}
