#' Run StateDiscovery
#'
#' \code{StateDiscovery.R} Discovers states from refined ct-specific gep
#'
#' @param Statescope Statescope obj from StatescopeRefinement.
#' @param Ncores number of cores to use for paralellization.
#' @param max_clusters maximum allowed states per cell type.
#' @param n_iter Number of initial cNMF restarts.
#' @param n_final_iter Number of final cNMF restarts.
#' @param min_cophenetic Minimum cophenetic coefficient to determine K.
#'
#' @return updated Statescope S4 object with statescores per celltype added
#' @import reticulate basilisk
#' @export
#'
#' @examples
#' #' ## Load data
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
#' ## Perform Deconvolution with BLADE, refine gene expression estimates
#' Statescope <- BLADE_deconvolution(signature, pseudobulk, selected_genes, 1L)
#' Statescope <- Refinement(Statescope, signature, pseudobulk, 1L)
#'
#' ## Discover states
#' Statescope <- StateDiscovery(Statescope)
#'
#' ## Look at output
#' statescores(Statescope)
#'
StateDiscovery <- function(Statescope, max_clusters = 10L, n_iter = 10L,
    n_final_iter = 100L, min_cophenetic = 0.9, Ncores = 1L) {
    ## start basilisk & run StateDiscovery
    setBasiliskShared(FALSE)
    proc <- basiliskStart(statescope)
    Statescope <- basiliskRun(proc, fun <- function(Statescope,
    max_clusters, n_iter, n_final_iter, min_cophenetic, Ncores) {
        ## source cNMF code
        reticulate::source_python(system.file("python/cNMF_functions.py",
            package = "StatescopeR"))
        reticulate::source_python(system.file("python/cNMF_helper_functions.py",
            package = "StatescopeR"))

        ## perform state discovery per cell type
        statescores <- list()
        for (ct in names(ct_specific_gep(Statescope))) {
            ## get ct_specific_gep for state clustering
            data_scaled <- as.matrix(assay(ct_specific_gep(Statescope)[[ct]],
                                                                "weighted_gep"))
            ## Run initial NMF runs for k selection
            nclust <- select_k(data_scaled, max_clusters, n_iter, Ncores,
                                min_cophenetic)
            ## Run final model
            final_cNMF_result <- cNMF(data_scaled, as.integer(nclust),
                                        n_final_iter, Ncores)
            final_H <- DataFrame(t(final_cNMF_result[[1]]$H))
            rownames(final_H) <- colnames(fractions(Statescope))
            ## Add result to ct lists
            statescores[ct] <- final_H}
        ## Add statescores to Statescope obj
        statescores(Statescope) <- statescores
        Statescope
    },
    Statescope = Statescope, max_clusters = max_clusters, n_iter = n_iter,
    n_final_iter = n_final_iter, min_cophenetic = min_cophenetic,
    Ncores = Ncores)

    ## stop basilisk
    basiliskStop(proc)

    return(Statescope)
}
