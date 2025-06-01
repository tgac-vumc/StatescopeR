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
#' @importFrom Matrix rowSums
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
#' ## Perform Deconvolution with BLADE, refine gene expression estimates
#' Statescope <- BLADE_deconvolution(
#'     signature, pseudobulk, selected_genes,
#'     prior, 1L
#' )
#' Statescope <- Refinement(Statescope, signature, pseudobulk, 1L)
#'
#' ## Discover states
#' Statescope <- StateDiscovery(Statescope)
#'
#' ## Look at output
#' statescores(Statescope)
#' stateloadings(Statescope)
#'
StateDiscovery <- function(Statescope, max_clusters = 10L, n_iter = 10L,
        n_final_iter = 100L, min_cophenetic = 0.9, Ncores = 1L) {
    ## start basilisk & run StateDiscovery
    setBasiliskShared(FALSE)
    proc <- basiliskStart(statescope)
    Statescope <- basiliskRun(proc, fun <- function(Statescope, max_clusters,
                                                    n_iter, n_final_iter,
                                                    min_cophenetic, Ncores) {
        ## source cNMF code
        reticulate::source_python(system.file("python/cNMF_functions.py",
            package = "StatescopeR"))
        reticulate::source_python(system.file("python/cNMF_helper_functions.py",
            package = "StatescopeR"))

        ## perform state discovery per cell type
        statescores <- list()
        stateloadings <- list()
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
            final_H <- t(final_cNMF_result[[1]]$H)
            final_H_norm <- DataFrame(final_H / rowSums(final_H)) # Sum to 1
            rownames(final_H_norm) <- colnames(fractions(Statescope))

            final_W <- DataFrame(final_cNMF_result[[1]]$W)
            rownames(final_W) <- rownames(ct_specific_gep(Statescope)[[1]])
            ## Add result to ct lists
            statescores[ct] <- final_H_norm
            stateloadings[ct] <- final_W
        }
        ## Add statescores to Statescope obj
        statescores(Statescope) <- statescores
        stateloadings(Statescope) <- stateloadings
        Statescope},
    Statescope = Statescope, max_clusters = max_clusters, n_iter = n_iter,
    n_final_iter = n_final_iter, min_cophenetic = min_cophenetic,
    Ncores = Ncores)
    ## stop basilisk
    basiliskStop(proc)

    return(Statescope)
}
