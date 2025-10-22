#' Run StateDiscovery
#'
#' \code{StateDiscovery.R} Discovers states from refined ct-specific gep
#'
#' @param Statescope SummarizedExperiment object from Statescope Refinement.
#' @param k number of cluster to choose, default is NA for automatic selection
#' @param Ncores number of cores to use for paralellization.
#' @param max_clusters maximum allowed states per cell type.
#' @param n_iter Number of initial cNMF restarts.
#' @param n_final_iter Number of final cNMF restarts.
#' @param min_cophenetic Minimum cophenetic coefficient to determine K.
#'
#' @return SummarizedExperiment object with statescores per celltype added
#' @import reticulate basilisk
#' @importFrom Matrix rowSums
#' @importFrom S4Vectors DataFrame metadata<-
#' @importFrom SummarizedExperiment assay
#' @importFrom methods is
#' @export
#'
#' @examples
#' ## Load Refined Statescope object
#' load(system.file("extdata", "example_Statescope_Refined.RData",
#'     package = "StatescopeR"
#' ))
#'
#' ## Discover states
#' Statescope <- StateDiscovery(Statescope, k = 2L, Ncores = 2L)
#'
#' ## Look at statescores and stateloadings
#' S4Vectors::metadata(Statescope)$statescores
#' S4Vectors::metadata(Statescope)$stateloadings
#'
StateDiscovery <- function(Statescope, k = NA, max_clusters = 10L, n_iter = 10L,
    n_final_iter = 100L, min_cophenetic = 0.9, Ncores = 1L) {
    if (!is(Statescope, 'SummarizedExperiment')){ ## Check Statescope input
        stop('Statescope is not a SummarizedExperiment object')}
    ## start basilisk & run StateDiscovery
    proc <- basiliskStart(statescope)
    Statescope <- basiliskRun(proc, fun <- function(Statescope, max_clusters,
    n_iter, n_final_iter, min_cophenetic, Ncores) {
        ## source cNMF code
        reticulate::source_python(system.file("python/cNMF_functions.py",
            package = "StatescopeR"))
        reticulate::source_python(system.file("python/cNMF_helper_functions.py",
            package = "StatescopeR"))

        ## perform state discovery per cell type
        statescores <- list()
        stateloadings <- list()
        for (ct in names(S4Vectors::metadata(Statescope)$ct_specific_gep)) {
            ## get ct_specific_gep for state clustering
            data_scaled <- as.matrix(SummarizedExperiment::assay(
                S4Vectors::metadata(Statescope)$ct_specific_gep[[ct]],
                "weighted_gep"))
            ## Run initial NMF runs for k selection
            if (is.na(k)) {nclust <- select_k(data_scaled, max_clusters, n_iter,
                    Ncores, min_cophenetic)
            } else {nclust <- k} # or use preselected k
            ## Run final model
            final_cNMF_result <- cNMF(data_scaled, as.integer(nclust),
                                    n_final_iter, Ncores)
            final_H <- t(final_cNMF_result[[1]]$H)
            final_H_norm <- S4Vectors::DataFrame(final_H / rowSums(final_H))
            rownames(final_H_norm) <- colnames(Statescope)

            final_W <- S4Vectors::DataFrame(final_cNMF_result[[1]]$W)
            rownames(final_W) <- rownames(S4Vectors::metadata(
                Statescope)$ct_specific_gep[[1]])
            ## Add result to ct lists
            statescores[ct] <- final_H_norm
            stateloadings[ct] <- final_W}
        ## Add statescores to Statescope obj
        S4Vectors::metadata(Statescope)$statescores <- statescores
        S4Vectors::metadata(Statescope)$stateloadings <- stateloadings
        Statescope},Statescope = Statescope, max_clusters = max_clusters,
    n_iter= n_iter, n_final_iter = n_final_iter, min_cophenetic = min_cophenetic
    , Ncores = Ncores)
    ## stop basilisk
    basiliskStop(proc)

    return(Statescope)
}
