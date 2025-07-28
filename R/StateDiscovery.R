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
#' ## Load Refined Statescope object
#' load(system.file('extdata', 'example_Statescope_Refined.RData',
#' package = 'StatescopeR'))
#'
#' ## Discover states
#' Statescope <- StateDiscovery(Statescope, Ncores = 2L, max_clusters = 4L)
#'
#' ## Look at statescores and stateloadings
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
