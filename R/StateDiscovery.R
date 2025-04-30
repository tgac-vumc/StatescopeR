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
#' @return updated Statescope S4 object with states per celltype added
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
#' Statescope <- StateDiscovery(Statescope, 1L)
#'
#' ## Look at output
#' states(Statescope)
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
        states <- list()
        for (ct in names(ct_specific_gep(Statescope))) {
            ## get ct_specific_gep for state clustering
            data_scaled <- as.matrix(assay(ct_specific_gep(Statescope)[[ct]],
                                                                "weighted_gep"))
            ## Run initial NMF runs for k selection
            data_dict <- list()
            for (k in seq(2, max_clusters)) {
                cNMF_result <- cNMF(data_scaled, max_clusters, n_iter, Ncores)
                H <- cNMF_result[[1]]$H
                cluster_assignment <- list()
                for (i in seq(ncol(H))) {
                    cluster_assignment <- append(cluster_assignment,
                                                which(H[, i] == max(H[, i])))}
                data_dict[k] <- SimpleList("model" = cNMF_result[[1]],
                                            "cophcor" = cNMF_result[[2]],
                                            "consensus" = cNMF_result[[3]],
                                    "cluster_assignments" = cluster_assignment)}
            ## 1.2 Choose k
            ks <- c()
            cophcors <- c()
            for (k in seq(2, max_clusters)) {
                ks <- append(ks, k)
                cophcors <- append(cophcors, data_dict[[k]]$cophcor)}
            nclust <- find_threshold(cophcors, ks, min_cophenetic)
            drop <- biggest_drop(cophcors)
            if (!nclust) {
                nclust <- drop}
            ## 1.3 Run final model
            final_cNMF_result <- cNMF(data_scaled, as.integer(nclust),
                                        n_final_iter, Ncores)
            final_H <- DataFrame(t(final_cNMF_result[[1]]$H))
            rownames(final_H) <- colnames(fractions(Statescope))
            ## Add result to states
            states[ct] <- final_H}
        ## Add states to Statescope obj
        states(Statescope) <- states
        Statescope
    },
    Statescope = Statescope, max_clusters = max_clusters, n_iter = n_iter,
    n_final_iter = n_final_iter, min_cophenetic = min_cophenetic,
    Ncores = Ncores)

    ## stop basilisk
    basiliskStop(proc)

    return(Statescope)
}
