#' Run Statescope
#'
#' \code{Statescope.R} Discovers states in purified resuls
#'
#' @param Statescope Statescope obj from BLADE_deconvolution
#' @param Ncores number of cores to use for paralellization
#' @param max_clusters ...
#' @param n_iter ...
#' @param n_final_iter ...
#' @param min_cophenetic ...
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
#' ## Perform Deconvolution with BLADE
#' Statescope <- BLADE_deconvolution(signature, pseudobulk, selected_genes, 1L)
#' Statescope <- BLADE_purification(Statescope, signature, pseudobulk, 1L)
#' Statescope <- Statescope(Statescope, 1L)
#'
#' ## Look at output
#' states(Statescope)
#'
Statescope <- function(Statescope,
    max_clusters = 10L,
    n_iter = 10L,
    n_final_iter = 100L,
    min_cophenetic = 0.9,
    Ncores = 1L) {
    ## init list for states
    states <- list()
    ## get celltypes
    cts <- names(ct_specific_gep(Statescope))

    ## start basilisk
    setBasiliskShared(FALSE)
    proc <- basiliskStart(statescope)

    ## perform state discovery within Basilisk
    states <- basiliskRun(proc, fun <- function(states, cts, Statescope,
    max_clusters, n_iter, n_final_iterL, min_cophenetic,
    Ncores) {
        ## source cNMF code
        reticulate::source_python(system.file("python/cNMF_functions.py",
            package = "StatescopeR"
        ))
        reticulate::source_python(system.file(
            "python/cNMF_helper_functions.py",
            package = "StatescopeR"
        ))

        ## perform state discovery within Basilisk
        for (ct in cts) {
            ## get ct_specific_gep for state clustering
            data_scaled <- as.matrix(assay(
                ct_specific_gep(Statescope)[[ct]], "weighted_gep"
            ))
            #---------------------------------------------------------------
            # 1.1 Run initial NMF runs for k selection
            #---------------------------------------------------------------
            data_dict <- list()
            for (k in seq(2, max_clusters)) {
                cNMF_result <- cNMF(
                    data_scaled, max_clusters, n_iter,
                    Ncores
                )
                H <- cNMF_result[[1]]$H
                cluster_assignment <- list()
                for (i in seq(ncol(H))) {
                    cluster_assignment <- append(
                        cluster_assignment,
                        which(H[, i] == max(H[, i]))
                    )
                }
                data_dict[k] <- SimpleList(
                    "model" = cNMF_result[[1]],
                    "cophcor" = cNMF_result[[2]],
                    "consensus" = cNMF_result[[3]],
                    "cluster_assignments" = cluster_assignment
                )
            }
            #---------------------------------------------------------------
            # 1.2 Choose k
            #---------------------------------------------------------------
            ## Extract ks and cophcors
            ks <- c()
            cophcors <- c()
            for (k in seq(2, max_clusters)) {
                ks <- append(ks, k)
                cophcors <- append(cophcors, data_dict[[k]]$cophcor)
            }

            nclust <- find_threshold(cophcors, ks, min_cophenetic)
            drop <- biggest_drop(cophcors)
            if (!nclust) {
                nclust <- drop
            }
            #---------------------------------------------------------------
            # 1.3 Run final model
            #---------------------------------------------------------------
            final_cNMF_result <- cNMF(
                data_scaled, as.integer(nclust),
                n_final_iter, Ncores
            )
            final_H <- DataFrame(t(final_cNMF_result[[1]]$H))
            rownames(final_H) <- colnames(fractions(Statescope))

            ## Add result to states
            states[ct] <- final_H
        }
        states
    },
    states = states, cts = cts, Statescope = Statescope,
    max_clusters = max_clusters, n_iter = n_iter,
    n_final_iter = n_final_iter, min_cophenetic = min_cophenetic,
    Ncores = Ncores
    )

    ## stop basilisk
    basiliskStop(proc)

    ## Add States to Statescope object
    states(Statescope) <- states

    return(Statescope)
}
