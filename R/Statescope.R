#' Run Statescope
#'
#' \code{Statescope.R} Discovers states in purified resuls
#'
#' @param Statescope Statescope obj from BLADE_deconvolution
#' @param cores number of cores to use for paralellization
#'
#' @return updated Statescope S4 object with states per celltype added
#' @import reticulate basilisk
#' @export
#'
#' @examples
#' Statescope = Statescope(Statescope, 1L)
Statescope <- function(Statescope,
                       max_clusters = 10L,
                       n_iter = 10L,
                       n_final_iter = 100L,
                       min_cophenetic = 0.9,
                       Ncores = 1L) {
    ## init list for states
    states = list()
    ## get celltypes
    cts = names(Statescope@ct_specific_gep)

    ## start basilisk
    proc <- basiliskStart(statescope)

    ## perform state discovery within Basilisk
    states <- basiliskRun(proc, fun = function(states, cts, Statescope,
                       max_clusters, n_iter, n_final_iterL, min_cophenetic,
                       Ncores) {
            ## source cNMF code
            reticulate::source_python(system.file('python/cNMF_functions.py',
                                                  package ='StatescopeR'))
            reticulate::source_python(system.file(
                'python/cNMF_helper_functions.py',package ='StatescopeR'))

            ## perform state discovery within Basilisk
            for (ct in cts) {
                ## get ct_specific_gep for state clustering
                data_scaled = as.matrix(assay(Statescope@ct_specific_gep[[ct]],
                                              'weighted_gep'))
                #-----------------------------------------------------------------------
                # 1.1 Run initial NMF runs for k selection
                #-----------------------------------------------------------------------
                data_dict = list()
                for (k in seq(2, max_clusters)) {
                    cNMF_result = cNMF(data_scaled, max_clusters, n_iter,
                                       Ncores)
                    H = cNMF_result[[1]]$H
                    cluster_assignment = list()
                    for (i in seq(ncol(H))) {
                        cluster_assignment = append(cluster_assignment,
                                                    which(H[, i] == max(H[, i])))
                    }
                    data_dict[k] = SimpleList(
                        'model' = cNMF_result[[1]],
                        'cophcor' = cNMF_result[[2]],
                        'consensus' = cNMF_result[[3]],
                        'cluster_assignments' = cluster_assignment
                    )

                }
                #-----------------------------------------------------------------------
                # 1.2 Choose k
                #-----------------------------------------------------------------------
                ## Extract ks and cophcors
                ks = c()
                cophcors = c()
                for (k in seq(2, max_clusters)) {
                    ks = append(ks, k)
                    cophcors = append(cophcors, data_dict[[k]]$cophcor)
                }

                nclust = find_threshold(cophcors, ks, min_cophenetic)
                drop = biggest_drop(cophcors)
                if (!nclust) {
                    nclust = drop
                }
                #-----------------------------------------------------------------------
                # 1.3 Run final model
                #-----------------------------------------------------------------------
                final_cNMF_result = cNMF(data_scaled, as.integer(nclust),
                                         n_final_iter, Ncores)
                final_H = DataFrame(t(final_cNMF_result[[1]]$H))
                rownames(final_H) = colnames(Statescope@fractions)

                ## Add result to states
                states[ct] = final_H
            }
            states
        }, states = states, cts = cts, Statescope = Statescope,
        max_clusters = max_clusters, n_iter = n_iter,
        n_final_iter = n_final_iter, min_cophenetic = min_cophenetic,
        Ncores = Ncores)

    ## stop basilisk
    basiliskStop(proc)

    ## Add States to Statescope object
    Statescope@states = states

    return(Statescope)

}
