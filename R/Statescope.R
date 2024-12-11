#' Run Statescope
#'
#' \code{Statescope.R} Discovers states in purified resuls
#'
#' @param Statescope Statescope obj from BLADE_deconvolution
#' @param cores number of cores to use for paralellization
#'
#' @return updated Statescope S4 object with states per celltype added
#' @export
#'
#' @examples
#' Statescope = Statescope(Statescope, 1L)
Statescope <- function(Statescope, cores = 1L) {
    ## get celltypes
    cts = names(Statescope@ct_specific_gep)
    cts = 'gamma cell'
    for (ct in cts){
        ## get ct_specific_gep for state clustering
        data_scaled = as.matrix(assay(Statescope@ct_specific_gep[[ct]],
                                      'weighted_gep'))
        data_scaled = data_scaled[selected_genes,]
        #-----------------------------------------------------------------------
        # 1.1 Run initial NMF runs for k selection
        #-----------------------------------------------------------------------
        data_dict = list()
        for (k in range(2,3)){
            source_python('Python/cNMF_functions.py')
            source_python('Python/cNMF_helper_functions.py')
            cNMF_result = cNMF(data_scaled, 2L,10L, 5L)
            H = cNMF_result[[1]]$H
            cluster_assignment = list()
            for (i in seq(ncol(H))){
                cluster_assignment = append(cluster_assignment,
                                            which(H[,i] == max(H[,i])))
            }
            data_dict[k] = SimpleList('model' = cNMF_result[[1]],
                                      'cophcor' = cNMF_result[[2]],
                                      'consensus' = cNMF_result[[3]],
                                      'cluster_assignments' =
                                          cluster_assignment)

        }
        #-----------------------------------------------------------------------
        # 1.2 Choose k
        #-----------------------------------------------------------------------


    }






  return(Statescope)

}

