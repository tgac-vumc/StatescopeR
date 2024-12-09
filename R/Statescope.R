#' Run Statescope
#'
#' \code{Statescope.R} Discovers states in purified resuls
#'
#' @param Statescope Statescope obj from BLADE_deconvolution
#' @param cores number of cores to use for paralellization
#'
#' @return updated Statescope S4 object with states per celltype added
#' @import cNMF_functions
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
        data_scaled = as.matrix(assay(Statescope@ct_specific_gep[[ct]],'weighted_gep'))

        ## run initial cNMF runs
        cNMF(data_scaled, 2L,100L, 1L)

    }






  return(Statescope)

}

