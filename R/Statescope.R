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
  ## get final BLADE object from result ready for Python code
  BLADE_obj = list('final_obj' = Statescope@BLADE_output[[1]],
                   'outs' = Statescope@BLADE_output[[4]])



  return(Statescope)

}

