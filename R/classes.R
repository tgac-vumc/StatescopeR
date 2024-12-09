## Define classes for use in package

#' S4 class to store Statescope results
#
#' @slot BLADE_output list with all BLADE output
#' @slot fractions predicted fraction per sample per cell type
#' @slot ct_specific_gep predicted cell type-specific gene expression profiles
#' per sample
#' @import S4Vectors
#' @export
setClass("Statescope",
         slots = c(
           BLADE_output = "list",
           fractions = 'DataFrame',
           ct_specific_gep = 'SimpleList'),
         prototype = list(
           BLADE_output = list(NA),
           fractions = DataFrame(NA),
           ct_specific_gep = SimpleList(NA)
         )
         )
