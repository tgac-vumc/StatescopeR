## Define classes for use in package

#' S4 class to store Statescope results
#
#' @slot BLADE_output list with all BLADE output
#' @slot fractions predicted fraction per sample per cell type
#' @slot ct_specific_gep predicted cell type-specific gene expression profiles
#' @slot states cell type-specific states based on ct_specific_gep
#' per sample
#' @importFrom S4Vectors DataFrame
#' @export
setClass("Statescope",
         slots = c(
           BLADE_output = "list",
           fractions = 'DataFrame',
           ct_specific_gep = 'list',
           states = 'list'),
         prototype = list(
           BLADE_output = list(NA),
           fractions = S4Vectors::DataFrame(NA),
           ct_specific_gep = list(NA),
           states = list(NA)
         )
         )
