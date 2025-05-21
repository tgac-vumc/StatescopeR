## Define classes for use in package

#' S4 class to store Statescope results
#
#' @slot BLADE_output list with all BLADE output
#' @slot fractions predicted fraction per sample per cell type
#' @slot ct_specific_gep predicted cell type-specific gene expression profiles
#' @slot statescores cell type-specific state scores based on ct_specific_gep
#' per sample
#' @importFrom S4Vectors DataFrame
#' @importFrom methods setClass setGeneric setMethod
#' @importFrom utils globalVariables
#' @export
methods::setClass("Statescope",
    slots = c(
        BLADE_output = "list",
        fractions = "DataFrame",
        ct_specific_gep = "list",
        statescores = "list"
    ),
    prototype = list(
        BLADE_output = list(NA),
        fractions = S4Vectors::DataFrame(NA),
        ct_specific_gep = list(NA),
        statescores = list(NA)
    )
)


## Make Python and dplyr (.) functions global for check
utils::globalVariables(c(
    "Framework_Iterative", "Purify_AllGenes", "biggest_drop",
    "cNMF", "find_threshold", "select_k", "."
))
