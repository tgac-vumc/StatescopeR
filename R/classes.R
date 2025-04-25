## Define classes for use in package

#' S4 class to store Statescope results
#
#' @slot BLADE_output list with all BLADE output
#' @slot fractions predicted fraction per sample per cell type
#' @slot ct_specific_gep predicted cell type-specific gene expression profiles
#' @slot states cell type-specific states based on ct_specific_gep
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
        states = "list"
    ),
    prototype = list(
        BLADE_output = list(NA),
        fractions = S4Vectors::DataFrame(NA),
        ct_specific_gep = list(NA),
        states = list(NA)
    )
)
## Set Accessor & Setter functions
## BLADE_output
#' @export
methods::setGeneric(
    "BLADE_output",
    function(x) standardGeneric("BLADE_output")
)
methods::setGeneric(
    "BLADE_output<-",
    function(x, value) standardGeneric("BLADE_output<-")
)

methods::setMethod("BLADE_output", "Statescope", function(x) x@BLADE_output)

methods::setMethod("BLADE_output<-", "Statescope", function(x, value) {
    x@BLADE_output <- value
    methods::validObject(x)
    x
})

## fractions
#' @export
methods::setGeneric(
    "fractions",
    function(x) standardGeneric("fractions")
)
methods::setGeneric(
    "fractions<-",
    function(x, value) standardGeneric("fractions<-")
)

methods::setMethod("fractions", "Statescope", function(x) x@fractions)
methods::setMethod("fractions<-", "Statescope", function(x, value) {
    x@fractions <- value
    methods::validObject(x)
    x
})

## ct_specific_gep
#' @export
methods::setGeneric(
    "ct_specific_gep",
    function(x) standardGeneric("ct_specific_gep")
)
methods::setGeneric(
    "ct_specific_gep<-",
    function(x, value) standardGeneric("ct_specific_gep<-")
)

methods::setMethod(
    "ct_specific_gep", "Statescope",
    function(x) x@ct_specific_gep
)
methods::setMethod("ct_specific_gep<-", "Statescope", function(x, value) {
    x@ct_specific_gep <- value
    methods::validObject(x)
    x
})

## states
#' @export
methods::setGeneric(
    "states",
    function(x) standardGeneric("states")
)
methods::setGeneric("states<-", function(x, value) standardGeneric("states<-"))

methods::setMethod("states", "Statescope", function(x) x@states)
methods::setMethod("states<-", "Statescope", function(x, value) {
    x@states <- value
    methods::validObject(x)
    x
})

## Make Python and dplyr (.) functions global for check
utils::globalVariables(c(
    "Framework_Iterative", "Purify_AllGenes",
    "biggest_drop", "cNMF", "find_threshold", "."
))
