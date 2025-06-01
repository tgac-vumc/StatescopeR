# Defines methods for new generics for Statescope class

#' @export
setMethod("BLADE_output", "Statescope", function(x) {
    output <- x@BLADE_output
    return(output)
})

#' @importFrom methods validObject
#' @export
setReplaceMethod("BLADE_output", "Statescope", function(x, value) {
    x@BLADE_output <- value
    validObject(x)
    return(x)
})


#' @export
setMethod("fractions", "Statescope", function(x) {
    output <- x@fractions
    return(output)
})

#' @importFrom methods validObject
#' @export
setReplaceMethod("fractions", "Statescope", function(x, value) {
    x@fractions <- value
    validObject(x)
    return(x)
})


#' @export
setMethod("ct_specific_gep", "Statescope", function(x) {
    output <- x@ct_specific_gep
    return(output)
})

#' @importFrom methods validObject
#' @export
setReplaceMethod("ct_specific_gep", "Statescope", function(x, value) {
    x@ct_specific_gep <- value
    validObject(x)
    return(x)
})


#' @export
setMethod("statescores", "Statescope", function(x) {
    output <- x@statescores
    return(output)
})

#' @importFrom methods validObject
#' @export
setReplaceMethod("statescores", "Statescope", function(x, value) {
    x@statescores <- value
    validObject(x)
    return(x)
})


#' @export
setMethod("stateloadings", "Statescope", function(x) {
    output <- x@stateloadings
    return(output)
})

#' @importFrom methods validObject
#' @export
setReplaceMethod("stateloadings", "Statescope", function(x, value) {
    x@stateloadings <- value
    validObject(x)
    return(x)
})
