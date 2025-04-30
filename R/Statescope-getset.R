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
    return(fractions)
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
setMethod("states", "Statescope", function(x) {
    output <- x@states
    return(output)
})

#' @importFrom methods validObject
#' @export
setReplaceMethod("states", "Statescope", function(x, value) {
    x@states <- value
    validObject(x)
    return(x)
})
