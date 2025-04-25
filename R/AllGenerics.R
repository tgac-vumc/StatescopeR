## Getter & setter methods for Statescope class

#' @export
setGeneric("BLADE_output", function(x) standardGeneric("BLADE_output"))

#' @export
setGeneric("BLADE_output<-", function(x, value) standardGeneric("BLADE_output<-"))

#' @export
setGeneric("fractions", function(x) standardGeneric("fractions"))

#' @export
setGeneric("fractions<-", function(x, value) standardGeneric("fractions<-"))

#' @export
setGeneric("ct_specific_gep", function(x) standardGeneric("ct_specific_gep"))

#' @export
setGeneric("ct_specific_gep<-", function(x, value) standardGeneric("ct_specific_gep<-"))

#' @export
setGeneric("states", function(x) standardGeneric("states"))

#' @export
setGeneric("states<-", function(x, value) standardGeneric("states<-"))

## subset operators

#' @export
getGeneric("[")

#' @export
getGeneric("[[")
