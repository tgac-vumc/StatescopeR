\name{Getter/setter methods}
\alias{BLADE_output}
\alias{fractions}
\alias{ct_specific_gep}
\alias{statescores}
\alias{stateloadings}

\alias{BLADE_output,Statescope-method}
\alias{fractions,Statescope-method}
\alias{ct_specific_gep,Statescope-method}
\alias{statescores,Statescope-method}
\alias{stateloadings,Statescope-method}

\alias{BLADE_output<-}
\alias{fractions<-}
\alias{ct_specific_gep<-}
\alias{statescores<-}
\alias{stateloadings<-}

\alias{BLADE_output<-,Statescope-method}
\alias{fractions<-,Statescope-method}
\alias{ct_specific_gep<-,Statescope-method}
\alias{statescores<-,Statescope-method}
\alias{stateloadings<-,Statescope-method}

\title{Statescope getters/setters}
\description{Getter/setter methods for the Statescope class.}

\usage{
\S4method{BLADE_output}{Statescope}(x)

\S4method{BLADE_output}{Statescope}(x) <- value

\S4method{fractions}{Statescope}(x)

\S4method{fractions}{Statescope}(x) <- value

\S4method{ct_specific_gep}{Statescope}(x)

\S4method{ct_specific_gep}{Statescope}(x) <- value

\S4method{statescores}{Statescope}(x)

\S4method{statescores}{Statescope}(x) <- value

\S4method{stateloadings}{Statescope}(x)

\S4method{stateloadings}{Statescope}(x) <- value
}

\arguments{
\item{x}{A Statescope object.}
\item{value}{An appropriate value to assign to the relevant slot.}
}

\details{
Any \code{value} to assign to \code{BLADE_output}, \code{ct_specific_gep}, \code{statescores} and \code{stateloadings} should be lists,
while \code{fractions} should be a DataFrame - see \linkS4class{Statescope} for details.

}

\value{
For the getter methods \code{BLADE_output}, \code{fractions}, \code{ct_specific_gep}, \code{statescores}  and \code{stateloadings}, the value of the slot with the same name is returned.
For the corresponding setter methods, a Statescope object is returned with modifications to the named slot.
}

\seealso{
\linkS4class{Statescope}
}

\examples{
Statescope <- new("Statescope",
                BLADE_output = list(NA)
            )

BLADE_output(Statescope)
BLADE_output(Statescope) <- list(NA)

fractions(Statescope)
fractions(Statescope) <- S4Vectors::DataFrame(NA)

ct_specific_gep(Statescope)
ct_specific_gep(Statescope) <- list(NA)

statescores(Statescope)
statescores(Statescope) <- list(NA)

stateloadings(Statescope)
stateloadings(Statescope) <- list(NA)
}
