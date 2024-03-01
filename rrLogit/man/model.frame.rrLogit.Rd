\name{model.frame.rrLogit}
\alias{model.frame.rrLogit}
\alias{model.matrix.rrLogit}
\alias{formula.rrLogit}
\alias{terms.rrLogit}

\title{
Model-Related Objects from an \code{rrLogit} Fit
}

\description{
These methods, which are typically invoked using the
generic function names \code{\link{model.frame}}, \code{\link{model.matrix}},
\code{\link{formula}} and \code{\link{terms}}, extract various objects
from a model that has been fit with \code{\link{rrLogit}}.
}

\usage{
\method{model.frame}{rrLogit}(formula, ...)
\method{model.matrix}{rrLogit}(object, ...)
\method{formula}{rrLogit}(x, ...)
\method{terms}{rrLogit}(x, ...)
}

\arguments{
  \item{formula}{an object of class \code{"rrLogit"} returned by
  the function \code{rrLogit}}
  
  \item{object}{an object of class \code{"rrLogit"} returned by
  the function \code{rrLogit}}
  
  \item{x}{an object of class \code{"rrLogit"} returned by
  the function \code{rrLogit}}
  
  \item{\dots}{additional arguments passed to or from other methods}
}

\value{
\code{model.frame} returns the model frame, a data frame that contains
all the variables appearing in the model formula, along with
frequencies or weights, if present.

\code{model.matrix} returns the model (design) matrix used to fit the
model.

\code{formula} returns the model formula.

\code{terms} returns the model's \code{\link[=terms.object]{terms object}}.
}

\references{
Chambers, J. M. and Hastie, T. J. (1992) Statistical models. Chapter 2
of \emph{Statistical Models in S}, eds. J. M. Chambers and T. J. Hastie,
Wadsworth & Brooks/Cole.

Chambers, J. M.(1992) Data for models. Chapter 3 of \emph{Statistical
Models in S}, eds. J. M. Chambers and T. J. Hastie, Wadsworth & Brooks/Cole.
}

\author{Joe Schafer \email{Joseph.L.Schafer@census.gov} }

\seealso{
\code{\link{rrLogit}}
}

\examples{
myFit <- rrLogit( Food ~ Lake + Sex + Size, freq=Freq, data=alligatorNarrow)
model.frame(myFit)
model.matrix(myFit)
formula(myFit)
terms(myFit)
}
