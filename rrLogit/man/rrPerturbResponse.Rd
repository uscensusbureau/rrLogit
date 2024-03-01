\name{rrPerturbResponse}
\alias{rrPerturbResponse}

\title{
Inject Noise into a Categorical Variable
}

\description{
Distorts a categorical (factor) variable by injecting random noise
according to a given perturbation matrix
}

\usage{
rrPerturbResponse( y, pertMat )
}

\arguments{
  \item{y}{factor variable to be distorted with random noise. This
  variable is assumed to contain microdata, with each element
  corresponding to a single individual or observational unit.}

  \item{pertMat}{perturbation matrix}
}

\details{
  \code{pertMat} must be a square stochastic matrix of order
  \code{nlevels(y)}. A stochastic matrix is a matrix whose elements
  are probabilities and whose columns sum to one.
}

\value{a noise-infused factor variable having the same length, same
\code{levels} and other attributes as \code{y}
}

\author{Joe Schafer \email{Joseph.L.Schafer@census.gov} }

\seealso{
\code{\link{rrPertMat}}
}

\examples{
pertMat <- rrPertMat( privLoss=1.0, nLevels=5 )
alligatorMicro$noisyFood <- rrPerturbResponse( alligatorMicro$Food,
   pertMat )
}
