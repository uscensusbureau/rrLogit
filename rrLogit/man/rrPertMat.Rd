\name{rrPertMat}
\alias{rrPertMat}

\title{
Perturbation Matrix for Noise Infusion
}

\description{
Computes the perturbation or distortion matrix for infusing a
categorical variable (factor) with random noise while maintaining a
given level of privacy loss
}

\usage{
rrPertMat( privLoss=Inf, nLevels=2L, method=c("WWH") )
}

\arguments{
  \item{privLoss}{privacy-loss parameter}

  \item{nLevels}{number of response categories}

  \item{method}{currently the only option is "WWH" for the method of
  Wang, Wu and Hu (2016); others may be added in the future} 
}

\details{
  A perturbation or distortion matrix describes a mechanism for
  infusing random noise into a categorical variable. It is a full-rank
  square matrix of probabilities of order \eqn{C}, where \eqn{C} is
  the number of response categories, with each column summing to
  one. Element \eqn{(r,s)} of this matrix is the conditional
  probability that the noise-infused response is equal to \eqn{r},
  given that the true response is \eqn{s}.

  Wang, Wu and Hu (2016) present a method for noise infusion within a
  framework of differential privacy. Their perturbation matrix is symmetric,
  with diagonal elements \eqn{p} and off-diagonal elements
  \eqn{(1-p)/(C-1)}, where \deqn{p=( 1 + (C-1) e^\epsilon)^{-1}}{p = 1
  / ( 1 + (C-1) * exp( -\epsilon) ) } and \eqn{\epsilon > 0} is the
  privacy-loss parameter.  
}

\value{a square numeric matrix of order \code{nLevels}
}

\references{
Wang, Y., Wu, X. and Hu, D. (2016) Using randomized response for
differential privacy preserving data collection. \emph{Proceedings of
EDBT/ICDT Workshops}, \bold{1558},
\url{http://ceur-ws.org/Vol-1558/paper35.pdf}.
}

\author{Joe Schafer \email{Joseph.L.Schafer@census.gov} }

\seealso{
\code{\link{rrPerturbResponse}}
}

\examples{
rrPertMat(1.0, 5)
}
