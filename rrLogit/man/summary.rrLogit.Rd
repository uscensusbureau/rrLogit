\name{summary.rrLogit}
\alias{summary.rrLogit}
\alias{print.summary.rrLogit}
\alias{print.rrLogit}


\title{
Summarize and Print Results from an \code{rrLogit} Fit
}

\description{
Summarizes the results from a model that has been fit with
\code{\link{rrLogit}}
}

\usage{
\method{summary}{rrLogit}(object, showCoef=FALSE, digits=4L,
   dispersion=1.0, ...)
\method{print}{summary.rrLogit}(x, digits=x$digits, ...)
\method{print}{rrLogit}(x, ...)
}

\arguments{
  \item{object}{an object of class \code{"rrLogit"} returned by
  the function \code{rrLogit}}

  \item{showCoef}{if \code{TRUE}, the table of coefficients will be
  displayed when the \code{print} method is invoked.}

  \item{digits}{for printing}

  \item{dispersion}{dispersion parameter used to adjust standard errors
  to account for overdispersion or (more rarely) underdispersion. Ignored if
  \code{object} was generated with \code{method="MCMC"} or
  \code{method="approxBayes"}.}

  \item{\dots}{additional arguments passed to or from other methods}

  \item{x}{an object of class \code{"rrLogit"} or \code{"summary.rrLogit"}}
}

\value{
for the \code{summary} method, an object of class
\code{"summary.rrLogit"} with at least the following components:
   \item{coefArray }{a table displaying estimated coefficients,
   standard errors, z-statistics and p-values. (To obtain the estimated
   coefficients and their covariance matrix, use
   \code{\link[=coef.rrLogit]{coef}}.)}
   \item{iter }{number of iterations of the fitting method
   (\code{"EM"}, \code{"NR"} or \code{"FS"}) performed}
   \item{converged }{TRUE if the estimation procedure converged,
   otherwise FALSE.  To adjust the convergence criterion, see
   \code{\link{rrLogitControl}}.}
   \item{boundary}{TRUE if the estimated parameters lie on or near a
   boundary of the parameter space, otherwise FALSE. To adjust the
   criterion for judging proximity to a boundary, see
   \code{\link{rrLogitControl}}.}
   \item{loglik}{final value of the loglikelihood function}
   \item{logprior}{final value of the log-prior density function, if
   the model was fit using \code{prior="DAP"}, or \code{0} if
   \code{prior="none"}}
}

\note{
If \code{obj} is an object of class \code{"rrLogit"}, then
\code{print(obj)} is equivalent to \code{print(summary(obj))}.
}

\author{Joe Schafer \email{Joseph.L.Schafer@census.gov} }

\seealso{
\code{\link{rrLogit}}
}

\examples{
myFit <- rrLogit( Food ~ Lake + Sex + Size, freq=Freq, data=alligatorNarrow)
summary(myFit)
}
