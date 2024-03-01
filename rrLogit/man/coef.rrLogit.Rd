\name{coef.rrLogit}
\alias{coef.rrLogit}

\title{
Model Coefficients from an \code{rrLogit} Fit
}

\description{
\code{coef.rrLogit}, which is typically invoked using the generic name
\code{coef} or its alias
\code{coefficients}, extracts the estimated coefficients from a model
that has been fit with \code{\link{rrLogit}}.
}

\usage{
\method{coef}{rrLogit}(object, covMat = FALSE, meanSeries = TRUE,
  completeData = FALSE, ...)
}

\arguments{
  \item{object}{an object of class \code{"rrLogit"} returned by
  the function \code{rrLogit}}
  
  \item{covMat}{if \code{TRUE}, an estimated covariance matrix
  is returned along with the estimated coefficients; see DETAILS}

  \item{meanSeries}{applies when \code{object} is the result from a
  simulation run. If \code{TRUE}, results will be based on from a
  running average of simulated coefficients over all iterations after
  the burn-in period. If \code{FALSE}, results will be based only on
  the simulated parameter values at the end of the run}

  \item{completeData}{applies when \code{object} was produced with
  \code{method="EM"}. If \code{TRUE}, the covariance matrix returned
  when \code{covMat=TRUE} is based on the expected complete-data
  loglikelihood, and it estimates the variability in the model
  coefficients that would have resulted if the true response variable
  were known (i.e., if no noise had been added); see DETAILS}

  \item{\dots}{additional arguments passed to or from other methods}
}

\details{
If \code{object} was produced with \code{method} equal to \code{"EM"},
\code{"NR"}, or \code{"FS"}, the covariance matrix is computed as the
inverse of (minus one times) the actual second-derivative matrix of the
loglikelihood function (or, if the model was fit using
\code{prior="DAP"}, the log-posterior density function), evaluated at
the estimated parameters.

If \code{completeData=TRUE}, then the covariance matrix is based not
on the actual loglikelihood function, but on the expected
complete-data loglikelihood function being maximized in the final
M-step of EM. By comparing this covariance matrix to the one resulting
from \code{completeData=FALSE}, we can compute width inflation factors
and fractions of missing information for the model parameters due to
noise infusion; see Chen and Savalei (2021).

If \code{object} was produced with \code{method} equal to \code{"MCMC"}
or \code{"approxBayes"}, the covariance matrix is computed as the
sample covariance matrix of the simulated coefficients over all
iterations after the burn-in period.
}

\value{
If \code{covMat=FALSE}, a matrix of estimated coefficients with
\eqn{p} rows and \eqn{C} columns, where \eqn{p} is the number of
covariates (i.e., the number of columns appearing in the
\code{\link[=model.matrix.rrLogit]{model.matrix}}) and \eqn{C} is the
number of response categories.  The column of coefficients correspording to the
baseline category is present but filled with zeros.

If \code{covMat=TRUE}, a list with two components:
  \item{coefVec }{a vector of length \eqn{p(C-1)} obtained by stacking
  the columns of the coefficient matrix, with the column of zeros removed} 
  \item{covMat }{the estimated covariance matrix for \code{coefVec}}
}

\note{
To obtain a table of estimated coefficients with standard
errors, z-statistics and p-values, use the
\code{\link[=summary.rrLogit]{summary}} method.
}

\references{
Chen, L. and Savalei, V. (2021) Three sample estimates of fraction of
missing information from full information maximum likelihood,
\emph{Frontiers in Psychology}, \bold{12}, 667802.
}

\author{Joe Schafer \email{Joseph.L.Schafer@census.gov} }

\seealso{
\code{\link{rrLogit}},
\code{\link[=summary.rrLogit]{summary}}
}

\examples{
fitA <- rrLogit( cbind(Fish,Inv,Rept,Bird,Other) ~ Lake + Sex + Size,
   data=alligatorWide)
coef(fitA)
}
