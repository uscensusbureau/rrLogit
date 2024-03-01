\name{coefSeries}
\alias{coefSeries}
\alias{impList}
\alias{loglik}
\alias{logP}
\alias{margProbSeries}
\alias{minus2logPSeries}
\alias{df.residual.rrLogit}

\title{Extract Information from an \code{rrLogit} Fit}

\description{This group of functions will extract various
summaries from a model fit by \code{\link{rrLogit}}.
}

\usage{

coefSeries(obj, coda = TRUE)

\method{df.residual}{rrLogit}(object, ...)

impList(obj)

loglik(obj)

logP(obj)

margProbSeries(obj, noisy = FALSE, coda = TRUE)

minus2logPSeries(obj, startValShift = TRUE, coda = TRUE)

}


\arguments{

\item{obj}{an object resulting from a call to \code{\link{rrLogit}}}

\item{coda}{if \code{TRUE}, the series from an MCMC run is returned as
an \code{mcmc} object for plotting and diagnostic analysis with the
\code{coda} package. If \code{FALSE}, a one-dimensional series is
returned as a numeric vector, and a multidimensional series is
returned as a numeric matrix with rows corresponding to iterations and
columns corresponding to elements of the multidimensional quantities
being monitored.}

\item{object}{an object resulting from a call to \code{\link{rrLogit}}}

\item{noisy}{if \code{FALSE}, marginal probabilities correspond to
categories of the true (noise-free) categorical response variable; if
\code{TRUE}, for the observed noise-infused response.}

\item{startValShift}{the function \code{get.minus2logPSeries} extracts a
saved series from an MCMC run containing  the values of (minus 2 times)
the log-posterior density function. If \code{startValShift} is true,
the series is shifted by (minus 2 times) the log-posterior density at
the starting value, if the starting value appears to be a mode.}

\item{\ldots}{additional arguments passed to or from other methods}
}

\details{

The series objects returned by \code{coefSeries} and
\code{minus2logPSeries} omit results from the burn-in period, if any,
and may also be thinned. The default behavior is no burn-in period and
no thinning. The burn-in period and thinning interval are set by
control parameters provided to \code{rrLogit} through the the
\code{control} argument, which are handled by the function
\code{\link{rrLogitControl}}; the relevant control parameters are
\code{control$burnMCMC} and \code{control$thinMCMC}.

The number of imputations returned by \code{impList} is determined by
the control parameters \code{control$iterMCMC},
\code{control$imputeEvery}, \code{control$iterApproxBayes}, and
\code{control$imputeApproxBayes}. 

Other useful information from a model fit can be extracted
with the \code{summary} method for an \code{rrLogit} object, and with the
functions listed below in SEE ALSO.

}

\value{

\code{coefSeries} returns a series of logistic coefficients created during a
run of \code{rrLogit} with \code{method="MCMC"} or
\code{method="approxBayes"}, either
as a numeric matrix or as an \code{mcmc} object used by the
\code{coda} package.

\code{df.residual} returns the degrees of freedom, the number of free
parameters separating the current model from its saturated counterpart.

\code{loglik} and \code{get.logP} return the value of the
loglikelihood function or log-posterior density at the final value of
the model parameters held in the \code{rrLogit} object.

\code{impList} returns a list holding multiple imputations of the true
(noise-free) response created during a run with \code{method="MCMC"}
or \code{method="approxBayes"}. If the dataset originally supplied to
\code{rrLogit} represents microdata, each component of the list is a
single factor variable whose length is equal to the number of rows in
the dataset. Otherwise, each component is a matrix of frequencies,
with rows corresponding to the rows of the original dataset and
columns corresponding to response categories.

\code{margProbSeries} returns a series of marginal probabilities from
the iterations of MCMC for the categories of the true (noise-free)
response if \code{noisy=FALSE}, or the observed response if
\code{noisy=TRUE}.

\code{minus2logPSeries} returns a series of (minus 2 times) the
log-posterior density values from the iterations of MCMC, either as a
numeric vector or as an \code{mcmc} object used by the \code{coda}
package. 

}

\author{Joe Schafer \email{Joseph.L.Schafer@census.gov} }

\seealso{
\code{\link[=anova.rrLogit]{anova}},
\code{\link[=coef.rrLogit]{coef}},
\code{\link[=fitted.rrLogit]{fitted}},
\code{\link[=formula.rrLogit]{formula}},
\code{\link[=impute.rrLogit]{impute}},
\code{\link[=model.frame.rrLogit]{model.frame}},
\code{\link[=model.matrix.rrLogit]{model.matrix}},
\code{\link[=predict.rrLogit]{predict}},
\code{\link[=residuals.rrLogit]{residuals}},
\code{\link[=summary.rrLogit]{summary}},
\code{\link[=terms.rrLogit]{terms}}
}
