\name{residuals.rrLogit}
\alias{residuals.rrLogit}

\title{Extract Model Residuals}

\description{
\code{residuals.rrLogit}, which is typically invoked using the generic
name \code{residuals} or its alias \code{resid}, extracts residuals
from a model that has been fit with \code{\link{rrLogit}}. 
}

\usage{
\method{residuals}{rrLogit}(object, type=c("pearson", "response"),
    dispersion = 1.0, covPatt = FALSE, include.predvars=FALSE, digits=NULL, ...)
}

\arguments{
  \item{object}{an object produced by \code{rrLogit}}

  \item{type}{type of residuals: \code{"pearson"} for Pearson
  residuals, \code{"response"} for raw differences between observed
  frequencies and their estimated means.}

  \item{dispersion}{dispersion parameter for adjusting Pearson
  residuals to account for overdispersion or (more rarely)
  underdispersion. Ignored if
  \code{object} was generated with \code{method="MCMC"} or
  \code{method="approxBayes"}.}

  \item{covPatt}{if \code{TRUE}, the returned matrix of residuals will
  be given an attribute \code{"covPatt"} that maps the rows of the
  original dataset to the covariate patterns.}

  \item{include.predvars}{if \code{TRUE}, all predictors variables
  appearing the model will be included in the result.}

  \item{digits}{number of decimal places to which fitted values will be
  rounded; default is no rounding.}

  \item{\dots}{additional arguments passed to or from other methods}
}

\details{
Residuals for an \code{rrLogit} model are computed within covariate
patterns, the distinct combinations of values of model predictors
occurring in a dataset.  For more discussion of covariate patterns, see
\code{\link[=fitted.rrLogit]{fitted}}.

Within a covariate pattern, residuals are computed with respect to the
noisy response variable \eqn{Y^*} for each response category
\eqn{1,\ldots,C}. A raw residual (\code{type="response"}) is the
observed frequency for a given response category minus the estimated
expected value for that frequency (its fitted mean) under the model. A
Pearson residual (\code{type="pearson"}) is the raw residual divided
by the square root of the estimated expected value, with \code{0/0}
taken to be \code{0}. The Pearson goodness-of-fit statistic can be
computed by summing the squared Pearson residuals over the covariate
patterns and response categories.

If \code{object} was produced with \code{method} equal to \code{"EM"},
\code{"NR"}, or \code{"FS"}, estimated expected values are computed at
the final parameter estimates. If code{object} was produced with
\code{method} equal to \code{"MCMC"} or \code{"approxBayes"},
estimated expected values are averaged over all iterations after the
burn-in period.

If a dispersion parameter is provided, Pearson residuals are divided
by the square root of this parameter to adjust for overdispersion or
underdispersion.

}

\note{
When an \code{rrLogit} model is saturated, raw residuals are not
necessarily all equal to zero, because the observed sample proportions
within a covariate pattern might fall outside of the parameter space
implied by the perturbation matrix. When this happens, some Pearson
residuals may be undefined.
}

\value{
If \code{include.predvars=FALSE},
a matrix of residuals having one row for each covariate
pattern and one column for each response category.

If \code{include.predvars=TRUE}, the result will be a data frame that
includes all of the predictor variables appearing in
the model along with the residuals,
which makes it easier to identify the covariate patterns.

If \code{covPatt=TRUE}, the result will have an attribute 
\code{"covPatt"} that maps the rows of the original dataset to the covariate
patterns.
}

\author{Joe Schafer \email{Joseph.L.Schafer@census.gov} }

\seealso{
\code{\link{rrLogit}},
\code{\link[=fitted.rrLogit]{fitted}}
}

\examples{
modelFit <- rrLogit( Food ~ Lake + Sex + Size, freq=Freq,
    data=alligatorNarrow, saturated=TRUE, privLoss=1)
resids <- residuals(modelFit)
}
