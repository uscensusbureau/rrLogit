\name{anova.rrLogit}
\alias{anova.rrLogit}

\title{Comparing the Fit of Two or More Models}

\description{
Compares the fit of two or more \code{\link{rrLogit}} objects
}

\usage{
\method{anova}{rrLogit}(object, ..., method = c("lrt", "logP", "AIC", "BIC"),
   pval = FALSE, pvalDigits = 4L, showRank = NULL )
}

\arguments{
  \item{object}{an object produced by \code{rrLogit}}

  \item{\dots}{additional \code{rrLogit} objects}

  \item{method}{criterion for model comparison: \code{"lrt"} uses -2
  times the loglikelihood function for a standard likelihood-ratio test;
  \code{"logP"} uses -2 times the penalized loglikelihood or
  log-posterior density; \code{"AIC"} uses Akaike information criterion;
  and \code{"BIC"} uses Bayesian information criterion.}

  \item{pval}{if \code{TRUE} then p-values will be computed if method
  is \code{"lrt"} or \code{"logP"}.}

  \item{pvalDigits}{digits for rounding of p-values}

  \item{showRank}{if \code{TRUE}, models will be ranked from best to
  worst (with rank \code{1} being the best) according to the fit measure
  specified by \code{method}.  Defaults to \code{TRUE} if \code{method} is
  \code{"AIC"} or \code{"BIC"}}
}

\details{
The p-values reported for the \code{"lrt"} and \code{"logP"}
methods use a standard chi-squared approximation, with degrees of
freedom equal to the difference in the number of parameters for the
models being compared.  This approximation is valid only if the models
being compared are properly nested and ordered, with the simplest
model appearing first in the argument list.  The chi-squared
approximation can be poor if the degrees of freedom for the comparison
is large, and if the model corresponding to the null hypothesis (i.e.,
the smaller one) has fitted values that are too small.

The residual degrees of freedom are the difference between the number
of free parameters in a saturated model minus the number of free
parameters in the current model.  The saturated model estimates
response probabilities independently within each non-empty covariate
pattern. No adjustments are made for estimates on a boundary of the
parameter space.

For \code{"BIC"}, the sample size is taken to be the total number of
observations or total frequency in the data supplied by the user,
which does not include fictitious prior observations added when the
model is fit using \code{prior="DAP"}.
}

\value{an object of class \code{c("anova","data.frame")}}

\author{Joe Schafer \email{Joseph.L.Schafer@census.gov} }

\seealso{
\code{\link{rrLogit}}
}

\examples{
M1 <- rrLogit( Food ~ Lake + Sex + Size, freq=Freq,
    data=alligatorNarrow)
M2 <- rrLogit( Food ~ Lake + Sex + Size, freq=Freq,
    data=alligatorNarrow, saturated=TRUE)
anova(M1, M2, pval=TRUE)
}
