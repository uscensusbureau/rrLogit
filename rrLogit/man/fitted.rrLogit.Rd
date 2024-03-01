\name{fitted.rrLogit}
\alias{fitted.rrLogit}

\title{Extract Model Fitted Values}

\description{
\code{fitted.rrLogit}, which is typically invoked using the generic
name \code{fitted} or its alias \code{fitted.values}, extracts the
fitted values from a model that has been fit with
\code{\link{rrLogit}}. Fitted values are computed within covariate
patterns, the unique combinations of predictor values appearing in the
dataset.
}

\usage{
\method{fitted}{rrLogit}(object, type=c("prob", "link", "mean"),
   noisy=FALSE, covPatt=FALSE, include.predvars=FALSE, digits=NULL, ...)
}

\arguments{
  \item{object}{an object produced by \code{rrLogit}}

  \item{type}{type of fitted values to be extracted: \code{"prob"}
  for probabilities, \code{"link"} for linear predictors (logits),
  \code{"mean"} for expected frequencies}

  \item{noisy}{if \code{TRUE}, fitted values will describe
  the observed, noise-infused response variable; if \code{FALSE},
  to the unobserved ``true'' or noise-free response}

  \item{covPatt}{if \code{TRUE}, the returned matrix of fitted values will
  be given an attribute \code{"covPatt"} that maps the rows of the
  original dataset to the covariate patterns.}

  \item{include.predvars}{if \code{TRUE}, all predictors variables
  appearing the model will be included in the result.}

  \item{digits}{number of decimal places to which fitted values will be
  rounded; default is no rounding.}

  \item{\dots}{additional arguments passed to or from other methods}
}

\details{
A baseline-category logistic regression model describes
relationships between a categorical variable \eqn{Y}, taking values
\eqn{1,\dots,C}, and a set of predictors encoded into a covariate
vector \eqn{x}.  The model is \deqn{log \left( \frac{P(Y=y)}{P(Y=b)}
\right) = x' \beta_y}{ log( P(Y=y) / P(Y=b) ) = x' \beta_y} for
\eqn{y=1,\ldots,C}, where \eqn{\beta_1,\ldots,\beta_C} are vectors of
coefficients to be estimated, \eqn{b} is a response category
selected to serve as the baseline, and the coefficients associated with
the baseline category \eqn{\beta_b} are constrained to be
zero. \code{rrLogit} is designed to fit this model in situations where
\eqn{Y} is withheld from the data analyst and replaced by a perturbed
or noise-infused version \eqn{Y^*}.  Parameter estimates
computed by \code{rrLogit} based on the noisy response \eqn{Y^*} are
still consistent, but in general, their standard errors will be larger than
they would have been if the true response \eqn{Y} had been used.

\emph{Covariate patterns.}  Before fitting a model, \code{rrLogit}
processes the supplied dataset, grouping together rows that are
identical with respect to the model predictors, and then computing
frequencies for the response categories \eqn{1,\ldots,C} within each
group. This grouping of rows into covariate patterns happens whether
the data are supplied to \code{rrLogit} in wide format, narrow format,
or as microdata. Fitted values are returned as a matrix with
one row per covariate pattern.  This matrix has an attribute named
\code{"covPatt"}, a vector of integers whose length is equal to the
number of rows in the original dataset, which maps the rows to the
covariate patterns.

\emph{Fitted value types.}  Fitted values may be requested either for
\eqn{Y} (\code{noisy=FALSE}) or for \eqn{Y^*} (\code{noisy=TRUE}).  When
\code{type="prob"}, fitted values for \eqn{Y} are computed as
\deqn{P(Y=c) = \frac{\exp (x'\beta_c)}{\sum_{d=1}^C\exp (x'\beta_d)}}{
P(Y=c) = exp(x'\beta_c) / \sum_d exp(x'\beta_d)}
for \eqn{c=1,\ldots,C}, and fitted values for \eqn{Y^*} are obtained by
pre-multiplying the vector of fitted values for \eqn{Y} by the
perturbation matrix.  When \code{type="mean"}, the fitted values are
obtained by multiplying the fitted probabilities for \eqn{Y} or
\eqn{Y^*} in each covariate pattern by the total frequency for the
covariate pattern.  When \code{type="link"}, fitted values for
\eqn{Y} are computed as \eqn{x'\beta_c} for \eqn{c=1,\ldots,C}, and
fitted values for \eqn{Y^*} are undefined, because the relationship
between \eqn{P(Y^*=c)} and \eqn{x} does not have a
logistic-linear representation.

\emph{Fitting method.} If \code{object} was produced with
\code{method} equal to \code{"EM"}, \code{"NR"}, or \code{"FS"},
fitted values are computed at the final parameter estimates. If
code{object} was produced with \code{method} equal to \code{"MCMC"} or
\code{"approxBayes"}, fitted values are averaged over all iterations
after the burn-in period.
}

\value{
If \code{include.predvars=FALSE}, a matrix of fitted values having one row
for each unique covariate pattern  
appearing in the data and one column for each response category.

If \code{include.predvars=TRUE}, the result will be a data frame that
includes all of the predictor variables appearing in
the model along with the fitted values,
which makes it easier to interpret the covariate patterns.

If \code{covPatt=TRUE}, the result will have an attribute 
"covPatt" that maps the rows of the original dataset to the covariate patterns.
}

\author{Joe Schafer \email{Joseph.L.Schafer@census.gov} }

\note{
If the model was fit with the option \code{saturated=TRUE}, then some
rows of the matrix returned by \code{fitted} may be
filled with \code{NA}. These rows correspond to covariate patterns in
the original dataset that are empty in the sense that all
reported frequencies within the pattern are zero.

To get fitted values that correspond to actual rows of the dataset,
rather than its unique covariate patterns, use
\code{\link[=predict.rrLogit]{predict}}.
}

\seealso{
\code{\link{rrLogit}},
\code{\link[=predict.rrLogit]{predict}},
\code{\link[=residuals.rrLogit]{residuals}}
}

\examples{
modelFit <- rrLogit( Food ~ Lake + Sex + Size, freq=Freq,
    data=alligatorNarrow)
fittedVal <- fitted(modelFit)
}
