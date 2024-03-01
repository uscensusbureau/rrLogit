\name{predict.rrLogit}
\alias{predict.rrLogit}

\title{Predictions from \code{rrLogit} Model Fit}

\description{
\code{predict.rrLogit}, which is typically invoked using the generic
name \code{predict}, computes predicted values for a multinomial
response variable based on the estimated regression function contained
in an \code{\link{rrLogit}} object. The values of predictor variables
may come from the data frame used to fit the model or another data
frame. Predictions from \code{predict} are similar to fitted
values from \code{\link[=fitted.rrLogit]{fitted}}, except that they
are computed for the rows of the given data frame rather than its
unique covariate patterns. 
}

\usage{
\method{predict}{rrLogit}(object, newdata=NULL, freq=NULL,
   type=c("prob", "link", "mean"), noisy=FALSE, se.fit=FALSE,
   na.action=na.pass, completeData=FALSE, ...)
}

\arguments{
  \item{object}{an object produced by \code{rrLogit}}

  \item{newdata}{an optional data frame where \code{predict} looks for
  predictor variables to generate the predictions. If \code{NULL},
  predictions will be based on the data used to fit the model
  contained in \code{object}.}

  \item{freq}{an optional numeric variable in \code{newdata}
  containing frequencies. If \code{newdata=NULL}, \code{freq} is
  ignored and frequencies are taken from the data frame used to fit
  the model.}

  \item{type}{type of predictions: \code{"prob"}
  for probabilities, \code{"link"} for linear predictors (logits),
  \code{"mean"} for expected frequencies}

  \item{noisy}{if \code{TRUE}, predictions will correspond
  the noise-infused response variable; if \code{FALSE},
  to the unobserved ``true'' or noise-free response}

  \item{se.fit}{if \code{TRUE}, standard errors for the predictions
  will be computed, along with an array of covariance matrices; see
  DETAILS}

  \item{na.action}{action to be taken if predictors contain missing
  values. For now, \code{NA} values are not allowed, but this may change
  in future versions.}

  \item{completeData}{applies when \code{object} was produced with
  \code{method="EM"}. If \code{TRUE}, then the standard errors
  returned when \code{se.fit=TRUE} will be based on the expected 
  complete-data loglikelihood, and will approximate the standard
  errors that would have been seen if the response variable had no
  added noise. These standard errors will be larger. Dividing the
  standard errors resulting from \code{completeData=TRUE} by those
  resulting from \code{completeData=FALSE} will give the width
  inflation factors, a diagnostic that describes the impact of
  noise infusion on measures of uncertainty; see Chen and Savalei (2021).}

  \item{\dots}{additional arguments passed to or from other methods}
}

\details{
For an explanation of the types of predicted values and how they are
computed, see \code{\link[=fitted.rrLogit]{fitted}}.

If \code{se.fit=TRUE}, standard errors and covariance matrices are
computed by the delta method, using the fact that each row of
predicted values is a known function of the logistic
coefficients. Standard errors are not available if the model was fit
with \code{saturated=TRUE}.

\emph{Fitting method.} If \code{object} was produced with
\code{method} equal to \code{"EM"}, \code{"NR"}, or \code{"FS"},
predicted values are computed at the final parameter estimates. If
code{object} was produced with \code{method} equal to \code{"MCMC"} or
\code{"approxBayes"}, predicted values are based on simulated
parameters averaged over all iterations after the burn-in period, and
standard errors are based on the sample covariance matrix of these
simulated parameters.
}

\value{
if \code{se.fit=FALSE}, a matrix of predicted values having one row for
each row in the dataset (the original data frame or \code{newdata})
and one column for each response category

if \code{se.fit=TRUE}, a list with three components:
   \item{fit}{the matrix of predicted values}
   \item{se.fit}{a matrix of the same dimensions as \code{fit} containing
   standard errors for the predicted values}
   \item{cov.fit.array}{an array of estimated covariance matrices for
   the predicted values. Dimensions of this array are \code{c(n,C,C)},
   where \code{n} is the number or rows in the dataset (the original
   data frame or \code{newdata}) and \code{C} is the number of
   response categories. Layer \code{cov.fit.array[r,,]} contains the
   estimated covariance matrix for \code{fit[r,]}. Each of these
   covariance matrices is rank-deficient, because the response
   probabilities for each row of data are constrained to sum to one.}
}

\references{
Chen, L. and Savalei, V. (2021) Three sample estimates of fraction of
missing information from full information maximum likelihood,
\emph{Frontiers in Psychology}, \bold{12}, 667802.
}

\author{Joe Schafer \email{Joseph.L.Schafer@census.gov} }

\note{
If the model was fit with the option \code{saturated=TRUE}, then some
rows of the matrix of predicted values may be filled with
\code{NA}. These rows correspond to covariate patterns in the original
dataset that were empty in the sense that all reported frequencies
within the pattern were zero, or (if \code{newdata} is provided) to
covariate patterns that were not found in the original dataset.
}

\seealso{
\code{\link{rrLogit}},
\code{\link[=fitted.rrLogit]{fitted}},
}

\examples{
modelFit <- rrLogit( Food ~ Lake + Sex + Size, freq=Freq,
    data=alligatorNarrow)
fittedVal <- predict(modelFit)
}
