\name{impute.rrLogit}
\alias{impute.rrLogit}
\alias{impute}

\title{Impute True Responses from an \code{rrLogit} Model Fit}

\description{
\code{impute.rrLogit}, which can be invoked using the generic
name \code{impute}, imputes values of the ``true'' (noise-free)
response variable from their posterior predictive distribution given
predictors, observed (noisy) responses, and parameter values held in 
\code{\link{rrLogit}} object. The values of predictors and the
observed response may come from the data frame used to fit the model
or another data frame.
}

\usage{
% the generic function
impute(object, ...)

\method{impute}{rrLogit}(object, newdata=NULL, freq=NULL,
   type=c("random", "condMean"), impVarName=NULL,
   meanSeries=FALSE, na.action=na.pass, ...)
}

\arguments{
  \item{object}{an object produced by \code{rrLogit}}

  \item{newdata}{an optional data frame where \code{impute} looks for
  predictor variables and observed responses to generate the imputed
  values. If \code{NULL}, imputations will be based on the data used
  to fit the model contained in \code{object}.}

  \item{freq}{an optional numeric variable in \code{newdata}
  containing frequencies. If \code{newdata=NULL}, \code{freq} is
  ignored and frequencies are taken from the data frame used to fit
  the model.}

  \item{type}{type of imputed values: \code{"random"} for a random 
  draw from the multinomial predictive distribution of the unseen true
  (noise-free) responses given the observed responses, predictors, and
  parameters held in \code{object}; \code{"condMean"} for the expected
  values from this distribution.}

  \item{impVarName}{name to give the imputed response variable. If
  \code{NULL}, then the imputed response variable is set to the name
  of the observed noisy response variable, with the added suffix
  \code{".imp"}.}

  \item{meanSeries}{applies when \code{obj} contains results from a
     simulation run. If \code{TRUE}, then the imputations are based on
     a running mean of cell probabilities over all iterations 
     after the burn-in period. If \code{FALSE}, then the imputations
     are based only on the cell probabilities from the final
     iteration, and (assuming the run was sufficiently long, if MCMC)
     can be regarded as a single draw from their posterior distribution;
     see NOTE.}

  \item{na.action}{action to be taken if predictors or responses in
  \code{newdata} contain missing values. For now, \code{NA} values are
  not allowed, but this may change in future versions.}

  \item{\dots}{additional arguments passed to or from other methods}
}

\details{
If the original dataset used to fit the model held in \code{object}
was provided in wide format, then the response frequencies in
\code{newdata} should also be arranged in wide format. Otherwise,
\code{newdata} should be arranged in narrow format with
frequencies specified through \code{freq}, or as microdata with
\code{freq} omitted. For descriptions of these data formats, see
\code{\link{rrLogit}}.
}

\note{
When this function is used within a process for multiple imputation,
\code{meanSeries} should 
be set to \code{FALSE}, otherwise the imputations will not correctly reflect
uncertainty about model parameters. 
}

\value{
if \code{type="random"} and the dataset (original data frame or
\code{newdata}) contains microdata, the imputed values are expresses as
a single factor variable; otherwise, the imputed values are expressed as a
matrix having one column for each response category.

If \code{newdata=NULL}, the imputed values are returned as a single
factor variable or a matrix. If a new data frame was provided through
\code{newdata}, that entire data frame is returned with the imputed
values appended as additional variables using \code{rbind()}.
}

\author{Joe Schafer \email{Joseph.L.Schafer@census.gov} }

\note{
If the model was fit with the option \code{saturated=TRUE}, then some
imputed values may be \code{NA}. These missing values correspond to
covariate patterns in the original 
dataset that were empty in the sense that all reported frequencies
within the pattern were zero, or (if \code{newdata} is provided) to
covariate patterns that were not present in the original dataset.
}

\seealso{
\code{\link{rrLogit}},
\code{\link[=predict.rrLogit]{predict}},
}

\examples{
# infuse microdata with noise using randomized-response mechanism,
# setting privacy-loss parameter to 2.0
myPertMat <- rrPertMat( 2.0, 5 )
set.seed(789)
alligatorMicro$noisyFood <- rrPerturbResponse( alligatorMicro$Food,
    myPertMat )

# fit main-effects model to noise-infused response variable,
# stabilizing the fit with a data-augmentation prior
fit <- rrLogit( noisyFood ~ Lake + Sex + Size, data=alligatorMicro,
    pertMat=myPertMat, prior="DAP")
# impute the true (noise-free) responses from their predictive
# distribution
imp <- impute(fit)
}
