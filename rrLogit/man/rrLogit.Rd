\name{rrLogit}
\alias{rrLogit}
\alias{rrLogit.formula}
\alias{rrLogit.rrLogit}
\alias{rrLogit.default}

\title{
Multinomial Logistic Regression with Randomized Response
}

\description{
Fits a baseline-category logit model to a categorical response
variable that has been infused with random noise. Typical applications
include surveys on sensitive topics that use randomized response, and
person characteristics that have been deliberately distorted by the
post-randomization method for statistical disclosure avoidance.  The
noise-infusion mechanism is assumed to be known and is supplied by the
user in the form of a perturbation matrix.  If no such matrix is
given, the response is assumed to have no noise added, and the
procedure reverts to ordinary baseline-category logistic regression.
Data may be supplied in wide format, in narrow format with
frequencies, or as microdata. Methods include mode-finding procedures
for parameter estimation (EM, Newton-Raphson, Fisher scoring) and
Markov chain Monte Carlo (MCMC) for Bayesian posterior
simulation. With MCMC, a user can generate multiple imputations of the
unseen true response for use in subsequent analyses.
}

\usage{

% the generic function
rrLogit(obj, ...)

\method{rrLogit}{formula}(obj, data, baseline = 1L, freq, weight,
    pertMat = NULL, privLoss = NULL, saturated = FALSE,
    prior = c("none", "DAP"), priorFreqTot = NULL, priorAlloc = NULL,
    method = c("EM", "NR", "FS", "MCMC", "approxBayes"),
    startVal = NULL, control = list(), ...)

\method{rrLogit}{rrLogit}(obj, method = obj$method, 
    saturated = obj$saturated, control = NULL, startVal = NULL, ...)
   
}

\arguments{
  \item{obj}{an object used to select a method: either a formula that
  specifies a model to be fit, or the result from a previous call to
  rrLogit. If a formula, the syntax of the formula depends on whether
  the data are supplied in wide format or narrow format; see DETAILS.}

  \item{data}{an optional data frame, list or environment (or object
  coercible to a data frame by \code{as.data.frame}) containing the variables
  in the model.  If not found in \code{data}, the variables are taken from
  \code{environment(obj)}, typically the environment from which
  \code{rrLogit} is called.}

  \item{baseline}{response category to serve as the baseline (the
  denominator) in the baseline-category logits, specified either as an
  integer or character string; see DETAILS.}

  \item{freq}{an optional numeric variable for holding 
  frequencies when the data are supplied in narrow format and the
  observations are grouped.  If data are in narrow format and
  \code{freq} is omitted, the observations are assumed to represent
  microdata, and all frequencies are set to one.}

  \item{weight}{an optional numeric variable containing survey weights, which
  are used when computing pseudo-maximum likelihood (PML) estimates
  from survey data.  If \code{weight} is given, then the supplied data
  are interpreted as microdata, with each row having a frequency of
  one.}

  \item{pertMat}{perturbation or distortion matrix describing the
  noise-infusion mechanism.  Must be a full-rank square matrix of
  probabilities of order \eqn{C}, where \eqn{C} is the number of
  response categories, and each column must sum to one.  Element
  \eqn{(r,s)} of this matrix is the conditional probability that
  the noise-infused response is equal to \eqn{r}, given that the true
  response is \eqn{s}.  If \code{NULL}, then \code{pertMat} is taken to
  be the identity matrix.}

  \item{privLoss}{privacy-loss parameter which provides an alternative
  way to specify \code{pertMat} under the differentially-private
  randomized-response
  mechanism described by Wang, Wu and Hu (2016).  Produces a symmetric
  perturbation matrix with diagonal elements \eqn{p} and off-diagonal
  elements \eqn{(1-p)/(C-1)}, where
  \deqn{p=( 1 + (C-1) e^\epsilon)^{-1}}{p = 1 / ( 1 + (C-1) * exp(
  -\epsilon) ) } and \eqn{\epsilon} is the privacy-loss parameter.}

  \item{saturated}{if \code{TRUE}, fits a saturated model with
  response probabilities estimated independently within each covariate
  pattern; see DETAILS.}

  \item{prior}{if \code{"none"}, then no prior distribution is
  used.  If \code{"DAP"}, then a data-augmentation prior as described by
  Clogg et al. (1991) is applied with a total prior frequency of
  \code{priorFreqTot}, divided equally among the unique covariate
  patterns and allocated to the response categories in proportions
  specified by \code{priorAlloc}; see DETAILS.}

  \item{priorFreqTot}{total prior frequency when \code{prior="DAP"},
  representing the total number of imaginary prior observations added
  to the dataset when fitting the model.  Ignored when \code{prior="none"}.}

  \item{priorAlloc}{proportions for allocating imaginary prior
  observations to response categories when \code{prior="DAP"}.  This
  represents a prior guess for the marginal proportions of the true
  (noise-free) response variable in the population.  If \code{NULL},
  then \code{priorAlloc} is obtained by estimating the marginal
  proportions for the observed (noise-infused) response variable in
  the sample, and then pre-multiplying those proportions by the
  inverse of \code{pertMat}.  Ignored when \code{prior="none"}.}

  \item{method}{procedure for fitting the model: \code{"EM"} for
  an expectation-maximization algorithm; \code{"NR"} for
  Newton-Raphson; \code{"FS"} for Fisher scoring; \code{"MCMC"} for
  Markov chain Monte Carlo; and \code{"approxBayes"} for approximate
  Bayes. See DETAILS.}

  \item{startVal}{optional starting values for the model fitting
  procedure.  If \code{saturated=FALSE}, \code{startVal} should be
  a matrix of coefficients with rows corresponding to model predictors
  and columns corresponding to response categories, with zeros in the
  column corresponding to the baseline.  If \code{saturated=TRUE},
  \code{startVal} should be a matrix of probabilities with rows
  corresponding to distinct covariate patterns and columns
  corresponding to response categories. This argument is rarely or
  never used; in the vast majority of situations, it should be kept at
  \code{NULL}. If \code{startVal=NULL} and \code{obj} is a model
  formula, the function finds reasonable starting values by its
  own default procedure. If \code{startVal=NULL} and \code{obj} is the
  result from a previous call to \code{rrLogit}, then the parameter
  values stored in that object are taken as starting values. See DETAILS.}

  \item{control}{a named list containing control parameters which are
  passed to \code{\link{rrLogitControl}}.  Control parameters determine the
  maximum number of iterations, criteria for judging convergence and
  proximity to a boundary, and so on.  Control parameters that are
  not provided through this argument are set to default values, which
  can be viewed by calling \code{rrLogitControl()}.}

  \item{\dots}{additional arguments passed to or from other methods}

}

\details{
Traditional baseline-category logistic regression (Agresti, 2013) relates
a categorical response variable \eqn{Y} with categories
\eqn{y=1,\ldots,C} to a set of predictors.  The model is
\deqn{log \left( \frac{P(Y=y)}{P(Y=b)} \right) = x' \beta_y}{
log( P(Y=y) / P(Y=b) ) = x' \beta_y} for \eqn{y=1,\ldots,C},
where \eqn{x} is a vector of unmodeled predictors,
\eqn{\beta_1,\ldots,\beta_C} are vectors of unknown coefficients to be
estimated, and \eqn{b} is a response category selected to serve as the
baseline.  The coefficients associated with the baseline category
\eqn{\beta_b} are constrained to be zero.  As with other functions for
regression modeling such as \code{\link{lm}} and \code{\link{glm}},
the vector \eqn{x} is not directly supplied by the user, but is
computed automatically from the predictor variables appearing in the
model formula. Maximum-likelihood (ML) estimates for the unknown
coefficients are typically computed by an iterative Newton-Raphson (NR)
procedure. 

The function \code{rrLogit} is designed to fit baseline-category
logit models in situations where \eqn{Y} is withheld from the data
analyst and replaced by a perturbed or noise-infused version
\eqn{Y^*}. This happens when survey data on sensitive topics are
collected using randomized response (RR), or when key items that could be
used to identify individuals are deliberately distorted by the
post-randomization method (PRAM) (van den Hout and van der Heijden,
2007). The stochastic relationship between \eqn{Y} and \eqn{Y^*} is described
by a perturbation or distortion matrix which is provided to the data
analyst. With noise infusion, the interpretation of the model
has not changed; the parameters still describe relationships between
predictors and the ``true'' response \eqn{Y}, but the model-fitting is
carried out using the noisy response \eqn{Y^*}.  Parameter estimates
computed by \code{rrLogit} based on \eqn{Y^*} are still consistent, but
in general, standard errors will be larger than they would have been
without the addition of noise.

The loglikelihood function based on \eqn{Y^*} is more complicated and
often less well behaved than the loglikelihood based on \eqn{Y}; for
example, it may not be concave everywhere, causing NR to fail unless
the starting values are already near the solution. Fisher
scoring (FS), which is closely related to NR, does not fare much better.
A more reliable procedure is an expectation-maximization (EM) algorithm that
treats the noise-free response \eqn{Y} as ``missing data,'' repeatedly
maximizing a simpler and better behaved loglikelihood function
where the unknown frequencies for \eqn{Y} have been replaced by their expected
values.  EM is the default procedure used by \code{rrLogit}, but other
methods are available; see \emph{Methods} below.

\emph{Model formulas for data in wide format.}  Data may be provided
to \code{rrLogit} in wide format, narrow format, or as microdata.  In
wide format, each row of the dataset represents a group of
observational units whose predictor variables are identical, and
frequencies for the response categories are stored as numeric
variables, with one variable per response category.  With wide-format
data, the left-hand side of the model formula should be the names of
the frequency variables, separated by commas and enclosed by
\code{cbind()}; the right-hand side should be the predictors,
separated by the standard regression-formula operators \code{+},
\code{*} and \code{:}.  For example, suppose that there are \eqn{C=3}
response categories with frequency variables named \code{f1},
\code{f2}, and \code{f3}.  The formula \code{cbind(f1,f2,f3) ~ x1 +
x2*x3} specifies a model with main effects for \code{x1}, \code{x2},
\code{x3} and the \code{x2:x3} interaction.

\emph{Model formulas for data in narrow format and microdata.}  In
narrow format, each row of the dataset represents a group of
observations that are identical with respect to the predictor
variables and the response.  The response is a factor
variable with \eqn{C} levels.  With narrow-format data, the response
variable should appear on the left-hand side of the model formula, as
in \code{Y ~ x1 + x2*x3}.  If frequencies are present, the name of the
numeric variable containing the frequencies should be provided as the
argument \code{freq}.  If \code{freq} is not provided, the dataset is
assumed to contain microdata, with each row representing a single
observational unit (i.e., all frequencies are assumed to be one).

\emph{Baseline category.}  The baseline category is determined by the
argument \code{baseline} and may be supplied either as an integer or
as a character string.  For example, suppose the data are in wide
format and the model formula is \code{cbind(f1,f2,f3) ~ x1 + x2*x3};
in this case, the last category may be selected as the baseline by
setting \code{baseline=3} or \code{baseline="f3"}.  Or suppose the data
are in narrow or microdata format, the model formula is \code{Y ~ x1 +
x2*x3}, and \code{Y} is a factor whose \code{levels} attribute is
\code{c("red","green","blue")}; in this case, the last category may
be selected as the baseline by setting \code{baseline=3} or
\code{baseline="blue"}.

\emph{Survey weights.}  With microdata, the user has the option of
providing survey weights.  Weights play the same role as frequencies
for computing parameter estimates, but special procedures must be
applied afterward to compute standard errors in a manner consistent
with the survey design.  Survey weights, if present, should be
supplied through the argument \code{weight}, and all frequencies are
then assumed to be one.  Before model fitting, the weights are
automatically scaled to sum to the overall sample size, which has no
effect on the estimates.

\emph{Saturated model.}  To assess the quality of a model, it is often
helpful to compare its fit to that of a saturated model.  The
saturated model estimates the probabilities \eqn{P(Y=1), P(Y=2),
\ldots, P(Y=C)} independently within each of the covariate patterns
(unique combinations of values of the predictor variables) occurring
in the dataset.  To fit a saturated model, supply the argument
\code{saturated=TRUE}.  The parameters of a saturated model are
expressed as a matrix of probabilities with one row for each covariate
pattern and one column for each response category.  Each row of this
matrix sums to one.  If a covariate pattern happens to be empty in the
sense that it has a total frequency or total survey weight of zero,
the corresponding row of the matrix will be filled with \code{NA}s.
To extract the estimated parameters from a saturated model, use
\code{\link[=fitted.rrLogit]{fitted}}.

\emph{Prior distributions.}  By default, \code{rrLogit} computes
maximum-likelihood (ML) estimates for model parameters.  In some
cases, ML estimation can be unstable.  This may happen when data are
sparse, when predictors are highly intercorrelated, or under
conditions analogous to complete or quasi-complete separation in
ordinary logistic models (Albert and Anderson, 1984).  To stabilize
the model fit in these situations, it may be helpful to apply a
data-augmentation prior (DAP) using the argument \code{prior="DAP"}.
The DAP is a penalty function that smooths the fitted values for all
covariate patterns toward a common set of probabilities.  It is
functionally equivalent to augmenting the observed data by inserting
fictitious fractional observations into each covariate pattern.  The
overall strength of the DAP is determined by \code{priorFreqTot}, a
positive number indicating the total number of prior observations to
be added to the dataset.  By default (i.e., if \code{prior="DAP"} and
\code{priorFreqTot=NULL}), the total number of prior observations is
set to the number of unknown free parameters in the model, as
suggested by Clogg et al. (1991).  This total frequency is distributed
equally among the non-empty covariate patterns appearing in the
dataset.  Within each covariate pattern, the frequency is then
allocated to categories of the ``true'' (noise-free) response in
proportions given by \code{priorAlloc}.  If \code{priorAlloc=NULL},
these proportions are estimated from the data by computing the
marginal proportions for the noisy response from the observed data and
pre-multiplying that vector by the inverse of \code{pertMat}. If any
element of the resulting vector is negative, an error message is
given, in which case the user should supply a \code{priorAlloc} vector
whose elements are positive and sum to one. A typical choice for
\code{priorAlloc} is \code{rep(1/C,C)}, where \code{C} is the number
of response categories.

\emph{Methods.}  The default estimation procedure \code{method="EM"} is
an expectation-maximization (EM) algorithm that regards the true
(noise-free) response frequencies for each covariate pattern as
missing data.  At each cycle, EM computes the expected true frequencies
given the observed (noisy) response frequencies under the current
estimate of the parameters, then re-estimates the parameters by
maximizing the complete-data loglikelihood function based on these
expected true frequencies.  \code{method="NR"} (Newton-Raphson)
operates on the loglikelihood function based on the observed
frequencies, iteratively updating the parameter estimates using the
loglikelihood's first and second derivatives.  \code{method="FS"}
(Fisher scoring) is similar to \code{"NR"} but substitutes expected
second derivatives for the actual ones.  Regardless of which method is
used, convergence is judged not by examining changes in the estimated
coefficients, but by comparing the estimated response probabilities
(fitted values) within covariate patterns from one iteration to the
next, stopping when all fitted values change by less than a given
threshold.  That convergence criterion, and other control parameters,
are set by the argument \code{control} and the function
\code{\link{rrLogitControl}}.  When \code{saturated=TRUE}, the
saturated model is fit by applying an EM algorithm independently within
each non-empty covariate pattern.

\emph{Starting values.}  If \code{rrLogit} is called with a model
formula as its first argument and no starting values are supplied
(\code{startVal = NULL}), then default starting values are used.
Default starting values for the coefficients of a non-saturated model
(\code{saturated=FALSE}) are obtained by setting each row of the
fitted response-probability matrix to the estimated marginal
proportions, applying the logit transformation and regressing on the
model matrix.  For a saturated model (\code{saturated=TRUE}), default
starting values for the response probabilities within each covariate
pattern are set to uniform values \eqn{1/C}, where \eqn{C} is the
number of response categories.

\emph{Applying the \code{rrLogit} function to an \code{rrlogit}
object.} If \code{rrLogit} is called with an \code{rrLogit} object as
its first argument, then the data, model and prior distribution will
be taken from that object, and (unless \code{startVal} was supplied),
starting values will be set to the parameter values stored in
that object.

\emph{Markov chain Monte Carlo (MCMC).} With \code{method="MCMC"}, the
\code{rrLogit} function simulates random draws from the posterior
distribution of the model parameters under a uniform prior (if
\code{prior="none"}) or a data-augmentation prior (if
\code{prior="DAP"}). In typical applications, \code{method="MCMC"}
would be applied to an \code{rrLogit} object that holds the results
from \code{method="EM"}, \code{method="NR"}, or
\code{method="FS"}. The MCMC method may be used on
its own, or in tandem with \code{\link[=impute.rrLogit]{impute}},
to create multiple imputations (MIs) of the unseen true 
response for subsequent analyses.

When \code{saturated=FALSE}, two versions of MCMC are available. The
default procedure is a random-walk Metropolis (RWM) algorithm (Gelman
\emph{et al}., 2013) whose proposal is a multivariate
t-distribution centered at the current parameters, and whose scale matrix
is equal to the Hessian-based covariance matrix for the parameters
multiplied by a scale factor. An alternative to RWM is a
data-augmentation (DA) procedure that resembles EM, alternately
drawing values of the noise-free response given the parameters and then the
parameters given the noise-free response. When \code{saturated=TRUE},
only DA is available. Constants that select the MCMC method and
determine its behavior are set through the \code{control} argument,
which passes them to \code{\link{rrLogitControl}}.

\emph{Convergence diagnostics for MCMC.} When using
\code{method="MCMC"}, summaries of the
output stream are stored in the resulting \code{rrLogit} object. These
may be extracted using the functions documented
\link[=coefSeries]{here} and then examined using functions in the 
\code{coda} package.

\emph{Multiple imputation.} Multiple imputations of the unseen true
response may be generated in two different ways. Perhaps the easiest
way is to do a long run of MCMC, storing the imputations along the way
using \code{control=list(imputeEvery=...)}, and then extracting the
imputations with \code{\link{impList}}. The other way is to run
mutiple sequential or parallel chains and then apply
\code{\link[=impute.rrLogit]{impute}} to the result of each one.


\emph{Approximate Bayesian inference.} \code{method="approxBayes"} is
experimental; it does not perform well and should not be used.

}

\value{an object of class \code{c("rrLogit","list")} containing the
results of a model fit.  Components of this list may be examined
directly, but it is recommended that users extract information from
the object using the methods specifically written for that purpose,
including \code{\link[=summary.rrLogit]{summary}} and others listed
below under SEE ALSO.
}

\references{
Agresti, A. (2013) \emph{Categorical Data Analysis, Third
Edition}. Hoboken, NJ: John Wiley & Sons.

Albert, A, and Anderson, J.L. (1984) On the existence of
maximum-likelihood estimates in logistic regression models,
\emph{Biometrika}, \bold{71}, 1--10.

Clogg, C.C., Rubin, D. B., Schenker, N., Schultz, B. and Weidman,
L. (1991) Multiple imputation of industry and occupation codes in
census public-use samples using Bayesian logistic regression, 
\emph{Journal of the American Statistical Association}, \bold{79},
762--771.

Gelman, A., Stern, H.S., Carlin, J.B., Dunson, D.B., Vehtari, A. and
Rubin, D.B. (2013) \emph{Bayesian Data Analysis, Third Edition}. Boca
Raton, FL: Chapman & Hall/CRC Press.

van den Hout, A., van der Heijden, P.G. (2007) Randomized response,
statistical disclosure control and misclassification: a
review, \emph{International Statistical Review}, \bold{70}, 269-–288.

Wang, Y., Wu, X. and Hu, D. (2016) Using randomized response for
differential privacy preserving data collection. \emph{Proceedings of
EDBT/ICDT Workshops}, \bold{1558},
\url{http://ceur-ws.org/Vol-1558/paper35.pdf}.
}

\author{Joe Schafer \email{Joseph.L.Schafer@census.gov} }


\seealso{
\code{\link[=anova.rrLogit]{anova}},
\code{\link[=coef.rrLogit]{coef}},
\code{\link{coefSeries}},
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

\examples{
#######################
# model fitting to a response without random noise, using the alligator data

# fit main-effects model to alligator data in wide format, show coefficients
fitA <- rrLogit( cbind(Fish,Inv,Rept,Bird,Other) ~ Lake + Sex + Size,
   data=alligatorWide)
coef(fitA)

# fit same model to data in narrow format, verifying that
# estimated coefficients are the same
fitB <- rrLogit( Food ~ Lake + Sex + Size, freq=Freq, data=alligatorNarrow )
all.equal( coef(fitA), coef(fitB) )

# repeat for data in micro format
fitC <- rrLogit( Food ~ Lake + Sex + Size, data=alligatorMicro )
all.equal( coef(fitA), coef(fitC) )

#######################
# model fitting to a response containing random noise, using alligator data
# (purely for demonstration, as noise infusion makes no actual sense here)

# infuse alligator microdata with noise using randomized-response mechanism,
# setting privacy-loss parameter to 2.0
myPertMat <- rrPertMat( 2.0, 5 )
set.seed(789)
alligatorMicro$noisyFood <- rrPerturbResponse( alligatorMicro$Food,
    myPertMat )

# fit main-effects model to noise-infused response variable
fitD <- rrLogit( noisyFood ~ Lake + Sex + Size, data=alligatorMicro,
    pertMat=myPertMat)

# stabilize the fit using a data-augmentation prior
fitE <- rrLogit( noisyFood ~ Lake + Sex + Size, data=alligatorMicro,
    pertMat=myPertMat, prior="DAP")

#######################
# simple example involving two survey items with randomized response
# and no predictors: rule violations in unemployment benefits,
# published by van den Hout and van der Heijden (2007) and analyzed
# in the rrLogit package vignette 

# marginal analysis of first item
df.1 <- data.frame( Yes=120, No=292 )
Tmat <- matrix( c(.8, .2, .2, .8), 2, 2 )
fit.1 <- rrLogit( cbind(Yes, No) ~ 1, data=df.1, pertMat=Tmat )   
predict(fit.1)
predict(fit.1, se.fit=TRUE)$se.fit

# marginal analysis of responses to second item
df.2 <- data.frame( Yes=171, No=241 )
fit.2 <- rrLogit( cbind(Yes, No) ~ 1, data=df.2, pertMat=Tmat )   
predict(fit.2)
predict(fit.2, se.fit=TRUE)$se.fit

# joint analysis of both items
df.12 <- data.frame( Yes.Yes=68, Yes.No=52, No.Yes=103, No.No=189 )
fit.12 <- rrLogit( cbind(Yes.Yes, Yes.No, No.Yes, No.No) ~ 1,
   data=df.12, pertMat=kronecker(Tmat, Tmat) )   
predict(fit.12)
predict(fit.12, se.fit=TRUE)$se.fit

}