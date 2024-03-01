\name{rrLogitControl}
\alias{rrLogitControl}

\title{
Control Parameters for \code{rrLogit}
}

\description{
The \code{\link{rrLogit}} function fits baseline-category logistic
regression models to categorical response variables that have been
distorted by random noise. Its model-fitting
procedures are governed by parameters in a \code{rrLogitControl} object
created by the auxiliary function documented here. This function was designed
for internal use by \code{rrLogit}; the main reason why a user might
invoke \code{rrLogitControl} directly would be to display the control
parameters and their default values.
}

\usage{
rrLogitControl( iterMaxNR = 50L, iterMaxFS = 200L, iterMaxEM = 1000L,
   iterMaxMstep = 25L, critConverge = 1e-08, critBoundary = 1e-08,
   iterApproxBayes = 1L, imputeApproxBayes = FALSE,
   iterMCMC = 5000L, burnMCMC = 0L, thinMCMC = 1L, imputeEvery = 0L,
   typeMCMC = c("RWM","DA"), tuneDA = c(100,.2,.5),
   tuneRWM = c(1000,.2), stuckLimit = 25L, startValJitter = 0 )
}

\arguments{
  \item{iterMaxNR}{maximum number of Newton-Raphson iterations performed when
  \code{method = "NR"}}

  \item{iterMaxFS}{maximum number of Fisher scoring iterations performed when
  \code{method = "FS"}}

  \item{iterMaxEM}{maximum number of EM iterations performed when
  \code{method = "EM"}}

  \item{iterMaxMstep}{maximum number of iterations of Newton-Raphson
  performed during an M-step of EM}

  \item{critConverge}{convergence criterion for stopping rule}

  \item{critBoundary}{criterion for testing whether any estimated
  response probabilities are close to a boundary, in which case a
  warning is given}

  \item{iterApproxBayes}{number of simulated coefficient
  vectors to be drawn from their approximate posterior distribution
  when \code{method="approxBayes"}.}

  \item{imputeApproxBayes}{if \code{TRUE} then, for each draw of the
  logistic coefficients from their approximate posterior distribution,
  the true (noise-free) responses will be imputed as a random draw
  from their conditional posterior distribution.}

  \item{iterMCMC}{number of iterations of Markov chain Monte Carlo
  after the burn-in period when \code{method="MCMC"}.}

  \item{burnMCMC}{number of iterations of Markov chain Monte Carlo
  performed as a burn-in period, for which the results are
  discarded. The total number of iterations performed is
  \code{burnMCMC+iterMCMC}.} 

  \item{thinMCMC}{thinning interval for saving the results from MCMC
  as a series.} 

  \item{imputeEvery}{imputation interval for saving imputed
  frequencies for the complete-data table. If \code{0}, then no
  imputations are saved.} 

  \item{typeMCMC}{either \code{"DA"} (data augmentation) or
  \code{"RWM"} (random-walk Metropolis); see DETAILS.}

  \item{tuneDA}{tuning parameters for data augmentation MCMC; see DETAILS.}

  \item{tuneRWM}{tuning parameter for random-walk Metropolis MCMC; see
  DETAILS.}

  \item{stuckLimit}{criterion for deciding if the MCMC algorithm has
  gotten stuck.}

  \item{startValJitter}{standard deviation for Gaussian random noise added to
  starting values. If \code{rrLogit} is called with
  \code{saturated=FALSE}, the logistic coefficients are perturbed by
  this amount; if \code{saturated=TRUE}, the log-probabilities for
  each covariate pattern are perturbed by this amount and renormalized to sum
  to one.}
}

\details{
When  \code{method="NR"} or \code{method="FS"}, the \code{rrLogit}
function performs Newton-Raphson or Fisher scoring, stopping when the
absolute changes in fitted probabilities for all covariate patterns from one
iteration to the next fall below the threshold for convergence
(\code{critConverge}) or when the maximum number of iterations
(\code{iterMaxNR} or \code{iterMaxFS}) has been reached, whichever
comes first.

The EM algorithm is doubly iterative. The main iteration of EM
consists of an expectation or E-step followed by a maximization or
M-step. Each M-step is also iterative, carried out by a Newton-Raphson
procedure with its own stopping rule. The number of iterations
performed within an M-step is governed by \code{critConverge} and
\code{iterMaxMstep}, and the number of main iterations is governed by
\code{critConverge} and \code{iterEM}. When \code{method="EM"}, the
number of iterations performed reported by
\code{\link[=summary.rrLogit]{print.summary}} refers to main
iterations.

When \code{rrLogit} is called with \code{method="MCMC"}, the algorithm
that is run depends on \code{typeMCMC} and on whether the model is fit
with \code{saturated=TRUE}.
\itemize{
\item If \code{saturated=FALSE} and
\code{typeMCMC="DA"}, then the algorithm is a data-augmentation
procedure that resembles EM. At each cycle, the true (noise-free)
responses are imputed by taking a random draw from their posterior
predictive distribution given the observed (noisy) responses,
predictor variables, and the current values of the parameters; the
parameters are then updated by taking one step of a
Metropolis-Hastings algorithm that mimics Newton-Raphson and
conditions on the imputed true responses.  The proposal distribution
is multivariate-t and can be adjusted by tuning constants in 
\code{tuneDA}, a numeric vector containing the degrees of
freedom, step size and scale factor.
\item If \code{saturated=FALSE} and \code{typeMCMC="RWM"},
the true responses are not imputed, and the logistic
coefficients are updated by a step of random-walk
Metropolis. The proposal is mutivariate-t and can be adjusted by
tuning constants in \code{tuneRWM}, a numeric vector containing the
degrees of freedom and scale factor.
\item If \code{saturated=TRUE}, then the algorithm is a
data-augmentation procedure that requires no tuning.}
}

\value{a list of control parameters for internal use by the
function \code{rrLogit}.
}

\author{Joe Schafer \email{Joseph.L.Schafer@census.gov} }

\seealso{
\code{\link{rrLogit}}
}

\examples{
# display all control parameters and their default values
rrLogitControl()
}
