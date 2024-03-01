\name{alligatorMicro}

\docType{data}

\alias{alligatorMicro}

\title{Alligator Food Choice Data, Microdata Format}

\description{
This dataset, originally reported by Delany and Moore (1987) and
published by Agresti (2013), comes from a study of the primary food
choices of alligators in four Florida lakes. The same data arranged in
other formats are found in \code{\link{alligatorWide}} and
\code{\link{alligatorNarrow}}. In this 
microdata version, each row corresponds to one alligator, and
the response variable (\code{Food}) is a factor with five
levels.
}

\usage{alligatorMicro}

\format{a data frame with 219 rows and four factor variables:

\describe{    

    \item{\code{Lake}}{lake where the alligator
    was captured, with levels \code{"Hancock"}, \code{"Oklawaha"},
    \code{"Trafford"}, and \code{"George"}}

    \item{\code{Sex}}{sex of the alligator (\code{"M"} or \code{"F"})}

    \item{\code{Size}}{size of the alligator, classified as
    \code{"small"} or \code{"large"}}

    \item{\code{Food}}{primary food choice, with
    levels \code{"Fish"}, \code{"Inv"}, \code{"Rept"}, \code{"Bird"},
    and \code{"Other"}.}
}
}

\source{
Delany, M.F. and Moore, C.T. (1987) American alligator food habits in
Florida. Unpublished manuscript.

Agresti, A. (2013) \emph{Categorical Data Analysis, Third
Edition}. Hoboken, NJ: John Wiley & Sons.
}

\keyword{ datasets }

\seealso{
\code{\link{alligatorWide}},
\code{\link{alligatorNarrow}}
}
