\name{alligatorNarrow}

\docType{data}

\alias{alligatorNarrow}

\title{Alligator Food Choice Data, Narrow Format}

\description{
This dataset, originally reported by Delany and Moore (1987) and
published by Agresti (2013), comes from a study of the primary food
choices of alligators in four Florida lakes. The same data arranged in
other formats are found in \code{\link{alligatorWide}} and
\code{\link{alligatorMicro}}. In this 
narrow-format version, the response variable (\code{Food}) is a factor with
five levels, and the frequencies are given in a single column.
}

\usage{alligatorNarrow}

\format{a data frame with 80 rows and five variables:

\describe{    

    \item{\code{Lake}}{factor indicating the lake where the alligator
    was captured, with levels \code{"Hancock"}, \code{"Oklawaha"},
    \code{"Trafford"}, and \code{"George"}}

    \item{\code{Sex}}{sex of the alligator (\code{"M"} or \code{"F"})}

    \item{\code{Size}}{size of the alligator, classified as
    \code{"small"} or \code{"large"}}

    \item{\code{Food}}{factor indicating the primary food choice, with
    levels \code{"Fish"}, \code{"Inv"}, \code{"Rept"}, \code{"Bird"},
    and \code{"Other"}}

    \item{\code{Freq}}{frequency}
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
\code{\link{alligatorMicro}}
}
