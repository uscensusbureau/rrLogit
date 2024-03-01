\name{alligatorWide}

\docType{data}

\alias{alligatorWide}

\title{Alligator Food Choice Data, Wide Format}

\description{
This dataset, originally reported by Delany and Moore (1987) and
published by Agresti (2013), comes from a study of the primary food
choices of alligators in four Florida lakes. Researchers classified
the the stomach contents of 219 captured alligators into five
categories: Fish (the most common primary food choice), Invertebrate
(snails, insects, crayfish, etc.), Reptile (turtles, alligators),
Bird, and Other (amphibians, plants, household pets, stones, and other
debris). These data are arranged in wide format, with frequencies for
each of the five food choices appearing in a separate
column. The same dataset arranged in other formats are provided as
\code{\link{alligatorNarrow}} and \code{\link{alligatorMicro}}.
}

\usage{alligatorWide}

\format{a data frame with 16 rows and eight variables:

\describe{    

    \item{\code{Lake}}{factor indicating the lake where the alligator
    was captured, with levels \code{"Hancock"}, \code{"Oklawaha"},
    \code{"Trafford"}, and \code{"George"}}

    \item{\code{Sex}}{sex of the alligator (\code{"M"} or \code{"F"})}

    \item{\code{Size}}{size of the alligator, classified as
    \code{"small"} or \code{"large"}}

    \item{\code{Fish}}{number of alligators with primary food choice Fish}

    \item{\code{Inv}}{number of alligators with primary food choice
    Invertebrate}

    \item{\code{Rept}}{number of alligators with primary food choice Reptile}

    \item{\code{Bird}}{number of alligators with primary food choice Bird}

    \item{\code{Other}}{number of alligators with primary food choice Other}

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
\code{\link{alligatorNarrow}},
\code{\link{alligatorMicro}}
}
