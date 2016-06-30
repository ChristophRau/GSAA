\name{matchkeep}
\alias{matchkeep}
\title{ matchkeep }
\description{
  This function looks at a row/column, determines if a given element has a match in a master list, and keeps the row/column if it does.}

\arguments{
  \item{x}{row/column vector}
  \item{names_to_keep}{master list to match to.}
  \item{element}{which element to attempt to match with}
  }

\value{
  x if match, NULL if no match
}

\usage{
  matchkeep(x,names_to_keep, element)
}
\author{ Christoph Rau }
\note{
}
\seealso{ \code{\link{GSAA}},\code{\link{calculateIncidenceTable}}}