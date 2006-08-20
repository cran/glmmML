\encoding{UTF-8}
\name{ghq}
\alias{ghq}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Gauss-Hermite}
\description{
Calculates the zeros and weights needed for Gauss-Hermite quadrature. 
}
\usage{
ghq(n.points = 1, modified = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{n.points}{Number of points.}
  \item{modified}{Multiply by exp(zeros**2)? Default is TRUE.}
}
\details{

}
\value{
 A list vith components
  \item{zeros}{The zeros (abscissas).}
  \item{weights}{The weights}
}
\references{Gauss-hermite}
\author{Jianming Jin, Univ. of Illinois, Urbana-Campaign}
\note{
The code is modified to suit the purpose of glmmML, with the permission
of professor Jin.}
\seealso{\code{\link{glmmML}}}
\examples{
ghq(15, FALSE)
}
\keyword{math}% at least one, from doc/KEYWORDS

