\encoding{UTF-8}
\name{print.glmmboot}
\alias{print.glmmboot}
\title{Prints a 'glmmML' object.}
\description{
A glmmboot object is the output of \code{glmmboot}.
}
\usage{
print.glmmboot(x, digits = max(3, getOption("digits") - 3), na.print = "", ...)
}
%- maybe also `usage' for other objects documented here.
\arguments{
  \item{x}{The glmmboot object}
  \item{digits}{Number of printed digits.}
  \item{na.print}{How to print NAs}
  \item{\dots}{Additional parameters, which are ignored.}
}
\details{
Nothing in particular.
}
\value{
  A short summary of the object is printed.
}
%\references{Broström (2003)}
\author{Göran Broström}
\note{This is the only summary method available for the moment.}

\seealso{\code{\link{glmmboot}}}

%\examples{
%}
\keyword{print}% at least one, from doc/KEYWORDS
