\name{print.glmmML}
\alias{print.glmmML}
\title{Prints a 'glmmML' object.}
\description{
A glmmML object is the output of \code{glmmML}.
}
\usage{
print.glmmML(x, digits = max(3, getOption("digits") - 3), na.print = "", ...)
}
%- maybe also `usage' for other objects documented here.
\arguments{
  \item{x}{The glmmML object}
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
\references{Broström (2003)}
\author{Göran Broström}
\note{This is the only summary method available for the moment.}

\seealso{\code{\link{glmmML}}}

\examples{
## The function is currently defined as
function (x, digits = max(3, getOption("digits") - 3), na.print = "", 
    ...) 
{
  cat("\nCall: ", deparse(x$call), "\n\n")
  cat("Coefficients:\n")
  print.default(format(x$coefficients, digits = digits), print.gap = 2, 
                quote = FALSE)
  cat("Standard errors:\n")
  print.default(format(x$sd, digits = digits), print.gap = 2, 
                quote = FALSE)
  if(x$mixed){
    cat("Standard deviation in mixing distribution: ", x$sigma, "\n")
    cat("Std. Error:                                ", x$sigma.sd, "\n")
  }
  cat("\nDegrees of Freedom:", 
      x$df.residual, "Residual\n")
  cat("Residual Deviance:",
      format(signif(x$deviance, digits)),
      "\tAIC:",
      format(signif(x$aic, digits)), "\n")
  invisible(x)
  }
}
\keyword{print}% at least one, from doc/KEYWORDS
