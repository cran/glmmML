\name{glmmML}
\alias{glmmML}

\title{Generalized Linear Models with random intercept}
\description{
Fits GLMs with random intercept by Maximum Likelihood and numerical
integration via Gauss-Hermite quadrature.
}
\usage{
glmmML(formula, data = list(), cluster = NULL, family = binomial,
start.coef = NULL, start.sigma = NULL, offset = NULL, method = "vmmin",
control = glm.control(epsilon = 1e-08, maxit = 100, trace = FALSE),
n.points = 16) 
}
%- maybe also `usage' for other objects documented here.
\arguments{
  \item{formula}{ a symbolic description of the model to be fit. The details of
          model specification are given below.}
  \item{data}{an optional data frame containing the variables in the model.
           By default the variables are taken from
          `environment(formula)', typically the environment from which
          `glmmML' is called.
}
  \item{cluster}{Factor indicating which items are correlated}
  \item{family}{Currently, the only valid values are \code{binomial} and
  \code{poisson}. The binomial family allows for the \code{logit} and
  \code{cloglog} links, but can only be represented as binary data.}
  \item{start.coef}{starting values for the parameters in the linear predictor.
 Defaults to zero.}
  \item{start.sigma}{starting value for the mixing standard
    deviation. Defaults to 0.5.}
  \item{offset}{this can be used to specify an a priori known component to be
          included in the linear predictor during fitting.}
  \item{method}{the method to be used in fitting the model. The default (and
          presently only) method `vmmin' uses the BFGS method in the
	  'optim' function.}
  \item{control}{Controls the convergence criteria. See
  \code{\link{glm.control}} for details.} 
  \item{n.points}{Number of points in the Gauss-hermite quadrature. If
    n.points == 1, an ordinary glm is fitted.}
}
\details{
After the 'vmmin' function has converged, an ordinary Newton-Raphson
procedure finishes the maximization. As a by-product, the
variance-covariance is estimated. 
}
\value{
  The return value is a list, an object of class 'glmmML'.
}
\references{Broström (2003). Generalized linear models with random
  intercepts. \url{http://www.stat.umu.se/forskning/glmmML.pdf}
 }
\author{Göran Broström}
\note{This is a very preliminary version. The optimization may fail with
the default value of \code{start.sigma}. In that case, try different
start values for sigma.}

\seealso{\code{\link{optim}},  \code{\link[repeated]{glmm}} in Lindsey's
\code{repeated} package, and \code{\link[MASS]{glmmPQL}} in \code{MASS}.}

\examples{
x <- cbind(rep(1, 14), rnorm(14))
y <- rbinom(14, prob = 0.5, size = 1)
id <- rep(1:7, 2)
dat <- data.frame(y = y, x1 = x[, 2])
glmmML(y ~ x1, data = dat, cluster = id)
}
\keyword{regression}% at least one, from doc/KEYWORDS
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
