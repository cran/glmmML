\name{glmmML.fit}
\alias{glmmML.fit}

\title{Generalized Linear Model with random intercept}
\description{
This function is called by \code{glmmML}, but it can also be called
directly by the user.
}
\usage{
glmmML.fit(X, Y, weights = rep(1, NROW(Y)),
start.coef = NULL, start.sigma = NULL,
fix.sigma = FALSE, mixed = TRUE,
cluster = NULL, offset = rep(0, nobs), family = binomial(), n.points =
15, control = list(epsilon = 1.e-8, maxit = 200, trace = FALSE),
method = 0,
intercept = TRUE, boot = 0, prior = 0) 
}
%- maybe also `usage' for other objects documented here.
\arguments{
  \item{X}{Design matrix of covariates.}
  \item{Y}{Response vector. Or two-column matrix.}
  \item{weights}{Case weights. Defaults to one.}
  \item{start.coef}{Starting values for the coefficients. }
  \item{start.sigma}{Starting value for the mixing standard deviation.}
  \item{fix.sigma}{Should sigma be fixed at start.sigma?}
  \item{mixed}{Logical. If FALSE, an ordinary glm is fitted.}
  \item{cluster}{The clustring variable.}
  \item{offset}{The offset in the model.}
  \item{family}{Family of distributions. Defaults to binomial with logit
  link. Other possibilities are binomial with cloglog link and poisson
  with log link.}
  \item{n.points}{Number of points in the Gauss-hermite quadrature.}
  \item{control}{Control of the iterations. See \code{\link{glm.control}}}
  \item{method}{Which optimizer? 0 = "lbfgsb" (default), 1 = "vmmin".}
  \item{intercept}{Logical. If TRUE, an intercept is fitted.}
  \item{boot}{Integer. If > 0, bootstrapping with \code{boot}
    replicates.}
  \item{prior}{Which prior distribution? 0 for "gaussian", 1 for "logistic".}
}
\details{
"vmmin" is followed by some Newton-Raphson steps, until convergence. As
a by-product we get the estimated variance-covariance matrix.
}
\value{
  A list. For details, see the code, and \code{glmmML}.
}
\references{Broström (2003)}
\author{Göran Broström}
\note{A preliminary version, with high potential for bugs. However, when
it works, it is very fast, compared to \code{glmm} and \code{glmmPQL}}

\seealso{\code{\link{glmmML}}, \code{\link[MASS]{glmmPQL}}, and \code{\link[repeated]{glmm}}}

\examples{
x <- cbind(rep(1, 14), rnorm(14))
y <- rbinom(14, prob = 0.5, size = 1)
id <- rep(1:7, 2)

glmmML.fit(x, y, cluster = id, mixed = TRUE, method = 1)


}
\keyword{regression}% at least one, from doc/KEYWORDS
