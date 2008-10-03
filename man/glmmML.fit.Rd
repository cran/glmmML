\encoding{latin1}
\name{glmmML.fit}
\alias{glmmML.fit}

\title{Generalized Linear Model with random intercept}
\description{
This function is called by \code{glmmML}, but it can also be called
directly by the user.
}
\usage{
glmmML.fit(X, Y, weights = rep(1, NROW(Y)), cluster.weights = rep(1, NROW(Y)),
start.coef = NULL, start.sigma = NULL,
fix.sigma = FALSE,
cluster = NULL, offset = rep(0, nobs), family = binomial(),
method = 1, n.points = 1,
control = list(epsilon = 1.e-8, maxit = 200, trace = FALSE),
intercept = TRUE, boot = 0, prior = 0) 
}
%- maybe also `usage' for other objects documented here.
\arguments{
  \item{X}{Design matrix of covariates.}
  \item{Y}{Response vector. Or two-column matrix.}
  \item{weights}{Case weights. Defaults to one.}
  \item{cluster.weights}{Cluster weights. Defaults to one.}
  \item{start.coef}{Starting values for the coefficients. }
  \item{start.sigma}{Starting value for the mixing standard deviation.}
  \item{fix.sigma}{Should sigma be fixed at start.sigma?}
  \item{cluster}{The clustering variable.}
  \item{offset}{The offset in the model.}
  \item{family}{Family of distributions. Defaults to binomial with logit
  link. Other possibilities are binomial with cloglog link and poisson
  with log link.}
\item{method}{Laplace (1) or Gauss-hermite (0)?}
\item{n.points}{Number of points in the Gauss-Hermite
    quadrature. Default is \code{n.points = 1}, which is equivalent to
    Laplace approximation.}
  \item{control}{Control of the iterations. See \code{\link{glm.control}}.}
  \item{intercept}{Logical. If TRUE, an intercept is fitted.}
  \item{boot}{Integer. If > 0, bootstrapping with \code{boot}
    replicates.}
  \item{prior}{Which prior distribution? 0 for "gaussian", 1 for
    "logistic", 2 for "cauchy".}
}
\details{
  In the optimisation, "vmmin" (in C code) is used.
  }
\value{
  A list. For details, see the code, and \code{glmmML}.
}
\references{\enc{Broström}{Bridgewater} (2003)}
\author{\enc{Göran Broström}{George Bridgewater}}
%\note{The profiling method is very fast, compared to \code{glmm} and
%  \code{glmmPQL}} 

\seealso{\code{\link{glmmML}}, \code{\link[MASS]{glmmPQL}}, and
  \code{\link[Matrix]{lmer}}.} 

\examples{
x <- cbind(rep(1, 14), rnorm(14))
y <- rbinom(14, prob = 0.5, size = 1)
id <- rep(1:7, 2)

glmmML.fit(x, y, cluster = id)


}
\keyword{regression}% at least one, from doc/KEYWORDS
