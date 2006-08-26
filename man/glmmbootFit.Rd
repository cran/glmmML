\encoding{UTF-8}
\name{glmmbootFit}
\alias{glmmbootFit}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Generalized Linear Models with fixed effects grouping}
\description{
'glmmbootFit' is the workhorse in the function \code{glmmboot}. It is
suitable to call instead of 'glmmboot', e.g. in simulations.
}
\usage{
glmmbootFit(X, Y, weights = rep(1, NROW(Y)),
start.coef = NULL, cluster = rep(1, length(Y)),
offset = rep(0, length(Y)), family = binomial(),
control = list(epsilon = 1.e-8, maxit = 200, trace
= FALSE), boot = 0)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{The design matrix (n * p).}
  \item{Y}{The response vector of length n.}
  \item{weights}{Case weights.}
  \item{start.coef}{start values for the parameters in the linear
    predictor (except the intercept).}
  \item{cluster}{Factor indicating which items are correlated.}
  \item{offset}{this can be used to specify an a priori known component to be
    included in the linear predictor during fitting.}
  \item{family}{Currently, the only valid values are \code{binomial} and
  \code{poisson}. The binomial family allows for the \code{logit} and
  \code{cloglog} links.}
  \item{control}{A list. Controls the convergence criteria. See
  \code{\link{glm.control}} for details.} 
  \item{boot}{number of bootstrap replicates. If equal to zero, no test
    of significance of the grouping factor is performed.}
}
\details{

}
\value{
A list with components
  \item{coefficients}{Estimated regression coefficients (note: No intercept).}
  \item{logLik}{The maximised log likelihood.}
  \item{cluster.null.deviance}{deviance from a moddel without cluster.}
  \item{frail}{The estimated cluster effects.}
  \item{bootLog}{The maximised bootstrap log likelihood values. A vector
    of length \code{boot}.}
  \item{bootP}{The bootstrap p value.}
  \item{variance}{The variance-covariance matrix of the fixed effects
    (no intercept).}
  \item{sd}{The standard errors of the \code{coefficients}.}
  \item{boot_rep}{The number of bootstrap replicates.}
}
%\references{}
\author{\enc{Göran Broström}{Goran Brostrom}}
\note{
A profiling approach is used to estimate the cluster effects. 
}
\seealso{\code{\link{glmmboot}}}
\examples{
## Not run
x <- matrix(rnorm(1000), ncol = 1)
id <- rep(1:100, rep(10, 100))
y <- rbinom(1000, size = 1, prob = 0.4)
fit <- glmmbootFit(x, y, cluster = id, boot = 2000)
summary(fit)
## End(Not run)
## Should show no effects.
}
\keyword{regression}% at least one, from doc/KEYWORDS
\keyword{nonlinear}% __ONLY ONE__ keyword per line
