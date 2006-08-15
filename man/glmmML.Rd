\name{glmmML}
\alias{glmmML}

\title{Generalized Linear Models with random intercept}
\description{
Fits GLMs with random intercept by Maximum Likelihood and numerical
integration via Gauss-Hermite quadrature.
}
\usage{
glmmML(formula, family = binomial, data, cluster, weights, subset, na.action, 
offset, prior = c("gaussian", "logistic"), start.coef = NULL, start.sigma
= NULL, fix.sigma = FALSE, 
control = list(epsilon = 1e-08, maxit = 200, trace = FALSE),
method = c("Laplace", "ghq"), n.points = 1, boot = 0) 
}
%- maybe also `usage' for other objects documented here.
\arguments{
  \item{formula}{ a symbolic description of the model to be fit. The details of
          model specification are given below.}
 \item{family}{Currently, the only valid values are \code{binomial} and
  \code{poisson}. The binomial family allows for the \code{logit} and
  \code{cloglog} links, but can only be represented as binary data.}
  \item{data}{an optional data frame containing the variables in the model.
           By default the variables are taken from
          `environment(formula)', typically the environment from which
          `glmmML' is called.
}
\item{cluster}{Factor indicating which items are correlated.}
\item{weights}{Case weights. Defaults to one.}
  \item{subset}{an optional vector specifying a subset of observations
    to be used in the fitting process.}
  \item{na.action}{See glm.}
  \item{start.coef}{starting values for the parameters in the linear predictor.
 Defaults to zero.}
  \item{start.sigma}{starting value for the mixing standard
    deviation. Defaults to 0.5.}
  \item{fix.sigma}{Should sigma be fixed at start.sigma?}
  \item{offset}{this can be used to specify an a priori known component to be
          included in the linear predictor during fitting.}
\item{prior}{Which "prior" distribution (for the random effects)?
Possible choices are "gaussian" (default) and "logistic".} 
  \item{control}{Controls the convergence criteria. See
	  \code{\link{glm.control}} for details.}
  \item{method}{There are two choices "Laplace" (default) and "ghq"
    (gauss-hermite). This parameter is in fact unused; the chice is
    presently made by n.points. }
  \item{n.points}{Number of points in the Gauss-Hermite quadrature. If
    n.points == 1, the Gauss-Hermite is the same as Laplace approximation.}
  \item{boot}{Do you want a bootstrap estimate of cluster effect? The default
    is \emph{No} (\code{boot = 0}). If you want to say yes, enter a
    positive integer here. It should be equal to the number of bootstrap
    samples you want to draw. A recomended absolute \emph{minimum value} is
    \code{boot = 2000}.}
}
\details{
The integrals in the log likelihood function are evaluated by the
Laplace approximation (default) or Gauss-Hermite quadrature. The latter
is now fully adaptive; however, only approximate estimates of variances
are available for the Gauss-hermite (n.points > 1) method. Wil be
changed soon. 
}
\value{
  The return value is a list, an object of class 'glmmML'. The components are:
  \item{boot}{No. of boot replicates}
  \item{converged}{Logical}
  \item{coefficients}{Estimated regression coefficients}
  \item{coef.sd}{Their standard errors}
  \item{sigma}{The estimated random effects' standard deviation}
  \item{sigma.sd}{Its standard error}
  \item{variance}{The estimated variance-covariance matrix. The last
    column/row corresponds to the log of the standard
    deviation of the random effects (log(sigma))} 
  \item{aic}{AIC}
  \item{bootP}{Bootstrap p value from testing the null hypothesis of no
    random effect (sigma = 0)} 
  \item{deviance}{Deviance}
  \item{mixed}{Logical}
  \item{df.residual}{Degrees of freedom}
  \item{cluster.null.deviance}{Deviance from a glm with no clustering}
  \item{cluster.null.df}{Its degrees of freedom}
  \item{posterior.modes}{Estimated posterior modes of the random effects}
  \item{posterior.means}{Estimated posterior means of the random effects}
  \item{terms}{The terms object}
  \item{info}{From hessian inversion. Should be 0. If not, no variances
    could be estimated. You could try fixing sigma at the estimated
    value and rerun.}
  \item{prior}{Which prior was used?}
  \item{call}{The function call}
}
\references{Broström (2003). Generalized linear models with random
  intercepts. \url{http://www.stat.umu.se/forskning/reports/glmmML.pdf}
 }
\author{Göran Broström}
\note{The optimization may not converge with
the default value of \code{start.sigma}. In that case, try different
start values for sigma. If still no convergence, consider the
possibility to fix the value of sigma at several values and study the
profile likelihood.}

\seealso{\code{\link{glmmboot}}, \code{\link{optim}},
  \code{\link[repeated]{glmm}} in Lindsey's 
\code{repeated} package, \code{\link[Matrix]{lmer}} in \code{Matrix}and
\code{\link[MASS]{glmmPQL}} in \code{MASS}.} 

\examples{
id <- factor(rep(1:20, rep(5, 20)))
y <- rbinom(100, prob = rep(runif(20), rep(5, 20)), size = 1)
x <- rnorm(100)
dat <- data.frame(y = y, x = x, id = id)
glmmML(y ~ x, data = dat, cluster = id)
}
\keyword{regression}% at least one, from doc/KEYWORDS
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
