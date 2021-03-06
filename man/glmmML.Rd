\encoding{latin1}
\name{glmmML}
\alias{glmmML}

\title{Generalized Linear Models with random intercept}
\description{
Fits GLMs with random intercept by Maximum Likelihood and numerical
integration via Gauss-Hermite quadrature.
}
\usage{
glmmML(formula, family = binomial, data, cluster, weights,
cluster.weights, subset, na.action, 
offset, contrasts = NULL, prior = c("gaussian", "logistic", "cauchy"),
start.coef = NULL, start.sigma = NULL, fix.sigma = FALSE, x = FALSE, 
control = list(epsilon = 1e-08, maxit = 200, trace = FALSE),
method = c("Laplace", "ghq"), n.points = 8, boot = 0) 
}
%- maybe also `usage' for other objects documented here.
\arguments{
  \item{formula}{ a symbolic description of the model to be fit. The details of
          model specification are given below.}
 \item{family}{Currently, the only valid values are \code{binomial} and
  \code{poisson}. The binomial family allows for the \code{logit} and
  \code{cloglog} links.}
  \item{data}{an optional data frame containing the variables in the model.
           By default the variables are taken from
          `environment(formula)', typically the environment from which
          `glmmML' is called.
}
\item{cluster}{Factor indicating which items are correlated.}
\item{weights}{Case weights. Defaults to one.}
\item{cluster.weights}{Cluster weights. Defaults to one.}
  \item{subset}{an optional vector specifying a subset of observations
    to be used in the fitting process.}
  \item{na.action}{See glm.}
  \item{start.coef}{starting values for the parameters in the linear predictor.
 Defaults to zero.}
  \item{start.sigma}{starting value for the mixing standard
    deviation. Defaults to 0.5.}
  \item{fix.sigma}{Should sigma be fixed at start.sigma?}
  \item{x}{If TRUE, the design matrix is returned (as x).}
  \item{offset}{this can be used to specify an a priori known component to be
    included in the linear predictor during fitting.}
  \item{contrasts}{an optional list. See the 'contrasts.arg' of 'model.matrix.default'.}
\item{prior}{Which "prior" distribution (for the random effects)?
Possible choices are "gaussian" (default), "logistic", and "cauchy".} 
  \item{control}{Controls the convergence criteria. See
	  \code{\link{glm.control}} for details.}
  \item{method}{There are two choices "Laplace" (default) and "ghq"
    (Gauss-Hermite). }
  \item{n.points}{Number of points in the Gauss-Hermite quadrature. If
    n.points == 1, the Gauss-Hermite is the same as Laplace
    approximation. If \code{method} is set to "Laplace", this parameter
    is ignored.}
  \item{boot}{Do you want a bootstrap estimate of cluster effect? The default
    is \emph{No} (\code{boot = 0}). If you want to say yes, enter a
    positive integer here. It should be equal to the number of bootstrap
    samples you want to draw. A recomended absolute \emph{minimum value} is
    \code{boot = 2000}.}
}
\details{The integrals in the log likelihood function are evaluated by the
Laplace approximation (default) or Gauss-Hermite quadrature. The latter
is now fully adaptive; however, only approximate estimates of variances
are available for the Gauss-Hermite (n.points > 1) method.

For the binomial families, the response can be a two-column matrix, see
the help page for glm for details.}
\value{
  The return value is a list, an object of class 'glmmML'. The components are:
  \item{boot}{No. of boot replicates}
  \item{converged}{Logical}
  \item{coefficients}{Estimated regression coefficients}
  \item{coef.sd}{Their standard errors}
  \item{sigma}{The estimated random effects' standard deviation}
  \item{sigma.sd}{Its standard error}
  \item{variance}{The estimated variance-covariance matrix. The last
    column/row corresponds to the standard
    deviation of the random effects (\code{sigma})} 
  \item{aic}{AIC}
  \item{bootP}{Bootstrap p value from testing the null hypothesis of no
    random effect (sigma = 0)} 
  \item{deviance}{Deviance}
  \item{mixed}{Logical}
  \item{df.residual}{Degrees of freedom}
  \item{cluster.null.deviance}{Deviance from a glm with no
    clustering. Subtracting \code{deviance} gives a test statistic for
    the null hypothesis of no clustering. Its asymptotic distribution is
  a symmetric mixture a constant at zero and a chi-squared distribution
  with one df. The printed p-value is based on this.}
  \item{cluster.null.df}{Its degrees of freedom}
  \item{posterior.modes}{Estimated posterior modes of the random effects}
  \item{terms}{The terms object}
  \item{info}{From hessian inversion. Should be 0. If not, no variances
    could be estimated. You could try fixing sigma at the estimated
    value and rerun.}
  \item{prior}{Which prior was used?}
  \item{call}{The function call}
  \item{x}{The design matrix if asked for, otherwise not present}
}
\references{Brostr\"om, G. and Holmberg, H. (2011). Generalized linear models with
  clustered data: Fixed and random effects models. Computational
  Statistics and Data Analysis 55:3123-3134.
 }
\author{\enc{G\"oran Brostr\"om}{George Bridgewater}}
\note{The optimization may not converge with
the default value of \code{start.sigma}. In that case, try different
start values for sigma. If still no convergence, consider the
possibility to fix the value of sigma at several values and study the
profile likelihood.}

\seealso{\code{\link{glmmboot}}, \code{\link{glm}}, \code{\link{optim}},
  \code{\link[lme4]{lmer}} in \code{Matrix}and
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
