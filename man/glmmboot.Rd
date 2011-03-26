\encoding{UTF-8}
\name{glmmboot}
\alias{glmmboot}
\title{Generalized Linear Models with fixed effects grouping}
\description{
  Fits grouped GLMs with fixed group effects. The significance of the
  grouping is tested by simulation, with a bootstrap approach.
}
\usage{
glmmboot(formula, family = binomial, data, cluster, weights, subset, na.action,
offset, start.coef = NULL,
control = list(epsilon = 1e-08, maxit = 200, trace = FALSE), boot = 0)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{formula}{a symbolic description of the model to be fit. The details of
          model specification are given below.}
  \item{family}{Currently, the only valid values are \code{binomial} and
  \code{poisson}. The binomial family allows for the \code{logit} and
  \code{cloglog} links.}
  \item{data}{an optional data frame containing the variables in the model.
           By default the variables are taken from
          `environment(formula)', typically the environment from which
          `glmmML' is called.}
\item{cluster}{Factor indicating which items are correlated.}
\item{weights}{Case weights.}
  \item{subset}{an optional vector specifying a subset of observations
    to be used in the fitting process.}
  \item{na.action}{See glm.}
  \item{offset}{this can be used to specify an a priori known component to be
    included in the linear predictor during fitting.}
  \item{start.coef}{starting values for the parameters in the linear predictor.
 Defaults to zero.}
  \item{control}{Controls the convergence criteria. See
  \code{\link{glm.control}} for details.} 
  \item{boot}{number of bootstrap replicates. If equal to zero, no test
    of significance of the grouping factor is performed.}
}
\details{
The simulation is performed by 
  simulating new response vectors from the fitted probabilities without
  clustering, and comparing the maximized log likelihoods. The
maximizations are performed by profiling out the grouping factor. It is
a very fast procedure, compared to \code{\link{glm}}, when the grouping
factor has many levels.}
\value{
  The return value is a list, an object of class 'glmmboot'.
  \item{coefficients}{Estimated regression coefficients}
  \item{logLik}{the max log likelihood}
  \item{cluster.null.deviance}{Deviance without the clustering}
  \item{frail}{The estimated cluster effects}
  \item{bootLog}{The logLik values from the bootstrap samples}
  \item{bootP}{Bootstrap p value}
  \item{variance}{Variance covariance matrix}
  \item{sd}{Standard error of regression parameters}
  \item{boot_rep}{No. of bootstrap replicates}
  \item{mixed}{Logical}
  \item{deviance}{Deviance}
  \item{df.residual}{Its degrees of freedom}
  \item{aic}{AIC}
  \item{boot}{Logical}
  \item{call}{The function call}
}

%\references{ ~put references to the literature/web site here ~ }
\author{\enc{Göran Broström}{George Bridgewater}}

\note{There is no overall intercept for this model; each cluster has its
own intercept. See \code{frail}}

\seealso{\code{link{glmmML}}, \code{\link{optim}},
  \code{\link[lme4]{lmer}} in \code{Matrix}, and 
\code{\link[MASS]{glmmPQL}} in \code{MASS}.}

\examples{
## Not run:
id <- factor(rep(1:20, rep(5, 20)))
y <- rbinom(100, prob = rep(runif(20), rep(5, 20)), size = 1)
x <- rnorm(100)
dat <- data.frame(y = y, x = x, id = id)
res <- glmmboot(y ~ x, cluster = id, data = dat, boot = 5000)
## End(Not run)
##system.time(res.glm <- glm(y ~ x + id, family = binomial))
}
\keyword{regression}% at least one, from doc/KEYWORDS
\keyword{nonlinear}
