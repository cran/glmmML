\name{glmmboot}
\alias{glmmboot}
\title{Generalized Linear Models with fixed effects grouping}
\description{
  Fits grouped GLMs with fixed group effects. The significance of the
  grouping is tested by simulation, with a bootstrap approach.
}
\usage{
glmmboot(formula, family = binomial, data, cluster, subset, na.action,
offset, start.coef = NULL,
control = glm.control(epsilon = 1e-08, maxit = 100, trace = FALSE), boot = 0)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{formula}{a symbolic description of the model to be fit. The details of
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
The simulation is performed by making random permutations of the
grouping factor and comparing the maximized loglikelihoods. The
maximizations are performed by profiling out the grouping factor. It is
a very fast procedure, compared to \code{\link{glm}}, when the grouping
factor has many levels. 
}
\value{
The return value is a list, an object of class 'glmmboot'.
}
\references{ ~put references to the literature/web site here ~ }
\author{Göran Broström}
\note{This is a preliminary version and not well tested.}

\seealso{\code{link{glmmML}}, \code{\link{optim}},
  \code{\link[repeated]{glmm}} in Lindsey's 
\code{repeated} package, \code{\link[lme4]{GLMM}} in \code{lme4}and
\code{\link[MASS]{glmmPQL}} in \code{MASS}.}
 
\examples{
id <- factor(rep(1:20, rep(5, 20)))
y <- rbinom(100, prob = rep(runif(20), rep(5, 20)), size = 1)
x <- rnorm(100)
dat <- data.frame(y = y, x = x, id = id)
res <- glmmboot(y ~ x, cluster = id, data = dat, boot = 5000)
##system.time(res.glm <- glm(y ~ x + id, family = binomial))
}
\keyword{regression}% at least one, from doc/KEYWORDS
\keyword{nonlinear}
