glmmML <- function(formula,
                   data = list(),
                   cluster = NULL,
                   family = binomial,
                   start.coef = NULL,
                   start.sigma = NULL,
                   offset = NULL,
                   method = "vmmin",
                   control = glm.control(epsilon = 1.e-8,
                     maxit = 100, trace = FALSE),
                   n.points = 16){

  method <- (method == "vmmin") ## 1 if vmmin, 0 otherwise
  if (!method) stop("Use default method (the only available at present)")
  cl <- match.call()

  if (is.character(family)) 
    family <- get(family)
  if (is.function(family)) 
    family <- family()
  if (is.null(family$family)) {
    print(family)
    stop("`family' not recognized")
  }

  if (missing(data))
    data <- environment(formula)
  
  mf <- match.call(expand.dots = FALSE)
                     # get a copy of the call; result: a list.

  mf$family <- mf$start.coef <- mf$start.sigma <- NULL
  mf$control <- mf$maxit <- NULL
  mf$n.points <- mf$method <- mf$start.coef <- mf$start.sigma <- NULL
  mf[[1]] <- as.name("model.frame") # turn into a call to model.frame
  mf <- eval(mf, environment(formula)) # run model.frame

  # Pick out the parts.
  mt <-  attr(mf, "terms")

  
  xvars <- as.character(attr(mt, "variables"))[-1]
  if ((yvar <- attr(mt, "response")) > 0) 
    xvars <- xvars[-yvar]
  xlev <- if (length(xvars) > 0) {
    xlev <- lapply(mf[xvars], levels)
    xlev[!sapply(xlev, is.null)]
  }
  
  X <- if (!is.empty.model(mt)) 
    model.matrix(mt, mf, contrasts)

  p <- NCOL(X)
 
  Y <- model.response(mf, "numeric")
  offset <- model.offset(mf)

  if (NCOL(Y) >  1) stop("Response must be univariate")
  
  if (!is.null(offset) && length(offset) != NROW(Y)) 
    stop(paste("Number of offsets is", length(offset), ", should equal", 
               NROW(Y), "(number of observations)"))

  mixed <- ( !is.null(cluster) ) && ( n.points >= 2 )

  fit <- glmmML.fit(X, Y,
                    start.coef,
                    start.sigma,
                    mixed,
                    cluster,
                    offset,
                    family,
                    n.points,
                    control,
                    method,
                    intercept = ( attr(mt, "intercept") > 0) )
                    
  if (!fit$convergence) return(list(convergence = fit$convergence))
  bdim <- p + 1
  res <- list()

  res$convergence <- as.logical(fit$convergence)
  res$aic <- -2 * fit$loglik + 2 * (p + as.integer(mixed))
  res$variance <- fit$coef.variance
  if (mixed){
    res$sigma <- fit$sigma
    res$sigma.sd <- sqrt(fit$sigma.vari)
  }else{
    res$sigma = 0
    res$sigma.sd = 0
  }
  res$coefficients <- fit$beta
  names(res$coefficients) <- c(colnames(X))
  res$deviance <- fit$deviance
##   options(show.error.messages = FALSE)
##  vari <- try(solve(-res$hessian))
#  if(is.numeric(vari)){
#    se <- sqrt(diag(vari))
#  }else{
#    se <- rep(NA, p + 1)
#  }
  res$call <- cl
  res$df.residual <- fit$df.residual
  res$sd <- sqrt(diag(res$variance))
  names(res$sd) <- names(res$coefficients)
  res$mixed <- mixed
  if (mixed){
    res$frail <- fit$frail
  }
  class(res) <- "glmmML"
  res
}

print.glmmML <- function(x,
                         digits = max(3, getOption("digits") - 3),
                         na.print = "",
                         ...){ 
  
  cat("\nCall: ", deparse(x$call), "\n\n")
  savedig <- options(digits = digits)
  on.exit(options(savedig))
  coef <- x$coefficients
  se <- x$sd
  tmp <- cbind(coef,
               se,
               coef/se,
               signif(1 - pchisq((coef/se)^2, 1), digits - 1)
               )
  dimnames(tmp) <- list(names(coef),
                        c("coef", "se(coef)", "z", "Pr(>|z|)")
                        )
  cat("\n")
  prmatrix(tmp)

  if(x$mixed){
    cat("\nStandard deviation in mixing distribution: ", x$sigma,  "\n")
    cat("Std. Error:                                ", x$sigma.sd, "\n")
  }
  cat("\nResidual deviance:",
      format(signif(x$deviance, digits)), "on",
      x$df.residual, "degrees of freedom", 
      "\tAIC:",
      format(signif(x$aic, digits)), "\n")
}
