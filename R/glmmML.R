glmmML <- function(formula,
                   family = binomial,
                   data,
                   cluster,
                   subset,
                   na.action,
                   offset,
                   start.coef = NULL,
                   start.sigma = NULL,
                   control = list(epsilon = 1.e-8,
                       maxit = 200, trace = FALSE),
                   n.points = 16,
                   boot = 0){

    if (is.list(control)) {
        if (is.null(control$epsilon))
          control$epsilon <- 1e-08
        if (is.null(control$maxiter))
          control$maxit <- 200
        if (is.null(control$trace))
          control$trace <- FALSE
    }
    else {
        stop("control must be a list")
    }

    method <- 1 ## Always vmmin! 1 if vmmin, 0 otherwise
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
    ## get a copy of the call; result: a list.
    
    mf$family <- mf$start.coef <- mf$start.sigma <- NULL
    mf$control <- mf$maxit <- mf$boot <- NULL
    mf$n.points <- mf$method <- mf$start.coef <- mf$start.sigma <- NULL
    mf[[1]] <- as.name("model.frame") # turn into a call to model.frame
    mf <- eval(mf, environment(formula)) # run model.frame
    
    ## Pick out the parts.
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
 
    cluster <- mf$"(cluster)"

    no.cluster <- (missing(cluster) || is.null(cluster) ||
                   (length(unique(cluster)) <= 1))
    if (no.cluster){
        warning("No (or constant) 'cluster'; consider using 'glm'")
        return(NULL)
    }


    ##    return(clus)
    
    if (NCOL(Y) >  1) stop("Response must be univariate")
    
    if (!is.null(offset) && length(offset) != NROW(Y)) 
        stop(paste("Number of offsets is", length(offset), ", should equal", 
                   NROW(Y), "(number of observations)"))
    
    mixed <- TRUE # Not used any more, here for compatility.

    if (n.points <= 0) n.points <- 1 # Should give 'Laplace'(?)
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
                      intercept = ( attr(mt, "intercept") > 0),
                      boot
                      )
    
    if (!fit$convergence)
      warning("'vmmin' did not converge. Increase 'maxit'?")
##    if (fit$info) return(list(info = fit$info,
##                              convergence = fit$convergence,
##                              sigma = fit$sigma,
##                              coefficients = fit$beta,
##                              deviance = fit$deviance)
##                         )
    bdim <- p + 1
    res <- list()
    res$boot <- boot
    res$converged <- as.logical(fit$convergence)
    res$coefficients <- fit$beta
    res$coef.sd <- fit$beta.sd
    res$sigma <- fit$sigma
    res$sigma.sd <- fit$sigma.sd
    if (fit$cluster.null.deviance <= fit$deviance){
          res$sigma = 0
          res$sigma.sd = NA
      }
    res$variance <- fit$variance
    res$aic <- fit$aic
    names(res$coef.sd) <- names(res$coefficients)
    
    res$bootP <- fit$bootP
    res$deviance <- fit$deviance
    res$mixed <- TRUE
    res$df.residual <- fit$df.residual
    res$cluster.null.deviance <- fit$cluster.null.deviance
    res$cluster.null.df <- fit$cluster.null.df
    if (mixed){
        res$posterior.modes <- fit$post.mode
        res$posterior.means <- fit$post.mean
    }
    res$call <- cl
    names(res$coefficients) <- c(colnames(X))
    class(res) <- "glmmML"
    res
}

