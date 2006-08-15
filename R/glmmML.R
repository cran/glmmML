glmmML <- function(formula,
                   family = binomial,
                   data,
                   cluster,
                   weights,
                   subset,
                   na.action,
                   offset,
                   prior = c("gaussian", "logistic"),
                   start.coef = NULL,
                   start.sigma = NULL,
                   fix.sigma = FALSE,
                   control = list(epsilon = 1.e-8,
                       maxit = 200, trace = FALSE),
                   method = c("Laplace", "ghq"),
                   n.points = 1,
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

    ## don't use this for the moment!
    ## Not used:
    
    method <- as.numeric(method[1] == "Laplace")

    if (!prior[1] %in% c("gaussian", "logistic"))
      stop("Prior distribution not known")

    prior <- as.numeric(prior[1] == "logistic")
    ## 'gaussian' is the default
    
    ##if (!method) stop("Use default method (the only available at present)")
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
    
    mf$family <- mf$start.coef <- mf$start.sigma <- mf$fix.sigma <- NULL
    mf$weights <- mf$control <- mf$maxit <- mf$boot <- NULL
    mf$n.points <- mf$method <- mf$prior <- NULL
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
    
    ##if (NCOL(Y) >  1) stop("Response must be univariate")
    
    if (!is.null(offset) && length(offset) != NROW(Y)) 
        stop(paste("Number of offsets is", length(offset), ", should equal", 
                   NROW(Y), "(number of observations)"))
    
    
    if (missing(weights)) weights <- rep.int(1, NROW(Y))
    if (any(weights < 0)) stop("negative weights not allowed")
    mixed <- TRUE # Not used any more, here for compatility.

    if (n.points <= 0) n.points <- 1 # Should give 'Laplace'(?)
    fit <- glmmML.fit(X, Y,
                      weights,
                      start.coef,
                      start.sigma,
                      fix.sigma,
                      mixed,
                      cluster,
                      offset,
                      family,
                      n.points,
                      control,
                      method,
                      intercept = ( attr(mt, "intercept") > 0),
                      boot,
                      prior # gaussian by default
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
    res$sigma <- abs(fit$sigma) # Note 
    res$sigma.sd <- fit$sigma.sd
    ## For the time being: Show the attained max!!!!!!!!!!!!!
    ##if (fit$cluster.null.deviance <= fit$deviance){
    ##      res$sigma = 0
    ##      res$sigma.sd = NA
    ##  }
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
    res$prior <- if (prior) "logistic" else "gaussian"
    res$terms = mt
    res$info <- fit$info # From inverting the hessian! Should be zero.
    res$call <- cl
    names(res$coefficients) <- c(colnames(X))
    class(res) <- "glmmML"
    res
}

