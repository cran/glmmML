glmmboot <- function(formula,
                     family = binomial,
                     data,
                     cluster,
                     subset,
                     na.action,
                     offset,
                     conditional = FALSE,
                     start.coef = NULL,
                     control = list(epsilon = 1.e-8,
                     maxit = 200, trace = FALSE),
                     boot = 0){

    if (is.list(control)) {
        if (is.null(control$epsilon))
          control$epsilon <- 1e-08
        if (is.null(control$maxit))
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
    mf$control <- mf$maxit <- mf$boot <- mf$conditional <- NULL
    mf$n.points <- mf$method <- mf$start.coef <- NULL
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

    if (NCOL(Y) >  1) stop("Response must be univariate")
    
    if (!is.null(offset) && length(offset) != NROW(Y)) 
        stop(paste("Number of offsets is", length(offset), ", should equal", 
                   NROW(Y), "(number of observations)"))
    
    ## Remove eventual intercept from X.
    ## Taken care of thru separate intercepts for each 'cluster'.
    ## Weck: NOW,  31 jan 2005, we change that! ##
    ## Weck: NOW,  6 jun 2006, we move it to glmmbootFit! ##
    
    ##if (!is.na(coli <- match("(Intercept)", colnames(X))))
    ##    X <- X[, -coli, drop = FALSE]

    fortran <- TRUE
    res <- glmmbootFit(X, Y,
                       start.coef,
                       cluster,
                       offset,
                       family,
                       conditional,
                       control,
                       boot,
                       method,
                       fortran)
    
    res$mixed <- FALSE # ??????????????????
    res$deviance <- -2 * res$logLik
    nvars <- NCOL(X) - 1 + length(unique(cluster))
    res$df.residual <- length(Y) - nvars
    res$aic <- res$deviance + 2 * nvars
    res$boot <- TRUE
    res$call <- cl
    if (!is.null(res$coefficients))
      names(res$coefficients) <- c(colnames(X))[-1]
    class(res) <- "glmmboot"
    res
}

