glmmboot <- function(formula,
                   family = binomial,
                   data,
                   cluster,
                   subset,
                   na.action,
                   offset,
                   start.coef = NULL,
                   control = glm.control(epsilon = 1.e-8,
                       maxit = 100, trace = FALSE),
                   boot = 0){

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
    
    if (NCOL(Y) >  1) stop("Response must be univariate")
    
    if (!is.null(offset) && length(offset) != NROW(Y)) 
        stop(paste("Number of offsets is", length(offset), ", should equal", 
                   NROW(Y), "(number of observations)"))
    
    ## Remove eventual intercept from X.
    ## Taken care of thru separate intercepts for each 'cluster'.
    if (!is.na(coli <- match("(Intercept)", colnames(X))))
        X <- X[, -coli, drop = FALSE]

    fortran <- TRUE
    res <- glmmbootFit(X, Y,
                       start.coef,
                       cluster,
                       offset,
                       family,
                       control,
                       method,
                       boot,
                       fortran)
    
    res$mixed <- FALSE
    res$deviance <- -2 * res$logLik
    nvars <- NCOL(X) + length(unique(cluster))
    res$df.residual <- length(Y) - nvars
    res$aic <- res$deviance + 2 * nvars
    res$boot <- TRUE
    res$call <- cl
    names(res$coefficients) <- c(colnames(X))
    class(res) <- "glmmboot"
    res
}

