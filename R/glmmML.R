glmmML <- function(formula,
                   family = binomial,
                   data,
                   cluster,
                   subset,
                   na.action,
                   offset,
                   start.coef = NULL,
                   start.sigma = NULL,
                   control = glm.control(epsilon = 1.e-8,
                       maxit = 100, trace = FALSE),
                   n.points = 16,
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
    ##    return(clus)
    
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
                      intercept = ( attr(mt, "intercept") > 0),
                      boot
                      )
    
    if (!fit$convergence) return(list(convergence = fit$convergence))
    if (fit$info) return(list(info = fit$info,
                              convergence = fit$convergence,
                              sigma = fit$sigma,
                              coefficients = fit$beta,
                              deviance = fit$deviance)
                         )
    bdim <- p + 1
    res <- list()

    res$boot.rep<- boot
    res$boot <- (boot > 0)
    res$convergence <- as.logical(fit$convergence)
    res$aic <- -2 * fit$loglik + 2 * (p + as.integer(mixed))
    res$variance <- fit$coef.variance
    if (mixed){
        res$sigma <- fit$sigma
        res$sigma.sd <- fit$sigma.vari
    }else{
        res$sigma = 0
        res$sigma.sd = 0
    }
    res$bootP <- fit$bootP
    res$coefficients <- fit$beta
    res$deviance <- fit$deviance
    ##   options(show.error.messages = FALSE)
    ##  vari <- try(solve(-res$hessian))
    ##  if(is.numeric(vari)){
    ##    se <- sqrt(diag(vari))
    ##  }else{
    ##    se <- rep(NA, p + 1)
    ##  }
    res$df.residual <- fit$df.residual
    res$sd <- sqrt(diag(res$variance))
    names(res$sd) <- names(res$coefficients)
    res$mixed <- mixed
    if (mixed){
        res$frail <- fit$frail
    }
    res$call <- cl
    names(res$coefficients) <- c(colnames(X))
    class(res) <- "glmmML"
    res
}

