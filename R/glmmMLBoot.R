glmmbootFit <- function (X, Y, 
                         start.coef = NULL, 
                         cluster = rep(1, length(Y)),                        
                         offset = rep(0, length(Y)),
                         family = binomial(),
                         conditional = FALSE,
                         control = glm.control(),
                         method,
                         boot,
                         fortran = TRUE){

    X <- as.matrix(X)

    predicted <- glm.fit(X, Y,
                         start = start.coef,
                         offset = offset,
                         family = family,
                         control = control,
                         intercept = TRUE
                         )$fitted.values

    if (is.null(offset)) offset <- rep(0, length(Y))
    p <- ncol(X)
    if (is.null(start.coef)){
        start.coef <- numeric(p) # Start values equal to zero,
        ## if (FALSE){
        if (family$family == "binomial"){
            start.coef[1] <- log(mean(Y) / (1 - mean(Y)))
        }else if (family$family == "poisson"){
            start.coef[1] <- log(mean(Y))
        }else{ ## this is a proviso!!
            start.coef[1] <- mean(Y)
        }
    ##}
    }else{                   
        if (length(start.coef) != p) stop("beta.start has wrong length")
    }

    ord <- order(cluster)

    Y <- Y[ord]
    X <- X[ord, ,drop = FALSE]
    cluster <- cluster[ord]

    if (family$family == "binomial"){
        if (family$link == "logit"){
            fam <- 0
        }else if (family$link == "cloglog"){
            fam <- 1
        }else{
            stop("Unknown link function; only 'logit' and 'cloglog' implemented")
        }
    }else if (family$family == "poisson"){
        fam <- 2
    }else{
        stop("Unknown family; only 'binomial' and 'poisson' implemented")
    }
    
    famSize <- as.vector(table(cluster))
    nFam <- length(famSize)
  
    ## cat("nFam = ", nFam, "\n")

    if (fortran){
        if (p >= 1){
            means <- colMeans(X)
            X <- scale(X, center = TRUE, scale = FALSE)
    
            fit <- .C("glmm_boot",
                      as.integer(fam),
                      as.integer(method),
                      as.integer(p),
                      as.double(start.coef),
                      as.integer(cluster),
                      as.double(t(X)),       # Note! #
                      as.integer(Y),
                      as.double(offset),
                      as.integer(famSize),
                      as.integer(nFam),
                      as.integer(conditional),
                      as.double(control$epsilon),
                      as.integer(control$maxit),
                      as.integer(control$trace),
                      as.integer(boot),
                      beta = double(p),
                      predicted = as.double(predicted), # Watch up! #
                      loglik = double(1),
                      hessian = double(p * p),
                      frail = double(nFam),
                      bootP = double(1),
                      bootLog = double(boot),
                      convergence = integer(1)
                      )
            res <- list(coefficients = fit$beta,
                        logLik = fit$loglik,
                        frail = fit$frail,
                        bootLog = fit$bootLog,
                        bootP = fit$bootP)
            res$variance <- solve(-matrix(fit$hessian, nrow = p, ncol = p))
            res$sd <- sqrt(diag(res$variance))
            res$boot_rep <- boot

            return(res)
        }else{ # A null model:

            fit <- .C("glmm_boot0",
                      as.integer(fam),
                      as.integer(method),
                      ##as.integer(p),
                      ##as.double(start.coef),
                      as.integer(cluster),
                      ##as.double(t(X)),       ## Note! ##
                      as.integer(Y),
                      as.double(offset),
                      as.integer(famSize),
                      as.integer(nFam),
                      as.integer(conditional),
                      ##as.double(control$epsilon),
                      ##as.integer(control$maxit),
                      as.integer(control$trace),
                      as.integer(boot),
                      predicted = as.double(predicted),
                      ##beta = double(p),
                      loglik = double(1),
                      ##hessian = double(p * p),
                      frail = double(nFam),
                      bootP = double(1),
                      bootLog = double(boot),
                      convergence = integer(1)
                      )
            res <- list(coefficients = NULL,
                        logLik = fit$loglik,
                        frail = fit$frail,
                        bootLog = fit$bootLog,
                        bootP = fit$bootP)
            res$variance <- NULL
            res$sd <- NULL
            res$boot_rep <- boot

            return(res)
        }

    }else{
        
        fit <- snut(Y, X, cluster, offset)
        
        if (boot < 1) error("boot must be at least 1")
        logl <- numeric(boot + 1)
        logl[1] <- -fit$value
        
        for (i in 2:(boot + 1)){
            cat("****************** Replicate ", i-1, "\n")
            cluster <- sample(cluster)
            ord <- order(cluster)
            cluster <- cluster[ord]
            Y <- Y[ord]
            X <- X[ord, ,drop = FALSE]
            
            logl[i] <- -snut(Y, X, cluster, offset)$value
        }
        
        ##cat("logl = ", logl, "\n")
        p.value <- 1 - rank(logl, ties.method = "first")[1] / (boot + 1)  
        
        nvars <- length(unique(cluster))
        nvar <- length(fit$par)
        aic.model <- 2 * fit$value + 2 * nvars
        ## Note difference between nvar & nvars!
        n.fam <- nvar - p
        list(beta = fit$par[(n.fam + 1):nvar],
             loglik = -fit$value,
             ##coef.variance = vari[1:p, 1:p, drop = FALSE],
             frail = fit$par[1:n.fam],
             ##residuals = residuals,
             ##fitted.values = fit$mu, 
             ##family = family, 
             deviance = 2*fit$value,
             aic = aic.model, 
             ##null.deviance = nulldev,
             ##df.residual = NROW(Y) - NCOL(X) - as.integer(mixed),
             ##df.null = NROW(Y) - as.integer(intercept),
             ##y = y,
             ##convergence = fit$convergence
             frail.p = p.value
             )
    }
}
