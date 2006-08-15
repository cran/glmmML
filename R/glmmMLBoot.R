glmmbootFit <- function (X, Y, weights = rep(1, NROW(Y)), 
                         start.coef = NULL, 
                         cluster = rep(1, length(Y)),                        
                         offset = rep(0, length(Y)),
                         family = binomial(),
                         conditional = FALSE,
                         control = list(epsilon = 1.e-8, maxit = 200,
                           trace = FALSE), 
                         boot = 0,
                         method = 1,
                         fortran = TRUE){

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

    method <- 1
    fortran <- TRUE
    X <- as.matrix(X)

    nobs <- NROW(Y)
    if (family$family == "binomial"){ # This will be better later!
        ## From 'binomial':
        if (NCOL(Y) == 1) {
            if (is.factor(Y)) Y <- Y != levels(Y)[1]
            n <- rep.int(1, nobs)
            if (any(Y < 0 | Y > 1)) stop("Y values must be 0 <= Y <= 1")
            ##mustart <- (weights * Y + 0.5)/(weights + 1)
            m <- weights * Y
            if (any(abs(m - round(m)) > 0.001))
              warning("non-integer #successes in a binomial glm!")
        } else if (NCOL(Y) == 2) {
            if (any(abs(Y - round(Y)) > 0.001))
              warning("non-integer counts in a binomial glm!")
            n <- Y[, 1] + Y[, 2]
            Y <- ifelse(n == 0, 0, Y[, 1]/n)
            weights <- weights * n
            ##mustart <- (n * Y + 0.5)/(n + 1)
        } else
        stop("for the binomial family, Y must be a vector of 0 and 1's\n",
             "or a 2 column matrix where col 1 is no. successes and col 2 is no. failures")
        ## End of 'from binomial'.
    }
    coli <- match("(Intercept)", colnames(X))
    with.intercept <- !is.na(coli)

    if (is.null(offset)) offset <- rep(0, length(Y))
    glmFit <- glm.fit(X, Y, weights = weights,
                      start = start.coef,
                      offset = offset,
                      family = family,
                      control = control,
                      intercept = with.intercept,
                      )
    predicted <- glmFit$fitted.values
    cluster.null.deviance <- glmFit$deviance
    if (with.intercept)
      X <- X[, -coli, drop = FALSE]
    p <- ncol(X)
    if (is.null(start.coef)){
        start.coef <- numeric(p) # Start values equal to zero,
        if (FALSE){## Not sensible below:
            if (family$family == "binomial"){
                start.coef[1] <- log(mean(Y) / (1 - mean(Y)))
            }else if (family$family == "poisson"){
                start.coef[1] <- log(mean(Y))
            }else{ ## this is a proviso!!
                start.coef[1] <- mean(Y)
            }
        } ## End 'if (FALSE)'
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
                      as.double(weights),
                      as.double(t(X)),       # Note! #
                      as.double(Y),
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
                        cluster.null.deviance = cluster.null.deviance,
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
                      as.double(weights),       ## Note! ##
                      as.double(Y),
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
                        cluster.null.deviance = cluster.null.deviance,
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
            if (control$trace)
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
        aic.model <- 2 * fit$value + 2 * nvars ## STRANGE ?! (nvar?)
        ## Note difference between nvar & nvars!
        n.fam <- nvar - p
        list(beta = fit$par[(n.fam + 1):nvar],
             loglik = -fit$value,
             ##coef.variance = vari[1:p, 1:p, drop = FALSE],
             frail = fit$par[1:n.fam],
             ##residuals = residuals,
             ##fitted.values = fit$mu, 
             ##family = family,
             null.cluster.deviance = null.cluster.deviance,
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
