glmmML.fit <- function (X, Y, 
                        start.coef = NULL, 
                        start.sigma = NULL,
                        mixed = FALSE,
                        cluster = NULL,                        
                        offset = rep(0, nobs),
                        family = binomial(),
                        n.points = 16,
                        control = glm.control(),
                        method,
                        intercept = TRUE,
                        boot = 0){
  
  X <- as.matrix(X)
  conv <- FALSE
  nobs <- NROW(Y)
  p <- NCOL(X)
  nvars <- p + as.integer(mixed)
  
  if (is.null(offset)) 
    offset <- rep(0, nobs)
  variance <- family$variance
  dev.resids <- family$dev.resids
  aic <- family$aic
  linkinv <- family$linkinv
  mu.eta <- family$mu.eta
  
  if (!is.function(variance) || !is.function(linkinv)) 
    stop("illegal `family' argument")

  if (is.null(start.coef)){
    start.coef <- numeric(p) # Start values equal to zero,
    if (family$family == "binomial"){
      start.coef[1] <- log(mean(Y) / (1 - mean(Y)))
    }else if (family$family == "poisson"){
      start.coef[1] <- log(mean(Y))
    }else{ ## this is a proviso!!
      start.coef[1] <- mean(Y)
    }
                           
  }else{                   
    if (length(start.coef) != p) stop("beta.start has wrong length")
  }
  
  if (mixed) {
    if (is.null(start.sigma)){
      start.sigma <- 0.5 ## More sofisticated choice is = ?
    }else{                  
      if (length(start.sigma) != 1) stop("sigma.start has wrong length")
    }
  }else{
    if (length(start.coef) != p) stop("beta.start has wrong length")
    n.points <- 1
    if (is.null(start.sigma)) start.sigma <- 0
  }
  
  ord <- order(cluster)
  Y <- Y[ord]
  X <- X[ord, ,drop = FALSE]

  ## Center the covariates so we avoid (some) numeric problems: 
  if (intercept){
    if (p >= 2){
      means <- numeric(p-1)
      for (i in 2:p){
        means[i-1] <- mean(X[, i])
        X[, i] <- X[, i] - means[i-1]
      }
    }
  }else{
    means <- numeric(p)
    for (i in 1:p){
      means[i] <- mean(X[, i])
      X[, i] <- X[, i] - means[i]
    }
  }
  
  cluster <- cluster[ord]
  fam.size <- as.vector(table(cluster))
  n.fam <- length(fam.size)
  
  glmFit <- glm.fit(X, Y,
                    start = start.coef,
                    offset = offset,
                    family = family,
                    control = control##,
                    ##intercept = with.intercept,
                    )
  predicted <- glmFit$fitted.values
  cluster.null.deviance <- glmFit$deviance
  
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
              
  fit <- .C("glmm_ml",
            as.integer(fam),
            as.integer(method),
            as.integer(p),
            as.double(start.coef),
            as.integer(cluster),
            as.double(start.sigma),
            as.double(t(X)),       ### Note CAREFULLY (03-01-09)!!!
            as.integer(Y),
            as.double(offset),
            as.integer(fam.size),
            as.integer(n.fam),
            as.integer(n.points),
            as.double(control$epsilon),
            as.integer(control$maxit),
            as.integer(control$trace),
            as.integer(boot),
            as.double(predicted),
            beta = double(p),  ## Return values from here.
            sigma = double(1),
            loglik = double(1),
            variance = double((p + 1) * (p + 1)),
            frail = double(n.fam),
            mu = double(nobs),
            bootP = double(1),
            bootLog = double(boot), 
            convergence = integer(1),
            info = integer(1),
            PACKAGE = "glmmML"
            )  

  vari <- matrix(fit$variance, ncol = (p + 1))
  ## Correct the estimate of the intercept for the centering:
  if (intercept){
    if (p >= 2){
      fit$beta[1] <- fit$beta[1] - sum(fit$beta[2:p] * means)
      aa <- numeric(p)
      aa[1] <- 1.0
      for (i in 2:p){ 
        aa[i] <- -means[i-1]
        ## "Restore" X (to what use?!):
        X[, i] <- X[, i] + means[i-1]
      }
      vari[1, 1] <- aa %*% vari[1:p, 1:p] %*% aa
    }
  }else{
    for (i in 1:p){ ## "Restore" mm (to what use?!):
      X[, i] <- X[, i] + means[i]
    }
  }
  
  if (mixed){
    sigma.vari <- vari[(p + 1), (p + 1)] * fit$sigma * fit$sigma
  }else{
    sigma.vari <- NULL
  }
  
  aic.model <- -2 * fit$loglik + 2 * nvars
  if (boot){
      bootP <- fit$bootP
  }else{
      bootP <- NA
  }
  list(beta = fit$beta,
       sigma = fit$sigma,
       loglik = fit$loglik,
       coef.variance = vari[1:p, 1:p, drop = FALSE],
       sigma.variance = sigma.vari,
       bootP = bootP,
       frail = fit$frail,
       residuals = residuals,
       fitted.values = fit$mu, 
       family = family, 
       deviance = -2*fit$loglik,
       aic = aic.model, 
       #null.deviance = nulldev,
       df.residual = NROW(Y) - NCOL(X) - as.integer(mixed),
       df.null = NROW(Y) - as.integer(intercept),
       #y = y,
       convergence = fit$convergence,
       info = fit$info) # 'info' coming from 'nr_opt'!
}
