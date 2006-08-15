extractAIC.glmmML <- function(fit, scale, k, ...){
    c(length(fit$coefficients) + 1, fit$aic)
}
