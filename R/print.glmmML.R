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
    if(x$boot){
        cat("\n Bootstrap p-value for fixed mixing: ",
        x$bootP, paste("(", x$boot.rep, ")\n", sep = ""))
    }
    cat("\nResidual deviance:",
        format(signif(x$deviance, digits)), "on",
        x$df.residual, "degrees of freedom", 
        "\tAIC:",
        format(signif(x$aic, digits)), "\n")
}
