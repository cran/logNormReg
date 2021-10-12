print.lognlm <-
function (x, digits = max(3L, getOption("digits") - 3L), ...) {
    cat("\nCall:\n", paste(deparse(x$call), sep = "\n", collapse = "\n"), 
        "\n\n", sep = "")
    if (length(coef(x))) {
        cat("Coefficients:\n")
        print.default(format(coef(x), digits = digits), print.gap = 2L, quote = FALSE)
		cat("\nst.dev:\t", format(sqrt(x$s2),digits=digits),"\n")
#        cat("\nst.dev:\t")
#        print.default(format(sqrt(x$s2), digits = digits), print.gap = 3L, quote = FALSE)

#        cat("\ngradient at solution:\n")
#        print.default(format(drop(x$grad), digits = digits+3), print.gap = 5L, quote = FALSE)
        if(x$lik) {
            Fnobj<- "\nLog Likelihood: "
            } else {
                Fnobj<-if(length(x$weights)<=0) "\nSum of squared residuals (logs):" else "\nSum of (weighted) squared residuals (logs):" 
            }
        cat("\nOptimizer:", paste(x$opt,",",sep=""), " Converg code:" ,paste(x$convergence,",",sep=""), " Max abs grad:", signif(max(abs(x$grad)),2),
            Fnobj , paste(format(x$loglik,digits=digits+2)," ",sep=""), "\n")
        
    }
    else cat("No coefficients\n")
    cat("\n")
    invisible(x)
}
