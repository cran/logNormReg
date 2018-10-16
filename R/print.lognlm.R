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
        Fnobj<- if(x$lik) "\nLog Lik: " else "\nSum log Res:"
        cat("\nOptimizer:", paste(x$opt,";",sep=""), " Convergence code: " ,x$convergence, Fnobj , paste(format(x$loglik,digits=digits),";",sep=""), " Max grad (abs value):", max(abs(x$grad)),"\n")
        
    }
    else cat("No coefficients\n")
    cat("\n")
    invisible(x)
}
