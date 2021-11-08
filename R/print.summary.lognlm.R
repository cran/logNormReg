print.summary.lognlm <-
function(x, digits = max(3L, getOption("digits") - 3L),  signif.stars = getOption("show.signif.stars"), ...) {
    cat("\nCall:\n", paste(deparse(x$call), sep = "\n", collapse = "\n"), 
        "\n", sep = "")
#    cat("Deviance Residuals: \n")
#    if (x$df.residual > 5) {
#        x$deviance.resid <- setNames(quantile(x$deviance.resid, 
#            na.rm = TRUE), c("Min", "1Q", "Median", "3Q", "Max"))
#    }
#    xx <- zapsmall(x$deviance.resid, digits + 1L)
#    print.default(xx, digits = digits, na.print = "", print.gap = 2L)
    cat("\nCoefficients:\n")
    coefs <- x$coefficients
    printCoefmat(coefs, digits = digits, signif.stars = signif.stars, na.print = "NA", ...)
            
#    if (length(x$aliased) == 0L) {
#        cat("\nNo Coefficients\n")
#    }
#    else {
#        df <- if ("df" %in% names(x)) x[["df"]]
#        else NULL
#        if (!is.null(df) && (nsingular <- df[3L] - df[1L])) 
#        cat("\nCoefficients: (", nsingular, " not defined because of singularities)\n", sep = "")
#        else cat("\nCoefficients:\n")
#        coefs <- x$coefficients
#        if (!is.null(aliased <- x$aliased) && any(aliased)) {
#            cn <- names(aliased)
#            coefs <- matrix(NA, length(aliased), 4L, dimnames = list(cn, 
#                colnames(coefs)))
#            coefs[!aliased, ] <- x$coefficients
#        }
#        printCoefmat(coefs, digits = digits, signif.stars = signif.stars, 
#            na.print = "NA", ...)
#    }
#
  if(x$lik) {
      Fnobj<- "Log Likelihood:"
  } else {
      Fnobj<-if(length(x$weights)<=0) "Sum of squared Residuals (logs):" else "Sum of (weighted) squared residuals (logs):" 
  }
#  Fnobj<- paste(Fnobj, " (on",,"")  
  V<-x$cov
	se.sd<-if(nrow(V)==(nrow(coefs)+1)) sqrt(V[nrow(coefs)+1,nrow(coefs)+1]) else NA
    cat("\nStandard deviation estimate: ", format(x$sigma, digits=max(5L, digits)), 
        "(St.Err =", paste(format(se.sd, digits=max(4, digits)),")", sep="")) 
    cat("\n")
    if (nzchar(mess <- naprint(x$na.action))) cat("  (", mess, ")\n", sep = "")
    cat(Fnobj, format(x$loglik, digits = max(5L, digits + 1L)), " (on", x$df.residual ,"degrees of freedom)",
        "\npseudo-R2:", formatC(x$r.squared, digits = digits), " Adj pseudo-R2:", formatC(x$adj.r.squared, digits = digits)
        )
    cat("\n")
    invisible(x)
}


