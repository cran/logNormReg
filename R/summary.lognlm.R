summary.lognlm <-
function(object, ...){
          n<-length(object$y)
          coef.p <- object$coefficients
          p<-length(coef.p)
          V<-vcov(object, se=FALSE, ...)
          s.err <- sqrt(diag(V[1:p,1:p]))
          tvalue <- coef.p/s.err
          dn <- c("Estimate", "Std. Error")
          df.r<- n-p-1
          pvalue <- 2 * pt(-abs(tvalue), df.r) #pnorm
          coef.table <- cbind(coef.p, s.err, tvalue, pvalue)
          dimnames(coef.table) <- list(names(coef.p), c(dn, "t value", "Pr(>|t|)"))
#          ans<-list(coefficients = coef.table, sigma=sqrt(object$s2), cov=V, call=object$call, loglik=object$loglik)
          
          keep <- match(c("call", "terms", "contrasts", "df.residual", "iter", "na.action"), names(object), 0L)
		      ans <- c(object[keep], list(coefficients = coef.table, sigma=sqrt(object$s2), cov=V, loglik=object$loglik), lik=object$lik)

#          ans <- c(object[keep], list(deviance.resid = residuals(object, 
#              type = "deviance"), coefficients = coef.table, aliased = aliased, 
#              dispersion = dispersion, df = c(object$rank, df.r, df.f), 
#              cov.unscaled = covmat.unscaled, cov.scaled = covmat))

          class(ans) <- "summary.lognlm"
          return(ans)
          }
