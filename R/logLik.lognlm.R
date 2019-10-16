logLik.lognlm <-
  function(object, full=FALSE, ...){ 
    if(!full) {
      ll<- object$loglik 
    } else { 
      if(!object$lik) stop("full=TRUE is meaningless with the fit 'object' ")
      mu<- fitted(object)
      s2<- object$s2
      ll<- sum(dlnorm(object$y, log(mu)-s2/2, sqrt(s2), log=TRUE))
    }
    attr(ll, "df")<-length(object$coefficients) +1
    attr(ll, "nobs")<-length(object$residuals)
    ll
  }
