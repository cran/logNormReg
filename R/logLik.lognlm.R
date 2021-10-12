logLik.lognlm <-
  function(object, full=FALSE, ...){ 
    if(!full) {
      ll<- object$loglik 
    } else { 
      if(!object$lik) stop("full=TRUE is meaningless with the non Lik-based fit 'object' ")
      mu<- fitted(object)
      n<-length(mu)
      p<-length(object$coefficients)
      s2<- object$s2*(n-p)/n 
      ll<- sum(dlnorm(object$y, log(mu)-s2/2, sqrt(s2), log=TRUE))
    }
    attr(ll, "df")<-length(object$coefficients) +1
    attr(ll, "nobs")<-length(object$residuals)
    ll
  }
