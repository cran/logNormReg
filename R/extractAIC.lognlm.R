extractAIC.lognlm<-function(fit, scale = 0, k = 2, ...){
  n <- length(fit$residuals)
  edf <- length(fit$coefficients) + 1 #n - fit$df.residual
  if(k<=0) k <-log(n)
  aic <- -2* as.numeric(logLik.lognlm(fit,...)) + k*edf
  c(edf, aic)
}
