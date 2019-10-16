extractAIC.lognlm<-function(fit, scale = 0, k = 2, ...){
  n <- length(fit$residuals)
  edf <- length(fit$coefficients) + 1 #n - fit$df.residual
  aic <- -2* as.numeric(logLik.lognlm(fit)) + 2*edf
  c(edf, aic)
}
