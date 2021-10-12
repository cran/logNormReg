confint.lognlm<-function(object, parm, level = 0.95, type=c("wald", "gradient", "lrt"), ...){
  if(missing(parm)) parm<-2
  if(length(parm)>1) stop("'parm' should not be a vector") 
  if(is.character(parm)) {
    if(! parm%in%names(object$coefficients)) stop(" 'parm' does not appear to be in the model")
  } else { #if numeric..
    if(parm > length(object$coefficients)) stop("parm is larger than the number of linear coefficients")
    parm<- names(object$coefficients)[parm]
  }
  type<-match.arg(type)
  if(type=="lrt" && !object$lik) stop("type='lrt' is allowed only for likelihood-based fits")
  lgrad.ML<-function(par,X,y,offs=0,ww=1,mu){
	            n<-length(y)
              b <- par[-length(par)]
	            s <- par[length(par)]
	            if(missing(mu)) mu<- pmax(X %*% b, 0.0001)+offs
	            mz<- log(mu)-s^2/2
	            #deriv.b<- t(X) %*% ((log(y)-mz)/mu)/s^2
	            deriv.b<- crossprod(X, (log(y) - mz)*ww/mu)/s^2 #t(X) %*% ((log(y)-mz)/mu)/s^2 
	            #questa  e' rispetto a s2:
              #deriv.s<- -n/(2*s^2) - (s^2*sum((log(y)-mz))-sum((log(y)-mz)^2))/(2*s^4) 
              #rispetto a s (Attenzione su qualche dataset usando la deriv s (e non s2) la stima di s viene negativa!!!)
              deriv.s<- -n/s - (s^2*sum((log(y)-mz))-sum((log(y)-mz)^2))/(s^3) 
	            #--
	            r<-c(-deriv.b, -deriv.s)
	            r
  }
  lgrad.MD <- function(par,X,y,offs=0, ww=1, mu){
    b <- par                      
    if(missing(mu)) mu<- pmax(drop(X %*% b),.0001) +offs
    #deriv.b<- -2*t(X)%*%(((log(y)-log(mu))*ww)/mu)
    deriv.b<- -crossprod(X, ww*(log(y)-log(mu))/mu)
    r<-deriv.b
    r
  }
  lgrad<-if(object$lik) lgrad.ML else lgrad.MD
  
################################################################################
  ci.grad<-function(obj, term, conf.level=0.95, stat=c("gradient","lrt"), 
                    lim=c(-3,3), values=NULL, return.val=FALSE){
    #term: numeric or character	
    est<-coef(obj)
    p<-length(est)
    nomi<-names(est)
    if(is.character(term)) {
            if(! term%in%names(est)) stop(" 'term' does not appear to be in the model")
            term<-match(term, nomi)
    }
    X<-model.matrix(obj)
    y<-obj$y
    offset<-obj$offset
    if(is.null(offset)) offset<-0
    weights<-obj$weights
    if(is.null(weights)) weights<-1
    if(is.character(term)) term<-match(term, nomi)
    id.others<- setdiff(1:p,term) #id for other covariates
    b<-est[term]
    x<-X[,term]
    name.X <- names(est)[term]
    name.others <- names(est)[-term]
    newX<-X[,-term,drop=FALSE]
    if(is.null(values)){
            se<-vcov(obj, se=TRUE)[term]
      	    if(length(lim)<=1) lim<-c(-abs(lim), abs(lim))
      	    valori<- seq(b+min(lim)*se, b+max(lim)*se, l=20)
            } else {
            valori<-values
    }
    if(length(valori)<=1) return.val<-TRUE
    if(length(valori)<10 && !return.val) warning("profiled based on few values.. probably not accurate")
    r<-rep(NA, length(valori))
    for(i in 1:length(valori)){
      		b0<-valori[i]
      		off0<-x*b0 + offset
      		start0<- if(obj$lik) c(est[-term], sqrt(obj$s2)) else est[-term] 
      		obj0<-lognlm.fit(newX, y, par=start0, lik=obj$lik, opt=obj$opt, offset=off0, weights=weights)
      		if(stat=="gradient"){
               est[term]<-b0 #est[name.X]<-b0
      		     est[id.others]<-obj0$coefficients #est[name.others]<-obj0$coefficients[name.others]
      		     estGrad<-if(obj$lik) c(est, sqrt(obj0$s2)) else est
      		     r[i]<-abs((b-b0)*lgrad(estGrad, X, y, offs=offset, ww=weights)[term])
      	       } else {
      	         r[i]<- 2*(obj$loglik-obj0$loglik)
      	       }
        }
      	r<-sqrt(r)*sign(b-valori)
        if(return.val){
           ris<-cbind(value=valori,stat=r)
           return(ris)
           }
      	alpha<-1-conf.level
      	z1 <- qnorm(alpha/2)
      	z2 <- qnorm(conf.level + alpha/2)
      	if(min(r)>z1) warning("too narrow range in the right endpoint")
      	if(max(r)<z2) warning("too narrow range in the left endpoint")
      	f<-splinefun(r,valori)
        ci<-f(c(z2, z1))
        #
      	ci
      }
################################################################################
    ci.wald<-function(obj, term, conf.level=0.95,...) {
            est<-coef(obj)
            nomi<-names(est)
            if(is.character(term)) {
                if(! term%in%names(est)) stop(" 'term' does not appear to be in the model")
                term<-match(term, nomi)
                }
            alpha<-1-conf.level
            z1 <- qnorm(alpha/2)
            z2 <- qnorm(conf.level + alpha/2)
            
            b<-est[term]
            se<-vcov(obj, se=TRUE,...)[term]
            r<-c(b+z1*se, b+z2*se)
            names(r)<-NULL
            r
            }
################################################################################
    r<-switch(type,
         wald = ci.wald(object, parm, level,...),
         gradient = ci.grad(object, parm, level, stat="gradient", ...),
         lrt = ci.grad(object, parm, level, stat="lrt", ...)
         )
    r<-matrix(r, ncol=2)
    colnames(r)<- paste("CI","(",level*100,"%",")",c(".low",".upp"),sep="")
    rownames(r)<- parm
    r
   }	
	