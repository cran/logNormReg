lognlm.fit <-
function(X, y, par, lik=TRUE, opt=c("nlminb","optim"), offset=NULL, weights=NULL, 
         ...){ #method=c("BFGS", "Nelder-Mead")
#fitter function to estimate multiple linear regression with logNormal errors
#X,y: design matrix and response
#par: optional, starting values for 
    opt<-match.arg(opt)
    if(min(y)<=0) stop("Data can include only positive values")
    n<-length(y)
    p<-ncol(X)
    if(!missing(par)) {
      if(lik && length(par)!=(ncol(X)+1) ) stop("if 'par' is provided, it has to be length(par)= ncol(X)+1")
      if(!lik && length(par)!=ncol(X) ) stop("if 'par' is provided, it has to be length(par)= ncol(X)")
    }
    #if (!is.null(offset)) y <- y - offset
    if (is.null(offset)) offset<-0
    #y<-y-offset #Non funziona perche' poi le y possono essere negative..
    if(is.null(weights)) weights<-1
    if(lik){
    mioNR <- function(start, X, y, h=.5, tol=1e-5,offs=0){
           par0<-start
           eps<-10
           it<-0
           while(eps>tol){
                   par1<- par0 - h* drop(solve(hess(par0,X,y,offs), lgrad(par0,X,y,offs)))
                   eps<- max(abs((par1-par0)/par0))
                   par0<-par1
                   it<-it+1
                   cat(it, par0, "\n")
                   }
           par1
           }

    llik<-function(par,X,y,offs=0,ww=1){
              b <- par[-length(par)]
              s <- par[length(par)]
              mu<- pmax(X %*% b,.0001)+offs
              mz<- log(mu)-s^2/2
              #r<- -sum(dlnorm(y,mz,s,log=TRUE))
      		    li<- -.5*log(s^2) -((log(y)-mz)^2)*ww/(2*s^2) #tolto "-log(y)" e poi dovrebbe essere "log(2*pi*s^2)".. non serve..
              r <- -sum(li)
      		    r
              }
	lgrad<-function(par,X,y,offs=0,ww=1,mu){
	      #NB ww=1 l'ho messo perche' lgrad() sotto lo ha.. (il check lo impone..)
	            b <- par[-length(par)]
	            s <- par[length(par)]
	            if(missing(mu)) mu<- pmax(X %*% b, 0.0001)+offs
	            mz<- log(mu)-s^2/2
	            deriv.b<- crossprod(X, (log(y) - mz)*ww/mu)/s^2 #t(X) %*% ((log(y)-mz)/mu)/s^2   
	            #questa  e' rispetto a s2:
              #deriv.s2<- -n/(s^2) - (s^2*sum((log(y)-mz)*ww)-sum(((log(y)-mz)^2)*ww))/(s^4) 
              #r<-c(-deriv.b, -deriv.s2)
              #rispetto a s (Attenzione su qualche dataset usando la deriv s (e non s2) la stima di s viene negativa!!!)
	            # invece (ott 21) sul dataset simulato in ?lognlm usando la deriv.s2, la stima ML viene negativa!!! 
              deriv.s<- -n/s - (s^2*sum((log(y)-mz))-sum((log(y)-mz)^2))/(s^3) 
              r<-c(-deriv.b, -deriv.s)
              #--
	            #r<-grad(llik, c(b,s), X=X,y=y)
	            r
	            }
 hess<-function(par, X, y, offs=0, ww=1, mu){
 #hessiana della logLik 
	            n<-length(y)
              b <- par[-length(par)]
	            s <- par[length(par)]
	            if(missing(mu)) mu <- pmax(drop(X %*% b), 0.0001) + offs
	            mz<- log(mu)-s^2/2
	            w <- drop((mz - log(y) - 1)/mu^2)
              #i w possono essere negativi..
              #r <- crossprod(X*sqrt(w))/(s^2)
              H.bb<- crossprod(X, w*X)/(s^2) #t(X) %*%diag(w) %*%X/(s^2) #
              H.bs<- - (2/s^3) * drop(crossprod(X, (log(y) - log(mu))/mu))# t(X) %*% ((log(y) - log(mu))/mu))
              H.ss<-n/(s^2)- (3/s^4)*sum((log(y)-log(mu))^2)-n/4
              r<- -cbind(rbind(H.bb, H.bs), c(H.bs, H.ss))
              r
              }
    if(missing(par)) {
                #b<-solve(crossprod(X),crossprod(X,log(y)))
                b<-solve(crossprod(X),crossprod(X,y))
                s<- log(sqrt(sum((y-drop(X%*%b))^2)/n))
                #s<- sqrt(sum((log(y)-drop(X%*%b))^2)/n)
                  } else {
                b<-par[1:p]
                s<-par[p+1]
                }
    opz<-list(...)
#browser()
    if(opt=="optim"){
      opz$par<-c(b,s)
      opz$fn<- llik
      opz$gr<-lgrad
      opz$offs<- quote(offset)
      opz$X<-quote(X)
      opz$y<-quote(y)
      opz$ww <-quote(weights)
      opz$hessian<-TRUE
      if(is.null(opz$method)) opz$method<-"BFGS"
      o<-do.call(optim, opz )
      ll<-o$value
      #hh<- o$hessian #hess(o$par,X,y) #o$hessian e' calcolata numericamente..
      } else {
      # browser()
      # o<-mioNR(c(b,s),X,y)
      opz$start<- c(b, s)
      opz$objective<- quote(llik)
      opz$gradient<- quote(lgrad)
      opz$hessian<- quote(hess)
      opz$offs<- quote(offset)
      opz$X<-quote(X)
      opz$y<-quote(y)
      opz$ww <-quote(weights)
      o<-do.call(nlminb, opz)
      ll<-o$objective
      #hh<- hess(o$par,X,y)
      }
    
	  est.b<-o$par[1:p]
    est.s2<-o$par[(p+1)]^2
    est.s2<-est.s2*n/(n-p)
    fits<-pmax(drop(X%*%est.b)+offset,.0001)
    gradd<-lgrad(o$par, X=X, y=y, offs=offset, ww=weights, mu=fits)
    hh<-hess(o$par,X,y,offs=offset,ww=weights,mu=fits)
    Ehh<-NULL
    #browser()
    
    #end link=TRUE
    } else {
######================================ Min distance
    dmin <- function(par,X,y,offs=0, ww=1){
        b <- par
        mu<- pmax(drop(X %*% b),.0001)+offs
        r<- sum(((log(y)-log(mu))^2)*ww)
        r
    }
    lgrad <- function(par,X,y,offs=0, ww=1, mu){
            b <- par                      
            if(missing(mu)) mu<- pmax(drop(X %*% b),.0001) +offs
            #deriv.b<- -2*t(X)%*%(((log(y)-log(mu))*ww)/mu)
            deriv.b<- -crossprod(X, ww*(log(y)-log(mu))/mu)
            r<-deriv.b
            r
    }
    hess <- function(par, X, y, offs=0, ww=1, mu){
      b<-par
      if(missing(mu)) mu<- pmax(drop(X %*% b), 0.0001) + offs
      w <- drop(ww*(log(y) - log(mu)+ 1)/mu^2) #possono essere negativi..
      #r <- crossprod(X*sqrt(w))/(s^2)
      H.bb<- crossprod(X, w*X) #t(X) %*%diag(w) %*%X/(s^2) #
      H.bb
    }
    E.hess<- function(par, X, y, offs=0, ww=1, mu){ #E[H]
      b<-par
      if(missing(mu)) mu<- pmax(drop(X %*% b), 0.0001) + offs
      w <- drop(ww/mu^2) #possono essere negativi..
      H.bb<- crossprod(X, w*X) #t(X) %*%diag(w) %*%X/(s^2) #
      H.bb
    }
    if(missing(par)) {
        #b<-solve(crossprod(X),crossprod(X,log(y)))
        b<-solve(crossprod(X),crossprod(X, y))
    } else {
        b<-par[1:p]
    }
    opz<-list(...)
    if(opt=="optim"){
        
        opz$par<-b
        opz$fn<- dmin
        opz$gr<-lgrad
        opz$offs<- quote(offset)
        opz$ww<- quote(weights)
        opz$X<-X
        opz$y<-y
        opz$hessian<-TRUE
        if(is.null(opz$method)) opz$method<-"BFGS"
        o<-do.call(optim, opz )
        ll<- -o$value
    } else {
        opz$start<- b
        opz$objective<- quote(dmin)
        opz$gradient<- quote(lgrad)
        opz$hessian<- quote(hess)
        opz$offs<- quote(offset)
        opz$X<-quote(X)
        opz$y<-quote(y)
        opz$ww<- quote(weights)
        o<-do.call(nlminb, opz)
        ll<- -o$objective
    }
    est.b<-o$par[1:p]
    fits<-pmax(drop(X%*%est.b)+offset,.0001)
    
    gradd<-lgrad(o$par, X=X, y=y, offs=offset, ww=weights, mu=fits)
    est.s2<- sum(weights*(log(y) - log(fits))^2)/(n-p)
    hh <- hess(o$par, X=X, y=y, offs=offset, ww=weights, mu=fits)  #o$hessian e' numerica!
    Ehh<- E.hess(o$par, X=X, y=y, offs=offset, ww=weights, mu=fits)
    } #end min distance

    if(o$convergence!=0) warning("Unsuccessful convergence") 
    
    est<-c(est.b, est.s2)
    
    r<-list(coefficients=est.b, loglik= -ll, s2=est.s2, fitted.values=fits, residuals=y-fits, 
             grad=drop(gradd) , hessian=hh, Ehessian=Ehh, convergence=o$convergence)
    return(r)
}
