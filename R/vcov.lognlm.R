vcov.lognlm <-
function(object, emp=FALSE, exH=TRUE, se=FALSE, ...) {
    if(object$lik) {
      exH<-FALSE
      emp<-TRUE
    }
    #browser()
    invH<- if(exH) solve(object$Ehessian) else solve(object$hessian)
    cof<-object$coefficients
    y<-object$y
    n<-length(y)
    X<-model.matrix(object, model.frame(object))
    s2<-object$s2
    s<-sqrt(s2)
    mu<-object$fitted.values
    w<-if(is.null(object$weights)) 1 else object$weights
    if(any(w<=0) && !emp) stop("set 'emp=TRUE' with null weights")
    if(!object$lik){#se min dist
        if(emp){
              U<- X * ((log(y)-log(mu))*w/mu)
              InF<- crossprod(U)
        } else {
              s2<-s2/w
              InF<- crossprod(X, s2*((w/mu)^2)*X)
        }
    } else { #se ML
        if(emp){
              mz<- log(mu)-s2/2
              U<- X * ((log(y)-mz)*w/mu)/s2
              #Rispetto a s2
              #U.s<- -1/(s2) - (s2*((log(y)-mz))-((log(y)-mz)^2))/(s2^2) #
              #Rispetto a s:
              #U.s<- -1/s - w*(s^2*(log(y)-mz)- (log(y)-mz)^2)/(s^3) #oppure (e' lo stesso):
              U.s<- -1/s -s*w/4+ (w*(log(y)-log(mu))^2)/(s^3) 
              U<-cbind(U, U.s)
              InF<- crossprod(U)
      } else {
              InF<- (1/s2)*crossprod(X/mu) #questa esclude la componente di s!
      }
    }
    V<- invH%*%InF%*%invH
    nomi<- if(ncol(V)> length(cof)) c(names(cof), "s") else names(cof)
    colnames(V)<-rownames(V)<-nomi
    if(se) V<-sqrt(diag(V))
    V
}
