vcov.lognlm <-
function(object, sandw=FALSE, emp=TRUE, se=FALSE, ...) {
            invH<- solve(object$hessian)
            cof<-object$coefficients
            if(sandw){
                      if(!object$lik) stop("'sandw=TRUE' is not allowed for non likelihood-based fits")
                      y<-object$y
                      n<-length(y)
                      X<-model.matrix(object)
                      s2<-object$s2
                      s<-sqrt(s2)
                      mu<-object$fitted.values
                      if(emp){
                              mz<- log(mu)-s2/2
                              U<- X * ((log(y)-mz)/mu)/s2
                              #U.s<- -1/(2*s2) - (s2*((log(y)-mz))-((log(y)-mz)^2))/(2*s2^2) #questa e' ripsetto a s2 (e non s)
                              #Rispetto a s:
                              #U.s<- -1/s - (s^2*(log(y)-mz)- (log(y)-mz)^2)/(s^3) #oppure (e' lo stesso):
                              U.s<- -1/s -s/4+ ((log(y)-log(mu))^2)/(s^3) 
                              U<-cbind(U, U.s)
                              INF<-crossprod(U)
                              } else {
                              INF<- (1/s2)*crossprod(X/mu) #questa esclude la componente di s!
                              }
               V<-invH %*%INF %*% invH
                            } else {
               V<-invH
               }
               nomi<- if(ncol(V)> length(cof)) c(names(cof), "s") else names(cof)
               colnames(V)<-rownames(V)<-nomi
               if(se) V<-sqrt(diag(V))
               V
               }
