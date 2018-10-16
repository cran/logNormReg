lognlm <-
function(formula, data, subset, weights, na.action, y = TRUE, 
          start, model=TRUE, lik=TRUE, opt=c("nlminb","optim"), ...) {#method=c("BFGS", "Nelder-Mead")
#...: additional arguments to be passed to optim or nlmnib depending on 'opt', most often a control argument.
    opt<-match.arg(opt)
    call <- match.call()
    if (missing(data)) data <- environment(formula)
    mf <- match.call(expand.dots = FALSE)
    m <- match(c("formula", "data", "subset", "weights", "na.action", "offset"), names(mf), 0L)
    mf <- mf[c(1, m)]
    mf$drop.unused.levels <- TRUE
    mf[[1L]] <- as.name("model.frame")
    mf <- eval(mf, parent.frame())
    mt <- attr(mf, "terms")
    intercMt <- attr(mt, "intercept")
    Y <- model.response(mf, "any")
    if (length(dim(Y)) == 1L) {
        nm <- rownames(Y)
        dim(Y) <- NULL
        if (!is.null(nm)) 
            names(Y) <- nm
    }
    X <- if (!is.empty.model(mt)) model.matrix(mt, mf, contrasts) else stop("error in the design matrix")
    p<-ncol(X)
    if(!missing(start)) {
              if(lik && (length(start)!= (p+1))) stop("if 'lik=TRUE', length(start) should have ncol(X)+1 values ")
              if(!lik && (length(start)!= (p))) stop("if 'lik=FALSE',  length(start) should have ncol(X) values  ")          
                        } 
    attrContr <- attr(X, "contrasts")
    weights <- as.vector(model.weights(mf))
    offset <- as.vector(model.offset(mf))
    if (!is.null(offset)) {
        if (length(offset) != NROW(Y)) 
            stop(gettextf("number of offsets is %d, should equal %d (number of observations)", 
                length(offset), NROW(Y)), domain = NA)
    }

    n <- nrow(X)
    if (!is.null(weights) && !is.numeric(weights)) 
        stop("'weights' must be a numeric vector")
    if (!is.null(weights) && any(weights < 0)) 
        stop("negative weights not allowed")
    if(missing(start)) {
                #b<-drop(solve(crossprod(X),crossprod(X,log(Y+.001))))
                b<-drop(solve(crossprod(X),crossprod(X, Y)))
                if(colnames(X)[1]=="(Intercept)") b[1]<-max(b[1],median(Y)) #oppure mean(Y)
                #s<-sqrt(sum((log(Y+.001)- log(drop(X%*%b)))^2)/n)
                #s<-log(sqrt(sum((y-drop(X%*%b))^2)/n))
                s<-log(mad(Y))
                
            } else {
                b<-start[1:p]
                s<-start[p+1]
            }
    par0 <- if(lik) c(b,s) else b
    obj<-lognlm.fit(X=X,y=Y, par=par0, lik=lik, opt=opt, offset=offset, ...)
    names(obj$coefficients)<-colnames(X)
    names(obj$s2)<-""
    obj$call<-call
    if(y) obj$y<-Y
    class(obj)<-"lognlm"
    obj$na.action <- attr(mf, "na.action")
    obj$offset <- offset
    obj$contrasts <- attr(X, "contrasts")
    obj$xlevels <- .getXlevels(mt, mf)
    obj$terms <- mt
    obj$opt<-opt
    obj$lik<-lik
    if (model) obj$model <- mf
#    fit$formula <- formula
    obj
    }
