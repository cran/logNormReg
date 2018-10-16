\name{lognlm}
\alias{lognlm}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Multiple linear regression with log Normal errors
}
\description{
The function fits simple multiple linear regression models with log Normals erros. Two objectives as well as two optimizing functions can be used. 
}
\usage{
lognlm(formula, data, subset, weights, na.action, y = TRUE, start, model = TRUE, 
     lik = TRUE, opt = c("nlminb", "optim"), ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{formula}{
a standard R formula with response and explanatory variables (and possible offset) specifying the regression model being fitted.
}
  \item{data}{
an optional data frame, list or environment containing some or all the variables in the model.
}
  \item{subset}{
an optional vector specifying a subset of observations to be used in the fitting process.
}
  \item{weights}{
an optional vector of weights to be used in the fitting process. Currently unimplemented.
}
  \item{na.action}{
an optional vector of weights to be used in the fitting process. 
}
  \item{y}{
logical. If \code{TRUE} the response vector is returned as \code{y} in the fit object.
}
  \item{start}{
(optional) starting values of the parameter to be estimated. If \code{start} is missing they are computed via ordinary least squares with the intercept \eqn{\hat\beta_0}{beta0} replaced by \eqn{max(\hat\beta_0,median\{y_i\})}{max(beta0,median(yi))}. \cr 
If \code{lik=TRUE} (i.e. a log Normal model is fitted), \code{start} refers to the regression parameters \emph{and} the error standard deviation; if \code{lik=FALSE}, \code{start} does \emph{not} include the starting guess for the standard deviation.
}
  \item{model}{
logical. If \code{TRUE} the model frame is returned as \code{model} in the fit object.
}
  \item{lik}{
If \code{TRUE} the log Normal log likelihood is optimized, otherwise the sum of squared residuals based on the logs (see Details).
}
  \item{opt}{
the optimization function to be used. \code{nlminb} has shown to be more efficient, probably because it uses (unlike \code{optim}) also the hessian matrix (supplied in the code).
}
  \item{\dots}{
optional arguments passed on to the optimizing functions (\code{nlminb} or \code{optim}), (and therefore should be \emph{consistent} with that).
}
}
\details{
\code{lognlm} fits simple linear regression models with log Normal errors and identity link. Actually two objectives could be used.

If \code{lik=TRUE} the usual log Normal likelihood is optimized, otherwise estimation is based on minimization of the following loss function
\deqn{\sum_i(\log y_i-\log \mu_i )^2}{ sum(log yi- log mui)^2}
where \eqn{\mu_i=x_i^T\beta}{mui=xi' beta} is the mean function equal to the linear predictor (as an identity link is exploited).
}
\value{
A list with components
  \item{coefficients}{the regression parameters estimate.}
  \item{loglik}{The objective function value, namely the log Normal log likelihood or the sum of the squared `log residuals' (depending on \code{lik} option).}
  \item{s2}{the error variance estimate.}
  \item{fitted.values}{the fitted values.}
  \item{residuals}{the raw residuals, i.e. 'observed - fitted'.}
  \item{grad}{the gradient at solution.}
  \item{hessian}{the hessian matrix at solution.}
  \item{convergence}{the convergence code coming from the fitter function. \code{0} means succefull convergence. }
  \item{call}{the matched call.}
  \item{y}{the response vector (provided that \code{y=TRUE} has been set).}
  \item{opt}{the employed optimizer.}
  \item{lik}{logical, indicating if the fit comes from a log Normal likelihood approach.}
  \item{xlevels}{(only where relevant) a record of the levels of the factors used in fitting.}
	\item{terms}{the terms object used.}
  \item{contrasts}{(only where relevant) the contrasts used.}
  \item{model}{if requested, i.e. \code{model=TRUE} has been set (the default), the model frame used.}
  \item{offset}{the (possible) offset used.}
}
%\references{
%}
\author{
Vito M.R. Muggeo
}
%\note{
%%  ~~further notes~~
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
See also \code{\link{print.lognlm}} and \code{\link{summary.lognlm}} to display results.
}
\examples{

n=200
s=.4
set.seed(1234)      #just to get reproducible results..

x<-seq(.1,10,l=n) #covariate
mu<- 10+2*x  #linear regression function
y<-rlnorm(n, log(mu)-s^2/2, s) #data..

o0<-lm(log(y)~x) #the usual but WRONG model
o<- lognlm(y~x) #fit the 'right' model

plot(x,y)
lines(x, mu, lwd=2)
points(x, exp(fitted(o0)), col=2, type="l", lwd=2)
points(x, fitted(o), col=3, type="l", lwd=2)
legend("topleft", legend=c("true", "lm(log(y)~x)", "lognlm(y~x)"), 
    col=c(1,2,3), lwd=2)

#Sometimes people would estimate parameters by minimizing a least square objective 
# (i.e. by setting 'lik=FALSE', see Details), wherein data would come from 
# Y = mu * exp(eps) where eps~N(0,s)..
y1<-mu*exp(rnorm(n,0,1)) #data..
o1<-lognlm(y1~x, lik=FALSE) #set 'lik=FALSE', see Details


}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ regression }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
