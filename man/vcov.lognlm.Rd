\name{vcov.lognlm}
\alias{vcov.lognlm}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Covariance matrix for lognlm fits
}
\description{
Computes covariance matrix of parameter estimates from a \code{lognlm} fit.
}
\usage{
\method{vcov}{lognlm}(object, sandw = FALSE, emp = TRUE, se = FALSE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{
a fitted model object of class \code{"lognlm"} returned by \code{lognlm()}
}
  \item{sandw}{
logical, if \code{TRUE} the returned covariance matrix comes from the sandwich formula, otherwise it is just the (minus) hessian. Meaningless if \code{object} has been obtained with \code{lik=FALSE}.
}
  \item{emp}{
logical; if \code{TRUE} (and \code{sandw=TRUE} too), the `meat' (i.e the information matrix) is computed empirically by the outer product of the individual score contributions. 
}
  \item{se}{
logical; if \code{TRUE} the square root of the elements of the main diagonal are returned (rather than the whole matrix).
}
  \item{\dots}{
additional arguments.
}
}
\details{
If \code{object} has been obtained via \code{lognlm(.., lik=TRUE)} the returned covariance matrix (or standard errors only) refers to regression coefficients \emph{and} the response standard deviation. Otherwise (if \code{lik=FALSE} has been set), it includes entries relevant to regression coefficients. 
}
\value{
The variance-covariance matrix of the parameter estimates, if \code{se=FALSE}; otherwise the square root of the main diagonal entries.
}
%\references{
%% ~put references to the literature/web site here ~
%}
\author{
Vito Muggeo
}
%\note{
%%  ~~further notes~~
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{lognlm}}
}
\examples{
\donttest{
n=50
s=.3
set.seed(1515)      #just to get reproducible results..

#covariates
x<-seq(.1,10,l=n) 
z<-rnorm(n)

#response
mu<- 10+.5*x- z  #linear regression function
y<-rlnorm(n, log(mu)-s^2/2, s) #data..

o<- lognlm(y~x+z) #the model
vcov(o) #the full covariance matrix 
vcov(o, se=TRUE, sandw=TRUE) #st.errs only, based on the sandwich formula 

o<- lognlm(y~x+z, lik=FALSE) #min. distance estimation
vcov(o, se=TRUE) #'sandw=TRUE' is not allowed if the fit has been obtained with 'lik=FALSE'

 }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{regression}