\name{vcov.lognlm}
\alias{vcov.lognlm}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Covariance matrix for lognlm fits
}
\description{
Computes covariance matrix of parameter estimates from a \code{lognlm} fit via the sandwich formula.
}
\usage{
\method{vcov}{lognlm}(object, emp = FALSE, exH = TRUE, se = FALSE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{
a fitted model object of class \code{"lognlm"} returned by \code{lognlm()}
}
  \item{emp}{
logical; if \code{TRUE}, the `meat' (i.e the information matrix) is computed empirically by the outer product of the individual score contributions. 
}
  \item{exH}{logical; if \code{TRUE} the \emph{expected} (rather than the observed) hessian is used in the sandwich formula.
}
  \item{se}{
logical; if \code{TRUE} the square root of the elements of the main diagonal are returned (rather than the whole matrix).
}
  \item{\dots}{
additional arguments.
}
}
\details{
If \code{object} has been obtained via \code{lognlm(.., lik=TRUE)} the returned covariance matrix (or standard errors only) refers to regression coefficients \emph{and} the log response standard deviation. Otherwise (if \code{lik=FALSE} has been set), it includes entries relevant to regression coefficients only. The var-covariance matrix comes from the sandwich formula using \emph{expected} (if \code{exH=TRUE}) or the observed (if \code{exH=FALSE}) hessian at solution. Some simulations under correct model specification show that \code{emp=TRUE} and \code{exH=FALSE} lead to somewhat more unstable standard errors.
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
\note{
Currently for likelihood-based fits, \code{exH=FALSE} and \code{emp=TRUE} are always set.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{lognlm}}
}
\examples{
\donttest{
n=50
s=.3

#covariates
x<-seq(.1,10,l=n) 
z<-rnorm(n)
#response
mu<- 10+.5*x- z  #linear regression function
y<-rlnorm(n, log(mu)-s^2/2, s) #data..

o<- lognlm(y~x+z, lik=TRUE) #the model
vcov(o) #the full covariance matrix 
vcov(o, se=TRUE) #st.errs only 

 }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{regression}