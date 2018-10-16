\name{logLik.lognlm}
\alias{logLik.lognlm}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Log Likelihood for log Normal linear regression
}
\description{
The function returns the log-likelihood value of the log Normal linear regression model evaluated at the estimated coefficients
}
\usage{
\method{logLik}{lognlm}(object, full=FALSE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{
A \code{lognlm} fit returned by \code{lognlm()}
}

  \item{full}{
If \code{FALSE}, only the \emph{kernel} of the log likelihood is returned, otherwise the complete log likelihood (with terms depending on data only)
}
  \item{\dots}{
optional  arguments (nothing in this method).
}
}
\details{
If \code{object} has been obtained via \code{lognlm(.., lik=TRUE)}, \code{logLik.lognlm} returns the log likelihood (kernel or complete, depending on argument \code{full}), otherwise the sum of log residuals, \eqn{\sum((\log(y_i)-\log(\hat\mu_i))^2)}{sum((log(y)-log(mu))^2)} .
}
\value{
The log likelihood (or the sum of log residuals squared) of the model fit \code{object}
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
# o is the fit object, see ?lognlm
n=50
s=.4

#covariates
x<-seq(.1,10,l=n) 

#response
set.seed(1234)      #just to get reproducible results..
mu<- 10+.5*x  #linear regression function
y<-rlnorm(n, log(mu)-s^2/2, s) #data..

o<- lognlm(y~x) #the model

logLik(o) #the kernel log likelihood value
logLik(o, full=TRUE)
   }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{regression} % use one of  RShowDoc("KEYWORDS")
