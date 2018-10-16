\name{summary.lognlm}
\alias{summary.lognlm}
\alias{print.summary.lognlm}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Summarizing  model fits for log Normal regression
}
\description{
summary method for class \code{lognlm}.
}
\usage{
\method{summary}{lognlm}(object, ...)

\method{print}{summary.lognlm}(x, digits = max(3L, getOption("digits") - 3L), 
       signif.stars = getOption("show.signif.stars"), ...)

}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{ object of class \code{"lognreg"}.}
  \item{x}{a \code{summary.segmented} object produced by \code{summary.segmented()}.}
  \item{digits}{controls number of digits printed in output.}
  \item{signif.stars}{logical, should stars be printed on summary tables of coefficients?}
  \item{\dots}{ further arguments to be passed to \code{vcov}, for instance \code{sandw=TRUE}. }
}
\details{
These functions compute and print some useful information relevant to \code{"lognlm"} fits, including point estimates, 
standard errors and p-values.
}
\value{
A list (similar to one returned by \code{lognlm} with additional components, such as the estimate standard errors and corresponding p-values.
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
See also \code{\link{lognlm}} and \code{\link{vcov.lognlm}}
}
\examples{
\dontrun{
n=20
s=.2
set.seed(10)      #just to get reproducible results..

#covariates
x<-seq(.1,10,l=n) 
z<-rnorm(n)

#response
mu<- 10+.5*x- z  #linear regression function
y<-rlnorm(n, log(mu)-s^2/2, s) #data..

o<- lognlm(y~x+z) #the model
summary(o, sandw=TRUE)
}
}
\keyword{ regression }