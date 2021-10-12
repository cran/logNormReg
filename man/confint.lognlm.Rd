\name{confint.lognlm}
\alias{confint.lognlm}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Confidence intervals for the parameters in log normal regression
}
\description{
Computes confidence intervals (based on the Gradient, Wald or Likelihood Ratio sattistic) for the linear parameters in a fitted `lognreg' model.
}
\usage{
\method{confint}{lognlm}(object, parm, level=0.95, type = c("wald", "gradient", "lrt"), ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{
  a fitted \code{lognlm} object.
}
  \item{parm}{
  the parameter of interest. Numeric (covariate number) or character (covariate name). If missing \code{parm=2} is taken, i.e. the coefficient of the first covariate, provided the intercept is in the model.
}
  \item{level}{
  the required confidence level (default to 0.95).
}
  \item{type}{
Which statistics should be used? Currently \code{"wald"}, \code{"gradient"}, or \code{"lrt"}. Names can be abbreviated. If \code{object} has been obtained with \code{lik=FALSE}, only \code{type="wald"} or \code{"gradient"} is permitted.
}
  \item{\dots}{
When \code{type} is not \code{"wald"}, other optional arguments to be passed on the internal functions: \cr 
- \code{lim} to specify the range of the evaluation points (default to (-3,3) resulting in the interval \eqn{\hat\beta \pm 3\times \mathrm{SE}(\hat\beta)}{est +- 3*SE}); \cr  
- \code{values} to set explicitly the evaluation point(s); \cr 
- \code{return.val} to return (if \code{TRUE}) the evaluation points and the corresponding statistic values (useful to plot the profiled statistic). If the supplied \code{values} includes just one scalar, \code{return.val} is set to \code{TRUE}.
}
}
\details{
Confidence intervals are computed and returned. Currently the Wald, Gradient or Likelihood ratio statistic can be used. Based on some simulation experiments the simple Wald based CIs appears adeguate to guarantee the nominal coverage levels.  
}
\value{
The end-points of confidence intervals.
}
\references{
For a gentle and general introduction about the likelihood-based statistics (including the gradient) see

Muggeo V.M.R., Lovison G. (2014), The 'three plus one' likelihood-based test statistics: unified geometrical and graphical interpretations. \emph{The American Statistician}, \bold{68}, 302-306.
}
\author{
Vito Muggeo
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{lognlm}}
}
\examples{
\donttest{

n=50
s=.4
set.seed(1515)      #just to get reproducible results..

#covariates
x<-seq(.1,10,l=n) 
z<-rnorm(n)

#response
mu<- 10+.5*x- z  #linear regression function
y<-rlnorm(n, log(mu)-s^2/2, s) #data..

o<- lognlm(y~x+z, lik=TRUE) #ML estimation

confint(o, "x", type="g")
confint(o, "z", type="w") #same than confint.default(o) 
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ regression }
\keyword{ nonlinear }
