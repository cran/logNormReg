\name{lognlm.fit}
\alias{lognlm.fit}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
The fitter function for log Normal Linear Models
}
\description{
\code{lognlm.fit} is called by \code{lognlm} to fit log Normal linear regression models. Two optimizing functions can be used, \code{nlminb} and 
\code{optim}. This function is not meant to be called by the user directly.
}
\usage{
lognlm.fit(X, y, par, lik = TRUE, opt = c("nlminb", "optim"), 
    offset=NULL, weights=NULL,...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{
design matrix for standard linear terms.
}
  \item{y}{
vector of observations of length \code{n}.
}
  \item{par}{
starting values of parameters to be estimated.
}
  \item{lik}{
logical. See \code{lik} in \code{\link{lognlm}}
}
  \item{opt}{
the optimizing algorithm. Default to \code{nlminb}.
}
  \item{offset}{
a possible offset term.
}
  \item{weights}{
apossible weights to be used if \code{lik=FALSE}.
}
  \item{\dots}{
other arguments to be passed to the optimizer specified in \code{opt}.
}
}
\details{
See \code{\link{lognlm}} for more details on the arguments and returned objects.
}
\value{
A list of fit information
}
%\references{
%% ~put references to the literature/web site here ~
%}
\author{
Vito M.R. Muggeo
}
\note{
This function should usually not be used directly by the user.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{nlminb}}, \code{\link{optim}}, \code{\link{lognlm}}
}
\examples{
## See ?lognlm
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ regression }% use one of  RShowDoc("KEYWORDS")
