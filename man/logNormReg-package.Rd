\name{logNormReg-package}
\alias{logNormReg-package}
\alias{logNormReg}
\docType{package}
\title{
\packageTitle{logNormReg} 
%Multiple linear regression with log Normal errors
}
\description{
%%\packageDescription{logNormReg} performs simple multiple linear regression with log Normal errors and identity link.
 Functions to perform and to display results from simple multiple linear regression with log Normal errors and identity link. Standard errors of estimates are based on the hessian or alternatively on the sandwich formula.
}
\details{
%The DESCRIPTION file:
%\packageDESCRIPTION{logNormReg}
%\packageIndices{logNormReg}
%Package \code{logNormReg} fits simple regression models with log Normal errors
\tabular{ll}{
Package: \tab logNormReg\cr
Type: \tab Package\cr
Version: \tab 0.2-0\cr
Date: \tab 2018-09-25\cr
License: \tab GPL\cr
}

}
\author{
Vito M.R. Muggeo <vito.muggeo@unipa.it>
%\packageAuthor{logNormReg}
%Maintainer: \packageMaintainer{logNormReg}
}
\references{
Gustavsson, S., Fagerberg, B., Sallsten, G.,  Andersson, E. M. (2014). Regression Models for Log-Normal Data: Comparing Different Methods for Quantifying the Association between Abdominal Adiposity and Biomarkers of Inflammation and Insulin Resistance . \emph{International Journal of Environmental Research and Public Health}, \bold{11}, 3521--3539. 

Muggeo, V.M.R. (2018) A note on regression with log Normal errors: linear and piecewise linear modelling in R, doi: 10.13140/RG.2.2.18118.16965.

}
\keyword{ package }
\seealso{
See \code{\link{lognlm}} for the main function with a toy example. \cr  
See also package \code{\link[gamlss]{gamlss}} for more general regression models including log Normal errors.
}
\section{Acknowledgement}{This package was inspired by a fruitful discussion with 
Andrew Beet (Marine Policy Center, Woods Hole Oceanographic Institution, U.S.).}