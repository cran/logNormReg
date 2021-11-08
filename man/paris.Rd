\name{paris}
\alias{paris}
\docType{data}
\title{
PM2.5 and PM10 measurements in Paris in 2019
}
\description{
Hourly time series of PM2.5 and PM10 measurements in Paris in 2019
}
\usage{data("paris")}
\format{
  A data frame with 647 observations on the following 4 variables.
  \describe{
    \item{\code{utc}}{date and time of measurements}
    \item{\code{pm25}}{The PM2.5 measurements}
    \item{\code{pm10}}{The PM10 measurements}
    \item{\code{hours}}{numeric, the measurement hours}
  }
}
\details{
  Ambient particulate matter measurements (PM 2.5 and PM 10) measured by reference-grade instruments for Paris - Centre, 
  hourly measurements for July 2019. Non-physical measurements (zero values and outliers with PM2.5 greater than 2xPM10) were     removed, as mass of PM2.5 is a fraction of mass of PM10 by definition.
}
\source{
European Environmental Agency, through Open AQ \code{(https://openaq.org/#/location/4146)}. Thanks to Vito Ilacqua (EPA) for pointing out that.
}
%\references{
%%  ~~ possibly secondary sources and usages ~~
%}
\examples{

data(paris)
plot(pm10~pm25, data=paris)
o<-lm(pm10~pm25, data=paris) #negative intercept! it's meaningless..
o1<-lognlm(pm10 ~ pm25, data=paris, lik=TRUE)
abline(coef=coef(o),  col=2)
abline(coef=coef(o1), col=3)
}
\keyword{datasets}
