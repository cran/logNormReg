\name{palermo}
\alias{palermo}
\docType{data}
\title{
Air quality in Palermo (Italy), 1997-2001
}
\description{
Daily time series of some pollutants and meteorological variables in Palermo, 1997-2001
}
\usage{data("palermo")}
\format{
  A data frame with 1826 observations on the following 8 variables.
  \describe{
    \item{\code{day}}{day of month}
    \item{\code{month}}{month of year}
    \item{\code{year}}{year}
    \item{\code{so2}}{Sulfur dioxide}
    \item{\code{no2}}{Nitrogen dioxide}
    \item{\code{pm10}}{particular matter}
    \item{\code{temp}}{tempearture (degrees Celsius)}
    \item{\code{hum}}{air humidity (\%)}
  }
}
\details{
Data refer to air pollution, temperature and humidity registered in Palermo, (Sicily, Italy) in 1997-2001. Data are averages from eight monitoring stations in the city.
}
%\source{
%%  ~~ reference to a publication or URL from which the data were obtained ~~
%}
%\references{
%%  ~~ possibly secondary sources and usages ~~
%}
\examples{
\dontrun{
data(palermo)

o<-lognlm(pm10 ~ temp + hum, data=palermo)


}
}
\keyword{datasets}
