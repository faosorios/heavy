\name{life}
\docType{data}
\alias{life}

\title{Data on per-capita income and life expectancy}
\description{
  The data represent life expectancies and per capita income in 1974 dollars for
  101 countries in the early 1970's. The data were originally published in the
  New York Times (September, 28, 1975, p. E-3), they also appear on Leinhardt and
  Wasserman (1979) and Regression Analysis by Ashish Sen and Muni Srivastava.
}
\usage{data(life)}
\format{
  A data frame with 101 observations on the following 3 variables.
  \describe{
    \item{Country}{ character string, containing names of the countries.}
    \item{life}{ life expectancy, in years for 1974.}
    \item{income}{ income per capita in 1974 US dollars.}
  }
}
\source{
  Sen, A., and Srivastava, M. (1990).
  \emph{Regression Analysis: Theory, Methods and Applications}.
  Springer-Verlag, New York.

  Leinhardt, S., and Vasserman, S.S. (1979).
  Teaching regression: An exploratory approach.
  \emph{The American Statistician} \bold{33}, 196-203.
}
\keyword{datasets}
