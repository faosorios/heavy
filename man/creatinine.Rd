\name{creatinine}
\docType{data}
\alias{creatinine}

\title{Clinical trial on endogenous creatinine clearance}
\description{
  Data from a clinical trial on endogenous creatinine clearance of 34 male patients.
  Of the 34 male patients, two had no recorded the body weight, and four were missing
  the serum creatinine.
}
\usage{data(creatinine)}
\format{
  A data frame with 34 observations on the following 4 variables.
  \describe{
    \item{WT}{ body weight in kilograms.}
    \item{SC}{ serum creatinine in mg/deciliter.}
    \item{Age}{ age of the patients in years.}
    \item{CR}{ endogenous creatinine clearance.}
  }
}
\source{
  Liu, C., and Rubin, D.B. (1995).
  ML estimation of the t distribution using EM and its extensions, ECM and ECME.
  \emph{Statistica Sinica} \bold{5}, 19-39.  
}
\keyword{datasets}