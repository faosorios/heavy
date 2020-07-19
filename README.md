# Robust estimation using heavy-tailed distributions

[![CRAN status](http://www.r-pkg.org/badges/version/heavy)](https://cran.r-project.org/package=heavy)
[![CRAN RStudio mirror downloads](http://cranlogs.r-pkg.org/badges/heavy)](https://cran.r-project.org/package=heavy)

The HEAVY package contains routines to perform robust estimation considering heavy-tailed distributions. Currently, the package includes linear regression, linear mixed-effect models, Grubbs' model, multivariate location and scatter estimation, multivariate regression, penalized splines, random variate generation and some support functions.

## Features

* Provide basic functionality for modeling using scale mixtures of normal distributions in R, via a package.
* Calculations associated with parameter estimation are performed by calling routines in C and Fortran.
* Estimation in linear regression, linear mixed effects models, multivariate regression and penalized splines using the EM algorithm.
* Estimation of location and Scatter using multivariate heavy-tailed distributions.
* Implemented families: normal, Cauchy, Student-t, slash and contaminated normal.
* Estimation of the shape parameters for Student-t and slash distributions.
* Multivariate random number generation for the implemented families and the uniform distribution on the p-dimensional sphere.
* Print and summary methods and some sample databases.

## Providing Feedback

Please report any bugs/suggestions/improvements to Felipe Osorio, Universidad Tecnica Federico Santa Maria. If you find these routines useful or not then please let me know. Also, acknowledgement of the use of the routines is appreciated.

Author: Felipe Osorio.

Project webpage: http://heavy.mat.utfsm.cl

Lastest binaries and sources for heavy are availables from CRAN at: https://cran.r-project.org/package=heavy
