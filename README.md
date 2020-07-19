# Robust estimation using heavy-tailed distributions

[![CRAN status](http://www.r-pkg.org/badges/version/heavy)](https://cran.r-project.org/package=heavy)
[![CRAN RStudio mirror downloads](http://cranlogs.r-pkg.org/badges/heavy)](https://cran.r-project.org/package=heavy)

The HEAVY package contains routines to perform robust estimation considering heavy-tailed distributions. Currently, the package includes linear regression, linear mixed-effect models, Grubbs' model, multivariate location and scatter estimation, multivariate regression, penalized splines, random variate generation and some support functions.

## Features

* Provide basic functionality for modeling using scale mixtures of normal distributions in R, via a package.
* Calculations associated with parameter estimation are performed by calling routines in C and Fortran.
* Estimation in linear regression, linear mixed effects models, Grubbs' model, multivariate regression and penalized splines using the EM algorithm.
* Estimation of location and Scatter using multivariate heavy-tailed distributions.
* Implemented families: normal, Cauchy, Student-t, slash and contaminated normal.
* Estimation of the shape parameters for Student-t and slash distributions.
* Multivariate random number generation for the implemented families and the uniform distribution on the p-dimensional sphere.
* Print and summary methods and some sample databases.

## Providing Feedback

Please report any bugs/suggestions/improvements to [Felipe Osorio](mailto:felipe.osorios@usm.cl), [Universidad Tecnica Federico Santa Maria](http://www.usm.cl). If you find these routines useful or not then please let me know. Also, acknowledgement of the use of the routines is appreciated.

## Resources

Latest binaries and sources for HEAVY are availables from [CRAN package repository](https://cran.r-project.org/package=heavy)

* [heavy_0.38.196.tar.gz](https://cran.r-project.org/src/contrib/heavy_0.38.196.tar.gz) - Package sources
* [heavy_0.38.196.zip](https://cran.r-project.org/bin/windows/contrib/4.0/heavy_0.38.196.zip) - Windows binaries (R-release)
* [heavy_0.38.196.tgz](https://cran.r-project.org/bin/macosx/contrib/4.0/heavy_0.38.196.tgz) - Mac OS binaries (R-release)
* [heavy.pdf](https://cran.r-project.org/web/packages/heavy/heavy.pdf) - Reference Manual

## Installation instructions

To install this package, start R and enter:
```
install.packages("heavy")
```

Alternatively, you can download the source as a tarball or as a zip file. Unpack this file (thereby creating a directory named, heavy) and install the package source by executing (at the console prompt)
```
R CMD INSTALL heavy
```

Next, you can load the package by using the command: `library(heavy)`

## Disclaimer

The package is provided under the [GPL](https://www.r-project.org/Licenses/). HEAVY is under active development: new features are being added and old features are being improved. Although the developer will make efforts to preserve backward compatibility, we cannot absolutely guarantee backward compatibility.

Author: Felipe Osorio.

Project webpage: http://heavy.mat.utfsm.cl

Lastest binaries and sources for heavy are availables from CRAN at: https://cran.r-project.org/package=heavy
