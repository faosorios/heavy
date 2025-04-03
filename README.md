# heavy: Robust estimation using heavy-tailed distributions

[![CRAN status](http://www.r-pkg.org/badges/version/heavy)](https://cran.r-project.org/package=heavy)
[![CRAN RStudio mirror downloads](http://cranlogs.r-pkg.org/badges/heavy)](https://cran.r-project.org/package=heavy)

The `HEAVY` package contains routines to perform robust estimation considering heavy-tailed distributions. Currently, the package includes linear regression, linear mixed-effect models, Grubbs' model, multivariate location and scatter estimation, multivariate regression, penalized splines, random variate generation and some support functions.

## Features
* Provide basic functionality for modeling using scale mixtures of normal distributions in R, via a package.
* Calculations associated with parameter estimation are performed by calling routines in C and Fortran.
* Estimation in linear regression, linear mixed effects models, Grubbs' model, multivariate regression and penalized splines using the EM algorithm.
* Estimation of location and Scatter using multivariate heavy-tailed distributions.
* Implemented families: normal, Cauchy, Student-t, slash and contaminated normal.
* Estimation of the shape parameters for Student-t and slash distributions.
* Multivariate random number generation for the implemented families and the uniform distribution on the p-dimensional sphere.
* Print and summary methods and some sample databases.

## Reference Manual
* [heavy.pdf](https://cran.r-project.org/web/packages/heavy/heavy.pdf) - Reference Manual

## Resources
Latest binaries and sources can be found at the [CRAN package repository](https://cran.r-project.org/package=heavy)
* [heavy_0.38.196.tar.gz](https://cran.r-project.org/src/contrib/heavy_0.38.196.tar.gz) - Package sources
* [heavy_0.38.196.zip](https://cran.r-project.org/bin/windows/contrib/4.0/heavy_0.38.196.zip) - Windows binaries (R-release)
* [heavy_0.38.196.tgz](https://cran.r-project.org/bin/macosx/contrib/4.0/heavy_0.38.196.tgz) - Mac OS binaries (R-release)

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

### To cite the heavy package in publications use:

Osorio, F. (2019). heavy: Robust estimation using heavy-tailed distributions. R package version 0.38.196.\
URL: [CRAN.R-project.org/package=heavy](https://CRAN.R-project.org/package=heavy)

### Some papers using heavy:
* di San Miniato, M.L., Pagui, K.E.C. (2024). Reliable simulation of extremely-truncated log-concave distributions. [Journal of Statistical Computation and Simulation](https://doi.org/10.1080/00949655.2024.2406954) 94, 3933-3956.
* Davie, S., Minto, C., Officer, R., Lordan, C. (2015). Defining value per unit effort in mixed métier fisheries. [Fisheries Research](https://doi.org/10.1016/j.fishres.2014.12.007) 165, 1-10.
* Osorio, F. (2016). Influence diagnostics for robust P-splines using scale mixture of normal distributions. [Annals of the Institute of Statistical Mathematics](https://doi.org/10.1007/s10463-015-0506-0) 68, 589-619.
* Singer, J.M., Rocha, F.M.M., Nobre, J.S. (2016). Graphical tools for detecting departures from linear mixed model assumptions and some remedial measures. [International Statistical Review](https://doi.org/10.1111/insr.12178) 85, 290-324.

## Providing Feedback

Please report any bugs/suggestions/improvements to [Felipe Osorio](https://faosorios.github.io/). If you find these routines useful or not then please let me know. Also, acknowledgement of the use of the routines is appreciated.

## About the Author

Felipe Osorio is an applied statistician and creator of several R packages. Webpage: [https://faosorios.github.io/](https://faosorios.github.io/)
